//! Shared direct-child supervision with capped capture or staged-file stdout.
//! No shell construction, implicit executable discovery, publication or biological policy.

use std::fs::File;
use std::io::{Read, Write};
use std::process::{Command, Output, Stdio};
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use std::thread;
use std::time::{Duration, Instant};

#[derive(Clone, Copy)]
pub(crate) struct Limits {
    pub timeout: Duration,
    pub stdout_bytes: usize,
    pub stderr_bytes: usize,
    pub refuse_overflow: bool,
}

#[derive(Debug)]
pub(crate) struct Captured {
    pub output: Output,
    pub stdout_truncated: bool,
    pub stderr_truncated: bool,
    pub stdout_bytes: u64,
    pub stderr_bytes: u64,
}

struct Drained {
    kept: Vec<u8>,
    truncated: bool,
    bytes: u64,
}

fn drain(
    mut source: impl Read,
    mut file: Option<File>,
    cap: usize,
    overflow: Arc<AtomicBool>,
    failed: Arc<AtomicBool>,
    done: Arc<AtomicBool>,
) -> std::io::Result<Drained> {
    let result = (|| {
        let mut kept = vec![];
        let mut buffer = [0u8; 16 * 1024];
        let (mut retained, mut bytes, mut truncated) = (0usize, 0u64, false);
        loop {
            let count = source.read(&mut buffer)?;
            if count == 0 {
                break;
            }
            bytes = bytes
                .checked_add(count as u64)
                .ok_or_else(|| std::io::Error::other("process byte count overflow"))?;
            let take = cap.saturating_sub(retained).min(count);
            if let Some(writer) = &mut file {
                writer.write_all(&buffer[..take])?;
            } else {
                kept.extend_from_slice(&buffer[..take]);
            }
            retained += take;
            if take < count {
                truncated = true;
                overflow.store(true, Ordering::Release);
            }
        }
        if let Some(writer) = &mut file {
            writer.flush()?;
            writer.sync_all()?;
        }
        Ok(Drained {
            kept,
            truncated,
            bytes,
        })
    })();
    failed.store(result.is_err(), Ordering::Release);
    done.store(true, Ordering::Release);
    result
}

/// File stdout stays in fixed-size buffers and closes before returning. On error
/// the caller owns/discards its private staging file; nothing is published here.
pub(crate) fn run(
    command: &mut Command,
    label: &str,
    limits: Limits,
    stdout_file: Option<File>,
    should_cancel: &mut dyn FnMut() -> bool,
) -> Result<Captured, String> {
    if should_cancel() {
        return Err(format!("{label} cancelled before spawn"));
    }
    if limits.timeout.is_zero() || limits.stdout_bytes == 0 || limits.stderr_bytes == 0 {
        return Err(format!("{label} invalid process limits"));
    }
    command
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = command
        .spawn()
        .map_err(|e| format!("Could not start {label}: {e}"))?;
    let (Some(stdout), Some(stderr)) = (child.stdout.take(), child.stderr.take()) else {
        let _ = child.kill();
        let _ = child.wait();
        return Err(format!("Could not capture {label} pipes"));
    };
    let overflow = Arc::new(AtomicBool::new(false));
    let out_failed = Arc::new(AtomicBool::new(false));
    let err_failed = Arc::new(AtomicBool::new(false));
    let out_done = Arc::new(AtomicBool::new(false));
    let err_done = Arc::new(AtomicBool::new(false));
    let (a, b, c) = (overflow.clone(), out_failed.clone(), out_done.clone());
    let out_reader = match thread::Builder::new()
        .name("gentle-process-stdout".into())
        .spawn(move || drain(stdout, stdout_file, limits.stdout_bytes, a, b, c))
    {
        Ok(reader) => reader,
        Err(error) => {
            let _ = child.kill();
            let _ = child.wait();
            return Err(format!("Could not capture {label}: {error}"));
        }
    };
    let (a, b, c) = (overflow.clone(), err_failed.clone(), err_done.clone());
    let err_reader = match thread::Builder::new()
        .name("gentle-process-stderr".into())
        .spawn(move || drain(stderr, None, limits.stderr_bytes, a, b, c))
    {
        Ok(reader) => reader,
        Err(error) => {
            let _ = child.kill();
            let _ = child.wait();
            let _ = out_reader.join();
            return Err(format!("Could not capture {label}: {error}"));
        }
    };
    let started = Instant::now();
    let mut status = None;
    let failure = loop {
        if should_cancel() {
            break Some(format!("{label} cancelled during probe"));
        }
        if limits.refuse_overflow && overflow.load(Ordering::Acquire) {
            break Some(format!("{label} output byte limit exceeded"));
        }
        if out_failed.load(Ordering::Acquire) || err_failed.load(Ordering::Acquire) {
            break Some(format!("Could not capture {label} output"));
        }
        if status.is_none() {
            match child.try_wait() {
                Ok(value) => status = value,
                Err(error) => break Some(format!("Could not poll {label}: {error}")),
            }
        }
        if status.is_some() && out_done.load(Ordering::Acquire) && err_done.load(Ordering::Acquire)
        {
            break None;
        }
        if started.elapsed() >= limits.timeout {
            break Some(format!(
                "{label} exceeded the {} second timeout",
                limits.timeout.as_secs()
            ));
        }
        thread::sleep(Duration::from_millis(10));
    };
    if failure.is_some() {
        let _ = child.kill();
        let _ = child.wait();
    }
    let stdout = out_reader
        .join()
        .map_err(|_| format!("{label} stdout worker panicked"))?
        .map_err(|e| format!("{label} stdout read/write failed: {e}"))?;
    let stderr = err_reader
        .join()
        .map_err(|_| format!("{label} stderr worker panicked"))?
        .map_err(|e| format!("{label} stderr read failed: {e}"))?;
    if let Some(error) = failure {
        return Err(error);
    }
    if limits.refuse_overflow && (stdout.truncated || stderr.truncated) {
        return Err(format!("{label} output byte limit exceeded"));
    }
    Ok(Captured {
        output: Output {
            status: status.ok_or_else(|| format!("{label} has no exit status"))?,
            stdout: stdout.kept,
            stderr: stderr.kept,
        },
        stdout_truncated: stdout.truncated,
        stderr_truncated: stderr.truncated,
        stdout_bytes: stdout.bytes,
        stderr_bytes: stderr.bytes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stub() -> (tempfile::TempDir, std::path::PathBuf) {
        let root = tempfile::tempdir().unwrap();
        let name = if cfg!(windows) {
            "synthetic process's stub.exe"
        } else {
            "synthetic process's stub"
        };
        let binary = root.path().join(name);
        let source = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("test_files/supervised_process/stub.rs");
        let compiler = std::env::var_os("RUSTC").unwrap_or_else(|| "rustc".into());
        let status = Command::new(compiler)
            .arg("--edition=2024")
            .arg(source)
            .arg("-o")
            .arg(&binary)
            .status()
            .unwrap();
        assert!(status.success());
        (root, binary)
    }

    fn limits() -> Limits {
        Limits {
            timeout: Duration::from_secs(5),
            stdout_bytes: 64,
            stderr_bytes: 64,
            refuse_overflow: true,
        }
    }

    #[test]
    fn supervised_process_file_output_is_closed_and_counted_before_return() {
        let (root, binary) = stub();
        let path = root.path().join("staged output's file");
        let file = File::create(&path).unwrap();
        let out = run(
            Command::new(binary).arg("ok"),
            "synthetic",
            limits(),
            Some(file),
            &mut || false,
        )
        .unwrap();
        assert_eq!(out.stdout_bytes, 8);
        assert!(out.output.stdout.is_empty());
        assert_eq!(std::fs::read(&path).unwrap(), b"ABCDEFGH");
        std::fs::rename(&path, root.path().join("renamed")).unwrap();
    }

    #[test]
    fn supervised_process_caps_refuse_or_truncate_without_unbounded_capture() {
        let (_root, binary) = stub();
        let capped = run(
            Command::new(&binary).arg("flood"),
            "synthetic",
            limits(),
            None,
            &mut || false,
        )
        .unwrap_err();
        assert!(capped.contains("byte limit"));
        let out = run(
            Command::new(binary).arg("finite"),
            "synthetic",
            Limits {
                refuse_overflow: false,
                ..limits()
            },
            None,
            &mut || false,
        )
        .unwrap();
        assert_eq!(out.output.stdout.len(), 64);
        assert!(out.stdout_truncated);
        assert_eq!(out.stdout_bytes, 8192);
    }

    #[test]
    fn supervised_process_cancel_and_timeout_reap_before_closing_staging() {
        let (root, binary) = stub();
        let path = root.path().join("partial");
        let mut checks = 0;
        let error = run(
            Command::new(&binary).arg("wait"),
            "synthetic",
            limits(),
            Some(File::create(&path).unwrap()),
            &mut || {
                checks += 1;
                checks > 3
            },
        )
        .unwrap_err();
        assert!(error.contains("cancelled"));
        std::fs::remove_file(&path).unwrap();
        let error = run(
            Command::new(binary).arg("wait"),
            "synthetic",
            Limits {
                timeout: Duration::from_millis(50),
                ..limits()
            },
            None,
            &mut || false,
        )
        .unwrap_err();
        assert!(error.contains("timeout"));
    }

    #[test]
    fn supervised_process_pre_cancel_does_not_resolve_or_spawn_executable() {
        let root = tempfile::tempdir().unwrap();
        let error = run(
            &mut Command::new(root.path().join("absent")),
            "synthetic",
            limits(),
            None,
            &mut || true,
        )
        .unwrap_err();
        assert_eq!(error, "synthetic cancelled before spawn");
    }
}
