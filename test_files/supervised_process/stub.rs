//! Synthetic standard-library-only subprocess fixture; not a biological executable.

use std::io::Write;

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_default();
    let mut stdout = std::io::stdout().lock();
    match mode.as_str() {
        "--version" => {
            stdout
                .write_all(b"samtools 1.25 synthetic\nUsing htslib 1.25 synthetic\n")
                .unwrap();
        }
        "view" => {
            let args: Vec<_> = std::env::args().skip(2).collect();
            for key in ["REF_PATH", "REF_CACHE"] {
                let value = std::env::var(key).unwrap_or_default();
                if value.is_empty() || value.contains("://") {
                    eprintln!("synthetic forbidden reference fallback");
                    std::process::exit(3);
                }
            }
            let header = args.iter().any(|a| a == "-H");
            let position = if header {
                args.len() - 1
            } else {
                args.len() - 2
            };
            let path = std::path::Path::new(&args[position]);
            let bytes = std::fs::read(path).unwrap();
            if header {
                for line in bytes
                    .split_inclusive(|b| *b == b'\n')
                    .take_while(|line| line.starts_with(b"@"))
                {
                    stdout.write_all(line).unwrap();
                }
            } else {
                if args.last().unwrap().contains("wrong-index") {
                    eprintln!("[W::hts_idx_load3] index warning");
                    std::process::exit(1);
                }
                if path.file_name().unwrap().to_string_lossy().contains("wait") {
                    stdout.write_all(b"partial").unwrap();
                    stdout.flush().unwrap();
                    std::thread::sleep(std::time::Duration::from_secs(30));
                } else {
                    stdout.write_all(&bytes).unwrap();
                }
            }
        }
        "ok" => {
            stdout.write_all(b"ABCDEFGH").unwrap();
        }
        "finite" => {
            stdout.write_all(&[b'X'; 8192]).unwrap();
        }
        "flood" => while stdout.write_all(&[b'X'; 8192]).is_ok() {},
        "wait" => {
            stdout.write_all(b"partial").unwrap();
            stdout.flush().unwrap();
            std::thread::sleep(std::time::Duration::from_secs(30));
        }
        _ => std::process::exit(2),
    }
}
