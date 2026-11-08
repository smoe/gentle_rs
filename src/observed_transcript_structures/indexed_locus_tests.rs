//! Synthetic native executable/selection controls. Real samtools is an explicit opt-in only.

use super::*;
fn fixture() -> (tempfile::TempDir, Request) {
    let root = tempfile::tempdir().unwrap();
    let path = root.path().canonicalize().unwrap();
    let (_, template, _) = sam_intake::tests::fixture();
    let executable = path.join(if cfg!(windows) {
        "synthetic samtool's.exe"
    } else {
        "synthetic samtool's"
    });
    let compiler = std::env::var_os("RUSTC").unwrap_or_else(|| "rustc".into());
    assert!(
        Command::new(compiler)
            .arg("--edition=2024")
            .arg(
                Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/supervised_process/stub.rs")
            )
            .arg("-o")
            .arg(&executable)
            .status()
            .unwrap()
            .success()
    );
    let sam=b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ctg_a\tLN:1000\n@SQ\tSN:ctg_b\tLN:1000\nsynthetic\t0\tctg_a\t101\t60\t30M\t*\t0\t0\t*\t*\n";
    let alignment = path.join("input's spaced.bam");
    fs::write(&alignment, sam).unwrap();
    let index = path.join("input's spaced.bai");
    fs::write(&index, b"synthetic index, not a BAI").unwrap();
    let dictionary = path.join("dictionary.json");
    let value = ReferenceContigDictionary {
        schema: crate::engine::REFERENCE_CONTIG_DICTIONARY_SCHEMA.into(),
        reference_id: template.reference_id.clone(),
        assembly: template.reference.assembly.clone(),
        reference_fasta_sha256: template.reference.fasta_sha256.clone(),
        contigs: vec![
            crate::engine::ReferenceContigDictionaryEntry {
                name: "ctg_a".into(),
                length: 1000,
            },
            crate::engine::ReferenceContigDictionaryEntry {
                name: "ctg_b".into(),
                length: 1000,
            },
        ],
    };
    fs::write(&dictionary, serde_json::to_vec(&value).unwrap()).unwrap();
    let bound = |p: &Path| fingerprint(p, 1024 * 1024, &mut || Ok(())).unwrap();
    let alignment = bound(&alignment);
    let request = Request {
        schema: SCHEMA.into(),
        samtools: Some(executable.to_string_lossy().into_owned()),
        format: Format::Bam,
        index_for_alignment_sha256: alignment.sha256.clone(),
        alignment,
        index: bound(&index),
        dictionary: bound(&dictionary),
        reference: template.reference,
        reference_id: template.reference_id,
        reference_fasta: None,
        reference_fai: None,
        contig: "ctg_a".into(),
        interval: (100, 150),
        minimum_mapq: None,
        require_flags: None,
        exclude_flags: None,
        verify_full_alignment: false,
        timeout_ms: 2000,
        max_sam_bytes: 1024 * 1024,
        max_stderr_bytes: 1024,
        max_units: 100,
        text_limits: template.limits,
    };
    (root, request)
}

#[test]
fn indexed_locus_native_stub_binds_partial_units_without_claiming_full_hash() {
    let (_root, request) = fixture();
    let selected = select(&request, &mut || Ok(())).unwrap();
    assert!(!selected.receipt.alignment_hash_verified);
    assert!(!selected.receipt.search_complete && !selected.receipt.views_complete);
    let manifest: sam_intake::Manifest =
        serde_json::from_reader(File::open(&selected.manifest.path).unwrap()).unwrap();
    assert_eq!(manifest.units.len(), 1);
    assert!(!manifest.units[0].search.as_ref().unwrap().search_complete);
    let sam = selected.sam.path.clone();
    drop(selected);
    assert!(!Path::new(&sam).exists());
}

#[test]
fn indexed_locus_missing_wrong_and_stale_indexes_refuse() {
    let (_root, mut request) = fixture();
    let index = Path::new(&request.index.path).with_file_name("wrong-index.bai");
    fs::rename(&request.index.path, &index).unwrap();
    request.index.path = index.to_string_lossy().into_owned();
    assert!(
        select(&request, &mut || Ok(()))
            .err()
            .unwrap()
            .message
            .contains("tool_refused")
    );
    fs::remove_file(&index).unwrap();
    assert!(select(&request, &mut || Ok(())).is_err());
    fs::write(&index, b"synthetic index, not a BAI").unwrap();
    std::thread::sleep(Duration::from_millis(1100));
    let bytes = fs::read(&request.alignment.path).unwrap();
    fs::write(&request.alignment.path, bytes).unwrap();
    assert!(
        select(&request, &mut || Ok(()))
            .err()
            .unwrap()
            .message
            .contains("stale_index")
    );
}

#[test]
fn indexed_locus_caps_cancellation_and_tool_unavailable_fail_closed() {
    let (_root, mut request) = fixture();
    request.samtools = None;
    assert!(
        select(&request, &mut || Ok(()))
            .err()
            .unwrap()
            .message
            .contains("unavailable")
    );
    let (_root, mut request) = fixture();
    request.max_sam_bytes = 8;
    assert!(
        select(&request, &mut || Ok(()))
            .err()
            .unwrap()
            .message
            .contains("byte limit")
    );
    request.max_sam_bytes = 1024 * 1024;
    let alignment = Path::new(&request.alignment.path).with_file_name("wait.bam");
    fs::rename(&request.alignment.path, &alignment).unwrap();
    request.alignment.path = alignment.to_string_lossy().into_owned();
    let start = std::time::Instant::now();
    let mut ticks = 0;
    let error = select(&request, &mut || {
        ticks += 1;
        if ticks > 40 {
            Err(fail("cancelled"))
        } else {
            Ok(())
        }
    })
    .err()
    .unwrap();
    assert!(error.message.contains("cancel"));
    assert!(start.elapsed() < Duration::from_secs(10));
    // Handles are closed: native Windows can rename/delete all input files after refusal.
    fs::rename(&request.alignment.path, alignment.with_extension("closed")).unwrap();
}

#[test]
fn indexed_locus_explicit_full_hash_and_cram_reference_refusals_are_separate() {
    let (_root, mut request) = fixture();
    request.verify_full_alignment = true;
    request.alignment.sha256 = format!("sha256:{}", "0".repeat(64));
    request.index_for_alignment_sha256 = request.alignment.sha256.clone();
    assert!(
        select(&request, &mut || Ok(()))
            .err()
            .unwrap()
            .message
            .contains("hash_binding")
    );
    request.format = Format::Cram;
    assert!(
        select(&request, &mut || Ok(()))
            .err()
            .unwrap()
            .message
            .contains("reference_required")
    );
}

#[test]
fn indexed_locus_cram_forbids_header_fallback_and_binds_fasta_fai() {
    let (root, mut request) = fixture();
    request.format = Format::Cram;
    let fasta = root.path().canonicalize().unwrap().join("reference.fa");
    fs::write(
        &fasta,
        format!(
            ">ctg_a\n{}\n>ctg_b\n{}\n",
            "A".repeat(1000),
            "C".repeat(1000)
        ),
    )
    .unwrap();
    fs::write(
        fasta.with_extension("fa.fai"),
        b"ctg_a\t1000\t7\t1000\t1001\nctg_b\t1000\t1015\t1000\t1001\n",
    )
    .unwrap();
    request.reference_fasta = Some(fingerprint(&fasta, 1024, &mut || Ok(())).unwrap());
    request.reference_fai =
        Some(fingerprint(&fasta.with_extension("fa.fai"), 1024, &mut || Ok(())).unwrap());
    assert!(
        select(&request, &mut || Ok(()))
            .err()
            .unwrap()
            .message
            .contains("reference_binding")
    );
    request.reference.fasta_sha256 = request.reference_fasta.as_ref().unwrap().sha256.clone();
    let mut dictionary: ReferenceContigDictionary =
        serde_json::from_reader(File::open(&request.dictionary.path).unwrap()).unwrap();
    dictionary.reference_fasta_sha256 = request.reference.fasta_sha256.clone();
    fs::write(
        &request.dictionary.path,
        serde_json::to_vec(&dictionary).unwrap(),
    )
    .unwrap();
    request.dictionary = fingerprint(
        Path::new(&request.dictionary.path),
        1024 * 1024,
        &mut || Ok(()),
    )
    .unwrap();
    let mut bytes = fs::read_to_string(&request.alignment.path).unwrap();
    bytes = bytes.replace(
        "SN:ctg_a\tLN:1000",
        "SN:ctg_a\tLN:1000\tUR:https://invalid.example/reference.fa",
    );
    fs::write(&request.alignment.path, bytes).unwrap();
    request.alignment = fingerprint(Path::new(&request.alignment.path), 1024 * 1024, &mut || {
        Ok(())
    })
    .unwrap();
    request.index_for_alignment_sha256 = request.alignment.sha256.clone();
    fs::write(&request.index.path, b"synthetic index, not a BAI").unwrap();
    assert!(
        select(&request, &mut || Ok(()))
            .err()
            .unwrap()
            .message
            .contains("fallback_refused")
    );
}

#[test]
fn indexed_locus_optional_real_samtools_requires_explicit_configuration() {
    let Some(executable) = std::env::var_os("GENTLE_TEST_SAMTOOLS") else {
        eprintln!("NOT RUN: GENTLE_TEST_SAMTOOLS not configured");
        return;
    };
    let (root, mut request) = fixture();
    let path = root.path().canonicalize().unwrap();
    let sam = path.join("synthetic.sam");
    fs::rename(&request.alignment.path, &sam).unwrap();
    let bam = path.join("actual.bam");
    assert!(
        Command::new(&executable)
            .args(["view", "-b", "-o"])
            .arg(&bam)
            .arg(&sam)
            .status()
            .unwrap()
            .success()
    );
    assert!(
        Command::new(&executable)
            .arg("index")
            .arg(&bam)
            .status()
            .unwrap()
            .success()
    );
    request.samtools = Some(PathBuf::from(executable).to_str().unwrap().into());
    request.alignment = fingerprint(&bam, 1024 * 1024, &mut || Ok(())).unwrap();
    request.index =
        fingerprint(&bam.with_extension("bam.bai"), 1024 * 1024, &mut || Ok(())).unwrap();
    request.index_for_alignment_sha256 = request.alignment.sha256.clone();
    let selected = select(&request, &mut || Ok(())).unwrap();
    assert!(selected.sam.bytes > 0);
    drop(selected);
    let fasta = path.join("synthetic.fa");
    fs::write(
        &fasta,
        format!(
            ">ctg_a\n{}\n>ctg_b\n{}\n",
            "A".repeat(1000),
            "C".repeat(1000)
        ),
    )
    .unwrap();
    let executable = request.samtools.as_ref().unwrap();
    assert!(
        Command::new(executable)
            .arg("faidx")
            .arg(&fasta)
            .status()
            .unwrap()
            .success()
    );
    let cram = path.join("actual.cram");
    assert!(
        Command::new(executable)
            .args(["view", "-C", "-T"])
            .arg(&fasta)
            .arg("-o")
            .arg(&cram)
            .arg(&sam)
            .status()
            .unwrap()
            .success()
    );
    assert!(
        Command::new(executable)
            .arg("index")
            .arg(&cram)
            .status()
            .unwrap()
            .success()
    );
    request.format = Format::Cram;
    request.reference_fasta = Some(fingerprint(&fasta, 1024 * 1024, &mut || Ok(())).unwrap());
    request.reference_fai =
        Some(fingerprint(&fasta.with_extension("fa.fai"), 1024 * 1024, &mut || Ok(())).unwrap());
    request.reference.fasta_sha256 = request.reference_fasta.as_ref().unwrap().sha256.clone();
    let mut dictionary: ReferenceContigDictionary =
        serde_json::from_reader(File::open(&request.dictionary.path).unwrap()).unwrap();
    dictionary.reference_fasta_sha256 = request.reference.fasta_sha256.clone();
    fs::write(
        &request.dictionary.path,
        serde_json::to_vec(&dictionary).unwrap(),
    )
    .unwrap();
    request.dictionary = fingerprint(
        Path::new(&request.dictionary.path),
        1024 * 1024,
        &mut || Ok(()),
    )
    .unwrap();
    request.alignment = fingerprint(&cram, 1024 * 1024, &mut || Ok(())).unwrap();
    request.index = fingerprint(&cram.with_extension("cram.crai"), 1024 * 1024, &mut || {
        Ok(())
    })
    .unwrap();
    request.index_for_alignment_sha256 = request.alignment.sha256.clone();
    assert!(select(&request, &mut || Ok(())).is_ok());
}

#[cfg(unix)]
#[test]
fn indexed_locus_symlink_inputs_refuse_before_tool_execution() {
    let (root, mut request) = fixture();
    let link = root.path().canonicalize().unwrap().join("aliased.bam");
    std::os::unix::fs::symlink(&request.alignment.path, &link).unwrap();
    request.alignment.path = link.to_string_lossy().into_owned();
    assert!(
        select(&request, &mut || Ok(()))
            .err()
            .unwrap()
            .message
            .contains("symlink")
    );
}
