//! Synthetic prepared-reader fixtures: generated FASTA/FAI/manifest bytes only.
//! No downloads, BLAST tools, environment changes or measured reference data.
//! `fixture` recreates both FASTA and FAI with explicit LF/CRLF line endings.

use super::*;

pub(crate) struct ReaderFixture {
    pub(crate) catalog: GenomeCatalog,
    pub(crate) manifest: PathBuf,
    pub(crate) fasta: PathBuf,
    pub(crate) index: PathBuf,
    pub(crate) annotation: PathBuf,
    pub(crate) cache_dir: PathBuf,
}

pub(crate) fn fixture(root: &Path, newline: &str) -> ReaderFixture {
    let root = root.join("reference with 'quote' and \u{e4}");
    fs::create_dir_all(&root).unwrap();
    let fasta = root.join("sequence.fa");
    let annotation = root.join("annotation.gtf");
    let index = root.join("sequence.fa.fai");
    fs::write(
        &fasta,
        ">chr1\nACGTACGT\nACGTACGT\n>NC_000099.1\nTTTTAACC\n".replace('\n', newline),
    )
    .unwrap();
    fs::write(
        &annotation,
        "# synthetic placeholder; no annotation parsing required\n",
    )
    .unwrap();
    build_fasta_index_with_progress(&fasta, &index, |_, _| true).unwrap();
    if newline == "\r\n" {
        let text = fs::read_to_string(&index).unwrap();
        fs::write(&index, text.replace('\n', newline)).unwrap();
    }
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    fs::write(
        &catalog_path,
        serde_json::to_vec(&serde_json::json!({
            "ReaderFixture": {"sequence_local":fasta, "annotations_local":annotation,
                "cache_dir":cache_dir}
        }))
        .unwrap(),
    )
    .unwrap();
    let catalog = GenomeCatalog::from_json_file(catalog_path.to_str().unwrap()).unwrap();
    let install = catalog.install_dir(
        "ReaderFixture",
        catalog.entry("ReaderFixture").unwrap(),
        None,
    );
    fs::create_dir_all(&install).unwrap();
    let manifest = install.join("manifest.json");
    let record: GenomeInstallManifest = serde_json::from_value(serde_json::json!({
        "genome_id":"ReaderFixture", "sequence_source":fasta,
        "annotation_source":annotation, "sequence_path":fasta,
        "annotation_path":annotation, "fasta_index_path":index,
        "gene_index_path":null, "installed_at_unix_ms":0
    }))
    .unwrap();
    GenomeCatalog::write_manifest(&manifest, &record).unwrap();
    ReaderFixture {
        catalog,
        manifest,
        fasta,
        index,
        annotation,
        cache_dir,
    }
}

#[test]
fn prepared_sequence_reader_preserves_lf_crlf_windows_aliases_errors_and_handle_closure() {
    for newline in ["\n", "\r\n"] {
        let temp = tempfile::tempdir().unwrap();
        let f = fixture(temp.path(), newline);
        let mut reader = PreparedSequenceRegionReader::new(&f.catalog);
        for (subject, start, end) in [
            ("chr1", 1, 16),
            ("1", 4, 12),
            ("chr1", 0, 3),
            ("chr1", 9, 8),
            ("chr1", 1, 17),
            ("missing", 1, 4),
            ("NC_000099.1", 1, 8),
        ] {
            assert_eq!(
                reader.get_sequence_region_with_cache("ReaderFixture", subject, start, end, None),
                f.catalog.get_sequence_region_with_cache(
                    "ReaderFixture",
                    subject,
                    start,
                    end,
                    None
                ),
                "{subject}:{start}-{end}"
            );
        }
        assert_eq!(reader.index_load_count(), 1);
        reader.validate_unchanged().unwrap();
        // Keep the reader alive: no cached handle may prevent Windows rename/delete.
        for path in [&f.fasta, &f.index] {
            let moved = path.with_extension("moved");
            fs::rename(path, &moved).unwrap();
            fs::rename(&moved, path).unwrap();
        }
        reader.validate_unchanged().unwrap();
        fs::remove_file(&f.index).unwrap();
        fs::remove_file(&f.fasta).unwrap();
    }
}

#[test]
fn prepared_sequence_reader_preserves_compatible_resolution_and_ambiguity_refusal() {
    for ambiguous in [false, true] {
        let temp = tempfile::tempdir().unwrap();
        let mut f = fixture(temp.path(), "\n");
        let entry = f.catalog.entry("ReaderFixture").unwrap().clone();
        f.catalog
            .entries
            .insert("ReaderAlias".into(), entry.clone());
        if ambiguous {
            f.catalog.entries.insert("OtherPrepared".into(), entry);
            let install = f.catalog.install_dir(
                "OtherPrepared",
                f.catalog.entry("OtherPrepared").unwrap(),
                None,
            );
            fs::create_dir_all(&install).unwrap();
            let mut manifest = GenomeCatalog::load_manifest(&f.manifest).unwrap();
            manifest.genome_id = "OtherPrepared".into();
            GenomeCatalog::write_manifest(&install.join("manifest.json"), &manifest).unwrap();
        }
        let mut reader = PreparedSequenceRegionReader::new(&f.catalog);
        let actual = reader.get_sequence_region_with_cache("ReaderAlias", "1", 1, 4, None);
        assert_eq!(
            actual,
            f.catalog
                .get_sequence_region_with_cache("ReaderAlias", "1", 1, 4, None)
        );
        assert_eq!(actual.is_err(), ambiguous);
        assert!(
            f.catalog
                .get_exact_sequence_region_with_cache("ReaderAlias", "1", 1, 4, None)
                .is_err()
        );
        reader.validate_unchanged().unwrap();
    }
}

#[test]
fn prepared_sequence_reader_parses_once_and_never_reuses_another_interpretation() {
    let temp = tempfile::tempdir().unwrap();
    let f = fixture(temp.path(), "\n");
    let mut reader = PreparedSequenceRegionReader::new(&f.catalog);
    for _ in 0..200 {
        assert_eq!(
            reader
                .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 8, None)
                .unwrap(),
            "ACGTACGT"
        );
        assert_eq!(
            reader
                .get_sequence_region_with_cache("ReaderFixture", "NC_000099.1", 1, 8, None)
                .unwrap(),
            "TTTTAACC"
        );
    }
    assert_eq!(reader.index_load_count(), 1);
    reader.validate_unchanged().unwrap();
    let mut next = PreparedSequenceRegionReader::new(&f.catalog);
    assert_eq!(next.index_load_count(), 0);
    next.get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 4, None)
        .unwrap();
    assert_eq!(next.index_load_count(), 1);
    next.validate_unchanged().unwrap();
}

#[test]
fn prepared_sequence_reader_keeps_missing_corrupt_errors_and_retries_failed_loads() {
    for missing in [true, false] {
        let temp = tempfile::tempdir().unwrap();
        let f = fixture(temp.path(), "\n");
        let bytes = fs::read(&f.index).unwrap();
        if missing {
            fs::remove_file(&f.index).unwrap();
        } else {
            fs::write(&f.index, "malformed index\n").unwrap();
        }
        let mut reader = PreparedSequenceRegionReader::new(&f.catalog);
        assert_eq!(
            reader.get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 4, None),
            f.catalog
                .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 4, None)
        );
        assert_eq!(reader.index_load_count(), 0);
        fs::write(&f.index, bytes).unwrap();
        assert_eq!(
            reader
                .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 4, None)
                .unwrap(),
            "ACGT"
        );
        assert_eq!(reader.index_load_count(), 1);
        reader.validate_unchanged().unwrap();
    }
}

#[test]
fn prepared_sequence_reader_rechecks_manifest_annotation_and_sequence_existence() {
    for name in ["manifest", "annotation", "fasta"] {
        let temp = tempfile::tempdir().unwrap();
        let f = fixture(temp.path(), "\n");
        let mut reader = PreparedSequenceRegionReader::new(&f.catalog);
        reader
            .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 4, None)
            .unwrap();
        let path = match name {
            "manifest" => &f.manifest,
            "annotation" => &f.annotation,
            _ => &f.fasta,
        };
        fs::remove_file(path).unwrap();
        assert_eq!(
            reader.get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 4, None),
            f.catalog
                .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 4, None)
        );
        assert!(reader.validate_unchanged().is_err(), "{name}");
    }
}

#[test]
fn prepared_sequence_reader_final_hash_rejects_equal_length_and_restored_mtime() {
    let temp = tempfile::tempdir().unwrap();
    let f = fixture(temp.path(), "\n");
    let mut reader = PreparedSequenceRegionReader::new(&f.catalog);
    reader
        .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 4, None)
        .unwrap();
    let before = fs::metadata(&f.index).unwrap();
    let original = fs::read_to_string(&f.index).unwrap();
    let changed = original.replace("chr1\t16\t", "chr1\t15\t");
    fs::write(&f.index, changed).unwrap();
    File::options()
        .write(true)
        .open(&f.index)
        .unwrap()
        .set_times(fs::FileTimes::new().set_modified(before.modified().unwrap()))
        .unwrap();
    assert_eq!(fs::metadata(&f.index).unwrap().len(), before.len());
    assert_eq!(
        fs::metadata(&f.index).unwrap().modified().unwrap(),
        before.modified().unwrap()
    );
    assert!(
        reader
            .validate_unchanged()
            .unwrap_err()
            .contains("changed during interpretation")
    );
    fs::write(&f.index, original).unwrap();
    assert!(reader.validate_unchanged().is_err());
}

#[test]
fn prepared_sequence_reader_observed_index_change_cannot_be_laundered_by_restoring_it() {
    let temp = tempfile::tempdir().unwrap();
    let f = fixture(temp.path(), "\n");
    let mut reader = PreparedSequenceRegionReader::new(&f.catalog);
    reader
        .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 4, None)
        .unwrap();
    let bytes = fs::read(&f.index).unwrap();
    fs::write(&f.index, "chr1\t1\t6\t1\t2\n").unwrap();
    assert!(
        reader
            .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 4, None)
            .is_err()
    );
    fs::write(&f.index, bytes).unwrap();
    assert!(reader.validate_unchanged().is_err());
}

#[test]
fn prepared_sequence_reader_refuses_manifest_changes_after_the_last_window() {
    let temp = tempfile::tempdir().unwrap();
    let f = fixture(temp.path(), "\n");
    let mut reader = PreparedSequenceRegionReader::new(&f.catalog);
    reader
        .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 4, None)
        .unwrap();
    let mut manifest = fs::read(&f.manifest).unwrap();
    manifest.push(b'\n');
    fs::write(&f.manifest, manifest).unwrap();
    assert!(
        reader
            .validate_unchanged()
            .unwrap_err()
            .contains("changed during interpretation")
    );
}

#[test]
fn prepared_sequence_reader_keeps_cache_overrides_separate_and_reads_fasta_fresh() {
    let temp = tempfile::tempdir().unwrap();
    let a = fixture(&temp.path().join("first"), "\n");
    let b = fixture(&temp.path().join("second"), "\n");
    let fasta = fs::read_to_string(&b.fasta)
        .unwrap()
        .replace("ACGTACGT", "TTTTAAAA");
    fs::write(&b.fasta, fasta).unwrap();
    let mut reader = PreparedSequenceRegionReader::new(&a.catalog);
    assert_eq!(
        reader
            .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 8, None)
            .unwrap(),
        "ACGTACGT"
    );
    assert_eq!(
        reader
            .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 8, b.cache_dir.to_str())
            .unwrap(),
        "TTTTAAAA"
    );
    assert_eq!(
        reader
            .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 8, None)
            .unwrap(),
        "ACGTACGT"
    );
    assert_eq!(reader.index_load_count(), 2);
    let fasta = fs::read_to_string(&a.fasta)
        .unwrap()
        .replace("ACGTACGT", "CCCCGGGG");
    fs::write(&a.fasta, fasta).unwrap();
    assert_eq!(
        reader
            .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 8, None)
            .unwrap(),
        a.catalog
            .get_sequence_region_with_cache("ReaderFixture", "chr1", 1, 8, None)
            .unwrap()
    );
    reader.validate_unchanged().unwrap();
}
