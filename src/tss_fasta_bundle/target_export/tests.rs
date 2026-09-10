//! Hand-crafted synthetic target-export v1 inputs, recreated by `Fixture::new`.
//! No biological source or external file is required by the ordinary reader tests.
//! The ignored acceptance test only reads user-materialized files from the stated
//! Git object; it never downloads, extracts, modifies, or commits those inputs.

use super::*;
use serde_json::{Value, json};
use tempfile::TempDir;

const GENOME: &str = "Synthetic Genome Release 9";
const ASSEMBLY: &str = "SYN9";
const DATASET: &str = "synthetic_promoterome_evidence_v1";

struct Fixture {
    temp: TempDir,
    root: PathBuf,
    manifest: Value,
    selection: Value,
}

impl Fixture {
    fn new() -> Self {
        let temp = tempfile::tempdir().unwrap();
        let root = temp.path().join("bundle");
        fs::create_dir(&root).unwrap();
        let mut members = Vec::new();
        let mut regions = Vec::new();
        for (id, gene, strand, tss, start, end, sequence) in [
            ("p.plus", "SYN_PLUS", "+", 100, 97, 101, "ACGTN"),
            ("p.minus", "SYN_MINUS", "-", 200, 199, 203, "TGNAC"),
        ] {
            let digest = sha256_hex_bytes(sequence.as_bytes());
            let transcripts = vec![format!("{id}.tx1"), format!("{id}.tx2")];
            let filename = format!("{id}.fasta");
            let header = format!(
                ">{gene}|gene_id={gene}.id|promoter_id={id}|assembly={ASSEMBLY}|chromosome=synthetic_chr|strand={strand}|tss_1based={tss}|genomic_1based={start}-{end}|window=minus3_plus1|orientation={ORIENTATION}|transcripts={}|sequence_sha256={digest}\n{sequence}\n",
                transcripts.join(",")
            );
            fs::write(root.join(&filename), header).unwrap();
            members.push(json!({
                "filename":filename, "gene_id":format!("{gene}.id"), "gene_symbol":gene,
                "record_count":1, "sha256":"",
                "records":[{
                    "promoter_id":id, "gene_id":format!("{gene}.id"),
                    "chromosome":"synthetic_chr", "strand":strand, "tss_1based":tss,
                    "genomic_start_1based":start, "genomic_end_1based":end,
                    "sequence_length_bp":5, "sequence_sha256":digest,
                    "transcript_ids":transcripts,
                }],
            }));
            let selected = if strand == "+" {
                vec![transcripts[0].clone()]
            } else {
                transcripts
            };
            regions.push(json!({
                "promoterome_id":id, "input_kind":"transcript_promoter",
                "region_id":format!("unrelated-report-name-{id}"), "gene_query":"not-a-join-key",
                "transcript_ids":selected,
                "genome_extraction": {
                    "gene_id":format!("{gene}.id"), "gene_name":gene, "genome_id":GENOME,
                    "chromosome":"synthetic_chr", "strand":strand, "tss_1based":tss,
                    "start_1based":1, "end_1based":999, "transcript_ids":selected,
                    "promoter_upstream_bp":700, "promoter_downstream_bp":298,
                },
                "sequence_length_bp":999, "sequence_sha256":format!("sha256:{}", "e".repeat(64)),
                "sequence_orientation":"biological_5prime_to_3prime",
                "cutrun_support": {
                    "factor":"SyntheticFactor", "criterion":"Synthetic experimental window mean exceeds matched control",
                    "non_claim":"Synthetic descriptive evidence only; no significance or activity claim",
                    "evidence":[{"mean":-109.91133953680091}],
                },
            }));
        }
        let mut fixture = Self {
            temp,
            root,
            manifest: json!({
                "schema":SCHEMA, "assembly_id":ASSEMBLY,
                "upstream_bp":3, "downstream_bp":1, "sequence_orientation":ORIENTATION,
                "files":members, "sha256sums_sha256":"",
                "source": {
                    "dataset_id":DATASET, "genome_id":GENOME,
                    "promoter_fasta_sha256":"a".repeat(64),
                    "promoter_transcripts_sha256":"b".repeat(64),
                    "promoter_windows_sha256":"c".repeat(64),
                    "promoterome_receipt_sha256":"d".repeat(64),
                },
                "source_revision":"1".repeat(40), "producer_sha256":"2".repeat(64),
                "record_policy":"One synthetic record per distinct physical TSS",
                "non_claims":["Synthetic fixture, not genomic evidence"],
            }),
            selection: json!({"schema":INTEGRATED_SELECTION_SCHEMA, "regions":regions}),
        };
        fixture.rebind();
        fixture.save_selection();
        fixture
    }

    fn request(&self) -> ComputeTssProfilesRequest {
        ComputeTssProfilesRequest {
            manifest: self
                .root
                .join("manifest.json")
                .to_string_lossy()
                .into_owned(),
            panel: "unused-by-reader.json".into(),
            fasta: Vec::new(),
            selection: Some(self.selection_path().to_string_lossy().into_owned()),
            expected_genome_id: GENOME.into(),
            expected_assembly: Some(ASSEMBLY.into()),
            expected_annotation_release: None,
            expected_dataset_id: Some(DATASET.into()),
        }
    }

    fn selection_path(&self) -> PathBuf {
        self.temp.path().join("selected.json")
    }

    fn save_manifest(&self) {
        fs::write(
            self.root.join("manifest.json"),
            serde_json::to_vec_pretty(&self.manifest).unwrap(),
        )
        .unwrap();
    }

    fn save_selection(&self) {
        fs::write(
            self.selection_path(),
            serde_json::to_vec_pretty(&self.selection).unwrap(),
        )
        .unwrap();
    }

    fn rebind(&mut self) {
        let mut sums = String::new();
        for file in self.manifest["files"].as_array_mut().unwrap() {
            let name = file["filename"].as_str().unwrap().to_owned();
            let digest = sha256_hex_bytes(&fs::read(self.root.join(&name)).unwrap());
            file["sha256"] = json!(digest);
            sums.push_str(&format!("{digest}  {name}\n"));
        }
        fs::write(self.root.join(CHECKSUMS_NAME), &sums).unwrap();
        self.manifest["sha256sums_sha256"] = json!(sha256_hex_bytes(sums.as_bytes()));
        self.save_manifest();
    }

    fn replace_fasta(&self, name: &str, before: &str, after: &str) {
        let path = self.root.join(name);
        let text = fs::read_to_string(&path).unwrap();
        assert!(text.contains(before), "missing replacement {before}");
        fs::write(path, text.replace(before, after)).unwrap();
    }
}

fn invalid(request: &ComputeTssProfilesRequest, diagnostic: &str) {
    let error = read_bundle(request).unwrap_err();
    assert!(
        format!("{error:?}").contains(diagnostic),
        "expected {diagnostic}: {error:?}"
    );
}

#[test]
fn target_roundtrip_preserves_nondefault_plus_minus_and_actual_source_hashes() {
    let fixture = Fixture::new();
    let before = fs::read(fixture.root.join("p.plus.fasta")).unwrap();
    let bundle = read_bundle(&fixture.request()).unwrap();
    assert_eq!(bundle.records.len(), 2);
    assert_eq!(bundle.records[0].0.promoter_id, "p.plus");
    assert_eq!(bundle.records[1].0.promoter_id, "p.minus");
    for (record, sequence, selected) in &bundle.records {
        assert_eq!(sequence.len(), 5);
        assert_eq!(record.geometry.upstream_bp, 3);
        assert_eq!(record.geometry.downstream_bp, 1);
        assert_eq!(
            record.geometry.genomic_at(3),
            Some(record.geometry.tss_1based)
        );
        assert_eq!(record.transcripts.len(), 2);
        assert!(*selected);
    }
    assert_eq!(bundle.records[1].1, "TGNAC"); // Already transcript-oriented; no second reverse complement.
    assert_eq!(bundle.reference.annotation_release, None);
    assert_eq!(bundle.source.schema, SCHEMA);
    assert_eq!(bundle.source.dataset_id.as_deref(), Some(DATASET));
    assert_eq!(bundle.source.source_revision, Some("1".repeat(40)));
    assert_eq!(bundle.source.producer_sha256, Some("2".repeat(64)));
    assert_eq!(
        bundle.source.manifest_sha256,
        sha256_hex_bytes(&fs::read(fixture.root.join("manifest.json")).unwrap())
    );
    assert_eq!(
        bundle
            .inputs
            .iter()
            .filter(|b| b.role == "declared_source")
            .count(),
        4
    );
    let selected = bundle
        .inputs
        .iter()
        .find(|b| b.role == "selection")
        .unwrap();
    assert_eq!(selected.name, "selected.json");
    assert_eq!(
        selected.sha256,
        sha256_hex_bytes(&fs::read(fixture.selection_path()).unwrap())
    );
    assert_eq!(bundle.selection_evidence.len(), 2);
    assert_eq!(fs::read(fixture.root.join("p.plus.fasta")).unwrap(), before);
}

#[test]
fn target_requires_exact_reference_dataset_and_declared_annotation() {
    let fixture = Fixture::new();
    for value in [
        "Synthetic Genome Release 8",
        "Release 9",
        "9",
        "SYN9",
        "Synthetic Genome Release 9 ",
    ] {
        let mut request = fixture.request();
        request.expected_genome_id = value.into();
        invalid(&request, "exactly match");
    }
    for case in 0..6 {
        let mut request = fixture.request();
        match case {
            0 => request.expected_assembly = None,
            1 => request.expected_dataset_id = None,
            2 => request.expected_assembly = Some("syn9".into()),
            3 => request.expected_dataset_id = Some(format!("{DATASET}x")),
            4 => request.expected_dataset_id = Some("evidence_v1".into()),
            _ => request.expected_annotation_release = Some("9".into()),
        }
        assert!(read_bundle(&request).is_err(), "case {case}");
    }
    for (field, value) in [
        ("genome_id", "Release 9"),
        ("dataset_id", "synthetic_promoterome_evidence_v2"),
    ] {
        let mut fixture = Fixture::new();
        fixture.manifest["source"][field] = json!(value);
        fixture.save_manifest();
        invalid(&fixture.request(), "exactly match");
    }
}

#[test]
fn target_checksum_chain_is_manifest_to_sums_to_fasta_not_circular() {
    for case in 0..6 {
        let mut fixture = Fixture::new();
        match case {
            0 => {
                fixture.manifest["sha256sums_sha256"] = json!("0".repeat(64));
                fixture.save_manifest();
            }
            1 => {
                let path = fixture.root.join(CHECKSUMS_NAME);
                let mut sums = fs::read_to_string(&path).unwrap();
                sums.push('\n');
                fs::write(path, sums).unwrap();
            }
            2 => fixture.replace_fasta("p.plus.fasta", "\nACGTN\n", "\nTCGTN\n"),
            3 => {
                fixture.replace_fasta("p.plus.fasta", "\nACGTN\n", "\nTCGTN\n");
                fixture.manifest["files"][0]["sha256"] = json!(sha256_hex_bytes(
                    &fs::read(fixture.root.join("p.plus.fasta")).unwrap()
                ));
                fixture.save_manifest(); // Manifest agrees with FASTA, checksum list does not.
            }
            4 => {
                fixture.replace_fasta("p.plus.fasta", "\nACGTN\n", "\nTCGTN\n");
                let previous = fixture.manifest["files"][0]["sha256"].clone();
                fixture.rebind();
                fixture.manifest["files"][0]["sha256"] = previous;
                fixture.save_manifest(); // Checksum list agrees with FASTA, manifest does not.
            }
            _ => {
                let path = fixture.root.join(CHECKSUMS_NAME);
                let mut sums = fs::read_to_string(&path).unwrap();
                sums.push_str(&format!("{}  manifest.json\n", "0".repeat(64)));
                fs::write(path, &sums).unwrap();
                fixture.manifest["sha256sums_sha256"] = json!(sha256_hex_bytes(sums.as_bytes()));
                fixture.save_manifest();
            }
        }
        invalid(
            &fixture.request(),
            if case < 2 {
                "sha256sums_sha256"
            } else {
                "SHA256SUMS"
            },
        );
    }
}

#[test]
fn target_manifest_geometry_and_lengths_are_recomputed_for_both_strands() {
    for file in 0..2 {
        for field in [
            "tss_1based",
            "genomic_start_1based",
            "genomic_end_1based",
            "sequence_length_bp",
        ] {
            let mut fixture = Fixture::new();
            let old = fixture.manifest["files"][file]["records"][0][field]
                .as_u64()
                .unwrap();
            fixture.manifest["files"][file]["records"][0][field] = json!(old + 1);
            fixture.save_manifest();
            assert!(read_bundle(&fixture.request()).is_err(), "{file}/{field}");
        }
    }
    let mut fixture = Fixture::new();
    fixture.manifest["upstream_bp"] = json!(1);
    fixture.manifest["downstream_bp"] = json!(3);
    fixture.save_manifest();
    invalid(&fixture.request(), "geometry");
}

#[test]
fn target_headers_use_only_the_versioned_grammar_and_exact_keys() {
    for (before, after) in [
        (">SYN_PLUS|", ">gene_symbol=SYN_PLUS|"),
        ("genomic_1based=97-101", "genomic_1based=97..101"),
        ("window=minus3_plus1", "window=-3..+1"),
        ("window=minus3_plus1", "window=minus03_plus1"),
        (
            "gene_id=SYN_PLUS.id",
            "gene_id=SYN_PLUS.id|gene_id=SYN_PLUS.id",
        ),
        ("gene_id=SYN_PLUS.id", "gene_id=SYN_PLUS.id|unknown=x"),
        ("gene_id=SYN_PLUS.id|", ""),
        ("tss_1based=100", "tss_1based=101"),
        ("assembly=SYN9", "assembly=syn9"),
        ("transcript_5prime_to_3prime", "genomic_5prime_to_3prime"),
        (
            "transcripts=p.plus.tx1,p.plus.tx2",
            "transcripts=p.plus.tx1",
        ),
        (
            "transcripts=p.plus.tx1,p.plus.tx2",
            "transcripts=p.plus.tx1,p.plus.tx1",
        ),
        ("sequence_sha256=", "sequence_sha256=sha256:"),
    ] {
        let mut fixture = Fixture::new();
        fixture.replace_fasta("p.plus.fasta", before, after);
        fixture.rebind();
        assert!(
            read_bundle(&fixture.request()).is_err(),
            "{before} -> {after}"
        );
    }
}

#[test]
fn target_sequences_normalize_only_ascii_whitespace_and_acgtn_case() {
    let mut fixture = Fixture::new();
    fixture.replace_fasta("p.plus.fasta", "\nACGTN\n", "\na c\tg\nTn\n");
    fixture.replace_fasta(
        "p.plus.fasta",
        "p.plus.tx1,p.plus.tx2",
        "p.plus.tx2,p.plus.tx1",
    );
    fixture.rebind();
    assert_eq!(
        read_bundle(&fixture.request()).unwrap().records[0].1,
        "ACGTN"
    );
    for sequence in ["ACGT", "ACGTNN", "ACGTU", "ACGTR", "TCGTN", "ACG\u{a0}TN"] {
        let mut fixture = Fixture::new();
        fixture.replace_fasta("p.plus.fasta", "\nACGTN\n", &format!("\n{sequence}\n"));
        fixture.rebind();
        assert!(read_bundle(&fixture.request()).is_err(), "{sequence}");
    }
}

#[test]
fn target_normalized_digest_cannot_be_forged_in_manifest_and_header_together() {
    let mut fixture = Fixture::new();
    let digest = fixture.manifest["files"][0]["records"][0]["sequence_sha256"]
        .as_str()
        .unwrap()
        .to_owned();
    fixture.replace_fasta("p.plus.fasta", &digest, &"0".repeat(64));
    fixture.manifest["files"][0]["records"][0]["sequence_sha256"] = json!("0".repeat(64));
    fixture.rebind();
    invalid(
        &fixture.request(),
        "Normalized sequence SHA-256 mismatch for p.plus",
    );
}

#[test]
fn target_duplicate_ids_physical_tsss_and_transcripts_fail_without_merging() {
    for case in 0..3 {
        let mut fixture = Fixture::new();
        match case {
            0 => fixture.manifest["files"][1]["records"][0]["promoter_id"] = json!("p.plus"),
            1 => {
                for field in [
                    "chromosome",
                    "strand",
                    "tss_1based",
                    "genomic_start_1based",
                    "genomic_end_1based",
                ] {
                    fixture.manifest["files"][1]["records"][0][field] =
                        fixture.manifest["files"][0]["records"][0][field].clone();
                }
            }
            _ => {
                fixture.manifest["files"][0]["records"][0]["transcript_ids"] =
                    json!(["p.plus.tx1", "p.plus.tx1"])
            }
        }
        fixture.save_manifest();
        invalid(
            &fixture.request(),
            [
                "Duplicate promoter",
                "multi-gene memberships",
                "Duplicate transcript",
            ][case],
        );
    }
}

#[test]
fn target_identical_sequences_at_distinct_loci_are_retained_with_warning() {
    let mut fixture = Fixture::new();
    let digest = fixture.manifest["files"][1]["records"][0]["sequence_sha256"]
        .as_str()
        .unwrap()
        .to_owned();
    fixture.replace_fasta("p.minus.fasta", &digest, &sha256_hex_bytes(b"ACGTN"));
    fixture.replace_fasta("p.minus.fasta", "\nTGNAC\n", "\nACGTN\n");
    fixture.manifest["files"][1]["records"][0]["sequence_sha256"] =
        json!(sha256_hex_bytes(b"ACGTN"));
    fixture.rebind();
    let bundle = read_bundle(&fixture.request()).unwrap();
    assert_eq!(bundle.records.len(), 2);
    assert_eq!(bundle.records[0].1, bundle.records[1].1);
    assert_ne!(bundle.records[0].0.geometry, bundle.records[1].0.geometry);
    assert_eq!(
        bundle
            .records
            .iter()
            .map(|r| r.0.transcripts.len())
            .sum::<usize>(),
        4
    );
    assert!(
        bundle
            .warnings
            .iter()
            .any(|w| w.contains("Repeated sequence")
                && w.contains("both positional identities retained"))
    );
}

#[test]
fn target_manifest_file_memberships_counts_and_gene_ids_are_exact() {
    for case in 0..4 {
        let mut fixture = Fixture::new();
        match case {
            0 => fixture.manifest["files"][0]["record_count"] = json!(2),
            1 => fixture.manifest["files"][0]["records"][0]["gene_id"] = json!("other-gene"),
            2 => {
                let first = fs::read(fixture.root.join("p.plus.fasta")).unwrap();
                let second = fs::read(fixture.root.join("p.minus.fasta")).unwrap();
                fs::write(fixture.root.join("p.plus.fasta"), second).unwrap();
                fs::write(fixture.root.join("p.minus.fasta"), first).unwrap();
                fixture.rebind();
            }
            _ => {
                let path = fixture.root.join("p.plus.fasta");
                let mut bytes = fs::read(&path).unwrap();
                bytes.extend(bytes.clone());
                fs::write(path, bytes).unwrap();
                fixture.rebind();
            }
        }
        fixture.save_manifest();
        assert!(read_bundle(&fixture.request()).is_err(), "case {case}");
    }
}

#[test]
fn target_selection_joins_only_promoterome_id_then_checks_all_identity_fields() {
    for (field, value) in [
        ("gene_id", json!("other")),
        ("gene_name", json!("syn_plus")),
        ("genome_id", json!("Release 9")),
        ("chromosome", json!("another_chr")),
        ("strand", json!("-")),
        ("tss_1based", json!(101)),
    ] {
        let mut fixture = Fixture::new();
        fixture.selection["regions"][0]["genome_extraction"][field] = value;
        fixture.save_selection();
        invalid(&fixture.request(), "identity mismatch for promoter p.plus");
    }
    for id in ["SYN_PLUS", "SYN_PLUS.id", "p.plus.tx1", "missing"] {
        let mut fixture = Fixture::new();
        fixture.selection["regions"][0]["promoterome_id"] = json!(id);
        fixture.save_selection();
        invalid(&fixture.request(), "promoterome_id");
    }
    let mut fixture = Fixture::new();
    fixture.selection["regions"][1] = fixture.selection["regions"][0].clone();
    fixture.save_selection();
    invalid(&fixture.request(), "Duplicate selected promoterome_id");
}

#[test]
fn target_selection_requires_subset_and_rejects_duplicates_in_both_arrays() {
    for extraction in [false, true] {
        for transcripts in [
            json!([]),
            json!(["p.plus.tx1", "p.plus.tx1"]),
            json!(["unknown"]),
        ] {
            let mut fixture = Fixture::new();
            let row = &mut fixture.selection["regions"][0];
            if extraction {
                row["genome_extraction"]["transcript_ids"] = transcripts;
            } else {
                row["transcript_ids"] = transcripts;
            }
            fixture.save_selection();
            assert!(read_bundle(&fixture.request()).is_err());
        }
    }
}

#[test]
fn target_selection_evidence_windows_are_not_display_geometry_or_digest_joins() {
    let mut fixture = Fixture::new();
    fixture.selection["regions"][0]["sequence_length_bp"] = json!(2201);
    fixture.selection["regions"][0]["sequence_sha256"] =
        json!(format!("sha256:{}", "f".repeat(64)));
    fixture.selection["regions"][0]["genome_extraction"]["start_1based"] = json!(10);
    fixture.selection["regions"][0]["genome_extraction"]["end_1based"] = json!(2210);
    fixture.selection["regions"][0]["cutrun_support"]["factor"] = json!("DifferentSyntheticFactor");
    fixture.save_selection();
    let before = fs::read(fixture.selection_path()).unwrap();
    let bundle = read_bundle(&fixture.request()).unwrap();
    let evidence = &bundle.selection_evidence["p.plus"];
    assert_eq!(evidence.factor.as_deref(), Some("DifferentSyntheticFactor"));
    assert_eq!(
        evidence.label,
        "Selected in integrated report \u{2014} DifferentSyntheticFactor CUT&RUN-supported TSS window"
    );
    assert_eq!(
        evidence.criterion,
        fixture.selection["regions"][0]["cutrun_support"]["criterion"]
            .as_str()
            .unwrap()
    );
    assert!(
        evidence
            .legend
            .contains("does not establish TSS usage, direct binding, or promoter activity")
    );
    assert!(
        evidence
            .legend
            .contains("Synthetic descriptive evidence only")
    );
    assert_eq!(bundle.records[0].1.len(), 5);
    assert_eq!(
        bundle
            .inputs
            .iter()
            .find(|b| b.role == "selection")
            .unwrap()
            .sha256,
        sha256_hex_bytes(&before)
    );
    assert_eq!(before, fs::read(fixture.selection_path()).unwrap());
}

#[test]
fn target_without_selection_does_not_infer_selected_gene_membership() {
    let fixture = Fixture::new();
    let mut request = fixture.request();
    request.selection = None;
    let bundle = read_bundle(&request).unwrap();
    assert!(bundle.records.iter().all(|r| !r.2));
    assert!(bundle.selection_evidence.is_empty());
}

#[test]
fn target_explicit_fasta_list_must_match_every_manifest_path() {
    let fixture = Fixture::new();
    let mut request = fixture.request();
    request.fasta = vec!["p.minus.fasta".into(), "p.plus.fasta".into()];
    assert!(read_bundle(&request).is_ok());
    request.fasta.pop();
    invalid(&request, "exactly match");
    request.fasta.push("p.minus.fasta".into());
    invalid(&request, "duplicate paths");
    request.fasta = vec!["../bundle/p.plus.fasta".into()];
    invalid(&request, "traversal");
}

#[test]
fn target_external_selection_requires_explicit_nontraversing_path_and_regular_file() {
    let fixture = Fixture::new();
    let mut request = fixture.request();
    request.selection = Some("../selected.json".into());
    invalid(&request, "traversal");
    request.selection = Some(
        fixture
            .temp
            .path()
            .join("bundle/../selected.json")
            .to_string_lossy()
            .into_owned(),
    );
    invalid(&request, "traversal");
    request.selection = Some(fixture.temp.path().to_string_lossy().into_owned());
    invalid(&request, "regular file");
    request.selection = Some(
        fixture
            .temp
            .path()
            .join("./selected.json")
            .to_string_lossy()
            .into_owned(),
    );
    invalid(&request, "traversal");
    fs::copy(
        fixture.selection_path(),
        fixture.root.join("selection.json"),
    )
    .unwrap();
    request.selection = Some("selection.json".into());
    assert!(read_bundle(&request).is_ok());
    request.selection = Some("p.plus.fasta".into());
    invalid(&request, "distinct");
}

#[test]
fn target_rejects_unsafe_manifest_paths_and_large_files_before_reading() {
    let mut fixture = Fixture::new();
    fixture.manifest["files"][0]["filename"] = json!("../outside.fasta");
    fixture.save_manifest();
    invalid(&fixture.request(), "traversal");
    let fixture = Fixture::new();
    File::options()
        .write(true)
        .open(fixture.root.join("p.plus.fasta"))
        .unwrap()
        .set_len(MAX_MEMBER_BYTES + 1)
        .unwrap();
    invalid(&fixture.request(), "byte limit");
    let fixture = Fixture::new();
    File::options()
        .write(true)
        .open(fixture.selection_path())
        .unwrap()
        .set_len(MAX_MANIFEST_BYTES + 1)
        .unwrap();
    invalid(&fixture.request(), "byte limit");
}

#[cfg(unix)]
#[test]
fn target_fasta_symlinks_cannot_escape_the_bundle() {
    let fixture = Fixture::new();
    let path = fixture.root.join("p.plus.fasta");
    let outside = fixture.temp.path().join("outside.fasta");
    fs::rename(&path, &outside).unwrap();
    std::os::unix::fs::symlink(outside, path).unwrap();
    invalid(&fixture.request(), "escapes its directory");
}

#[test]
fn target_schema_and_duplicate_json_keys_are_not_inferred_or_normalized() {
    let mut fixture = Fixture::new();
    fixture.manifest["schema"] = json!("gentle.target_tss_fasta_export.v2");
    fixture.save_manifest();
    invalid(&fixture.request(), "Unsupported TSS bundle schema");
    let fixture = Fixture::new();
    let path = fixture.root.join("manifest.json");
    let text = fs::read_to_string(&path).unwrap();
    fs::write(path, text.replacen('{', "{\"assembly_id\":\"other\",", 1)).unwrap();
    invalid(&fixture.request(), "Invalid bundle manifest JSON");
}

#[test]
#[ignore = "requires GENTLE_TSS_REAL_INPUTS pointing to the read-only materialized a106cbbd input tree"]
fn target_real_58_tsss_and_13_selected_read_only_acceptance() {
    let base = PathBuf::from(
        std::env::var_os("GENTLE_TSS_REAL_INPUTS")
            .expect("set GENTLE_TSS_REAL_INPUTS=/private/tmp/tss-glen-inputs"),
    )
    .join("docs/examples/regulatory_region_comparison");
    let request = ComputeTssProfilesRequest {
        manifest: base
            .join("target_tss_fastas/bundle/manifest.json")
            .to_string_lossy()
            .into_owned(),
        panel: base
            .join("five_target_tss_integrated/jaspar_30_track_panel.json")
            .to_string_lossy()
            .into_owned(),
        fasta: vec![],
        selection: Some(
            base.join("five_target_tss_integrated/selected_tss_candidate_regions.json")
                .to_string_lossy()
                .into_owned(),
        ),
        expected_genome_id: "Human GRCh38 Ensembl 116".into(),
        expected_assembly: Some("GRCh38".into()),
        expected_annotation_release: None,
        expected_dataset_id: Some("human_grch38_ensembl116_promoterome_2000_200_v1".into()),
    };
    let bundle = read_bundle(&request).unwrap();
    assert_eq!(bundle.records.len(), 58);
    assert_eq!(bundle.records.iter().filter(|r| r.2).count(), 13);
    let mut counts = BTreeMap::<String, (usize, usize, usize)>::new();
    for (record, sequence, selected) in &bundle.records {
        let row = counts.entry(record.gene_symbol.clone()).or_default();
        row.0 += 1;
        row.1 += record.transcripts.len();
        row.2 += usize::from(*selected);
        assert_eq!(sequence.len(), 701);
        assert_eq!(
            record.geometry.genomic_at(500),
            Some(record.geometry.tss_1based)
        );
    }
    assert_eq!(
        counts,
        BTreeMap::from([
            ("CD44".into(), (31, 109, 2)),
            ("TGFB1".into(), (6, 17, 2)),
            ("SERPINE1".into(), (4, 15, 3)),
            ("PATZ1".into(), (5, 13, 3)),
            ("TP73".into(), (12, 20, 3)),
        ])
    );
    assert!(bundle.reference.annotation_release.is_none());
    assert_eq!(bundle.selection_evidence.len(), 13);
    assert!(
        bundle
            .selection_evidence
            .values()
            .all(|e| e.factor.as_deref() == Some("TP73"))
    );
}
