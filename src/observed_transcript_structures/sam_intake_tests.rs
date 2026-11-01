//! Hand-crafted temporary SAM/JSON fixtures for private S2.2b file/unit admission.
//! Recreated by the `sam_intake::tests` filter; no biological records or external tools.

use super::*;
use crate::engine::{
    PRIMER_PAYLOAD_ANNOTATION_SCHEMA, PrimerPayloadAnnotation, PrimerPayloadTranscript,
    ReferenceContigDictionaryEntry,
};

const HEADER: &str =
    "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ctg_a\tLN:1000\n@SQ\tSN:ctg_b\tLN:1000\n@RG\tID:library\n";
const ROWS: &str = "same_name\t0\tctg_a\t101\t60\t10M10N20M\t*\t0\t0\t*\t*\tRG:Z:library\n\
same_name\t256\tctg_b\t201\t50\t30M\t*\t0\t0\t*\t*\tRG:Z:library\n\
same_name\t2048\tctg_a\t101\t60\t10M10N20M\t*\t0\t0\t*\t*\tRG:Z:library\n\
unmapped\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tRG:Z:library\n";

fn hash(c: char) -> String {
    format!("sha256:{}", c.to_string().repeat(64))
}
pub(crate) fn bound(path: &Path) -> PrimerPayloadFile {
    let raw = fs::read(path).unwrap();
    PrimerPayloadFile {
        path: path.display().to_string(),
        bytes: raw.len() as u64,
        sha256: crate::digest_utils::sha256_prefixed_bytes(&raw),
    }
}
pub(crate) fn manifest_file(r: &mut Request, m: &Manifest) {
    fs::write(&r.manifest.path, serde_json::to_vec_pretty(m).unwrap()).unwrap();
    r.manifest = bound(Path::new(&r.manifest.path));
}
pub(crate) fn complete() -> SearchDeclaration {
    SearchDeclaration {
        search_complete: true,
        views_complete: true,
        limitations: vec![],
        producer: "synthetic exhaustive declaration".into(),
        version: "1".into(),
        arguments: vec![],
        selection_scope: "all supplied competing mappings".into(),
    }
}
pub(crate) fn fixture() -> (tempfile::TempDir, Request, Manifest) {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("synthetic input.sam");
    fs::write(&sam, format!("{HEADER}{ROWS}")).unwrap();
    let dictionary = ReferenceContigDictionary {
        schema: REFERENCE_CONTIG_DICTIONARY_SCHEMA.into(),
        reference_id: "synthetic_reference".into(),
        assembly: "synthetic".into(),
        reference_fasta_sha256: hash('1'),
        contigs: ["ctg_a", "ctg_b"]
            .into_iter()
            .map(|name| ReferenceContigDictionaryEntry {
                name: name.into(),
                length: 1000,
            })
            .collect(),
    };
    let dictionary_path = dir.path().join("dictionary.json");
    fs::write(
        &dictionary_path,
        serde_json::to_vec_pretty(&dictionary).unwrap(),
    )
    .unwrap();
    let row = |ordinal, qname: &str| RowWitness {
        ordinal,
        qname: qname.into(),
        read_group: Some("library".into()),
    };
    let unit = |id: &str, rows: Vec<RowWitness>, length| Unit {
        unit_id: id.into(),
        unclipped_axis_length: length,
        rows,
        masks: None,
        ends: [vec![], vec![]],
        search: None,
    };
    let manifest = Manifest {
        schema: MANIFEST_SCHEMA.into(),
        sam_sha256: bound(&sam).sha256,
        units: vec![
            unit(
                "unit_a",
                vec![row(3, "same_name"), row(1, "same_name")],
                Some(30),
            ),
            unit("unit_b", vec![row(2, "same_name")], Some(30)),
            unit("unit_unmapped", vec![row(4, "unmapped")], None),
            unit("unit_missing", vec![], None),
        ],
    };
    let manifest_path = dir.path().join("query units.json");
    fs::write(
        &manifest_path,
        serde_json::to_vec_pretty(&manifest).unwrap(),
    )
    .unwrap();
    let request = Request {
        sam: bound(&sam),
        manifest: bound(&manifest_path),
        dictionary: bound(&dictionary_path),
        conversion: Conversion::DirectSam {
            producer: "hand-crafted fixture".into(),
            version: "1".into(),
            arguments: vec![],
            selection_scope: "explicit synthetic rows".into(),
        },
        reference_id: dictionary.reference_id,
        reference: ReferenceBinding {
            assembly: "synthetic".into(),
            fasta_sha256: hash('1'),
        },
        reference_gene: ReferenceGene {
            source_id: "source_a".into(),
            gene_id: "synthetic_gene".into(),
        },
        comparison_genes: vec![],
        annotation_sources: vec![AnnotationSource {
            source_id: "source_a".into(),
            provider: Some("synthetic".into()),
            annotation_sha256: hash('2'),
            annotation: PrimerPayloadAnnotation {
                schema: PRIMER_PAYLOAD_ANNOTATION_SCHEMA.into(),
                reference_id: "synthetic_reference".into(),
                assembly: "synthetic".into(),
                release: "synthetic_snapshot".into(),
                reference_fasta_sha256: hash('1'),
                coordinate_convention: "zero_based_half_open".into(),
                transcripts: vec![PrimerPayloadTranscript {
                    transcript_id: "synthetic_transcript".into(),
                    gene_ids: vec!["synthetic_gene".into()],
                    target_name: "ctg_a".into(),
                    target_length: 1000,
                    genomic_contig: Some("ctg_a".into()),
                    genomic_strand: Some('+'),
                    genomic_exons: vec![(100, 110), (120, 140)],
                }],
            },
        }],
        version_policy: sam_text::VersionPolicy::Strict16,
        unpadded_reference: true,
        required_flank_bases: 2,
        limits: sam_text::Limits {
            max_input_bytes: 1_000_000,
            max_line_bytes: 65536,
            max_lines: 100,
            max_rows: 100,
            max_dictionary_entries: 100,
            max_header_ids: 10,
            max_fields_per_line: 64,
            max_operations: 10000,
            max_sa_entries: 100,
            max_retained_bytes: 4 * 1024 * 1024,
            max_work: 100000,
        },
        max_units: 100,
        kernel_budgets: Budgets {
            max_queries: 100,
            max_records: 100,
            max_segments_per_query: 100,
            max_operations: 10000,
            max_skips_per_chain: 100,
            max_chain_families: 100,
            max_alternatives_per_record: 100,
            max_annotation_rows: 100,
            max_model_comparisons_per_family: 100,
            max_total_work: 1_000_000,
        },
    };
    (dir, request, manifest)
}
fn admitted(r: &Request) -> Prepared {
    prepare(r, &mut || Ok(())).unwrap()
}
fn refused(r: &Request) -> String {
    match prepare(r, &mut || Ok(())) {
        Ok(_) => panic!("unexpected admission"),
        Err(e) => e.message,
    }
}

#[test]
fn sam_intake_verifies_bytes_and_accounts_units_without_qname_inference() {
    let (_dir, r, m) = fixture();
    let prepared = admitted(&r);
    assert!(prepared.receipt.source_eof_verified);
    assert_eq!(prepared.receipt.coordinate_authority, "producer_declared");
    assert_eq!(prepared.receipt.identity_scope, "sam_declared_query_unit");
    assert_eq!(prepared.receipt.query_units, 4);
    assert_eq!(prepared.receipt.counts.alignment_rows, 4);
    assert_eq!(prepared.receipt.counts.secondary, 1);
    assert_eq!(prepared.receipt.counts.supplementary, 1);
    assert_eq!(prepared.receipt.counts.unmapped, 1);
    assert_eq!(
        prepared.input.units[0]
            .rows
            .iter()
            .map(|r| r.ordinal)
            .collect::<Vec<_>>(),
        [3, 1]
    );
    assert_eq!(
        prepared.input.units[0].rows[0].qname,
        prepared.input.units[1].rows[0].qname
    );
    assert_ne!(
        prepared.input.units[0].unit_id,
        prepared.input.units[1].unit_id
    );
    assert_eq!(
        serde_json::to_value(&prepared.receipt.query_unit_manifest).unwrap(),
        serde_json::to_value(&m).unwrap()
    );
    let analysis = sam::analyze(&prepared.input, &mut || Ok(())).unwrap();
    assert_eq!(analysis.summary.supports.len(), 4);
    assert!(analysis.summary.supports.iter().all(|s| s.definite_chains.is_empty() && s.reconciliation == Reconciliation::Unassessed));
    assert!(
        analysis
            .summary
            .observations
            .iter()
            .flat_map(|o| &o.flank_qualification)
            .all(|q| *q == Qualification::Unassessed)
    );
    assert!(
        analysis.units[2]
            .unavailable
            .contains(&sam::Unavailable::Unmapped)
    );
    assert!(
        analysis.units[3]
            .unavailable
            .contains(&sam::Unavailable::NotProvided)
    );
}

#[test]
fn sam_intake_refuses_missing_duplicate_out_of_range_and_wrong_identity_rows() {
    let (_dir, mut r, m) = fixture();
    for (case, expected) in [
        (0, "row_unassigned"),
        (1, "row_assignment"),
        (2, "row_assignment"),
        (3, "row_assignment"),
        (4, "row_witness"),
        (5, "row_witness"),
        (6, "unit_binding"),
        (7, "axis_binding"),
    ] {
        let mut next = m.clone();
        match case {
            0 => {
                next.units[0].rows.pop();
            }
            1 => next.units[1].rows[0] = next.units[0].rows[0].clone(),
            2 => next.units[0].rows[0].ordinal = 0,
            3 => next.units[0].rows[0].ordinal = 5,
            4 => next.units[0].rows[0].qname = "wrong".into(),
            5 => next.units[0].rows[0].read_group = None,
            6 => next.units[1].unit_id = next.units[0].unit_id.clone(),
            _ => next.units[0].unclipped_axis_length = Some(31),
        }
        manifest_file(&mut r, &next);
        assert_eq!(
            refused(&r),
            format!("observed_sam_intake_{expected}"),
            "case {case}"
        );
    }
}

#[test]
fn sam_intake_limits_sa_completeness_to_actual_rows_in_the_same_unit() {
    let (_dir, mut r, mut m) = fixture();
    let raw = format!("{HEADER}{ROWS}").replacen(
        "RG:Z:library\n",
        "RG:Z:library\tSA:Z:ctg_b,201,+,30M,50,0;\n",
        1,
    );
    fs::write(&r.sam.path, &raw).unwrap();
    r.sam = bound(Path::new(&r.sam.path));
    m.sam_sha256 = r.sam.sha256.clone();
    m.units[0].search = Some(complete());
    manifest_file(&mut r, &m);
    let p = admitted(&r);
    assert!(p.input.units[0].search_declared_complete);
    assert!(!p.input.units[0].views_declared_complete);
    assert_eq!(p.receipt.missing_sa_companions.len(), 1);
    assert_eq!(p.receipt.missing_sa_companions[0].row_ordinal, 1);
    let analysis = sam::analyze(&p.input, &mut || Ok(())).unwrap();
    assert!(
        analysis
            .summary
            .supports
            .iter()
            .all(|s| s.definite_chains.is_empty())
    );
    // Move the observed companion, not a fabricated SA segment, into the declared unit.
    let companion = m.units[1].rows.pop().unwrap();
    m.units[0].rows.push(companion);
    manifest_file(&mut r, &m);
    let p = admitted(&r);
    assert!(p.receipt.missing_sa_companions.is_empty());
    assert!(p.input.units[0].views_declared_complete);
    assert_eq!(p.input.units[0].rows.len(), 3);
}

#[test]
fn sam_intake_search_and_mask_declarations_remain_independent() {
    let (_dir, mut r, mut m) = fixture();
    m.units[0].search = Some(complete());
    manifest_file(&mut r, &m);
    let p = admitted(&r);
    let analysis = sam::analyze(&p.input, &mut || Ok(())).unwrap();
    let support = analysis
        .summary
        .supports
        .iter()
        .find(|s| s.key.unit_id == "unit_a")
        .unwrap();
    assert_eq!(support.definite_chains.len(), 1);
    assert!(
        analysis
            .summary
            .observations
            .iter()
            .filter(|o| o.key.unit_id == "unit_a")
            .all(|o| o.flank_qualification == [Qualification::Unassessed])
    );
    m.units[0].masks = Some(MaskWitness {
        coordinate_convention: "zero_based_half_open_unclipped_query".into(),
        unclipped_axis_length: 30,
        intervals: vec![],
        source_sha256: hash('3'),
        derivation: "synthetic explicit primer-free interval witness".into(),
    });
    manifest_file(&mut r, &m);
    let p = admitted(&r);
    let analysis = sam::analyze(&p.input, &mut || Ok(())).unwrap();
    assert!(
        analysis
            .summary
            .observations
            .iter()
            .filter(|o| o.key.unit_id == "unit_a")
            .all(|o| o.flank_qualification == [Qualification::Qualified])
    );
    m.units[0]
        .search
        .as_mut()
        .unwrap()
        .limitations
        .push("none".into());
    manifest_file(&mut r, &m);
    assert!(!admitted(&r).input.units[0].search_declared_complete);
    m.units[0].masks.as_mut().unwrap().unclipped_axis_length = 31;
    manifest_file(&mut r, &m);
    assert_eq!(refused(&r), "observed_sam_intake_mask_witness");
}

#[test]
fn sam_intake_lf_crlf_keep_physical_chains_but_not_raw_source_identity() {
    let (_dir, mut r, mut m) = fixture();
    let a = admitted(&r);
    fs::write(&r.sam.path, format!("{HEADER}{ROWS}").replace('\n', "\r\n")).unwrap();
    r.sam = bound(Path::new(&r.sam.path));
    m.sam_sha256 = r.sam.sha256.clone();
    manifest_file(&mut r, &m);
    let b = admitted(&r);
    assert_ne!(a.receipt.sam_sha256, b.receipt.sam_sha256);
    assert_eq!(a.receipt.counts, b.receipt.counts);
    let aa = sam::analyze(&a.input, &mut || Ok(())).unwrap();
    let bb = sam::analyze(&b.input, &mut || Ok(())).unwrap();
    assert_eq!(
        aa.summary.families[0].chain_id,
        bb.summary.families[0].chain_id
    );
    assert_ne!(
        aa.summary.observations[0].observation_id,
        bb.summary.observations[0].observation_id
    );
}

#[test]
fn sam_intake_refuses_tampered_sources_and_wrong_dictionary_or_conversion() {
    for file in ["sam", "manifest", "dictionary"] {
        let (_dir, r, _m) = fixture();
        let source = match file {
            "sam" => &r.sam,
            "manifest" => &r.manifest,
            _ => &r.dictionary,
        };
        let mut raw = fs::read(&source.path).unwrap();
        let i = if file == "sam" {
            // A valid changed QNAME reaches EOF/hash verification, unlike a broken @SQ.
            raw.windows(b"same_name".len())
                .position(|w| w == b"same_name")
                .unwrap()
        } else {
            raw.iter().position(|&b| b == b'a').unwrap()
        };
        raw[i] = if file == "sam" { b'S' } else { b'b' };
        fs::write(&source.path, raw).unwrap();
        assert_eq!(refused(&r), "observed_sam_intake_file_hash_or_change");
    }
    let (_dir, mut r, _m) = fixture();
    r.reference.assembly = "wrong".into();
    assert_eq!(refused(&r), "observed_sam_intake_dictionary_binding");
    r.reference.assembly = "synthetic".into();
    r.conversion = Conversion::BamConversion {
        bam_sha256: "unbound".into(),
        producer: "test".into(),
        version: "1".into(),
        arguments: vec![],
        selection_scope: "all".into(),
    };
    assert_eq!(refused(&r), "observed_sam_intake_conversion_provenance");
}

#[test]
fn sam_intake_admission_is_bounded_and_cancellation_is_atomic_and_read_only() {
    let (dir, r, _m) = fixture();
    let mut polls = 0;
    let result = prepare(&r, &mut || {
        polls += 1;
        Ok(())
    })
    .unwrap();
    for stop in 1..=polls {
        let mut count = 0;
        let error = prepare(&r, &mut || {
            count += 1;
            if count == stop {
                Err(bad("synthetic_cancel"))
            } else {
                Ok(())
            }
        });
        match error {
            Err(e) => assert_eq!(e.message, "synthetic_cancel"),
            Ok(_) => panic!("missed poll {stop}"),
        }
    }
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 3);
    for file in [&r.sam, &r.manifest, &r.dictionary] {
        assert_eq!(bound(Path::new(&file.path)).sha256, file.sha256);
    }
    let mut boundary = r.clone();
    boundary.limits.max_retained_bytes = result.receipt.retained_peak_bytes;
    admitted(&boundary);
    boundary.limits.max_retained_bytes -= 1;
    assert!(refused(&boundary).contains("memory_budget"));
    let mut boundary = r.clone();
    boundary.limits.max_input_bytes =
        (r.sam.bytes + r.manifest.bytes + r.dictionary.bytes) as usize;
    admitted(&boundary);
    boundary.limits.max_input_bytes -= 1;
    assert_eq!(refused(&boundary), "observed_sam_intake_admission");
    let mut boundary = r.clone();
    boundary.max_units = 3;
    assert_eq!(refused(&boundary), "observed_sam_intake_manifest_binding");
    for raw in [
        " ",
        "./input.sam",
        "sub/../input.sam",
        "input.sam/",
        "sub//input.sam",
    ] {
        let mut invalid = r.clone();
        invalid.sam.path = raw.into();
        assert_eq!(refused(&invalid), "observed_sam_intake_raw_path");
    }
}

#[test]
fn sam_intake_version_and_padding_policies_are_explicit_not_header_inferences() {
    let (_dir, mut r, mut m) = fixture();
    r.unpadded_reference = false;
    let prepared = admitted(&r);
    assert!(!prepared.receipt.unpadded_reference);
    let analysis = sam::analyze(&prepared.input, &mut || Ok(())).unwrap();
    assert!(analysis.summary.observations.is_empty());
    assert!(
        analysis.units[0]
            .unavailable
            .contains(&sam::Unavailable::PaddedReference)
    );
    r.unpadded_reference = true;
    fs::write(
        &r.sam.path,
        format!("{HEADER}{ROWS}").replace("VN:1.6", "VN:1.4"),
    )
    .unwrap();
    r.sam = bound(Path::new(&r.sam.path));
    m.sam_sha256 = r.sam.sha256.clone();
    manifest_file(&mut r, &m);
    assert_eq!(refused(&r), "observed_sam_text_header_version_or_order");
    r.version_policy = sam_text::VersionPolicy::Legacy10To16;
    let prepared = admitted(&r);
    assert_eq!(prepared.receipt.declared_version.as_deref(), Some("1.4"));
    assert_eq!(
        prepared.receipt.version_policy,
        sam_text::VersionPolicy::Legacy10To16
    );
}

#[test]
fn sam_intake_unanchored_competitor_withholds_a_complete_declared_unit() {
    let (_dir, mut r, mut m) = fixture();
    let raw = format!("{HEADER}{ROWS}").replace(
        "same_name\t2048\tctg_a\t101\t60\t10M10N20M",
        "same_name\t2048\tctg_a\t101\t60\t5I10N25M",
    );
    fs::write(&r.sam.path, raw).unwrap();
    r.sam = bound(Path::new(&r.sam.path));
    m.sam_sha256 = r.sam.sha256.clone();
    m.units[0].search = Some(complete());
    manifest_file(&mut r, &m);
    let p = admitted(&r);
    assert!(p.receipt.source_eof_verified);
    let analysis = sam::analyze(&p.input, &mut || Ok(())).unwrap();
    assert!(
        analysis.units[0]
            .unavailable
            .contains(&sam::Unavailable::UnanchoredSkip)
    );
    let support = analysis
        .summary
        .supports
        .iter()
        .find(|s| s.key.unit_id == "unit_a")
        .unwrap();
    assert_eq!(support.reconciliation, Reconciliation::Unassessed);
    assert!(support.hypotheses.is_empty());
    assert!(support.definite_chains.is_empty());
    assert!(support.possible_chains.is_empty());
    assert!(
        analysis
            .summary
            .observations
            .iter()
            .all(|o| o.key.unit_id == "unit_b")
    );
}

#[test]
fn sam_intake_charges_source_metadata_even_without_transcript_rows() {
    let (_dir, mut r, _m) = fixture();
    let mut source = r.annotation_sources[0].clone();
    source.annotation.transcripts.clear();
    r.annotation_sources = (0..64)
        .map(|i| {
            let mut copy = source.clone();
            copy.source_id = format!("source_{i}");
            copy
        })
        .collect();
    assert_eq!(refused(&r), "observed_sam_intake_memory_budget");
    r.annotation_sources.truncate(1);
    r.annotation_sources[0].annotation.release = "x".repeat(4097);
    assert_eq!(refused(&r), "observed_sam_intake_annotation_metadata");
}

#[test]
fn sam_intake_preserves_s21_reference_labels_without_treating_them_as_contig_names() {
    let (_dir, mut r, _m) = fixture();
    let mut dictionary: ReferenceContigDictionary =
        serde_json::from_slice(&fs::read(&r.dictionary.path).unwrap()).unwrap();
    dictionary.reference_id = "Synthetic reference label with spaces".into();
    r.reference_id = dictionary.reference_id.clone();
    r.annotation_sources[0].annotation.reference_id = dictionary.reference_id.clone();
    fs::write(
        &r.dictionary.path,
        serde_json::to_vec_pretty(&dictionary).unwrap(),
    )
    .unwrap();
    r.dictionary = bound(Path::new(&r.dictionary.path));
    let prepared = admitted(&r);
    assert_eq!(
        prepared.input.annotation_sources[0].annotation.reference_id,
        dictionary.reference_id
    );
    sam::analyze(&prepared.input, &mut || Ok(())).unwrap();
}

#[test]
fn sam_intake_unavailable_qname_is_accounted_without_inventing_read_identity() {
    let (_dir, mut r, mut m) = fixture();
    fs::write(
        &r.sam.path,
        format!("{HEADER}{ROWS}").replace("unmapped\t4", "*\t4"),
    )
    .unwrap();
    r.sam = bound(Path::new(&r.sam.path));
    m.sam_sha256 = r.sam.sha256.clone();
    m.units[2].rows[0].qname = "*".into();
    let unknown = m.units[2].rows.pop().unwrap();
    m.units[0].rows.push(unknown);
    m.units[0].search = Some(complete());
    manifest_file(&mut r, &m);
    let p = admitted(&r);
    assert_eq!(p.receipt.counts.alignment_rows, 4);
    assert_eq!(p.receipt.counts.unmapped, 1);
    assert_eq!(p.input.units[0].rows[2].qname, "*");
    let analysis = sam::analyze(&p.input, &mut || Ok(())).unwrap();
    assert!(
        analysis.units[0]
            .unavailable
            .contains(&sam::Unavailable::QnameNotProvided)
    );
    let support = analysis
        .summary
        .supports
        .iter()
        .find(|s| s.key.unit_id == "unit_a")
        .unwrap();
    assert!(support.possible_chains.is_empty());
    assert!(support.definite_chains.is_empty());
    assert_eq!(support.reconciliation, Reconciliation::Unassessed);
    assert!(
        analysis
            .summary
            .observations
            .iter()
            .all(|o| o.key.unit_id == "unit_b")
    );
}
