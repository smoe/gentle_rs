//! Exact-start inventory and atomic, approval-bound derivation of TSS collections.

use super::*;
use gb_io::seq::{Feature, Location};
use gentle_protocol::genomic_regions::GenomicRegionStrand;
use gentle_protocol::tss_profiles::TssStrand;
use gentle_protocol::tss_workspace::*;

const COLLECTIONS_KEY: &str = "tss_collections_v1";
const SNAPSHOT_ALGORITHM: &str = "gentle.tss_biological_snapshot.v1";

fn biological_snapshot(
    dna: &DNAsequence,
    anchor: &SequenceGenomeAnchorSummary,
) -> Result<String, EngineError> {
    hash(&(
        SNAPSHOT_ALGORITHM,
        dna.clone_seq_record(),
        dna.overhang(),
        (
            &anchor.genome_id,
            &anchor.chromosome,
            anchor.start_1based,
            anchor.end_1based,
            anchor.strand,
        ),
    ))
}

fn hash(value: &impl Serialize) -> Result<String, EngineError> {
    serde_json::to_vec(value)
        .map(|bytes| sha256_prefixed_bytes(&bytes))
        .map_err(|e| EngineError::internal(format!("TSS snapshot serialization: {e}")))
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;

    fn prepare_fixture(dna: &mut DNAsequence) {
        GentleEngine::prepare_sequence(dna);
        // Synthetic recognition rule guarantees many cached sites without a user REBASE catalog.
        *dna.restriction_enzymes_mut() =
            vec![serde_json::from_value(serde_json::json!({
            "name":"SyntheticTssCache", "sequence":"AACCGT", "note":null, "cut":2, "overlap":0
        })).unwrap()];
        dna.update_computed_features();
        assert!(dna.restriction_enzyme_groups().len() > 1);
    }

    // Entirely synthetic locus/annotations, recreated here. Not TP73 biological data.
    fn transcript(id: &str, gene: &str, start: i64, end: i64, reverse: bool) -> Feature {
        let location = Location::simple_range(start, end);
        Feature {
            kind: "mRNA".into(),
            location: if reverse {
                Location::Complement(Box::new(location))
            } else {
                location
            },
            qualifiers: vec![
                ("transcript_id".into(), Some(id.into())),
                ("gene".into(), Some(gene.into())),
                ("gene_id".into(), Some(format!("gene_{gene}"))),
                ("source".into(), Some("synthetic".into())),
            ],
        }
    }

    pub(crate) fn engine(anchor_reverse: bool) -> GentleEngine {
        let mut dna = DNAsequence::from_sequence(&"AACCGTGA".repeat(125)).unwrap();
        dna.features_mut().extend([
            transcript("tx1", "TOY", 300, 700, false),
            transcript("tx2", "TOY", 300, 800, false),
            transcript("tx3", "TOY", 400, 800, false),
            transcript("tx_other", "OTHER", 300, 900, false),
        ]);
        prepare_fixture(&mut dna);
        let mut state = ProjectState::default();
        state.sequences.insert("locus".into(), dna);
        state.metadata.insert(PROVENANCE_METADATA_KEY.into(), serde_json::json!({GENOME_EXTRACTIONS_METADATA_KEY: [{
            "seq_id":"locus", "genome_id":"GRCh38", "chromosome":"1", "start_1based":1001, "end_1based":2000, "anchor_strand":if anchor_reverse {"-"} else {"+"}
        }]}));
        GentleEngine::from_state(state)
    }

    fn request() -> TssInventoryRequest {
        TssInventoryRequest {
            seq_id: "locus".into(),
            gene_query: "TOY".into(),
            collection_id: "toy_tss".into(),
            upstream_bp: 50,
            downstream_bp: 20,
        }
    }

    pub(crate) fn approved(engine: &GentleEngine) -> TssMaterializeRequest {
        let preview = engine.inspect_tss_inventory(&request()).unwrap();
        TssMaterializeRequest {
            inventory: preview.request,
            expected_approval_sha256: preview.approval_sha256,
            selected_tss_ids: preview.rows.into_iter().map(|r| r.tss_id).collect(),
        }
    }

    #[test]
    fn tss_tutorial_starter_has_three_exact_starts_and_no_precompleted_task() {
        let example: crate::workflow_examples::WorkflowExample = serde_json::from_str(
            include_str!("../../../docs/examples/workflows/tss_collection_gui_starter.json"),
        )
        .unwrap();
        let dir = tempfile::tempdir().unwrap();
        let state = crate::workflow_examples::run_example_workflow_for_project_state(
            &example,
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")),
            dir.path(),
        )
        .unwrap();
        let mut e = GentleEngine::from_state(state);
        assert!(e.list_tss_collections().unwrap().collections.is_empty());
        let preview = e
            .inspect_tss_inventory(&TssInventoryRequest {
                seq_id: "tss_locus".into(),
                gene_query: "TOY".into(),
                collection_id: "tss_windows".into(),
                upstream_bp: 500,
                downstream_bp: 200,
            })
            .unwrap();
        assert_eq!(preview.rows.len(), 3);
        assert_eq!(
            preview
                .rows
                .iter()
                .map(|r| r.tss_local_0based)
                .collect::<Vec<_>>(),
            [600, 900, 1499]
        );
        assert_eq!(
            preview
                .rows
                .iter()
                .map(|r| r.local_strand.as_str())
                .collect::<Vec<_>>(),
            ["+", "+", "-"]
        );
        assert_eq!(preview.rows[0].transcript_ids, ["plus_a", "plus_b"]);
        let approval = TssMaterializeRequest {
            inventory: preview.request,
            expected_approval_sha256: preview.approval_sha256,
            selected_tss_ids: preview.rows.iter().map(|r| r.tss_id.clone()).collect(),
        };
        println!(
            "TSS_TUTORIAL_APPROVAL={}",
            serde_json::to_string(&approval).unwrap()
        );
        let out = e
            .apply(Operation::MaterializeTssWindows { request: approval })
            .unwrap();
        let collection = out.tss_collection.unwrap();
        assert_eq!(collection.members.len(), 3);
        for member in &collection.members {
            assert_eq!(e.state.sequences[&member.tss.output_seq_id].len(), 701);
        }
        let report = e.get_tss_collection("tss_windows").unwrap();
        println!(
            "TSS_TUTORIAL_MEMBERS={}",
            serde_json::to_string(
                &report
                    .members
                    .iter()
                    .map(|r| &r.tss.output_seq_id)
                    .collect::<Vec<_>>()
            )
            .unwrap()
        );
    }

    #[test]
    fn tss_collection_list_is_deterministic_read_only_and_never_validates_members() {
        let mut e = engine(false);
        assert!(e.list_tss_collections().unwrap().collections.is_empty());
        let request = approved(&e);
        let collection = e
            .apply(Operation::MaterializeTssWindows { request })
            .unwrap()
            .tss_collection
            .unwrap();
        // A real stale member must remain discoverable, never reported as validated.
        e.state
            .sequences
            .remove(&collection.members[0].tss.output_seq_id);
        let mut legacy = serde_json::to_value(collection.as_ref()).unwrap();
        legacy["collection_id"] = serde_json::json!("a_legacy");
        legacy["inventory"]["request"]["collection_id"] = serde_json::json!("a_legacy");
        legacy["inventory"]
            .as_object_mut()
            .unwrap()
            .remove("snapshot_algorithm");
        let registry = e
            .state
            .metadata
            .get_mut(COLLECTIONS_KEY)
            .unwrap()
            .as_object_mut()
            .unwrap();
        registry.insert("a_legacy".into(), legacy);
        registry.insert("b_invalid".into(), serde_json::json!({"broken":true}));
        let before = serde_json::to_value(e.snapshot()).unwrap();
        let revision = e.structural_revision();
        let listing = e
            .apply(Operation::ListTssCollections {})
            .unwrap()
            .tss_collection_list
            .unwrap();
        assert_eq!(
            listing
                .collections
                .iter()
                .map(|r| r.collection_id.as_str())
                .collect::<Vec<_>>(),
            ["a_legacy", "b_invalid", "toy_tss"]
        );
        assert_eq!(
            listing.collections[0].record_status,
            TssCollectionRecordStatus::Legacy
        );
        assert_eq!(
            listing.collections[1].record_status,
            TssCollectionRecordStatus::Invalid
        );
        assert_eq!(listing.collections[1].window_count, None);
        assert_eq!(
            listing.collections[2].record_status,
            TssCollectionRecordStatus::Readable
        );
        assert_eq!(listing.collections[2].window_count, Some(2));
        assert_eq!(
            listing.collections[2].source_seq_id.as_deref(),
            Some("locus")
        );
        assert!(
            listing
                .collections
                .iter()
                .all(|r| r.validation_status == TssCollectionValidationStatus::NotChecked)
        );
        assert!(e.get_tss_collection("toy_tss").is_err());
        assert_eq!(before, serde_json::to_value(e.snapshot()).unwrap());
        assert_eq!(revision, e.structural_revision());
        assert_eq!(
            serde_json::to_vec(listing.as_ref()).unwrap(),
            serde_json::to_vec(&e.list_tss_collections().unwrap()).unwrap()
        );
        e.state
            .metadata
            .insert(COLLECTIONS_KEY.into(), serde_json::json!([]));
        assert!(
            e.list_tss_collections()
                .unwrap_err()
                .message
                .contains("expected an object")
        );
    }

    #[test]
    fn tss_inventory_groups_exact_starts_and_keeps_gene_source_strand_separate() {
        let mut e = engine(false);
        let before = serde_json::to_value(e.snapshot()).unwrap();
        let report = e.inspect_tss_inventory(&request()).unwrap();
        assert_eq!(report.rows.len(), 2);
        assert_eq!(report.rows[0].transcript_ids, ["tx1", "tx2"]);
        assert_eq!(report.rows[0].genomic_tss.start_0based, 1300);
        assert_eq!(before, serde_json::to_value(e.snapshot()).unwrap());
        let dna = e.state.sequences.get_mut("locus").unwrap();
        dna.features_mut()
            .push(transcript("tx_reverse", "TOY", 100, 301, true));
        let mut other_source = transcript("tx_refseq", "TOY", 300, 700, false);
        other_source
            .qualifiers
            .retain(|(k, _)| k.as_ref() != "source");
        other_source
            .qualifiers
            .push(("source".into(), Some("other_annotation".into())));
        dna.features_mut().push(other_source);
        let report = e.inspect_tss_inventory(&request()).unwrap();
        assert_eq!(report.rows.len(), 4);
        assert_eq!(
            report
                .rows
                .iter()
                .filter(|r| r.genomic_tss.start_0based == 1300)
                .count(),
            3
        );
        assert!(
            report
                .rows
                .iter()
                .all(|r| !r.transcript_ids.contains(&"tx_other".to_owned()))
        );
    }

    #[test]
    fn tss_missing_flanks_are_not_clipped_and_distinct_starts_survive() {
        let mut e = engine(false);
        let mut req = request();
        req.upstream_bp = 450;
        let report = e.inspect_tss_inventory(&req).unwrap();
        assert_eq!(report.rows.len(), 2);
        assert!(
            report
                .rows
                .iter()
                .all(|r| r.availability == TssWindowAvailability::MissingFlanks)
        );
        let count = e.state.sequences.len();
        let apply = TssMaterializeRequest {
            inventory: req,
            expected_approval_sha256: report.approval_sha256,
            selected_tss_ids: report.rows.into_iter().map(|r| r.tss_id).collect(),
        };
        assert!(
            e.materialize_tss_windows(&apply)
                .unwrap_err()
                .message
                .contains("flanks")
        );
        assert_eq!(e.state.sequences.len(), count);
    }

    #[test]
    fn tss_missing_or_partial_annotation_fails_closed_for_imported_records() {
        let mut e = engine(false);
        e.state
            .sequences
            .get_mut("locus")
            .unwrap()
            .features_mut()
            .clear();
        assert!(
            e.inspect_tss_inventory(&request())
                .unwrap_err()
                .message
                .contains("not evidence")
        );
        let mut f = transcript("partial", "TOY", 10, 100, false);
        f.location = Location::Range(
            (10, gb_io::seq::Before(true)),
            (100, gb_io::seq::After(false)),
        );
        e.state
            .sequences
            .get_mut("locus")
            .unwrap()
            .features_mut()
            .push(f);
        let report = e.inspect_tss_inventory(&request()).unwrap();
        assert!(report.rows.is_empty());
        assert_eq!(
            report.excluded_transcripts[0].reason,
            TssTranscriptExclusionReason::UncertainFivePrimeEnd
        );
    }

    #[test]
    fn tss_process_worker() {
        let Ok(path) = std::env::var("GENTLE_TSS_REGRESSION_DIR") else {
            return;
        };
        let path = std::path::Path::new(&path);
        let state: ProjectState =
            serde_json::from_slice(&std::fs::read(path.join("state.json")).unwrap()).unwrap();
        let mut e = GentleEngine::from_state(state);
        match std::env::var("GENTLE_TSS_REGRESSION_PHASE")
            .unwrap()
            .as_str()
        {
            "preview" => {
                let command = crate::engine_shell::parse_shell_line(&format!(
                    "promoters tss-inventory '{}'",
                    serde_json::to_string(&request()).unwrap()
                ))
                .unwrap();
                let output = crate::engine_shell::execute_shell_command(&mut e, &command).unwrap();
                let preview: TssInventoryReport =
                    serde_json::from_value(output.output["result"]["tss_inventory"].clone())
                        .unwrap();
                let req = TssMaterializeRequest {
                    inventory: preview.request,
                    expected_approval_sha256: preview.approval_sha256,
                    selected_tss_ids: preview.rows.into_iter().map(|r| r.tss_id).collect(),
                };
                std::fs::write(
                    path.join("approval.json"),
                    serde_json::to_vec(&req).unwrap(),
                )
                .unwrap();
            }
            "materialize" => {
                let req: TssMaterializeRequest =
                    serde_json::from_slice(&std::fs::read(path.join("approval.json")).unwrap())
                        .unwrap();
                let command = crate::engine_shell::parse_shell_line(&format!(
                    "promoters tss-materialize '{}'",
                    serde_json::to_string(&req).unwrap()
                ))
                .unwrap();
                crate::engine_shell::execute_shell_command(&mut e, &command).unwrap();
                std::fs::write(
                    path.join("state.json"),
                    serde_json::to_vec(e.snapshot()).unwrap(),
                )
                .unwrap();
            }
            "reopen" => {
                let report = e.get_tss_collection("toy_tss").unwrap();
                for member in report.members {
                    e.apply(Operation::RecomputeFeatures {
                        seq_id: member.tss.output_seq_id,
                    })
                    .unwrap();
                }
                assert_eq!(e.get_tss_collection("toy_tss").unwrap().members.len(), 2);
            }
            _ => panic!("unknown regression phase"),
        }
    }

    #[test]
    fn tss_prepared_genbank_preview_materialize_and_reopen_across_processes() {
        let dir = tempfile::tempdir().unwrap();
        let mut e = engine(false);
        // Synthetic fixture from engine(), round-tripped through real GenBank import.
        let gb = dir.path().join("locus.gb");
        std::fs::write(&gb, e.state.sequences["locus"].to_genbank_string().unwrap()).unwrap();
        let mut imported = DNAsequence::from_genbank_file(gb.to_str().unwrap())
            .unwrap()
            .remove(0);
        prepare_fixture(&mut imported);
        assert!(
            !serde_json::to_value(&imported).unwrap()["restriction_enzyme_groups"]
                .as_array()
                .unwrap()
                .is_empty()
        );
        e.state.sequences.insert("locus".into(), imported);
        std::fs::write(
            dir.path().join("state.json"),
            serde_json::to_vec(e.snapshot()).unwrap(),
        )
        .unwrap();
        let worker = format!(
            "{}::tss_process_worker",
            module_path!().split_once("::").unwrap().1
        );
        for phase in ["preview", "materialize", "reopen"] {
            let output = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", &worker, "--nocapture"])
                .env("GENTLE_TSS_REGRESSION_DIR", dir.path())
                .env("GENTLE_TSS_REGRESSION_PHASE", phase)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{phase}: {} {}",
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
            assert!(
                String::from_utf8_lossy(&output.stdout).contains("1 passed"),
                "worker must actually run"
            );
        }
    }

    #[test]
    fn tss_cache_refresh_is_not_edit_but_anchor_and_record_edits_are() {
        let mut e = engine(false);
        let approved = approved(&e);
        e.apply(Operation::RecomputeFeatures {
            seq_id: "locus".into(),
        })
        .unwrap();
        let (report, _) = e.materialize_tss_windows(&approved).unwrap();
        for m in &report.members {
            e.apply(Operation::RecomputeFeatures {
                seq_id: m.tss.output_seq_id.clone(),
            })
            .unwrap();
        }
        assert!(e.get_tss_collection("toy_tss").is_ok());
        let id = &report.members[0].tss.output_seq_id;
        for change in ["bases", "topology", "anchor"] {
            let mut altered = GentleEngine::from_state(e.snapshot().clone());
            if change == "anchor" {
                let entries = altered
                    .state
                    .metadata
                    .get_mut(PROVENANCE_METADATA_KEY)
                    .unwrap()[GENOME_EXTRACTIONS_METADATA_KEY]
                    .as_array_mut()
                    .unwrap();
                let entry = entries.iter_mut().find(|v| v["seq_id"] == *id).unwrap();
                entry["anchor_strand"] = serde_json::json!("-");
            } else if change == "topology" {
                altered
                    .state
                    .sequences
                    .get_mut(id)
                    .unwrap()
                    .set_circular(true);
            } else {
                let mut seq = altered.state.sequences[id].clone_seq_record();
                seq.seq[0] = b'T';
                altered
                    .state
                    .sequences
                    .insert(id.clone(), DNAsequence::from_genbank_seq(seq));
            }
            assert!(
                altered
                    .get_tss_collection("toy_tss")
                    .unwrap_err()
                    .message
                    .contains("edited"),
                "{change}"
            );
        }
    }

    #[test]
    fn tss_forget_legacy_metadata_preserves_sequences_and_lineage_and_is_undoable() {
        let mut e = engine(false);
        let req = approved(&e);
        e.apply(Operation::MaterializeTssWindows { request: req })
            .unwrap();
        e.state.metadata.get_mut(COLLECTIONS_KEY).unwrap()["toy_tss"]["inventory"]
            .as_object_mut()
            .unwrap()
            .remove("snapshot_algorithm");
        assert!(
            e.get_tss_collection("toy_tss")
                .unwrap_err()
                .message
                .contains("Legacy")
        );
        let mut expected = serde_json::to_value(e.snapshot()).unwrap();
        e.apply(Operation::ForgetTssCollection {
            collection_id: "toy_tss".into(),
        })
        .unwrap();
        expected["metadata"][COLLECTIONS_KEY]
            .as_object_mut()
            .unwrap()
            .remove("toy_tss");
        assert_eq!(serde_json::to_value(e.snapshot()).unwrap(), expected);
        assert!(e.get_tss_collection("toy_tss").is_err());
        e.undo_last_operation().unwrap();
        assert!(e.state.metadata[COLLECTIONS_KEY].get("toy_tss").is_some());
    }

    #[test]
    fn tss_genome_bounds_reject_clipped_starts_on_both_strands() {
        for reverse in [false, true] {
            let mut e = engine(false);
            let record = GenomeTranscriptRecord {
                chromosome: "1".into(),
                transcript_id: "cropped".into(),
                gene_id: Some("gene_TOY".into()),
                gene_name: Some("TOY".into()),
                strand: Some(if reverse { '-' } else { '+' }),
                transcript_start_1based: if reverse { 1101 } else { 901 },
                transcript_end_1based: if reverse { 2100 } else { 1900 },
                exons_1based: if reverse {
                    vec![(1101, 1500), (2051, 2100)]
                } else {
                    vec![(901, 950), (1401, 1900)]
                },
                cds_1based: vec![],
            };
            let feature =
                GentleEngine::transcript_feature_from_genome_record(&record, 1001, 2000).unwrap();
            e.state
                .sequences
                .get_mut("locus")
                .unwrap()
                .features_mut()
                .push(feature);
            let report = e.inspect_tss_inventory(&request()).unwrap();
            assert_eq!(report.rows.len(), 2);
            assert_eq!(
                report.excluded_transcripts[0].reason,
                TssTranscriptExclusionReason::TruncatedFivePrimeEnd
            );
            assert_eq!(report.excluded_transcripts[0].transcript_id, "cropped");
        }
    }

    #[test]
    fn tss_extraction_marks_dropped_first_exon_and_clipped_start_not_three_prime() {
        for reverse in [false, true] {
            let mut e = engine(false);
            let dna = e.state.sequences.get_mut("locus").unwrap();
            dna.features_mut().clear();
            let mut joined = transcript("dropped_exon", "TOY", 100, 900, reverse);
            let location = Location::Join(vec![
                Location::simple_range(100, 200),
                Location::simple_range(400, 500),
                Location::simple_range(800, 900),
            ]);
            joined.location = if reverse {
                Location::Complement(Box::new(location))
            } else {
                location
            };
            dna.features_mut().extend([
                joined,
                transcript("clipped_start", "TOY", 100, 900, reverse),
                transcript(
                    "intact_start",
                    "TOY",
                    if reverse { 0 } else { 400 },
                    if reverse { 500 } else { 1000 },
                    reverse,
                ),
            ]);
            let (from, to) = if reverse { (50, 750) } else { (350, 950) };
            e.apply(Operation::ExtractRegion {
                input: "locus".into(),
                from,
                to,
                output_id: Some("cropped".into()),
            })
            .unwrap();
            // Extraction has no implicit genome anchor: bind the synthetic coordinates explicitly.
            let mut anchor = e.state.metadata[PROVENANCE_METADATA_KEY]
                [GENOME_EXTRACTIONS_METADATA_KEY][0]
                .clone();
            anchor["seq_id"] = serde_json::json!("cropped");
            anchor["start_1based"] = serde_json::json!(1001 + from);
            anchor["end_1based"] = serde_json::json!(1000 + to);
            e.state.metadata.get_mut(PROVENANCE_METADATA_KEY).unwrap()
                [GENOME_EXTRACTIONS_METADATA_KEY]
                .as_array_mut()
                .unwrap()
                .push(anchor);
            let mut req = request();
            req.seq_id = "cropped".into();
            let report = e.inspect_tss_inventory(&req).unwrap();
            assert_eq!(report.rows.len(), 1);
            assert_eq!(report.rows[0].transcript_ids, ["intact_start"]);
            assert_eq!(report.excluded_transcripts.len(), 2);
            assert!(
                report
                    .excluded_transcripts
                    .iter()
                    .all(|r| r.reason == TssTranscriptExclusionReason::TruncatedFivePrimeEnd)
            );
        }
    }

    #[test]
    fn tss_exact_start_allows_fuzzy_three_prime_and_unlinked_other_transcript() {
        for reverse in [false, true] {
            let mut e = engine(false);
            let dna = e.state.sequences.get_mut("locus").unwrap();
            dna.features_mut().clear();
            let mut f = transcript("exact_start", "TOY", 300, 700, reverse);
            let loc = Location::Range(
                (300, gb_io::seq::Before(reverse)),
                (700, gb_io::seq::After(!reverse)),
            );
            f.location = if reverse {
                Location::Complement(Box::new(loc))
            } else {
                loc
            };
            let mut unrelated = transcript("unlinked", "OTHER", 100, 200, false);
            unrelated
                .qualifiers
                .retain(|(k, _)| !matches!(k.as_ref(), "gene" | "gene_id"));
            dna.features_mut().extend([
                f,
                unrelated,
                transcript("other_gene", "OTHER", 250, 650, reverse),
            ]);
            let report = e.inspect_tss_inventory(&request()).unwrap();
            assert_eq!(report.rows.len(), 1);
            assert_eq!(
                report.unassigned_transcripts[0].reason,
                TssTranscriptExclusionReason::MissingGeneLink
            );
            assert!(report.excluded_transcripts.is_empty());
            let mut other = request();
            other.gene_query = "OTHER".into();
            let other_report = e.inspect_tss_inventory(&other).unwrap();
            assert_eq!(other_report.rows.len(), 1);
            assert!(other_report.excluded_transcripts.is_empty());
            assert_eq!(
                other_report.unassigned_transcripts,
                report.unassigned_transcripts
            );
        }
    }

    #[test]
    fn tss_compound_reverse_forms_share_exact_start_and_truncation_behavior() {
        let mut e = engine(false);
        let dna = e.state.sequences.get_mut("locus").unwrap();
        dna.features_mut().clear();
        let parts = vec![
            Location::simple_range(100, 200),
            Location::simple_range(400, 500),
            Location::simple_range(800, 900),
        ];
        for (id, location) in [
            (
                "outer",
                Location::Complement(Box::new(Location::Join(parts.clone()))),
            ),
            (
                "inner",
                Location::Join(
                    parts
                        .into_iter()
                        .rev()
                        .map(|p| Location::Complement(Box::new(p)))
                        .collect(),
                ),
            ),
        ] {
            let mut f = transcript(id, "TOY", 100, 900, true);
            f.location = location;
            dna.features_mut().push(f);
        }
        let report = e.inspect_tss_inventory(&request()).unwrap();
        assert_eq!(report.rows.len(), 1);
        assert_eq!(report.rows[0].tss_local_0based, 899);
        assert_eq!(report.rows[0].transcript_ids, ["inner", "outer"]);
        assert_eq!(report.rows[0].local_strand, "-");
        for (from, to, truncated) in [(50, 750, true), (50, 850, true), (350, 950, false)] {
            let cropped = e.state.sequences["locus"]
                .extract_region_preserving_features(from, to)
                .unwrap();
            assert_eq!(cropped.features().len(), 2);
            for f in cropped.features() {
                assert_eq!(
                    f.qualifiers
                        .iter()
                        .any(|(k, _)| k.as_ref() == "gentle_transcript_5prime_truncated"),
                    truncated
                );
            }
        }
        let req = approved(&e);
        let (collection, _) = e.materialize_tss_windows(&req).unwrap();
        assert_eq!(collection.members.len(), 1);
        assert!(e.get_tss_collection("toy_tss").is_ok());
    }

    #[test]
    fn tss_compound_fuzzy_and_mixed_strands_stay_fail_closed() {
        let mut e = engine(false);
        let dna = e.state.sequences.get_mut("locus").unwrap();
        dna.features_mut().clear();
        for (id, upper_fuzzy, mixed) in [
            ("fuzzy_three", false, false),
            ("fuzzy_five", true, false),
            ("mixed", false, true),
        ] {
            let mut f = transcript(id, "TOY", 100, 900, true);
            let first = Location::Range(
                (100, gb_io::seq::Before(true)),
                (200, gb_io::seq::After(false)),
            );
            f.location = Location::Join(vec![
                if mixed {
                    first
                } else {
                    Location::Complement(Box::new(first))
                },
                Location::Complement(Box::new(Location::Range(
                    (800, gb_io::seq::Before(false)),
                    (900, gb_io::seq::After(upper_fuzzy)),
                ))),
            ]);
            dna.features_mut().push(f);
        }
        let report = e.inspect_tss_inventory(&request()).unwrap();
        assert_eq!(report.rows.len(), 1);
        assert_eq!(report.rows[0].transcript_ids, ["fuzzy_three"]);
        assert_eq!(report.excluded_transcripts.len(), 2);
        assert!(
            report
                .excluded_transcripts
                .iter()
                .all(|x| x.reason == TssTranscriptExclusionReason::UncertainFivePrimeEnd)
        );
        let cropped = e.state.sequences["locus"]
            .extract_region_preserving_features(750, 950)
            .unwrap();
        let mixed = cropped
            .features()
            .iter()
            .find(|f| f.qualifier_values("transcript_id").any(|id| id == "mixed"))
            .unwrap();
        assert!(
            mixed
                .qualifiers
                .iter()
                .any(|(k, _)| k.as_ref() == "gentle_transcript_5prime_truncated")
        );
        let mut f = transcript("qualifier_reverse", "TOY", 100, 900, false);
        f.qualifiers.push(("strand".into(), Some("-".into())));
        let endpoint = crate::feature_location::transcript_five_prime_endpoint(&f).unwrap();
        assert!(endpoint.reverse && endpoint.exact);
        assert_eq!(endpoint.position, 899);
        let mut dna = DNAsequence::from_sequence(&"A".repeat(1000)).unwrap();
        dna.features_mut().push(f);
        let cropped = dna.extract_region_preserving_features(50, 750).unwrap();
        assert!(
            cropped.features()[0]
                .qualifiers
                .iter()
                .any(|(k, _)| k.as_ref() == "gentle_transcript_5prime_truncated")
        );
    }

    #[test]
    fn tss_mixed_gene_links_merge_only_unambiguous_same_source_and_strand() {
        for scenario in ["unique", "ambiguous", "other_source", "other_strand"] {
            let mut e = engine(false);
            let dna = e.state.sequences.get_mut("locus").unwrap();
            dna.features_mut().clear();
            let known = transcript("known", "TOY", 300, 700, false);
            let mut unknown = transcript("no_id", "TOY", 300, 700, false);
            unknown.qualifiers.retain(|(k, _)| k.as_ref() != "gene_id");
            if scenario == "other_source" {
                unknown.qualifiers.retain(|(k, _)| k.as_ref() != "source");
                unknown
                    .qualifiers
                    .push(("source".into(), Some("other".into())));
            }
            if scenario == "other_strand" {
                unknown.location = Location::Complement(Box::new(Location::simple_range(100, 301)));
            }
            dna.features_mut().extend([known, unknown]);
            if scenario == "ambiguous" {
                let mut conflict = transcript("conflict", "TOY", 300, 700, false);
                conflict.qualifiers.retain(|(k, _)| k.as_ref() != "gene_id");
                conflict
                    .qualifiers
                    .push(("gene_id".into(), Some("different_gene".into())));
                dna.features_mut().push(conflict);
            }
            let report = e.inspect_tss_inventory(&request()).unwrap();
            assert_eq!(
                report.rows.len(),
                match scenario {
                    "unique" => 1,
                    "ambiguous" => 3,
                    _ => 2,
                },
                "{scenario}"
            );
            if scenario == "unique" {
                assert_eq!(report.rows[0].transcript_ids, ["known", "no_id"]);
                let mut by_id = request();
                by_id.gene_query = "gene_TOY".into();
                assert_eq!(
                    e.inspect_tss_inventory(&by_id).unwrap().rows[0].transcript_ids,
                    ["known", "no_id"]
                );
            } else {
                assert!(
                    report
                        .rows
                        .iter()
                        .any(|r| r.gene_id.is_none() && r.transcript_ids == ["no_id"])
                );
            }
        }
    }

    #[test]
    fn tss_empty_unassigned_field_preserves_legacy_bound_json() {
        let report = engine(false).inspect_tss_inventory(&request()).unwrap();
        let value = serde_json::to_value(&report).unwrap();
        assert!(value.get("unassigned_transcripts").is_none());
        let decoded: TssInventoryReport = serde_json::from_value(value.clone()).unwrap();
        assert_eq!(serde_json::to_value(decoded).unwrap(), value);
    }

    #[test]
    fn tss_materialization_is_idempotent_persisted_and_scan_compatible() {
        let mut e = engine(false);
        let req = approved(&e);
        let result = e
            .apply(Operation::MaterializeTssWindows {
                request: req.clone(),
            })
            .unwrap();
        assert_eq!(result.created_seq_ids.len(), 2);
        let report = result.tss_collection.unwrap();
        for m in &report.members {
            assert!(
                e.state.sequences[&m.tss.output_seq_id]
                    .name()
                    .as_deref()
                    .unwrap()
                    .contains("TSS 1:")
            );
        }
        let again = e
            .apply(Operation::MaterializeTssWindows { request: req })
            .unwrap();
        assert!(again.created_seq_ids.is_empty());
        let restored = GentleEngine::from_state(
            serde_json::from_value(serde_json::to_value(e.snapshot()).unwrap()).unwrap(),
        );
        assert_eq!(
            restored
                .get_tss_collection("toy_tss")
                .unwrap()
                .members
                .len(),
            2
        );
        let scan_op = Operation::ScanTfbsHitsCollection {
            collection_subject: report.subject,
            member_bindings: vec![],
            motifs: vec!["AAC".into()],
            min_llr_bits: None,
            min_llr_quantile: None,
            per_tf_thresholds: vec![],
            max_hits_per_member: Some(10),
            path: None,
        };
        let scan = e.apply(scan_op.clone()).unwrap();
        assert!(scan.collection_tfbs_hit_scan.is_some());
        assert!(e.get_tss_collection("toy_tss").is_ok());
        e.state
            .sequences
            .get_mut(&report.members[0].tss.output_seq_id)
            .unwrap()
            .features_mut()
            .clear();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("must_not_be_written.json");
        let mut stale_op = scan_op;
        if let Operation::ScanTfbsHitsCollection { path: output, .. } = &mut stale_op {
            *output = Some(path.to_str().unwrap().into());
        }
        assert!(e.apply(stale_op).unwrap_err().message.contains("edited"));
        assert!(!path.exists());
    }

    #[test]
    fn tss_stale_preview_and_any_collision_prevent_partial_writes() {
        let mut e = engine(false);
        let req = approved(&e);
        e.state
            .sequences
            .get_mut("locus")
            .unwrap()
            .features_mut()
            .push(transcript("new", "TOY", 500, 900, false));
        assert!(
            e.materialize_tss_windows(&req)
                .unwrap_err()
                .message
                .contains("stale")
        );
        assert_eq!(e.state.sequences.len(), 1);
        let req = approved(&e);
        let preview = e.inspect_tss_inventory(&request()).unwrap();
        let occupied = preview.rows.last().unwrap().output_seq_id.clone();
        e.state
            .sequences
            .insert(occupied, DNAsequence::from_sequence("AAC").unwrap());
        assert!(
            e.materialize_tss_windows(&req)
                .unwrap_err()
                .message
                .contains("already exists")
        );
        assert_eq!(e.state.sequences.len(), 2);
        assert!(!e.state.metadata.contains_key(COLLECTIONS_KEY));
    }

    #[test]
    fn tss_negative_local_and_genomic_orientations_are_independent() {
        for anchor_reverse in [false, true] {
            for transcript_reverse in [false, true] {
                let mut e = engine(anchor_reverse);
                let dna = e.state.sequences.get_mut("locus").unwrap();
                dna.features_mut().clear();
                dna.features_mut()
                    .push(transcript("tx", "TOY", 200, 701, transcript_reverse));
                let req = approved(&e);
                let (report, _) = e.materialize_tss_windows(&req).unwrap();
                let m = &report.members[0];
                let input = &e.state.sequences["locus"];
                let s = m.tss.window_local_start_0based.unwrap();
                let end = m.tss.window_local_end_0based_exclusive.unwrap();
                let mut expected = input
                    .extract_region_preserving_features(s, end)
                    .unwrap()
                    .clone_seq_record();
                if transcript_reverse {
                    expected = expected.revcomp();
                }
                let output = &e.state.sequences[&m.tss.output_seq_id];
                assert_eq!(output.forward_bytes(), expected.seq);
                assert_eq!(output.len(), 71);
                assert_eq!(
                    m.tss.genomic_tss.strand == GenomicRegionStrand::Minus,
                    anchor_reverse != transcript_reverse
                );
                let local_tss = if transcript_reverse { 700 } else { 200 };
                assert_eq!(
                    m.tss.genomic_tss.start_0based + 1,
                    if anchor_reverse {
                        2000 - local_tss
                    } else {
                        1001 + local_tss
                    }
                );
                #[cfg(feature = "desktop-gui")]
                {
                    let view = crate::tss_sequence_view::TssSequenceView::from_dna(output).unwrap();
                    assert_eq!(view.geometry.upstream_bp, 50);
                }
            }
        }
    }

    #[test]
    fn tss_collection_rejects_edited_outputs_and_namespace_reassignment() {
        let mut e = engine(false);
        let req = approved(&e);
        let (report, _) = e.materialize_tss_windows(&req).unwrap();
        let mut narrower = req.clone();
        narrower.selected_tss_ids.pop();
        assert!(
            e.materialize_tss_windows(&narrower)
                .unwrap_err()
                .message
                .contains("namespace")
        );
        e.state
            .sequences
            .get_mut(&report.members[0].tss.output_seq_id)
            .unwrap()
            .features_mut()
            .clear();
        assert!(
            e.get_tss_collection("toy_tss")
                .unwrap_err()
                .message
                .contains("edited")
        );
    }
}

fn annotation_source(feature: &Feature) -> String {
    GentleEngine::first_nonempty_feature_qualifier(
        feature,
        &["annotation_source", "source", "annotation_release"],
    )
    .unwrap_or_else(|| "project_annotation".into())
}

fn feature(kind: &str, start: usize, end: usize, label: String, note: String) -> Feature {
    Feature {
        kind: kind.to_owned().into(),
        location: Location::simple_range(start as i64, end as i64),
        qualifiers: vec![("label".into(), Some(label)), ("note".into(), Some(note))],
    }
}

impl GentleEngine {
    /// Inspect only annotations on this anchored project sequence; never infer a TSS from a CDS.
    pub fn inspect_tss_inventory(
        &self,
        request: &TssInventoryRequest,
    ) -> Result<TssInventoryReport, EngineError> {
        if request.gene_query.trim().is_empty()
            || request.gene_query != request.gene_query.trim()
            || request.collection_id.is_empty()
            || request.collection_id.len() > 80
            || !request
                .collection_id
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || b"_-.".contains(&c))
        {
            return Err(EngineError::invalid_input(
                "TSS inventory requires a gene query and an ASCII collection ID (letters, digits, _, -, .; up to 80 characters)",
            ));
        }
        let window_len = request
            .upstream_bp
            .checked_add(request.downstream_bp)
            .and_then(|n| n.checked_add(1))
            .filter(|n| *n <= 2_000_000)
            .ok_or_else(|| EngineError::invalid_input("TSS windows are limited to 2 Mb"))?;
        let dna = self.state.sequences.get(&request.seq_id).ok_or_else(|| {
            EngineError::new(ErrorCode::NotFound, "TSS source sequence not found")
        })?;
        if dna.is_circular() || dna.len() > 20_000_000 || dna.features().len() > 100_000 {
            return Err(EngineError::invalid_input(
                "TSS inventory requires a linear locus up to 20 Mb and 100,000 features",
            ));
        }
        let anchor = self.sequence_genome_anchor_summary(&request.seq_id)?;
        let source_snapshot_sha256 = biological_snapshot(dna, &anchor)?;
        let mut groups: BTreeMap<String, TssInventoryRow> = BTreeMap::new();
        let mut excluded_transcripts = Vec::new();
        let mut unassigned_transcripts = Vec::new();
        let mut linkage_notes = Vec::new();
        // Only explicit, unambiguous aliases in this locus/source/strand may fill a missing ID.
        let mut gene_ids_by_label: BTreeMap<(String, bool, String), BTreeSet<String>> =
            BTreeMap::new();
        for f in dna.features() {
            if Self::construct_reasoning_role_from_feature(f) != Some(ConstructRole::Transcript) {
                continue;
            }
            if let (Some(label), Some(id), Some(endpoint)) = (
                Self::first_nonempty_feature_qualifier(f, &["gene", "gene_name", "locus_tag"]),
                Self::first_nonempty_feature_qualifier(f, &["gene_id", "locus_tag"]),
                crate::feature_location::transcript_five_prime_endpoint(f),
            ) {
                gene_ids_by_label
                    .entry((
                        annotation_source(f),
                        endpoint.reverse,
                        label.to_ascii_lowercase(),
                    ))
                    .or_default()
                    .insert(id);
            }
        }
        let mut transcript_count = 0;
        for (feature_id, f) in dna.features().iter().enumerate() {
            if Self::construct_reasoning_role_from_feature(f) != Some(ConstructRole::Transcript) {
                continue;
            }
            let transcript_id =
                Self::first_nonempty_feature_qualifier(f, &["transcript_id", "name", "label"])
                    .unwrap_or_else(|| format!("feature_{feature_id}"));
            let excluded = |reason, explanation: &str| TssExcludedTranscript {
                feature_id,
                transcript_id: transcript_id.clone(),
                reason,
                explanation: explanation.into(),
            };
            if Self::first_nonempty_feature_qualifier(
                f,
                &["gene", "gene_name", "gene_id", "locus_tag"],
            )
            .is_none()
            {
                unassigned_transcripts.push(excluded(
                    TssTranscriptExclusionReason::MissingGeneLink,
                    "Locus-level annotation with unavailable gene linkage; neither assigned to nor excluded from the requested gene",
                ));
                continue;
            }
            let (gene_label, mut gene_id) = Self::transcript_gene_metadata(dna, f, feature_id);
            let endpoint = crate::feature_location::transcript_five_prime_endpoint(f);
            let annotation_source = annotation_source(f);
            let mut inferred_link = false;
            if gene_id.is_none()
                && let (Some(label), Some(endpoint)) = (&gene_label, &endpoint)
                && let Some(ids) = gene_ids_by_label.get(&(
                    annotation_source.clone(),
                    endpoint.reverse,
                    label.to_ascii_lowercase(),
                ))
                && ids.len() == 1
            {
                gene_id = ids.first().cloned();
                inferred_link = true;
            }
            if !gene_label
                .as_deref()
                .is_some_and(|s| s.eq_ignore_ascii_case(&request.gene_query))
                && !gene_id
                    .as_deref()
                    .is_some_and(|id| id == request.gene_query)
                && !Self::feature_matches_identifier(
                    f,
                    feature_id,
                    &request.gene_query,
                    &["gene", "gene_id", "gene_name", "locus_tag"],
                )
            {
                continue;
            }
            transcript_count += 1;
            if transcript_count > 10_000 {
                return Err(EngineError::invalid_input(
                    "TSS inventory exceeds 10,000 matching transcript features",
                ));
            }
            let Some(endpoint) = endpoint.filter(|endpoint| endpoint.exact) else {
                excluded_transcripts.push(excluded(
                    TssTranscriptExclusionReason::UncertainFivePrimeEnd,
                    "Partial, fuzzy or unsupported transcript 5-prime end; exact TSS unavailable",
                ));
                continue;
            };
            if f.qualifiers
                .iter()
                .any(|(k, _)| k.as_ref() == "gentle_transcript_5prime_truncated")
            {
                excluded_transcripts.push(excluded(TssTranscriptExclusionReason::TruncatedFivePrimeEnd,
                    "Transcript 5-prime end was removed by extraction; surviving exon boundary is not a TSS"));
                continue;
            }
            let mut ranges = Vec::new();
            collect_location_ranges_usize(&f.location, &mut ranges);
            ranges.sort_unstable();
            if ranges.is_empty() || ranges.iter().any(|(s, e)| s >= e || *e > dna.len()) {
                return Err(EngineError::invalid_input(
                    "Transcript annotation lies outside the source sequence",
                ));
            }
            let reverse = endpoint.reverse;
            let tss = endpoint.position;
            let local_strand = if reverse {
                GenomicRegionStrand::Minus
            } else {
                GenomicRegionStrand::Plus
            };
            let (genomic_tss, _) = self.interval_and_projection_from_local(
                &request.seq_id,
                tss as u64,
                tss as u64 + 1,
                local_strand,
                None,
            )?;
            let genomic_start =
                Self::first_nonempty_feature_qualifier(f, &["genomic_start_1based"]);
            let genomic_end = Self::first_nonempty_feature_qualifier(f, &["genomic_end_1based"]);
            if genomic_start.is_some() || genomic_end.is_some() {
                let bounds = genomic_start
                    .as_deref()
                    .and_then(|s| s.parse::<u64>().ok())
                    .zip(genomic_end.as_deref().and_then(|s| s.parse::<u64>().ok()))
                    .filter(|(s, e)| *s > 0 && e >= s);
                let Some((s, e)) = bounds else {
                    excluded_transcripts.push(excluded(TssTranscriptExclusionReason::InvalidGenomicBounds,
                        "Transcript genomic bounds are incomplete or invalid; exact start cannot be checked"));
                    continue;
                };
                let expected = if genomic_tss.strand == GenomicRegionStrand::Minus {
                    e
                } else {
                    s
                };
                if genomic_tss.start_0based + 1 != expected {
                    excluded_transcripts.push(excluded(TssTranscriptExclusionReason::TruncatedFivePrimeEnd,
                        "Local 5-prime endpoint differs from annotated genomic transcript start; extend the locus before using this TSS"));
                    continue;
                }
            }
            if inferred_link {
                linkage_notes.push(format!("{transcript_id} (feature {feature_id}): missing gene_id resolved to {} through an unambiguous explicit gene-label link in the same source and strand; original annotations retained.", gene_id.as_deref().unwrap()));
            }
            let identity = hash(&(
                &request.seq_id,
                &gene_id,
                &gene_label,
                &annotation_source,
                &genomic_tss,
            ))?;
            let tss_id = format!("tss_{}", identity.trim_start_matches("sha256:"));
            if let Some(row) = groups.get_mut(&tss_id) {
                row.transcript_ids.push(transcript_id);
                row.transcript_feature_ids.push(feature_id);
                continue;
            }
            let bounds = gentle_engine::tss_window_geometry::window_bounds(
                tss as u64 + 1,
                if reverse {
                    TssStrand::Minus
                } else {
                    TssStrand::Plus
                },
                request.upstream_bp,
                request.downstream_bp,
            )
            .ok()
            .filter(|(_, end)| *end <= dna.len() as u64);
            groups.insert(tss_id.clone(), TssInventoryRow {
                output_seq_id: format!("{}_{}", request.collection_id, tss_id), tss_id,
                gene_id, gene_label, annotation_source, transcript_ids: vec![transcript_id], transcript_feature_ids: vec![feature_id],
                tss_local_0based: tss, local_strand: if reverse { "-" } else { "+" }.into(), genomic_tss,
                window_local_start_0based: bounds.map(|(s,_)| s as usize - 1),
                window_local_end_0based_exclusive: bounds.map(|(_,e)| e as usize),
                availability: if bounds.is_some() { TssWindowAvailability::Available } else { TssWindowAvailability::MissingFlanks },
                explanation: if bounds.is_some() { "Exact annotated transcript start; not experimentally established initiation" } else { "Requested flanks exceed the loaded locus; extend the anchored sequence and inspect again. No clipping or network retrieval performed" }.into(),
            });
        }
        if transcript_count == 0 {
            return Err(EngineError::new(
                ErrorCode::NotFound,
                "Transcript annotation unavailable for this gene in the loaded locus (including annotation formats without an indexed transcript view). This is not evidence that the gene has no TSSs; import annotated mRNA/transcript features first",
            ));
        }
        if groups.len() > 256
            || groups
                .len()
                .checked_mul(window_len)
                .is_none_or(|n| n > 32_000_000)
        {
            return Err(EngineError::invalid_input(
                "TSS inventory exceeds 256 starts or 32 Mb of derived windows; narrow the locus/flanks",
            ));
        }
        let mut rows = groups.into_values().collect::<Vec<_>>();
        rows.sort_by(|a, b| {
            a.genomic_tss
                .start_0based
                .cmp(&b.genomic_tss.start_0based)
                .then(a.tss_id.cmp(&b.tss_id))
        });
        for row in &mut rows {
            row.transcript_ids.sort();
            row.transcript_ids.dedup();
        }
        let mut report = TssInventoryReport {
            schema: "gentle.tss_inventory.v1".into(), snapshot_algorithm: SNAPSHOT_ALGORITHM.into(), request: request.clone(), source_snapshot_sha256,
            approval_sha256: String::new(), rows, excluded_transcripts, unassigned_transcripts,
            warnings: vec!["Scope: transcript features on this loaded, anchored project locus only; not all possible biological starts or a cross-source consensus. No TFBS or CUT&RUN analysis was run.".into()],
        };
        report.warnings.extend(linkage_notes);
        if !report.unassigned_transcripts.is_empty() {
            report.warnings.push(format!("{} locus-level transcript annotations have no gene link; see unassigned_transcripts. They are not counted as exclusions from this gene.", report.unassigned_transcripts.len()));
        }
        if !report.excluded_transcripts.is_empty() {
            report.warnings.push(format!(
                "{} transcript annotations excluded from exact-start candidates; inspect excluded_transcripts. Unavailable starts are not evidence of absent promoters.",
                report.excluded_transcripts.len()
            ));
        }
        report.approval_sha256 = hash(&report)?;
        Ok(report)
    }

    /// List stored registry metadata without hashing sequences or validating membership.
    /// Stale members therefore remain discoverable and must be checked with `get_tss_collection`.
    pub fn list_tss_collections(&self) -> Result<TssCollectionListReport, EngineError> {
        let mut collections = Vec::new();
        if let Some(value) = self.state.metadata.get(COLLECTIONS_KEY) {
            let registry = value.as_object().ok_or_else(|| {
                EngineError::invalid_input("Invalid TSS collection registry: expected an object; not an empty collection list")
            })?;
            for (id, value) in registry {
                let mut entry = TssCollectionListEntry {
                    collection_id: id.clone(),
                    source_seq_id: None,
                    gene_query: None,
                    window_count: None,
                    record_status: TssCollectionRecordStatus::Invalid,
                    validation_status: TssCollectionValidationStatus::NotChecked,
                    diagnostic: None,
                };
                match TssCollectionReport::deserialize(value) {
                    Ok(report)
                        if report.schema == "gentle.tss_collection.v1"
                            && report.collection_id == *id
                            && report.inventory.request.collection_id == *id =>
                    {
                        entry.source_seq_id = Some(report.inventory.request.seq_id);
                        entry.gene_query = Some(report.inventory.request.gene_query);
                        entry.window_count = Some(report.members.len());
                        if report.inventory.snapshot_algorithm == SNAPSHOT_ALGORITHM {
                            entry.record_status = TssCollectionRecordStatus::Readable;
                        } else {
                            entry.record_status = TssCollectionRecordStatus::Legacy;
                            entry.diagnostic = Some("Legacy snapshot algorithm; explicit re-derivation required. Sequences are retained.".into());
                        }
                    }
                    Ok(_) => entry.diagnostic = Some(
                        "Stored schema or collection identity does not match the registry entry"
                            .into(),
                    ),
                    Err(_) => {
                        entry.diagnostic = Some("Stored collection record cannot be decoded".into())
                    }
                }
                collections.push(entry);
            }
        }
        collections.sort_by(|a, b| a.collection_id.cmp(&b.collection_id));
        Ok(TssCollectionListReport {
            schema: "gentle.tss_collection_list.v1".into(),
            collections,
        })
    }

    /// Return a persisted collection only while its members still match their stored snapshots.
    pub fn get_tss_collection(
        &self,
        collection_id: &str,
    ) -> Result<TssCollectionReport, EngineError> {
        let value = self
            .state
            .metadata
            .get(COLLECTIONS_KEY)
            .and_then(|v| v.get(collection_id))
            .ok_or_else(|| {
                EngineError::new(
                    ErrorCode::NotFound,
                    format!("TSS collection '{collection_id}' not found"),
                )
            })?;
        let report: TssCollectionReport = serde_json::from_value(value.clone()).map_err(|e| {
            EngineError::invalid_input(format!("Invalid stored TSS collection: {e}"))
        })?;
        if report.inventory.snapshot_algorithm != SNAPSHOT_ALGORITHM {
            return Err(EngineError::invalid_input(
                "Legacy TSS collection uses cache-sensitive snapshots; inspect again under a new collection ID or forget its metadata (sequences are retained). No silent reapproval performed",
            ));
        }
        let expected_subject =
            gentle_protocol::collection_subjects::CollectionSubjectRef::TssCollection {
                collection_id: collection_id.into(),
            };
        let mut preview = report.inventory.clone();
        preview.approval_sha256.clear();
        let unique_ids = report
            .members
            .iter()
            .map(|m| &m.tss.output_seq_id)
            .collect::<BTreeSet<_>>();
        let members = report
            .members
            .iter()
            .map(|m| CollectionMemberRef {
                stable_member_id: m.tss.output_seq_id.clone(),
                seq_id: Some(m.tss.output_seq_id.clone()),
                parent_member_id: Some(report.inventory.request.seq_id.clone()),
                ..Default::default()
            })
            .collect::<Vec<_>>();
        let fingerprint =
            canonical_collection_membership_json(CollectionSubjectKind::ProjectSequences, &members);
        if report.schema != "gentle.tss_collection.v1"
            || report.collection_id != collection_id
            || report.collection_id != report.inventory.request.collection_id
            || report.subject != expected_subject
            || report.members.len() != unique_ids.len()
            || report.members.len() > 256
            || report.members.is_empty()
            || report.lifting_mode != CollectionLiftingMode::Derive
            || report.collection_membership_fingerprint_sha256 != sha256_prefixed_str(&fingerprint)
            || hash(&preview)? != report.inventory.approval_sha256
            || report
                .members
                .iter()
                .any(|m| !report.inventory.rows.contains(&m.tss))
        {
            return Err(EngineError::invalid_input(
                "TSS collection membership or preview binding is inconsistent",
            ));
        }
        for member in &report.members {
            let dna = self
                .state
                .sequences
                .get(&member.tss.output_seq_id)
                .ok_or_else(|| {
                    EngineError::invalid_input(
                        "TSS collection member is missing; no windows opened",
                    )
                })?;
            let anchor = self.sequence_genome_anchor_summary(&member.tss.output_seq_id)?;
            if biological_snapshot(dna, &anchor)? != member.record_snapshot_sha256 {
                return Err(EngineError::invalid_input(
                    "TSS collection member was edited; inspect the individual sequence instead of reusing stale TSS geometry",
                ));
            }
        }
        Ok(report)
    }

    /// Forget only the registry entry, even when legacy or stale; retain sequences and lineage.
    pub(super) fn forget_tss_collection(&mut self, collection_id: &str) -> Result<(), EngineError> {
        let collections = self
            .state
            .metadata
            .get_mut(COLLECTIONS_KEY)
            .and_then(serde_json::Value::as_object_mut)
            .ok_or_else(|| EngineError::new(ErrorCode::NotFound, "TSS collection not found"))?;
        if collections.remove(collection_id).is_none() {
            return Err(EngineError::new(
                ErrorCode::NotFound,
                "TSS collection not found",
            ));
        }
        Ok(())
    }

    pub(super) fn materialize_tss_windows(
        &mut self,
        request: &TssMaterializeRequest,
    ) -> Result<(TssCollectionReport, Vec<String>), EngineError> {
        let inventory = self.inspect_tss_inventory(&request.inventory)?;
        if inventory.approval_sha256 != request.expected_approval_sha256 {
            return Err(EngineError::invalid_input(
                "TSS preview is stale or parameters changed; inspect again before materialization",
            ));
        }
        let selected = request.selected_tss_ids.iter().collect::<BTreeSet<_>>();
        if selected.is_empty()
            || selected.len() != request.selected_tss_ids.len()
            || selected
                .iter()
                .any(|id| !inventory.rows.iter().any(|r| &r.tss_id == *id))
        {
            return Err(EngineError::invalid_input(
                "Select unique TSS IDs from the approved inventory",
            ));
        }
        let rows = inventory
            .rows
            .iter()
            .filter(|r| selected.contains(&r.tss_id))
            .cloned()
            .collect::<Vec<_>>();
        if rows
            .iter()
            .any(|r| r.availability != TssWindowAvailability::Available)
        {
            return Err(EngineError::invalid_input(
                "Selected TSS windows have unavailable flanks; nothing materialized",
            ));
        }
        if self
            .state
            .metadata
            .get(COLLECTIONS_KEY)
            .and_then(|v| v.get(&request.inventory.collection_id))
            .is_some()
        {
            let existing = self.get_tss_collection(&request.inventory.collection_id)?;
            if existing.inventory.approval_sha256 == inventory.approval_sha256
                && existing.members.iter().map(|m| &m.tss).eq(rows.iter())
            {
                return Ok((existing, Vec::new()));
            }
            return Err(EngineError::invalid_input(
                "TSS collection namespace already belongs to another derivation; choose a new collection ID",
            ));
        }
        if rows
            .iter()
            .any(|r| self.state.sequences.contains_key(&r.output_seq_id))
        {
            return Err(EngineError::invalid_input(
                "TSS output sequence ID already exists; no outputs overwritten",
            ));
        }
        let parent_anchor = self.latest_genome_anchor_for_seq(&request.inventory.seq_id)?;
        let parent = &self.state.sequences[&request.inventory.seq_id];
        let mut prepared = Vec::new();
        let mut members = Vec::new();
        let mut anchors = Vec::new();
        let mut feature_count = 0usize;
        for row in rows {
            let start = row.window_local_start_0based.unwrap();
            let end = row.window_local_end_0based_exclusive.unwrap();
            let extracted = parent
                .extract_region_preserving_features(start, end)
                .ok_or_else(|| {
                    EngineError::invalid_input("Cannot extract the complete TSS window")
                })?;
            let mut seq = extracted.clone_seq_record();
            if row.local_strand == "-" {
                seq = seq.revcomp();
            }
            let (window, _) = self.interval_and_projection_from_local(
                &request.inventory.seq_id,
                start as u64,
                end as u64,
                if row.local_strand == "-" {
                    GenomicRegionStrand::Minus
                } else {
                    GenomicRegionStrand::Plus
                },
                None,
            )?;
            let sequence_sha256 = sha256_prefixed_bytes(
                &seq.seq
                    .iter()
                    .map(u8::to_ascii_uppercase)
                    .collect::<Vec<_>>(),
            );
            let strand = if window.strand == GenomicRegionStrand::Minus {
                "-"
            } else {
                "+"
            };
            let gene = row
                .gene_label
                .as_deref()
                .or(row.gene_id.as_deref())
                .unwrap_or(&request.inventory.gene_query);
            seq.name = Some(format!(
                "{gene} TSS {}:{} ({strand})",
                window.reference.contig_name,
                row.genomic_tss.start_0based + 1
            ));
            seq.comments = vec![
                format!(
                    "GENtle promoter_id={}; sequence_sha256={}",
                    row.tss_id,
                    sequence_sha256.trim_start_matches("sha256:")
                ),
                format!(
                    "Reference={}; assembly={}; chromosome={}; genomic={}..{}; genomic_strand={strand}; local_axis=transcript_5prime_to_3prime; TSS_local_1based={}",
                    request.inventory.seq_id,
                    window.reference.assembly_name,
                    window.reference.contig_name,
                    window.start_0based + 1,
                    window.end_0based_exclusive,
                    request.inventory.upstream_bp + 1
                ),
                format!(
                    "Origin=project_annotation_derivation; source_snapshot={}; approval={}; transcripts={}; not a verified external report bundle",
                    inventory.source_snapshot_sha256,
                    inventory.approval_sha256,
                    row.transcript_ids.join(",")
                ),
            ];
            seq.features.retain(|f| f.kind.as_ref() != "source");
            for feature_id in &row.transcript_feature_ids {
                let transcript = &parent.features()[*feature_id];
                let id =
                    Self::first_nonempty_feature_qualifier(transcript, &["transcript_id", "label"])
                        .unwrap_or_else(|| format!("feature_{feature_id}"));
                let mut exons = Vec::new();
                collect_location_ranges_usize(&transcript.location, &mut exons);
                exons.sort_unstable();
                if row.local_strand == "-" {
                    exons.reverse();
                }
                for (i, (s, e)) in exons.into_iter().enumerate() {
                    let (s, e) = (s.max(start), e.min(end));
                    if s >= e {
                        continue;
                    }
                    let (local_s, local_e) = if row.local_strand == "-" {
                        (end - e, end - s)
                    } else {
                        (s - start, e - start)
                    };
                    seq.features.push(feature("exon", local_s, local_e, format!("{id} E{}", i+1), "Transcript exon projected into this TSS window; clipped to the displayed span when necessary".into()));
                    if seq.features.len() > 100_000 {
                        return Err(EngineError::invalid_input(
                            "TSS window exceeds 100,000 projected features; narrow the locus",
                        ));
                    }
                }
            }
            seq.features.insert(
                0,
                feature(
                    "source",
                    0,
                    end - start,
                    gene.into(),
                    "Project-derived annotated TSS window".into(),
                ),
            );
            seq.features.push(feature(
                "misc_feature",
                request.inventory.upstream_bp,
                request.inventory.upstream_bp + 1,
                "Annotated TSS candidate".into(),
                format!(
                    "Genomic {}; annotation-derived, not experimentally established initiation",
                    row.genomic_tss.start_0based + 1
                ),
            ));
            let dna = DNAsequence::from_genbank_seq(seq);
            feature_count += dna.features().len();
            if dna.features().len() > 100_000 || feature_count > 250_000 {
                return Err(EngineError::invalid_input(
                    "TSS materialization exceeds 100,000 features per window or 250,000 total; select fewer starts",
                ));
            }
            let record_snapshot_sha256 = biological_snapshot(
                &dna,
                &SequenceGenomeAnchorSummary {
                    seq_id: row.output_seq_id.clone(),
                    genome_id: parent_anchor.genome_id.clone(),
                    chromosome: window.reference.contig_name.clone(),
                    start_1based: window.start_0based as usize + 1,
                    end_1based: window.end_0based_exclusive as usize,
                    strand: Some(if strand == "-" { '-' } else { '+' }),
                    anchor_verified: parent_anchor.anchor_verified,
                },
            )?;
            anchors.push(GenomeExtractionProvenance {
                seq_id: row.output_seq_id.clone(),
                recorded_at_unix_ms: Self::now_unix_ms(),
                operation: "MaterializeTssWindows".into(),
                genome_id: parent_anchor.genome_id.clone(),
                catalog_path: parent_anchor.catalog_path.clone().unwrap_or_default(),
                cache_dir: parent_anchor.cache_dir.clone(),
                chromosome: Some(window.reference.contig_name.clone()),
                start_1based: Some(window.start_0based as usize + 1),
                end_1based: Some(window.end_0based_exclusive as usize),
                gene_query: Some(request.inventory.gene_query.clone()),
                occurrence: None,
                gene_extract_mode: Some("exact_annotated_tss".into()),
                transcript_id: None,
                tss_1based: Some(row.genomic_tss.start_0based as usize + 1),
                promoter_upstream_bp: Some(request.inventory.upstream_bp),
                promoter_downstream_bp: Some(request.inventory.downstream_bp),
                gene_id: row.gene_id.clone(),
                gene_name: row.gene_label.clone(),
                strand: Some(if strand == "-" { '-' } else { '+' }),
                anchor_strand: Some(if strand == "-" { '-' } else { '+' }),
                anchor_verified: parent_anchor.anchor_verified,
                sequence_source_type: Some("project_derivation".into()),
                annotation_source_type: Some("project_annotation".into()),
                sequence_source: Some(request.inventory.seq_id.clone()),
                annotation_source: Some(row.annotation_source.clone()),
                sequence_sha1: None,
                annotation_sha1: None,
            });
            prepared.push((row.output_seq_id.clone(), dna));
            members.push(TssCollectionMember {
                tss: row,
                sequence_sha256,
                record_snapshot_sha256,
            });
        }
        let seq_ids = prepared
            .iter()
            .map(|(id, _)| id.clone())
            .collect::<Vec<_>>();
        let collection_members = seq_ids
            .iter()
            .map(
                |seq_id| gentle_protocol::collection_subjects::CollectionMemberRef {
                    stable_member_id: seq_id.clone(),
                    seq_id: Some(seq_id.clone()),
                    parent_member_id: Some(request.inventory.seq_id.clone()),
                    ..Default::default()
                },
            )
            .collect::<Vec<_>>();
        let fingerprint =
            gentle_protocol::collection_subjects::canonical_collection_membership_json(
                gentle_protocol::collection_subjects::CollectionSubjectKind::ProjectSequences,
                &collection_members,
            );
        let report = TssCollectionReport {
            schema: "gentle.tss_collection.v1".into(),
            collection_id: request.inventory.collection_id.clone(),
            inventory,
            members,
            lifting_mode: gentle_protocol::collection_subjects::CollectionLiftingMode::Derive,
            collection_membership_fingerprint_sha256: sha256_prefixed_str(&fingerprint),
            subject: gentle_protocol::collection_subjects::CollectionSubjectRef::TssCollection {
                collection_id: request.inventory.collection_id.clone(),
            },
        };
        let value =
            serde_json::to_value(&report).map_err(|e| EngineError::internal(e.to_string()))?;
        // All extraction, hashes and collisions are validated before the first state change.
        let mut collections = self
            .state
            .metadata
            .get(COLLECTIONS_KEY)
            .cloned()
            .unwrap_or_else(|| serde_json::json!({}));
        collections
            .as_object_mut()
            .ok_or_else(|| EngineError::invalid_input("Invalid TSS collection registry"))?
            .insert(report.collection_id.clone(), value);
        let mut provenance = self
            .state
            .metadata
            .get(PROVENANCE_METADATA_KEY)
            .cloned()
            .unwrap_or_else(|| serde_json::json!({}));
        let entries = provenance
            .as_object_mut()
            .ok_or_else(|| EngineError::invalid_input("Invalid genome provenance registry"))?
            .entry(GENOME_EXTRACTIONS_METADATA_KEY)
            .or_insert_with(|| serde_json::json!([]))
            .as_array_mut()
            .ok_or_else(|| EngineError::invalid_input("Invalid genome extraction registry"))?;
        for anchor in anchors {
            entries.push(
                serde_json::to_value(anchor).map_err(|e| EngineError::internal(e.to_string()))?,
            );
        }
        self.state.sequences.extend(prepared);
        self.state
            .metadata
            .insert(PROVENANCE_METADATA_KEY.into(), provenance);
        self.state
            .metadata
            .insert(COLLECTIONS_KEY.into(), collections);
        Ok((report, seq_ids))
    }
}
