//! Synthetic cross-gene sequence sharing: portable annotations, never source-identity inference.

use super::tests::{request, shared_mature_sequence_sources, starter};
use super::*;

fn export(
    directory: &Path,
    historical: bool,
) -> (GentleEngine, RnaReadPanel, RnaReadPanelSummary, String) {
    let (engine, build, mature) = shared_mature_sequence_sources(directory, false);
    let mut panel = engine.build_rna_read_panel(&build).unwrap();
    if historical {
        panel.shared_transcript_sequences = None;
        seal_panel(&mut panel).unwrap();
    }
    let run = request(&panel, directory, &mature);
    engine
        .screen_rna_read_panel(run.clone(), &mut |_| true)
        .unwrap();
    let output = directory.join("summary");
    let summary = engine
        .export_rna_read_panel_summary(&ExportRnaReadPanelSummaryRequest {
            schema: RNA_PANEL_SUMMARY_REQUEST_SCHEMA.into(),
            report_paths: vec![
                Path::new(&run.output_dir)
                    .join("report.json")
                    .display()
                    .to_string(),
            ],
            output_dir: output.display().to_string(),
        })
        .unwrap();
    (
        engine,
        panel,
        summary,
        output.join("report.json").display().to_string(),
    )
}

#[test]
fn rna_panel_summary_sequence_sharing_reaches_json_tsv_and_bound_inspection_without_reassignment() {
    let td = tempfile::tempdir().unwrap();
    let (engine, panel, summary, path) = export(td.path(), false);
    let set = &summary.co_nomination_sets.as_ref().unwrap()[0];
    let sharing = set.sequence_sharing.as_ref().unwrap();
    assert_eq!(sharing.status, "recorded");
    assert_eq!(sharing.all_members_share_one_sequence, Some(true));
    assert_eq!(
        sharing.group_sha256s,
        vec![sha256_prefixed_bytes(
            &encode(&panel.shared_transcript_sequences.as_ref().unwrap()[0]).unwrap()
        )]
    );
    assert_eq!(summary.samples[0].reads, 1);
    assert_eq!(set.reads, 1);
    assert_eq!(set.member_ids, vec!["locus_0:0", "locus_1:0"]);
    assert!(
        summary
            .transcripts
            .iter()
            .all(|t| t.single_candidate_reads == 0 && t.co_nominated_reads == Some(1))
    );
    let (_, _, _, reopened) = engine.inspect_rna_read_panel(&path, None, 0).unwrap();
    assert_eq!(
        encode(&reopened.unwrap()).unwrap(),
        encode(&summary).unwrap()
    );
    engine
        .verify_rna_read_panel_summary_artifacts(&path, &summary)
        .unwrap();
    let mut table = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(td.path().join("summary/artifacts/co_nomination_sets.tsv"))
        .unwrap();
    let headers = table.headers().unwrap().clone();
    let row = table.records().next().unwrap().unwrap();
    let column = |name| headers.iter().position(|h| h == name).unwrap();
    assert_eq!(&row[column("sequence_sharing_status")], "recorded");
    assert_eq!(&row[column("all_members_share_one_sequence")], "true");
    assert_eq!(
        serde_json::from_str::<Vec<String>>(&row[column("shared_sequence_group_sha256s_json")])
            .unwrap(),
        sharing.group_sha256s
    );
    let mut tampered = summary.clone();
    tampered.co_nomination_sets.as_mut().unwrap()[0]
        .sequence_sharing
        .as_mut()
        .unwrap()
        .all_members_share_one_sequence = Some(false);
    tampered.content_sha256.clear();
    tampered.content_sha256 = sha256_prefixed_bytes(&encode(&tampered).unwrap());
    fs::write(&path, serde_json::to_vec_pretty(&tampered).unwrap()).unwrap();
    assert_eq!(
        engine
            .inspect_rna_read_panel(&path, None, 0)
            .unwrap_err()
            .message,
        "panel_summary_shared_sequence_binding"
    );
}

#[test]
fn rna_panel_summary_sequence_sharing_historical_unknown_is_not_recorded_empty() {
    let td = tempfile::tempdir().unwrap();
    let (engine, _, mut summary, path) = export(td.path(), true);
    assert!(summary.shared_transcript_sequences.is_none());
    let sharing = summary.co_nomination_sets.as_ref().unwrap()[0]
        .sequence_sharing
        .as_ref()
        .unwrap();
    assert_eq!(sharing.status, "unknown");
    assert!(sharing.group_sha256s.is_empty());
    assert_eq!(sharing.all_members_share_one_sequence, None);
    // A historical summary omits both additions; canonical bytes stay unchanged on reopen.
    for set in summary.co_nomination_sets.as_mut().unwrap() {
        set.sequence_sharing = None;
    }
    summary.content_sha256.clear();
    summary.content_sha256 = sha256_prefixed_bytes(&encode(&summary).unwrap());
    let bytes = encode(&summary).unwrap();
    let restored: RnaReadPanelSummary = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(encode(&restored).unwrap(), bytes);
    fs::write(&path, serde_json::to_vec_pretty(&restored).unwrap()).unwrap();
    engine.inspect_rna_read_panel(&path, None, 0).unwrap();
}

#[test]
fn rna_panel_summary_different_mature_sequences_preserve_genuine_cross_gene_ambiguity() {
    let td = tempfile::tempdir().unwrap();
    let (mut engine, build) = starter(2, td.path());
    let first = engine.state().sequences["source_0"].get_forward_string();
    let mut second = engine.state().sequences["source_1"].get_forward_string();
    second.replace_range(0..40, &first[..40]);
    let features = engine.state().sequences["source_1"].features().clone();
    let mut dna = DNAsequence::from_sequence(&second).unwrap();
    dna.features_mut().extend(features);
    engine.state_mut().sequences.insert("source_1".into(), dna);
    let panel = engine.build_rna_read_panel(&build).unwrap();
    assert!(
        panel
            .shared_transcript_sequences
            .as_ref()
            .unwrap()
            .is_empty()
    );
    let run = request(&panel, td.path(), &first[..40]);
    engine
        .screen_rna_read_panel(run.clone(), &mut |_| true)
        .unwrap();
    let summary = engine
        .export_rna_read_panel_summary(&ExportRnaReadPanelSummaryRequest {
            schema: RNA_PANEL_SUMMARY_REQUEST_SCHEMA.into(),
            report_paths: vec![
                Path::new(&run.output_dir)
                    .join("report.json")
                    .display()
                    .to_string(),
            ],
            output_dir: td.path().join("summary").display().to_string(),
        })
        .unwrap();
    let set = &summary.co_nomination_sets.as_ref().unwrap()[0];
    assert_eq!(set.class, "multi_locus_ambiguous");
    assert_eq!(set.member_ids.len(), 4);
    assert_eq!(set.reads, 1);
    let sharing = set.sequence_sharing.as_ref().unwrap();
    assert_eq!(sharing.status, "recorded");
    assert_eq!(sharing.all_members_share_one_sequence, Some(false));
    assert!(sharing.group_sha256s.is_empty());
}

#[test]
fn rna_panel_summary_sharing_groups_survive_zero_verified_files_and_reject_tampered_panel() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build, mature) = shared_mature_sequence_sources(td.path(), false);
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let run = request(&panel, td.path(), &mature);
    let report = engine
        .screen_rna_read_panel(run.clone(), &mut |_| false)
        .unwrap();
    assert!(report.files_completed.is_empty());
    let summary_request = ExportRnaReadPanelSummaryRequest {
        schema: RNA_PANEL_SUMMARY_REQUEST_SCHEMA.into(),
        report_paths: vec![
            Path::new(&run.output_dir)
                .join("report.json")
                .display()
                .to_string(),
        ],
        output_dir: td.path().join("summary").display().to_string(),
    };
    let summary = engine
        .export_rna_read_panel_summary(&summary_request)
        .unwrap();
    assert!(summary.transcripts.is_empty());
    assert_eq!(
        summary.shared_transcript_sequences.as_ref().unwrap().len(),
        1
    );
    engine
        .inspect_rna_read_panel(
            &td.path().join("summary/report.json").to_string_lossy(),
            None,
            0,
        )
        .unwrap();
    let mut tampered = panel;
    tampered.shared_transcript_sequences.as_mut().unwrap()[0]
        .member_ids
        .pop();
    seal_panel(&mut tampered).unwrap();
    fs::write(
        Path::new(&run.output_dir).join("artifacts/panel.json"),
        serde_json::to_vec_pretty(&tampered).unwrap(),
    )
    .unwrap();
    let mut again = summary_request;
    again.output_dir = td.path().join("another-summary").display().to_string();
    assert_eq!(
        engine
            .export_rna_read_panel_summary(&again)
            .unwrap_err()
            .message,
        "panel_shared_sequence_binding"
    );
}
