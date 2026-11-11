//! Invented two-exon nuclear locus, CDS and protein; generated inline in private tempdirs.
//! No patient, cell-line, gene-specific rule or downloaded biological fixture is used.

use super::*;
use crate::engine::{Engine, GentleEngine, Operation};
use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};
use std::{fs, path::Path};

fn bound(path: &Path) -> PrimerPayloadFile {
    let raw = fs::read(path).unwrap();
    PrimerPayloadFile {
        path: path.to_str().unwrap().into(),
        sha256: digest(&raw),
        bytes: raw.len() as u64,
    }
}
fn save(path: &Path, raw: &[u8]) -> PrimerPayloadFile {
    fs::write(path, raw).unwrap();
    bound(path)
}
fn replace_gtf(request: &mut CodingTranscriptModelsRequest, from: &str, to: &str) {
    let file = Path::new(&request.gtf.path);
    let text = fs::read_to_string(file).unwrap().replace(from, to);
    request.gtf = save(file, text.as_bytes());
}
fn fixture(
    reverse: bool,
    included_stop: bool,
    crlf: bool,
) -> (tempfile::TempDir, CodingTranscriptModelsRequest) {
    let root = tempfile::tempdir().unwrap();
    let dir = root
        .path()
        .canonicalize()
        .unwrap()
        .join("coding O'Brien space \u{03b1}");
    fs::create_dir(&dir).unwrap();
    let mut sequence = b"CCATGANNNNNNNNNNNNNNAACCCTAAGG".to_vec();
    assert_eq!(sequence.len(), 30);
    let strand = if reverse { '-' } else { '+' };
    if reverse {
        sequence = sequence
            .iter()
            .rev()
            .map(|b| crate::iupac_code::IupacCode::letter_complement(*b))
            .collect();
    }
    let range = |(a, z): (usize, usize)| if reverse { (230 - z, 230 - a) } else { (a, z) };
    let mut text = String::new();
    let rows = [
        ("transcript", (100, 130), "."),
        ("exon", (100, 106), "."),
        ("exon", (120, 130), "."),
        ("CDS", (102, 106), "0"),
        ("CDS", (120, if included_stop { 128 } else { 125 }), "2"),
        ("start_codon", (102, 105), "0"),
        ("stop_codon", (125, 128), "0"),
        ("five_prime_utr", (100, 102), "."),
        ("three_prime_utr", (128, 130), "."),
    ];
    for (kind, span, frame) in rows {
        let (a, z) = range(span);
        text.push_str(&format!("chrSynthetic\tsynthetic\t{kind}\t{}\t{z}\t.\t{strand}\t{frame}\tgene_id \"SYN_G\"; transcript_id \"SYN_T\"; transcript_version \"2\"; transcript_biotype \"protein_coding\"; tag \"mRNA_end_NF\";{}\n",a+1,
            if kind == "CDS" { " protein_id \"SYN_P\"; protein_version \"3\";" } else { "" }));
    }
    let newline = |s: String| if crlf { s.replace('\n', "\r\n") } else { s };
    let locus = newline(format!(
        ">chrSynthetic\n{}\n",
        String::from_utf8(sequence.clone()).unwrap()
    ));
    let reference_hash = digest(&sequence);
    let mut exons = vec![range((100, 106)), range((120, 130))];
    exons.sort_unstable();
    let annotation = serde_json::json!({"schema":crate::engine::PRIMER_PAYLOAD_ANNOTATION_SCHEMA,
        "reference_id":"synthetic-source", "assembly":"synthetic-assembly", "release":"synthetic-1",
        "reference_fasta_sha256":reference_hash, "coordinate_convention":"zero_based_half_open",
        "transcripts":[{"transcript_id":"SYN_T", "gene_ids":["SYN_G"], "target_name":"chrSynthetic",
            "target_length":1000, "genomic_contig":"chrSynthetic", "genomic_strand":strand, "genomic_exons":exons}]});
    let request = CodingTranscriptModelsRequest {
        schema: CODING_MODELS_REQUEST_SCHEMA.into(),
        gtf: save(&dir.join("annotation.gtf"), newline(text).as_bytes()),
        normalized_annotation: save(
            &dir.join("annotation.json"),
            newline(serde_json::to_string_pretty(&annotation).unwrap()).as_bytes(),
        ),
        reference: CdnaLocusReference {
            assembly: "synthetic-assembly".into(),
            contig: "chrSynthetic".into(),
            start_0based: 100,
            genomic_reference_sha256: reference_hash,
            locus_fasta: save(&dir.join("locus.fa"), locus.as_bytes()),
        },
        selector: CodingModelSelector {
            source_id: "synthetic-source".into(),
            gene_id: "SYN_G".into(),
        },
        provider: CodingModelProvider {
            name: "invented".into(),
            release: "synthetic-1".into(),
            cds_row_convention: if included_stop {
                CdsRowConvention::StopCodonIncludedInCdsRows
            } else {
                CdsRowConvention::StopCodonExcludedFromCdsRows
            },
        },
        reference_kind: CodingReferenceKind::Nuclear,
        translation_table: 1,
        protein_fasta: Some(save(
            &dir.join("protein.fa"),
            newline(">SYN_P.3\nMKP*\n".into()).as_bytes(),
        )),
        budgets: CodingModelBudgets::default(),
    };
    (root, request)
}
fn run(request: CodingTranscriptModelsRequest) -> CodingTranscriptModelsReport {
    derive_coding_transcript_models(request, &mut |_| true).unwrap()
}

#[test]
fn coding_models_split_codon_both_strands_conventions_and_line_endings() {
    for reverse in [false, true] {
        for included in [false, true] {
            for crlf in [false, true] {
                let (_root, request) = fixture(reverse, included, crlf);
                let report = run(request.clone());
                let model = &report.models[0];
                let coding = model.coding.as_ref().unwrap();
                assert_eq!(model.translation_status, "complete_verified");
                assert_eq!(coding.normalized_coding_cdna_0based_half_open, (2, 14));
                assert_eq!(
                    (coding.cdna_cds_start_1based, coding.cdna_cds_end_1based),
                    (3, 14)
                );
                assert_eq!(
                    coding.translated_protein_sha256.as_deref(),
                    Some(digest(b"MKP").as_str())
                );
                assert_eq!(coding.supplied_protein_comparison, "identical");
                assert_eq!(
                    model.transcript.versioned_accession.as_deref(),
                    Some("SYN_T.2")
                );
                assert_eq!(
                    coding
                        .protein
                        .as_ref()
                        .unwrap()
                        .versioned_accession
                        .as_deref(),
                    Some("SYN_P.3")
                );
                assert_eq!(model.transcript_completeness_tags, vec!["mRNA_end_NF"]);
                assert!(coding.cds_completeness_reasons.is_empty());
                assert_eq!(coding.start_codon.as_deref(), Some("ATG"));
                assert_eq!(coding.stop_codon.as_deref(), Some("TAA"));
                assert_eq!(hash(&run(request)).unwrap(), hash(&report).unwrap());
                let mut unhashed = report.clone();
                unhashed.content_sha256.clear();
                assert_eq!(hash(&unhashed).unwrap(), report.content_sha256);
            }
        }
    }
}

#[test]
fn coding_models_raw_hashes_differ_but_lf_crlf_semantics_match() {
    let (_a, lf) = fixture(false, false, false);
    let (_b, crlf) = fixture(false, false, true);
    assert_ne!(lf.gtf.sha256, crlf.gtf.sha256);
    let a = run(lf);
    let b = run(crlf);
    assert_eq!(hash(&a.models).unwrap(), hash(&b.models).unwrap());
    assert_eq!(a.sequence_sha256, b.sequence_sha256);
}

#[test]
fn coding_models_refuse_wrong_binding_geometry_and_stop_convention() {
    for fault in 0..8 {
        let (_root, mut request) = fixture(false, false, false);
        match fault {
            0 => request.gtf.sha256 = digest(b"wrong"),
            1 => request.gtf.bytes += 1,
            2 => request.provider.release = "wrong".into(),
            3 => request.reference.assembly = "wrong".into(),
            4 => replace_gtf(&mut request, "exon\t101\t106", "exon\t101\t107"),
            5 => replace_gtf(&mut request, "CDS\t121\t125", "CDS\t119\t125"),
            6 => request.provider.cds_row_convention = CdsRowConvention::StopCodonIncludedInCdsRows,
            7 => replace_gtf(&mut request, "CDS\t103\t106", "CDS\t103\t107"),
            _ => unreachable!(),
        }
        assert!(
            derive_coding_transcript_models(request, &mut |_| true).is_err(),
            "fault {fault}"
        );
    }
}

#[test]
fn coding_models_biological_irregularities_are_reported_not_repaired() {
    for (from, to, reason) in [
        (
            "tag \"mRNA_end_NF\";",
            "tag \"cds_start_NF\";",
            "incomplete_start",
        ),
        (
            "tag \"mRNA_end_NF\";",
            "tag \"cds_end_NF\";",
            "incomplete_end",
        ),
        ("+\t2\t", "+\t1\t", "phase_inconsistent"),
        ("stop_codon", "ignored_stop", "incomplete_end"),
        ("start_codon", "ignored_start", "incomplete_start"),
    ] {
        let (_root, mut request) = fixture(false, false, false);
        replace_gtf(&mut request, from, to);
        let report = run(request);
        let model = &report.models[0];
        assert!(
            model.translation_reasons.iter().any(|v| v == reason),
            "{reason}"
        );
        assert!(
            model
                .coding
                .as_ref()
                .unwrap()
                .translated_protein_sha256
                .is_none()
        );
    }
    for (from, to, reason) in [
        ("CCATGA", "CCCTGA", "start_not_atg"),
        ("AACCCTAA", "AATAATAA", "internal_stop"),
        ("CCATGA", "CCATGN", "ambiguous_sequence"),
    ] {
        let (_root, mut request) = fixture(false, false, false);
        let path = Path::new(&request.reference.locus_fasta.path);
        let text = fs::read_to_string(path).unwrap().replace(from, to);
        request.reference.locus_fasta = save(path, text.as_bytes());
        let report = run(request);
        assert!(
            report.models[0]
                .translation_reasons
                .iter()
                .any(|v| v == reason),
            "{reason}"
        );
    }
}

#[test]
fn coding_models_version_absence_conflict_and_protein_mismatch_remain_explicit() {
    let (_root, mut request) = fixture(false, false, false);
    replace_gtf(&mut request, " transcript_version \"2\";", "");
    let report = run(request);
    assert!(report.models[0].transcript.versioned_accession.is_none());
    assert_eq!(report.models[0].translation_status, "complete_verified");
    let (_root, mut request) = fixture(false, false, false);
    replace_gtf(
        &mut request,
        "transcript\t101\t130\t.\t+\t.\tgene_id \"SYN_G\"; transcript_id \"SYN_T\"; transcript_version \"2\";",
        "transcript\t101\t130\t.\t+\t.\tgene_id \"SYN_G\"; transcript_id \"SYN_T\"; transcript_version \"9\";",
    );
    assert_eq!(
        run(request).models[0].transcript.version_state,
        crate::observed_transcript_structures::annotation::VersionState::Inconsistent
    );
    let (_root, mut request) = fixture(false, false, false);
    let path = Path::new(&request.protein_fasta.as_ref().unwrap().path);
    request.protein_fasta = Some(save(path, b">SYN_P.3\nMKA\n"));
    let report = run(request);
    assert_eq!(
        report.models[0].translation_status,
        "supplied_protein_mismatch"
    );
    assert_eq!(
        report.models[0]
            .coding
            .as_ref()
            .unwrap()
            .supplied_protein_comparison,
        "differs"
    );
}

#[test]
fn coding_models_noncoding_and_undefined_cds_remain_distinct() {
    for (biotype, status) in [
        ("retained_intron", "noncoding"),
        ("protein_coding_CDS_not_defined", "coding_not_defined"),
    ] {
        let (_root, mut request) = fixture(false, false, false);
        let path = Path::new(&request.gtf.path);
        let text = fs::read_to_string(path)
            .unwrap()
            .lines()
            .filter(|l| !l.contains("\tCDS\t") && !l.contains("_codon\t"))
            .collect::<Vec<_>>()
            .join("\n")
            .replace("protein_coding", biotype);
        request.gtf = save(path, text.as_bytes());
        let report = run(request);
        assert_eq!(report.models[0].translation_status, status);
        assert!(report.models[0].coding.is_none());
    }
}

#[test]
fn coding_models_budgets_cancellation_and_non_nuclear_context_fail_closed() {
    for fault in 0..9 {
        let (_root, mut request) = fixture(false, false, false);
        match fault {
            0 => request.budgets.max_input_bytes = 1,
            1 => request.budgets.max_line_bytes = 1,
            2 => request.budgets.max_rows = 1,
            3 => request.budgets.max_segments = 1,
            4 => request.budgets.max_transcript_bases = 1,
            5 => request.budgets.max_work = 1,
            6 => request.budgets.max_output_bytes = 1,
            7 => request.translation_table = 2,
            8 => request.reference_kind = CodingReferenceKind::Mitochondrial,
            _ => unreachable!(),
        }
        assert!(
            derive_coding_transcript_models(request, &mut |_| true).is_err(),
            "fault {fault}"
        );
    }
    let (_root, request) = fixture(false, false, false);
    assert!(
        derive_coding_transcript_models(request, &mut |_| false)
            .unwrap_err()
            .message
            .contains("cancelled")
    );
}

#[test]
fn coding_models_operation_and_shared_shell_leave_state_history_and_input_bytes_unchanged() {
    let (_root, request) = fixture(false, false, false);
    let files = [
        &request.gtf,
        &request.normalized_annotation,
        &request.reference.locus_fasta,
        request.protein_fasta.as_ref().unwrap(),
    ];
    let before = files
        .iter()
        .map(|f| fs::read(&f.path).unwrap())
        .collect::<Vec<_>>();
    let mut engine = GentleEngine::new();
    let state = serde_json::to_value(engine.snapshot()).unwrap();
    let revision = engine.mutation_revision();
    let undo = engine.undo_available();
    let direct = engine
        .apply(Operation::DeriveCodingTranscriptModels {
            request: Box::new(request.clone()),
        })
        .unwrap();
    let journal = engine.journal_len();
    assert!(direct.coding_transcript_models.is_some());
    let line = format!(
        "rna-reads coding-models-derive {}",
        quote_shell_arg(&serde_json::to_string(&request).unwrap())
    );
    let command = parse_shell_line(&line).unwrap();
    let shell = execute_shell_command(&mut engine, &command).unwrap();
    assert!(!shell.state_changed);
    assert_eq!(
        serde_json::to_value(&direct.coding_transcript_models).unwrap(),
        shell.output["result"]["coding_transcript_models"]
    );
    assert_eq!(serde_json::to_value(engine.snapshot()).unwrap(), state);
    assert_eq!(engine.mutation_revision(), revision);
    assert_eq!(engine.undo_available(), undo);
    assert_eq!(engine.journal_len(), journal);
    for (file, raw) in files.iter().zip(before) {
        assert_eq!(fs::read(&file.path).unwrap(), raw);
    }
    assert!(parse_shell_line("rna-reads coding-models-derive '{}' extra").is_err());
}

#[test]
fn coding_models_mcp_requires_confirmation_before_io_and_keeps_project_absent() {
    let (root, request) = fixture(false, false, false);
    let expected = serde_json::to_value(run(request.clone())).unwrap();
    let state = root.path().join("absent MCP project.json");
    let operation = Operation::DeriveCodingTranscriptModels {
        request: Box::new(request.clone()),
    };
    let call = |args| {
        crate::mcp_server::mcp_tool_call_for_capability_surface_tests(
            state.to_str().unwrap(),
            "op",
            args,
        )
    };
    let mut invalid = request;
    invalid.gtf.path = root
        .path()
        .join("does not exist.gtf")
        .to_string_lossy()
        .into_owned();
    let denied = call(
        serde_json::json!({"operation":Operation::DeriveCodingTranscriptModels {request:Box::new(invalid)}}),
    );
    assert_eq!(denied["isError"], true);
    assert!(
        denied["content"][0]["text"]
            .as_str()
            .unwrap()
            .contains("without explicit confirm=true")
    );
    assert!(!state.exists());
    let allowed = call(serde_json::json!({"confirm":true,"operation":operation}));
    assert_ne!(allowed["isError"], true, "{allowed}");
    assert_eq!(
        allowed["structuredContent"]["result"]["coding_transcript_models"],
        expected
    );
    assert!(!state.exists());
}

#[cfg(feature = "lua-interface")]
#[test]
fn coding_models_lua_vm_matches_bound_service_without_state_mutation() {
    use mlua::LuaSerdeExt;
    let (_root, request) = fixture(false, false, false);
    let expected = serde_json::to_value(run(request.clone())).unwrap();
    let operation = Operation::DeriveCodingTranscriptModels {
        request: Box::new(request),
    };
    let state = crate::engine::ProjectState::default();
    let before = serde_json::to_value(&state).unwrap();
    let lua = crate::lua_interface::LuaInterface::new();
    lua.register_rust_functions().unwrap();
    lua.lua()
        .globals()
        .set("state_in", lua.lua().to_value(&state).unwrap())
        .unwrap();
    lua.lua()
        .globals()
        .set("operation", serde_json::to_string(&operation).unwrap())
        .unwrap();
    lua.lua()
        .load("result = apply_operation(state_in,operation)")
        .exec()
        .unwrap();
    let value: mlua::Value = lua.lua().globals().get("result").unwrap();
    let result: serde_json::Value = lua.lua().from_value(value).unwrap();
    assert_eq!(result["state"], before);
    assert_eq!(result["result"]["coding_transcript_models"], expected);
}

#[cfg(feature = "js-interface")]
#[test]
fn coding_models_js_vm_matches_bound_service_and_preserves_saved_project() {
    let (root, request) = fixture(false, false, false);
    let expected = run(request.clone());
    let operation = Operation::DeriveCodingTranscriptModels {
        request: Box::new(request),
    };
    let project = root.path().join("JS project O'Brien.json");
    crate::engine::ProjectState::default()
        .save_to_path(project.to_str().unwrap())
        .unwrap();
    let before = fs::read(&project).unwrap();
    let mut js = crate::js_interface::JavaScriptInterface::new();
    js.run_checked(format!(r#"
        function canon(v) {{
            if (Array.isArray(v)) return v.map(canon);
            if (v !== null && typeof v === 'object') return Object.fromEntries(Object.keys(v).sort().map(k=>[k,canon(v[k])]));
            return v;
        }}
        const input=load_project({path});
        const out=apply_operation(input,{op});
        if (JSON.stringify(canon(input)) !== JSON.stringify(canon(out.state))) throw Error('C1 changed state');
        if (JSON.stringify(canon(out.result.coding_transcript_models)) !== JSON.stringify(canon({expected}))) throw Error('C1 VM differs');
    "#,path=serde_json::to_string(&project.to_string_lossy()).unwrap(),
        op=serde_json::to_string(&operation).unwrap(),expected=serde_json::to_string(&expected).unwrap())).unwrap();
    assert_eq!(fs::read(&project).unwrap(), before);
}

#[test]
fn coding_models_non_triplet_cds_is_flagged_and_unknown_fields_are_refused() {
    let (_root, mut request) = fixture(false, false, false);
    replace_gtf(&mut request, "CDS\t121\t125", "CDS\t121\t124");
    replace_gtf(&mut request, "stop_codon\t126\t128", "stop_codon\t125\t127");
    replace_gtf(
        &mut request,
        "three_prime_utr\t129\t130",
        "three_prime_utr\t128\t130",
    );
    let report = run(request);
    assert!(
        report.models[0]
            .translation_reasons
            .contains(&"length_not_multiple_of_3".into())
    );
    let (_root, request) = fixture(false, false, false);
    let mut invalid = serde_json::to_value(request).unwrap();
    invalid["infer_provider_convention"] = serde_json::json!(true);
    assert!(serde_json::from_value::<CodingTranscriptModelsRequest>(invalid).is_err());
}
