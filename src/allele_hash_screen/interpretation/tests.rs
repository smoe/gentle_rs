//! Invented saved-SAM input from the existing comparison fixture plus synthetic CDS/protein.
//! Recreated in private tempdirs; no real gene, private reads or provider calls.

use super::*;
use crate::allele_hash_screen::{
    CdnaBasecallDocument, CdsRowConvention, CodingModelProvider, CodingModelSelector,
    CodingReferenceKind, CodingTranscriptModelsRequest,
};
use crate::engine::{Engine, GentleEngine, Operation, PrimerPayloadAnnotation};
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
fn save<T: Serialize>(path: &Path, value: &T) -> PrimerPayloadFile {
    fs::write(path, serde_json::to_vec_pretty(value).unwrap()).unwrap();
    bound(path)
}

pub(super) fn fixture() -> (tempfile::TempDir, Request) {
    let (root, mut comparison) = comparison::synthetic_comparison_fixture();
    let sequence = "CCATGTCCAAACCCGGGAAATAAGGGGGGG";
    assert_eq!(sequence.len(), 30);
    let locus = Path::new(&comparison.reference.locus_fasta.path);
    fs::write(
        locus,
        format!(
            ">ctg_a synthetic coding locus\n{}{}{}{}\n",
            &sequence[..10],
            "N".repeat(10),
            &sequence[10..],
            "N".repeat(10)
        ),
    )
    .unwrap();
    comparison.reference.locus_fasta = bound(locus);
    for (i, sample) in comparison.samples.iter_mut().enumerate() {
        let path = Path::new(&sample.basecalls.path);
        let mut calls: CdnaBasecallDocument =
            serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
        for call in &mut calls.reads {
            call.sequence = sequence.into();
            if i == 0 {
                call.sequence.replace_range(6..8, "AA");
            }
        }
        sample.basecalls = save(path, &calls);
    }
    let annotation: PrimerPayloadAnnotation =
        serde_json::from_slice(&fs::read(&comparison.normalized_annotation.path).unwrap()).unwrap();
    let mut gtf = String::new();
    for t in annotation
        .transcripts
        .iter()
        .filter(|t| t.gene_ids.contains(&comparison.gene_id))
    {
        let attributes = format!(
            "gene_id \"{}\"; transcript_id \"{}\"; transcript_version \"2\"; transcript_biotype \"protein_coding\"; tag \"mRNA_end_NF\";",
            comparison.gene_id, t.transcript_id
        );
        assert_eq!(t.genomic_exons, vec![(100, 110), (120, 140)]);
        for (feature, a, z, phase) in [
            ("transcript", 100, 140, "."),
            ("exon", 100, 110, "."),
            ("exon", 120, 140, "."),
            ("CDS", 102, 110, "0"),
            ("CDS", 120, 130, "1"),
            ("start_codon", 102, 105, "0"),
            ("stop_codon", 130, 133, "0"),
            ("five_prime_utr", 100, 102, "."),
            ("three_prime_utr", 133, 140, "."),
        ] {
            gtf.push_str(&format!(
                "ctg_a\tsynthetic\t{feature}\t{}\t{z}\t.\t+\t{phase}\t{attributes}{}\n",
                a + 1,
                if feature == "CDS" {
                    " protein_id \"synthetic_protein\"; protein_version \"3\";"
                } else {
                    ""
                }
            ));
        }
    }
    let directory = root.path().canonicalize().unwrap();
    let path = directory.join("C2 coding O'Brien.gtf");
    fs::write(&path, gtf).unwrap();
    let protein = directory.join("C2 protein.fa");
    fs::write(&protein, b">synthetic_protein.3\nMSKPGK\n").unwrap();
    let model_request = CodingTranscriptModelsRequest {
        schema: coding_models::CODING_MODELS_REQUEST_SCHEMA.into(),
        gtf: bound(&path),
        normalized_annotation: comparison.normalized_annotation.clone(),
        reference: comparison.reference.clone(),
        selector: CodingModelSelector {
            source_id: comparison.source_id.clone(),
            gene_id: comparison.gene_id.clone(),
        },
        provider: CodingModelProvider {
            name: "synthetic".into(),
            release: annotation.release,
            cds_row_convention: CdsRowConvention::StopCodonExcludedFromCdsRows,
        },
        reference_kind: CodingReferenceKind::Nuclear,
        translation_table: 1,
        protein_fasta: Some(bound(&protein)),
        budgets: Default::default(),
    };
    let models =
        coding_models::derive_coding_transcript_models(model_request, &mut |_| true).unwrap();
    assert!(
        models
            .models
            .iter()
            .all(|m| m.translation_status == "complete_verified")
    );
    let evidence = comparison::compare_cdna_files(comparison, &mut |_| true).unwrap();
    let files = GeneFiles {
        comparison_report: save(&directory.join("comparison report.json"), &evidence),
        coding_models: save(&directory.join("coding models.json"), &models),
    };
    (
        root,
        Request {
            schema: REQUEST_SCHEMA.into(),
            genes: vec![files],
            policy: Default::default(),
            linkage: Linkage::None,
            known_variants: None,
            budgets: Default::default(),
            format: Format::Json,
        },
    )
}

fn snapshot(root: &Path) -> BTreeMap<std::path::PathBuf, Vec<u8>> {
    fn visit(base: &Path, dir: &Path, out: &mut BTreeMap<std::path::PathBuf, Vec<u8>>) {
        for entry in fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                visit(base, &path, out);
            } else {
                out.insert(
                    path.strip_prefix(base).unwrap().into(),
                    fs::read(&path).unwrap(),
                );
            }
        }
    }
    let mut out = BTreeMap::new();
    visit(root, root, &mut out);
    out
}

#[test]
fn cdna_interpretation_file_service_reverifies_c1_and_keeps_every_source_byte() {
    let (root, mut request) = fixture();
    let before = snapshot(root.path());
    let plain = interpret_cdna_files(request.clone(), &mut |_| true).unwrap();
    request.linkage = Linkage::SameUnit;
    let linked = interpret_cdna_files(request.clone(), &mut |_| true).unwrap();
    assert!(linked.inputs[0].normalized_evidence_reverified);
    assert!(!plain.inputs[0].normalized_evidence_reverified);
    assert!(linked.genes[0].codon_combinations.iter().any(|c| {
        c.consequence
            .consequence_terms
            .contains(&"stop_gained".into())
    }));
    for (left, right) in plain.genes[0].sites.iter().zip(&linked.genes[0].sites) {
        assert_eq!(
            serde_json::to_value(&left.inherited_counts).unwrap(),
            serde_json::to_value(&right.inherited_counts).unwrap()
        );
    }
    assert_eq!(snapshot(root.path()), before);
    assert_eq!(
        interpret_cdna_files(request, &mut |_| true)
            .unwrap()
            .content_sha256,
        linked.content_sha256
    );
}

#[test]
fn cdna_interpretation_shell_projection_and_saved_json_regeneration_preserve_project() {
    let (root, request) = fixture();
    let before = snapshot(root.path());
    let report = interpret_cdna_files(request.clone(), &mut |_| true).unwrap();
    let decoded: Report = serde_json::from_slice(&serde_json::to_vec(&report).unwrap()).unwrap();
    for format in [Format::Json, Format::Tsv, Format::Summary] {
        assert_eq!(
            render(&report, format).unwrap(),
            render(&decoded, format).unwrap()
        );
    }
    let mut engine = GentleEngine::new();
    let state = serde_json::to_value(engine.snapshot()).unwrap();
    let journal = engine.journal_len();
    for (format, text) in [(Format::Tsv, "tsv"), (Format::Summary, "summary")] {
        let line = format!(
            "rna-reads allele-hash-screen --interpret-cdna {} --format {text}",
            quote_shell_arg(&serde_json::to_string(&request).unwrap())
        );
        let command = parse_shell_line(&line).unwrap();
        assert!(!command.is_state_mutating());
        let result = execute_shell_command(&mut engine, &command).unwrap();
        assert!(!result.state_changed);
        let mut changed = request.clone();
        changed.format = format;
        let expected = interpret_cdna_files(changed, &mut |_| true).unwrap();
        assert_eq!(
            result.output["cdna_interpretation_projection"]["text"],
            render(&expected, format).unwrap()
        );
    }
    assert_eq!(serde_json::to_value(engine.snapshot()).unwrap(), state);
    assert_eq!(engine.journal_len(), journal);
    assert_eq!(snapshot(root.path()), before);
}

#[test]
fn cdna_interpretation_report_model_and_raw_bindings_refuse_tampering_and_cancellation() {
    for target in 0..4 {
        let (_root, mut request) = fixture();
        match target {
            0 => request.genes[0].comparison_report.sha256 = digest(b"wrong"),
            1 => request.genes[0].coding_models.bytes += 1,
            2 => {
                let file = &request.genes[0].comparison_report;
                let mut report: comparison::CdnaFileComparisonReport =
                    serde_json::from_slice(&fs::read(&file.path).unwrap()).unwrap();
                report.comparison.gene = "different".into();
                request.genes[0].comparison_report = save(Path::new(&file.path), &report);
            }
            3 => {
                let file = &request.genes[0].coding_models;
                let mut report: coding_models::CodingTranscriptModelsReport =
                    serde_json::from_slice(&fs::read(&file.path).unwrap()).unwrap();
                report.models[0].genomic_exons_0based_half_open[0].1 -= 1;
                report.content_sha256.clear();
                report.content_sha256 = canonical(&report).unwrap();
                request.genes[0].coding_models = save(Path::new(&file.path), &report);
            }
            _ => unreachable!(),
        }
        assert!(
            interpret_cdna_files(request, &mut |_| true).is_err(),
            "target {target}"
        );
    }
    let (_root, request) = fixture();
    assert!(
        interpret_cdna_files(request, &mut |_| false)
            .unwrap_err()
            .message
            .contains("cancelled")
    );
}

#[test]
fn cdna_interpretation_shared_lowering_preserves_comparison_bytes_and_refuses_linkage_drift() {
    let (_root, mut request) = fixture();
    let file = &request.genes[0].comparison_report;
    let mut report: comparison::CdnaFileComparisonReport =
        serde_json::from_slice(&fs::read(&file.path).unwrap()).unwrap();
    let before = serde_json::to_vec(&report).unwrap();
    assert_eq!(
        serde_json::to_vec(
            &comparison::compare_cdna_files(report.request.clone(), &mut |_| true).unwrap()
        )
        .unwrap(),
        before
    );
    let lowered = comparison::lower_cdna_files(
        &report.request,
        &mut |_| true,
        true,
        2_000_000,
        report.request.max_input_bytes,
    )
    .unwrap();
    assert_eq!(
        struct_hash(&lowered.request).unwrap(),
        report.comparison.request_sha256
    );
    report.comparison.request_sha256 = digest(b"wrong normalized evidence");
    report.comparison.content_sha256.clear();
    report.comparison.content_sha256 = struct_hash(&report.comparison).unwrap();
    report.content_sha256.clear();
    report.content_sha256 = canonical(&report).unwrap();
    request.genes[0].comparison_report = save(Path::new(&file.path), &report);
    request.linkage = Linkage::SameUnit;
    assert!(
        interpret_cdna_files(request, &mut |_| true)
            .unwrap_err()
            .message
            .contains("linkage_input_drift")
    );
}

#[test]
fn cdna_interpretation_mcp_confirmation_precedes_file_io_and_keeps_state_absent() {
    let (root, request) = fixture();
    let state = root.path().join("absent interpretation project.json");
    let expected = interpret_cdna_files(request.clone(), &mut |_| true).unwrap();
    let operation = Operation::InterpretCdnaAlleles {
        request: Box::new(request),
    };
    let denied = crate::mcp_server::mcp_tool_call_for_capability_surface_tests(
        state.to_str().unwrap(),
        "op",
        serde_json::json!({"operation":operation}),
    );
    assert_eq!(denied["isError"], true);
    let allowed = crate::mcp_server::mcp_tool_call_for_capability_surface_tests(
        state.to_str().unwrap(),
        "op",
        serde_json::json!({"confirm":true,"operation":operation}),
    );
    assert_ne!(allowed["isError"], true, "{allowed}");
    assert_eq!(
        allowed["structuredContent"]["result"]["cdna_allele_interpretation"],
        serde_json::to_value(&expected).unwrap()
    );
    assert!(!state.exists());
}

#[test]
fn cdna_interpretation_duplicate_batch_cumulative_budget_and_saved_projection_refusals() {
    let (root, request) = fixture();
    let before = snapshot(root.path());
    let mut duplicate = request.clone();
    duplicate.genes.push(duplicate.genes[0].clone());
    assert!(
        interpret_cdna_files(duplicate, &mut |_| true)
            .unwrap_err()
            .message
            .contains("duplicate_gene")
    );
    let mut small = request.clone();
    small.budgets.max_input_bytes = small.genes[0].comparison_report.bytes;
    assert!(
        interpret_cdna_files(small, &mut |_| true)
            .unwrap_err()
            .message
            .contains("cumulative_input_budget")
    );
    let mut tiny = request.clone();
    tiny.budgets.max_output_bytes = 1;
    assert!(
        interpret_cdna_files(tiny, &mut |_| true)
            .unwrap_err()
            .message
            .contains("output_budget")
    );
    let report = interpret_cdna_files(request, &mut |_| true).unwrap();
    let mut tampered = report.clone();
    tampered.genes[0].sites[0].alternate_units += 1;
    assert!(
        render(&tampered, Format::Summary)
            .unwrap_err()
            .message
            .contains("projection_content_hash")
    );
    // Negative scientific-limit text is required, not rejected by a blind keyword lint.
    for text in [
        render(&report, Format::Json).unwrap(),
        render(&report, Format::Tsv).unwrap(),
        render(&report, Format::Summary).unwrap(),
    ] {
        assert!(
            text.contains("not DNA genotype"),
            "missing explicit negative nonclaim: {text}"
        );
        for forbidden in [
            "heterozygous",
            "homozygous",
            "clinical_significance",
            "pathogenicity_classification",
        ] {
            assert!(
                !text.contains(forbidden),
                "invented biological/clinical field: {forbidden}"
            );
        }
    }
    assert_eq!(snapshot(root.path()), before);
}

#[cfg(feature = "lua-interface")]
#[test]
fn cdna_interpretation_lua_vm_matches_bound_json_without_state_mutation() {
    use mlua::LuaSerdeExt;
    let (_root, request) = fixture();
    let expected = interpret_cdna_files(request.clone(), &mut |_| true).unwrap();
    let operation = Operation::InterpretCdnaAlleles {
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
    assert_eq!(
        result["result"]["cdna_allele_interpretation"],
        serde_json::to_value(&expected).unwrap()
    );
}

#[cfg(feature = "js-interface")]
#[test]
fn cdna_interpretation_js_vm_matches_bound_json_and_preserves_saved_state() {
    let (root, request) = fixture();
    let expected = interpret_cdna_files(request.clone(), &mut |_| true).unwrap();
    let operation = Operation::InterpretCdnaAlleles {
        request: Box::new(request),
    };
    let project = root.path().join("C2 JS O'Brien.json");
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
        const input=load_project({path});const out=apply_operation(input,{operation});
        if (JSON.stringify(canon(input)) !== JSON.stringify(canon(out.state))) throw Error('interpretation changed state');
        if (JSON.stringify(canon(out.result.cdna_allele_interpretation)) !== JSON.stringify(canon({expected}))) throw Error('interpretation differs');
    "#,path=serde_json::to_string(&project.to_string_lossy()).unwrap(),operation=serde_json::to_string(&operation).unwrap(),expected=serde_json::to_string(&expected).unwrap())).unwrap();
    assert_eq!(fs::read(&project).unwrap(), before);
}
