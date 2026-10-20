//! Hand-calculated synthetic integrals; legacy digest uses the original v1
//! serializer as an independent byte oracle. No experimental data is implied.

use super::*;
use crate::engine::{Engine, GentleEngine, Operation};

fn profile(id: &str, weight: f64) -> GelProfile {
    GelProfile {
        id: id.into(),
        label: id.into(),
        run_id: "synthetic".into(),
        lane_id: id.into(),
        source_sha256: "a".repeat(64),
        evidence_kind: GelEvidenceKind::VendorTrace,
        size_kind: GelSizeKind::LinearDnaBp,
        signal_units: "synthetic_area".into(),
        edges: vec![100.0, 300.0],
        weights: vec![Some(weight)],
        measure: "synthetic_bin_integral".into(),
        processing: vec!["hand-calculated synthetic integrals".into()],
        warnings: vec![],
    }
}

fn request(id: &str, left: &str, right: &str) -> GelAlgebraRequest {
    GelAlgebraRequest {
        output_id: id.into(),
        terms: vec![
            GelAlgebraTerm {
                profile_id: left.into(),
                coefficient: 1.0,
            },
            GelAlgebraTerm {
                profile_id: right.into(),
                coefficient: -1.0,
            },
        ],
        range: [100.0, 300.0],
    }
}

fn store() -> GelWorkbenchStore {
    let a = profile("a", 4.0);
    let b = profile("b", 2.0);
    let residual = super::super::algebra(&[&a, &b], &request("r", "a", "b")).unwrap();
    GelWorkbenchStore {
        profiles: [a, b, residual]
            .into_iter()
            .map(|p| (p.id.clone(), Arc::new(p)))
            .collect(),
        ..Default::default()
    }
}

fn engine(store: GelWorkbenchStore) -> GentleEngine {
    let original = GentleEngine::new();
    let mut value = serde_json::to_value(original.state()).unwrap();
    value["gel_images"]["workbench"] = serde_json::to_value(store).unwrap();
    GentleEngine::from_state(serde_json::from_value(value).unwrap())
}

fn export(
    engine: &mut GentleEngine,
    path: &std::path::Path,
    ids: Vec<String>,
) -> Result<crate::engine::OpResult, EngineError> {
    engine.apply(Operation::GelWorkbench {
        request: Box::new(GelWorkbenchRequest::Export(GelWorkbenchExportRequest {
            path: path.display().to_string(),
            analysis_ids: vec![],
            comparison_ids: vec![],
            algebra_ids: ids,
            comparison_display: None,
            title: "Synthetic signed algebra; not experimental enrichment".into(),
        })),
    })
}

#[test]
fn algebra_v1_keeps_legacy_digest_and_processing_bytes() {
    let s = store();
    let a = s.profiles["a"].as_ref();
    let b = s.profiles["b"].as_ref();
    let req = request("r", "a", "b");
    let legacy = sha256_hex_bytes(&serde_json::to_vec(&([a, b], &req)).unwrap());
    assert_eq!(s.profiles["r"].source_sha256, legacy);
    assert_eq!(
        legacy,
        "742d8517ddeff4b3f1e3d47555b22f9065e4f604f4fa5bb5a4629330d16b4a7d"
    );
    assert_eq!(
        s.profiles["r"].processing,
        vec![serde_json::to_string(&req).unwrap()]
    );
    let mut legacy_profile = a.clone();
    legacy_profile.id = "r".into();
    legacy_profile.label = "Signed algebra r".into();
    legacy_profile.source_sha256 = legacy;
    legacy_profile.evidence_kind = GelEvidenceKind::SignedResidual;
    legacy_profile.weights = vec![Some(2.0)];
    legacy_profile.processing = vec![serde_json::to_string(&req).unwrap()];
    legacy_profile.warnings = vec!["Signed linear combination with explicit coefficients, not a gene-specific enrichment measurement or validated mixture model. Negative residuals retained; no distribution/p-value.".into()];
    assert_eq!(*s.profiles["r"], legacy_profile);
    assert_eq!(
        inspect_algebra(&s, "r").unwrap().verification,
        GelAlgebraVerification::Verified
    );
}

#[test]
fn algebra_replay_retains_all_operand_metadata_and_signed_density_after_reopen() {
    let mut s = store();
    Arc::make_mut(s.profiles.get_mut("b").unwrap())
        .warnings
        .push("second-operand calibration warning".into());
    let p = super::super::algebra(
        &[&s.profiles["a"], &s.profiles["b"]],
        &request("r", "a", "b"),
    )
    .unwrap();
    s.profiles.insert("r".into(), Arc::new(p));
    let reopened: GelWorkbenchStore =
        serde_json::from_value(serde_json::to_value(s).unwrap()).unwrap();
    let p = inspect_algebra(&reopened, "r").unwrap();
    assert_eq!(p.verification, GelAlgebraVerification::Verified);
    assert_eq!(p.series.as_ref().unwrap().density, vec![Some(0.01)]);
    assert_eq!(p.series.as_ref().unwrap().weights, vec![Some(2.0)]);
    assert!(
        p.warnings
            .iter()
            .any(|w| w == "second-operand calibration warning")
    );
    assert_eq!(p.inputs.len(), 2);
}

#[test]
fn algebra_replay_detects_result_operand_and_receipt_tampering() {
    for target in ["r", "a"] {
        let mut s = store();
        Arc::make_mut(s.profiles.get_mut(target).unwrap()).weights[0] = Some(99.0);
        let p = inspect_algebra(&s, "r").unwrap();
        assert_eq!(p.verification, GelAlgebraVerification::Mismatch);
        assert!(p.series.is_none());
    }
    let mut s = store();
    Arc::make_mut(s.profiles.get_mut("r").unwrap()).processing[0] = "{}".into();
    assert_eq!(
        inspect_algebra(&s, "r").unwrap().verification,
        GelAlgebraVerification::Unverifiable
    );
    s.profiles.remove("a");
    assert!(inspect_algebra(&s, "missing").is_err());
    assert!(inspect_algebra(&s, "a").is_err());
}

#[test]
fn algebra_nested_dag_is_memoized_and_cycles_depth_and_missing_operands_are_bounded() {
    let mut s = store();
    let r = s.profiles["r"].clone();
    let nested = super::super::algebra(&[&r, &r], &request("nested", "r", "r")).unwrap();
    s.profiles.insert("nested".into(), Arc::new(nested));
    let p = inspect_algebra(&s, "nested").unwrap();
    assert_eq!(p.verification, GelAlgebraVerification::Verified);
    assert_eq!(p.series.unwrap().density, vec![Some(0.0)]);
    assert_eq!(
        p.inputs.len(),
        3,
        "shared operands are not dependency cycles"
    );
    let mut cycle = s.clone();
    Arc::make_mut(cycle.profiles.get_mut("r").unwrap()).processing =
        vec![serde_json::to_string(&request("r", "r", "r")).unwrap()];
    assert_eq!(
        inspect_algebra(&cycle, "r").unwrap().verification,
        GelAlgebraVerification::Mismatch
    );
    let mut deep = s.clone();
    let mut previous = "nested".to_string();
    for n in 0..18 {
        let id = format!("deep{n}");
        let source = deep.profiles[&previous].clone();
        let p = super::super::algebra(&[&source, &source], &request(&id, &previous, &previous))
            .unwrap();
        deep.profiles.insert(id.clone(), Arc::new(p));
        previous = id;
    }
    assert_eq!(
        inspect_algebra(&deep, &previous).unwrap().verification,
        GelAlgebraVerification::Unverifiable
    );
    s.profiles.remove("b");
    assert_eq!(
        inspect_algebra(&s, "nested").unwrap().verification,
        GelAlgebraVerification::Unverifiable
    );
}

#[test]
fn algebra_replay_work_and_distinct_profile_budgets_fail_closed() {
    let mut s = GelWorkbenchStore::default();
    let mut base = profile("base", 1.0);
    base.edges = (100..=20_100).map(f64::from).collect();
    base.weights = vec![Some(1.0); 20_000];
    s.profiles.insert(base.id.clone(), Arc::new(base));
    let mut previous = "base".to_string();
    for n in 0..5 {
        let id = format!("repeat{n}");
        let req = GelAlgebraRequest {
            output_id: id.clone(),
            terms: (0..16)
                .map(|_| GelAlgebraTerm {
                    profile_id: previous.clone(),
                    coefficient: 1.0,
                })
                .collect(),
            range: [100.0, 20_100.0],
        };
        let input = s.profiles[&previous].clone();
        let refs = vec![input.as_ref(); 16];
        let p = super::super::algebra(&refs, &req).unwrap();
        s.profiles.insert(id.clone(), Arc::new(p));
        previous = id;
    }
    let presentation = inspect_algebra(&s, &previous).unwrap();
    assert_eq!(
        presentation.verification,
        GelAlgebraVerification::Unverifiable
    );
    assert!(presentation.reason.unwrap().contains("bin visits"));
    assert!(presentation.series.is_none());

    let mut s = GelWorkbenchStore::default();
    let mut previous = Vec::new();
    for group in 0..5 {
        let mut req = request(&format!("group{group}"), "unused", "unused");
        req.terms.clear();
        for n in 0..13 {
            let id = format!("input{group}_{n}");
            s.profiles.insert(id.clone(), Arc::new(profile(&id, 1.0)));
            req.terms.push(GelAlgebraTerm {
                profile_id: id,
                coefficient: 1.0,
            });
        }
        let refs: Vec<_> = req
            .terms
            .iter()
            .map(|t| s.profiles[&t.profile_id].as_ref())
            .collect();
        let p = super::super::algebra(&refs, &req).unwrap();
        previous.push(req.output_id.clone());
        s.profiles.insert(req.output_id, Arc::new(p));
    }
    let mut req = request("many", "unused", "unused");
    req.terms = previous
        .iter()
        .map(|id| GelAlgebraTerm {
            profile_id: id.clone(),
            coefficient: 1.0,
        })
        .collect();
    let refs: Vec<_> = previous.iter().map(|id| s.profiles[id].as_ref()).collect();
    let p = super::super::algebra(&refs, &req).unwrap();
    s.profiles.insert("many".into(), Arc::new(p));
    let presentation = inspect_algebra(&s, "many").unwrap();
    assert_eq!(
        presentation.verification,
        GelAlgebraVerification::Unverifiable
    );
    assert!(presentation.reason.unwrap().contains("64 profiles"));
    assert!(presentation.series.is_none());
}

#[test]
fn algebra_density_is_flat_across_offset_bins_and_keeps_zero_coefficient_gaps() {
    let mut a = profile("a", 4.0);
    let mut b = profile("b", 2.0);
    a.edges = vec![100.0, 200.0, 300.0];
    a.weights = vec![Some(4.0), Some(4.0)];
    b.edges = vec![100.0, 150.0, 300.0];
    b.weights = vec![Some(1.0), Some(3.0)];
    let p = super::super::algebra(&[&a, &b], &request("r", "a", "b")).unwrap();
    let mut s = GelWorkbenchStore {
        profiles: [a, b, p]
            .into_iter()
            .map(|p| (p.id.clone(), Arc::new(p)))
            .collect(),
        ..Default::default()
    };
    assert_eq!(
        inspect_algebra(&s, "r").unwrap().series.unwrap().density,
        vec![Some(0.02), Some(0.02), Some(0.02)]
    );
    Arc::make_mut(s.profiles.get_mut("b").unwrap()).weights = vec![None, None];
    let mut req = request("zero", "a", "b");
    req.terms[1].coefficient = 0.0;
    let p = super::super::algebra(&[&s.profiles["a"], &s.profiles["b"]], &req).unwrap();
    s.profiles.insert("zero".into(), Arc::new(p));
    let p = inspect_algebra(&s, "zero").unwrap();
    assert!(p.series.unwrap().density.iter().all(Option::is_some));
    assert!(p.warnings.iter().any(|w| w.contains("Zero-coefficient")));
}

#[test]
fn algebra_all_missing_is_verified_unavailable_not_a_zero_distribution() {
    let mut a = profile("a", 0.0);
    a.weights[0] = None;
    let p = super::super::algebra(&[&a, &a], &request("r", "a", "a")).unwrap();
    let s = GelWorkbenchStore {
        profiles: [a, p]
            .into_iter()
            .map(|p| (p.id.clone(), Arc::new(p)))
            .collect(),
        ..Default::default()
    };
    let p = inspect_algebra(&s, "r").unwrap();
    assert_eq!(p.verification, GelAlgebraVerification::Verified);
    assert_eq!(p.series.unwrap().density, vec![None]);
    assert!(p.warnings.iter().any(|w| w.contains("No evaluated bins")));
}

#[test]
fn algebra_long_receipts_fail_explicitly_without_changing_v1_limits() {
    let profiles: Vec<_> = (0..16)
        .map(|n| profile(&format!("{n}{}", "x".repeat(220)), 1.0))
        .collect();
    let request = GelAlgebraRequest {
        output_id: "r".into(),
        range: [100.0, 300.0],
        terms: profiles
            .iter()
            .map(|p| GelAlgebraTerm {
                profile_id: p.id.clone(),
                coefficient: 1.0,
            })
            .collect(),
    };
    let refs: Vec<_> = profiles.iter().collect();
    assert!(
        super::super::algebra(&refs, &request)
            .unwrap_err()
            .message
            .contains("2048 bytes")
    );
}

#[test]
fn algebra_shared_inspection_is_read_only_and_json_export_is_saved_verified_and_no_clobber() {
    let mut engine = engine(store());
    let before = serde_json::to_value(engine.state()).unwrap();
    let command = crate::engine_shell::parse_shell_line(
        "gel-image workbench '{\"action\":\"inspect_algebra\",\"profile_id\":\"r\"}'",
    )
    .unwrap();
    assert!(!command.is_state_mutating());
    for (payload, valid) in [
        (
            serde_json::json!({"action":"algebra","output_id":"next","range":[100,300],"terms":[{"profile_id":"a","coefficient":1}]}),
            true,
        ),
        (
            serde_json::json!({"action":"export","path":"new.json","title":"Review before writing","analysis_ids":[],"algebra_ids":["r"]}),
            true,
        ),
        (
            serde_json::json!({"action":"import_profile","profile":profile("next", 1.0)}),
            true,
        ),
        (
            serde_json::json!({"action":"inspect_algebra","profile_id":"r","allow_unverified":true}),
            false,
        ),
    ] {
        let payload = serde_json::json!({"GelWorkbench":{"request":payload}}).to_string();
        assert_eq!(serde_json::from_str::<Operation>(&payload).is_ok(), valid);
        let command = crate::engine_shell::ShellCommand::Op { payload };
        assert!(
            command.is_state_mutating(),
            "Only the named valid inspection is exempt"
        );
    }
    let output = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
    assert!(!output.state_changed);
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("saved.json");
    export(&mut engine, &out, vec!["r".into()]).unwrap();
    let bytes = std::fs::read(&out).unwrap();
    let json: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(json["algebra"][0]["verification"], "verified");
    assert_eq!(json["algebra"][0]["series"]["density"][0], 0.01);
    assert_eq!(json["algebra"][0]["inputs"].as_array().unwrap().len(), 2);
    assert_eq!(before, serde_json::to_value(engine.state()).unwrap());
    assert!(export(&mut engine, &out, vec!["r".into()]).is_err());
    assert_eq!(std::fs::read(&out).unwrap(), bytes);
}

#[test]
fn algebra_import_and_failed_verification_never_create_a_report_or_mutate_state() {
    let mut e = engine(store());
    let mut imported = e.state().gel_images.workbench.profiles["r"]
        .as_ref()
        .clone();
    imported.id = "forged".into();
    let before = serde_json::to_value(e.state()).unwrap();
    assert!(
        e.apply(Operation::GelWorkbench {
            request: Box::new(GelWorkbenchRequest::ImportProfile { profile: imported })
        })
        .is_err()
    );
    assert_eq!(before, serde_json::to_value(e.state()).unwrap());
    for mismatching in [true, false] {
        let mut s = store();
        let p = Arc::make_mut(s.profiles.get_mut("r").unwrap());
        if mismatching {
            p.weights[0] = Some(99.0);
        } else {
            p.processing = vec!["Legacy imported residual".into()];
        }
        let mut e = engine(s);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("refused.pdf");
        assert!(export(&mut e, &path, vec!["r".into()]).is_err());
        assert!(!path.exists());
        let report = e
            .apply(Operation::GelWorkbench {
                request: Box::new(GelWorkbenchRequest::InspectAlgebra {
                    profile_id: "r".into(),
                }),
            })
            .unwrap()
            .gel_workbench
            .unwrap()
            .algebra
            .unwrap();
        assert_eq!(
            report.verification,
            if mismatching {
                GelAlgebraVerification::Mismatch
            } else {
                GelAlgebraVerification::Unverifiable
            }
        );
    }
}

#[test]
#[ignore = "Manual vector-PDF visual review; set GENTLE_GEL_ALGEBRA_TEST_REPORT to a new path"]
fn write_verified_algebra_pdf_for_visual_review() {
    let mut s = store();
    let mut a = s.profiles["a"].as_ref().clone();
    let mut b = s.profiles["b"].as_ref().clone();
    for p in [&mut a, &mut b] {
        p.evidence_kind = GelEvidenceKind::Synthetic;
        p.source_sha256 = sha256_hex_bytes(b"GENtle synthetic gel algebra v1");
        p.processing.push("Source hash binds the UTF-8 recipe label GENtle synthetic gel algebra v1, without a trailing newline; no experimental source".into());
    }
    a.edges = vec![100.0, 150.0, 200.0, 300.0];
    a.weights = vec![Some(2.0), None, Some(6.0)];
    b.edges = a.edges.clone();
    b.weights = vec![Some(3.0), Some(1.0), Some(3.0)];
    let residual = super::super::algebra(&[&a, &b], &request("r", "a", "b")).unwrap();
    s.profiles = [a, b, residual]
        .into_iter()
        .map(|p| (p.id.clone(), Arc::new(p)))
        .collect();
    let output =
        std::env::var("GENTLE_GEL_ALGEBRA_TEST_REPORT").expect("Explicit new PDF path required");
    let mut e = engine(s);
    export(&mut e, std::path::Path::new(&output), vec!["r".into()]).unwrap();
    let doc = lopdf::Document::load(&output).unwrap();
    assert!(doc.get_pages().len() >= 2);
}

#[test]
fn algebra_tutorial_commands_replay_through_the_shared_shell() {
    let prose = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("docs/tutorial/10-03_measured_gel_sizing_gui_cli.md"),
    )
    .unwrap();
    let lesson = prose
        .split("<!-- gel-algebra-lesson:start -->")
        .nth(1)
        .unwrap()
        .split("<!-- gel-algebra-lesson:end -->")
        .next()
        .unwrap();
    let mut e = GentleEngine::new();
    let dir = tempfile::tempdir().unwrap();
    let mut count = 0;
    for line in lesson.lines().filter(|l| l.starts_with("gentle_cli ")) {
        let shell = line
            .strip_prefix("gentle_cli --state /tmp/gentle-algebra-lesson/state.json ")
            .unwrap();
        let json = shell
            .strip_prefix("gel-image workbench '")
            .unwrap()
            .strip_suffix('\'')
            .unwrap();
        let mut request: GelWorkbenchRequest = serde_json::from_str(json).unwrap();
        if let GelWorkbenchRequest::ImportProfile { profile } = &request {
            assert_eq!(profile.evidence_kind, GelEvidenceKind::Synthetic);
            assert_eq!(
                profile.source_sha256,
                sha256_hex_bytes(b"GENtle synthetic gel algebra v1")
            );
        }
        if let GelWorkbenchRequest::Export(request) = &mut request {
            request.path = dir
                .path()
                .join(std::path::Path::new(&request.path).file_name().unwrap())
                .display()
                .to_string();
        }
        let line = format!(
            "gel-image workbench {}",
            crate::engine_shell::quote_shell_arg(&serde_json::to_string(&request).unwrap())
        );
        let command = crate::engine_shell::parse_shell_line(&line).unwrap();
        crate::engine_shell::execute_shell_command(&mut e, &command).unwrap();
        count += 1;
    }
    assert_eq!(count, 6);
    let p = inspect_algebra(&e.state().gel_images.workbench, "teaching_residual").unwrap();
    let s = p.series.unwrap();
    assert_eq!(s.weights, vec![Some(-1.0), None, Some(3.0)]);
    assert_eq!(s.density, vec![Some(-0.02), None, Some(0.03)]);
    assert!(dir.path().join("residual.pdf").exists());
    assert_eq!(
        lopdf::Document::load(dir.path().join("residual.pdf"))
            .unwrap()
            .get_pages()
            .len(),
        2
    );
}
