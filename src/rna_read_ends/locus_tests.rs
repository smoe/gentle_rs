//! Synthetic saved-file locus controls, including raw LF/CRLF bindings and original A/T letters.

use super::*;
pub(super) fn fixture(ending: &str) -> (tempfile::TempDir, locus::Build) {
    let (root, mut old) = tests::fixture();
    tests::attach_original_bases(&root, &mut old, ending);
    let sample = old.samples.remove(0);
    let annotation: PrimerPayloadAnnotation =
        serde_json::from_slice(&fs::read(&sample.normalized_annotation.path).unwrap()).unwrap();
    let transcript = &annotation.transcripts[0];
    let region = core::locus::Region {
        contig: transcript.genomic_contig.clone().unwrap(),
        contig_length: 1000,
        exon: transcript.genomic_exons[0],
        flank_left: 100,
        flank_right: 300,
        anchor_transcript_id: Some(transcript.transcript_id.clone()),
    };
    let request = locus::Build {
        schema: locus::BUILD_SCHEMA.into(),
        source_id: old.source_id,
        region,
        samples: vec![locus::Sample {
            sample_id: sample.sample_id,
            analysis: sample.analysis,
            annotation: Some(sample.normalized_annotation),
            model: Some(sample.model),
            original_bases: sample.original_bases,
        }],
        terminal_runs: oriented::Policy::default(),
        limits: old.limits,
        output_dir: root
            .path()
            .canonicalize()
            .unwrap()
            .join("locus run's α")
            .to_string_lossy()
            .into_owned(),
        max_output_bytes: 16 * 1024 * 1024,
        cancel_path: old.cancel_path,
    };
    (root, request)
}
fn command(request: locus::Build) -> Command {
    Command::Locus {
        command: Box::new(locus::Command::Build {
            request: Box::new(request),
        }),
    }
}
#[test]
fn read_end_locus_service_lf_crlf_original_runs_roundtrip_and_no_clobber() {
    for ending in ["\n", "\r\n"] {
        let (root, request) = fixture(ending);
        let original = request.samples[0].original_bases.clone().unwrap();
        let before = fs::read(&original.path).unwrap();
        let cmd = command(request);
        let out = execute(&GentleEngine::new(), cmd.clone(), &mut |_| true).unwrap();
        assert!(execute(&GentleEngine::new(), cmd, &mut |_| true).is_err());
        let report: PrimerPayloadFile = serde_json::from_value(out["report"].clone()).unwrap();
        let inspect = Command::Locus {
            command: Box::new(locus::Command::Inspect {
                report: report.clone(),
                view: core::locus::View::default(),
            }),
        };
        assert!(!inspect.requires_confirmation());
        let value = execute(&GentleEngine::new(), inspect, &mut |_| true).unwrap();
        assert!(!value["page"]["reads"].as_array().unwrap().is_empty());
        let json = serde_json::to_string(&value).unwrap();
        assert!(json.contains("\"base\":\"A\"") && json.contains("\"base\":\"T\""));
        assert!(json.contains(core::locus::RULE_VERSION));
        let export = Command::Locus {
            command: Box::new(locus::Command::Export {
                report,
                view: core::locus::View::default(),
                formats: vec![Format::Json, Format::Tsv, Format::Svg],
                output_dir: root
                    .path()
                    .canonicalize()
                    .unwrap()
                    .join("exports ' α")
                    .to_string_lossy()
                    .into_owned(),
                max_output_bytes: 16 * 1024 * 1024,
                cancel_path: None,
            }),
        };
        execute(&GentleEngine::new(), export.clone(), &mut |_| true).unwrap();
        assert!(execute(&GentleEngine::new(), export, &mut |_| true).is_err());
        assert_eq!(before, fs::read(original.path).unwrap());
    }
}
#[test]
fn read_end_locus_missing_annotation_preserves_genome_but_bad_hash_refuses() {
    let (_root, mut request) = fixture("\n");
    request.samples[0].annotation.as_mut().unwrap().sha256 = format!("sha256:{}", "f".repeat(64));
    assert!(
        execute(&GentleEngine::new(), command(request.clone()), &mut |_| {
            true
        })
        .is_err()
    );
    request.samples[0].annotation = None;
    let correct_model = request.samples[0].model.as_ref().unwrap().sha256.clone();
    request.samples[0].model.as_mut().unwrap().sha256 = format!("sha256:{}", "f".repeat(64));
    assert!(
        execute(&GentleEngine::new(), command(request.clone()), &mut |_| {
            true
        })
        .is_err()
    );
    request.samples[0].model.as_mut().unwrap().sha256 = correct_model;
    let out = execute(&GentleEngine::new(), command(request), &mut |_| true).unwrap();
    let report: PrimerPayloadFile = serde_json::from_value(out["report"].clone()).unwrap();
    let record: core::locus::Report = completed_record(&report, &mut || Ok(())).unwrap();
    assert!(!record.reads.is_empty() && record.models.is_empty());
    assert!(
        record
            .reads
            .iter()
            .all(|r| r.genomic_basecalls.iter().any(|e| !e.variants.is_empty()))
    );
}
#[test]
fn read_end_locus_cancel_is_request_bound_and_cannot_replace_inputs() {
    let (_root, mut request) = fixture("\n");
    request.cancel_path = Some(
        request.samples[0]
            .original_bases
            .as_ref()
            .unwrap()
            .path
            .clone(),
    );
    assert!(marker_path(&command(request.clone())).is_err());
    request.cancel_path = Some(
        Path::new(&request.samples[0].analysis.report_path)
            .parent()
            .unwrap()
            .join("new-control")
            .to_string_lossy()
            .into_owned(),
    );
    assert!(marker_path(&command(request.clone())).is_err());
    request.cancel_path = Some(format!("{}.cancel", request.output_dir));
    let cmd = command(request);
    execute(
        &GentleEngine::new(),
        Command::Cancel {
            command: Box::new(cmd.clone()),
        },
        &mut |_| true,
    )
    .unwrap();
    assert!(execute(&GentleEngine::new(), cmd, &mut |_| true).is_err());
}

#[test]
fn read_end_locus_declared_absent_annotation_and_cds_are_distinct_from_omitted_and_bad_hash() {
    for annotation_missing in [true, false] {
        let (root, mut request) = fixture("\n");
        let file = if annotation_missing {
            request.samples[0].annotation.as_mut().unwrap()
        } else {
            request.samples[0].model.as_mut().unwrap()
        };
        let original = file.clone();
        let before = fs::read(&original.path).unwrap();
        file.path = root
            .path()
            .canonicalize()
            .unwrap()
            .join("declared absent ' α.json")
            .to_string_lossy()
            .into_owned();
        let out = execute(&GentleEngine::new(), command(request.clone()), &mut |_| {
            true
        })
        .unwrap();
        let binding: PrimerPayloadFile = serde_json::from_value(out["report"].clone()).unwrap();
        let report: core::locus::Report = completed_record(&binding, &mut || Ok(())).unwrap();
        let prefix = if annotation_missing {
            "annotation_declared_not_found"
        } else {
            "cds_model_declared_not_found"
        };
        assert!(
            report.samples[0]
                .diagnostics
                .contains(&format!("{prefix}:{}", original.sha256))
        );
        assert!(!report.reads.is_empty());
        if annotation_missing {
            assert!(report.models.is_empty());
        } else {
            assert!(report.reads.iter().flat_map(|r| &r.assignments).all(|a| {
                a.paths
                    .iter()
                    .all(|p| matches!(&p.utr5, core::Metric::Missing { .. }))
            }));
        }
        request.output_dir.push_str("-bad-hash");
        let file = if annotation_missing {
            request.samples[0].annotation.as_mut().unwrap()
        } else {
            request.samples[0].model.as_mut().unwrap()
        };
        *file = original.clone();
        file.sha256 = format!("sha256:{}", "f".repeat(64));
        let err = execute(&GentleEngine::new(), command(request), &mut |_| true).unwrap_err();
        assert!(err.message.contains("file_hash_or_change"), "{err:?}");
        assert_eq!(fs::read(original.path).unwrap(), before);
    }
}

#[test]
fn read_end_locus_declared_absent_original_bases_still_refuse() {
    let (root, mut request) = fixture("\n");
    request.samples[0].original_bases.as_mut().unwrap().path = root
        .path()
        .canonicalize()
        .unwrap()
        .join("absent-originals.json")
        .to_string_lossy()
        .into_owned();
    let output = PathBuf::from(&request.output_dir);
    assert!(execute(&GentleEngine::new(), command(request), &mut |_| true).is_err());
    assert!(!output.exists());
}

#[test]
fn read_end_locus_json_only_export_survives_original_run_attachments_outside_window() {
    let (root, mut request) = fixture("\n");
    request.region.flank_left = 0;
    request.region.flank_right = 0;
    let out = execute(&GentleEngine::new(), command(request), &mut |_| true).unwrap();
    let report: PrimerPayloadFile = serde_json::from_value(out["report"].clone()).unwrap();
    let output = root
        .path()
        .canonicalize()
        .unwrap()
        .join("json only outside window");
    let cmd = Command::Locus {
        command: Box::new(locus::Command::Export {
            report,
            view: core::locus::View::default(),
            formats: vec![Format::Json],
            output_dir: output.to_string_lossy().into_owned(),
            max_output_bytes: 16 * 1024 * 1024,
            cancel_path: None,
        }),
    };
    execute(&GentleEngine::new(), cmd, &mut |_| true).unwrap();
    let page: core::locus::Page =
        serde_json::from_slice(&fs::read(output.join("artifacts/page.json")).unwrap()).unwrap();
    assert!(
        page.reads
            .iter()
            .flat_map(|r| &r.genomic_basecalls)
            .flat_map(|e| &e.variants)
            .flat_map(|v| &v.runs)
            .any(|r| r.reference_anchor.is_some_and(|a| a >= 110))
    );
}
