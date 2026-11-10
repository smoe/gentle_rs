//! Synthetic file/service/approval controls; no private reads, live tools or biological acceptance.

use super::*;
use crate::engine::{ObservedSamAnnotationBinding, Operation};

pub(super) fn fixture() -> (tempfile::TempDir, Build) {
    let (root, mut r) = GentleEngine::observed_sam_test_fixture();
    let path = root.path().canonicalize().unwrap();
    r.output_dir = path
        .join("saved S2's report")
        .to_string_lossy()
        .into_owned();
    let annotation = match &r.annotation_sources[0] {
        ObservedSamAnnotationBinding::Normalized { annotation, .. } => {
            let mut f = annotation.clone();
            f.path = Path::new(&f.path)
                .canonicalize()
                .unwrap()
                .to_string_lossy()
                .into_owned();
            f
        }
        _ => panic!("synthetic normalized fixture"),
    };
    let a: PrimerPayloadAnnotation =
        serde_json::from_reader(File::open(&annotation.path).unwrap()).unwrap();
    let model = core::ModelDocument {
        schema: core::MODEL_SCHEMA.into(),
        source_id: r.reference_gene.source_id.clone(),
        assembly: r.reference.assembly.clone(),
        reference_fasta_sha256: r.reference.fasta_sha256.clone(),
        annotation_sha256: annotation.sha256.clone(),
        coordinate_convention: "zero_based_half_open".into(),
        model_source: core::ModelSource {
            format: "gtf".into(),
            name: "synthetic model, not a biological annotation".into(),
            sha256: format!("sha256:{}", "f".repeat(64)),
        },
        transcripts: a
            .transcripts
            .iter()
            .map(|t| {
                let exons = t.genomic_exons.clone();
                let last = exons.last().unwrap();
                core::TranscriptModel {
                    transcript_id: t.transcript_id.clone(),
                    gene_id: r.reference_gene.gene_id.clone(),
                    contig: t.genomic_contig.clone().unwrap(),
                    strand: t.genomic_strand.unwrap(),
                    genomic_cds_0based_half_open: exons.clone(),
                    genomic_stop_codon_0based_half_open: vec![(last.1 - 3, last.1)],
                    genomic_exons_0based_half_open: exons,
                }
            })
            .collect(),
    };
    let model_path = path.join("synthetic model.json");
    fs::write(&model_path, json_bytes(&model).unwrap()).unwrap();
    let report = GentleEngine::new()
        .summarize_observed_transcript_structures_from_sam(r.clone(), &mut |_| true)
        .unwrap();
    let binding = GentleEngine::observed_sam_result_binding(&r, &report).unwrap();
    let build = Build {
        schema: "gentle.rna_read_end_build.v1".into(),
        source_id: r.reference_gene.source_id,
        gene_id: r.reference_gene.gene_id,
        transcript_id: model.transcripts[0].transcript_id.clone(),
        samples: vec![Sample {
            sample_id: "synthetic_sample".into(),
            analysis: binding,
            model: ingress::fingerprint(&model_path, CAP, &mut || Ok(())).unwrap(),
            normalized_annotation: annotation,
            composition: evidence::CompositionPolicy {
                orientation: evidence::Provenance::Undeclared,
                trimming: evidence::Provenance::Undeclared,
                adapter_sequences: vec![],
                primer_sequences: vec![],
                window_bp: 20,
                minimum_run_bp: 10,
                maximum_mismatches: 1,
            },
            tail: None,
            original_bases: None,
        }],
        exact_quantile_probabilities: vec![0.0, 0.5, 1.0],
        limits: core::Limits {
            max_units: 1000,
            max_observations: 10_000,
            ..core::Limits::default()
        },
        output_dir: path.join("read-end run's α").to_string_lossy().into_owned(),
        max_output_bytes: 8 * 1024 * 1024,
        cancel_path: Some(path.join("cancel.marker").to_string_lossy().into_owned()),
        terminal_runs: None,
    };
    (root, build)
}
pub(super) fn view() -> pileup::View {
    pileup::View {
        schema: pileup::VIEW_SCHEMA.into(),
        sort: pileup::SortKey::UnitKey,
        descending: false,
        filters: vec![],
        offset: 0,
        count: 10,
        max_sublanes: 100,
        width: 1400,
        scale: "genomic".into(),
        channels: pileup::Channels {
            utr5: true,
            utr3: true,
            read_length: true,
            clips: true,
            composition: true,
            tail: true,
        },
    }
}
pub(super) fn run(request: Build) -> PrimerPayloadFile {
    let out = Path::new(&request.output_dir).join("report.json");
    execute(
        &GentleEngine::new(),
        Command::Build {
            request: Box::new(request),
        },
        &mut |_| true,
    )
    .unwrap();
    ingress::fingerprint(&out, CAP, &mut || Ok(())).unwrap()
}

pub(super) fn attach_original_bases(root: &tempfile::TempDir, request: &mut Build, ending: &str) {
    let path = Path::new(&request.samples[0].analysis.report_path)
        .parent()
        .unwrap()
        .join("artifacts/analysis.json");
    let data: ObservedSamData = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    let analysis = data.analysis.unwrap();
    let bases = oriented::Document {
        schema: oriented::BASES_SCHEMA.into(),
        analysis_context_sha256: analysis.context_sha256,
        orientation: "original_sequencing_5_to_3".into(),
        provenance: "Synthetic original bases, not SAM reconstruction".into(),
        reads: analysis
            .units
            .iter()
            .filter_map(|u| {
                let length = u.unclipped_axis_length?;
                Some(oriented::OriginalRead {
                    unit: u.unit.clone(),
                    read_id: u.rows.first()?.qname.clone(),
                    sequence: format!("{}{}", "A".repeat(12), "T".repeat(length - 12)),
                })
            })
            .collect(),
    };
    let input = root
        .path()
        .canonicalize()
        .unwrap()
        .join("original read's bases α.json");
    fs::write(
        &input,
        serde_json::to_string_pretty(&bases)
            .unwrap()
            .replace('\n', ending),
    )
    .unwrap();
    request.samples[0].original_bases =
        Some(ingress::fingerprint(&input, CAP, &mut || Ok(())).unwrap());
}

#[test]
fn read_end_service_original_bases_are_bound_saved_and_no_clobber_under_lf_crlf() {
    for ending in ["\n", "\r\n"] {
        let (root, mut request) = fixture();
        attach_original_bases(&root, &mut request, ending);
        let input = request.samples[0]
            .original_bases
            .as_ref()
            .unwrap()
            .path
            .clone();
        let command = Command::Build {
            request: Box::new(request.clone()),
        };
        let receipt = run(request);
        let report = completed(&receipt, &mut || Ok(())).unwrap();
        let summary = report.orientation_summary.as_ref().unwrap();
        assert!(summary.units_with_bound_bases > 0);
        assert!(summary.units_with_terminal_a > 0 && summary.units_with_terminal_t > 0);
        let page = pileup::page(&report, &view()).unwrap();
        assert_eq!(page.total_gene_assigned, summary.gene_assigned_units);
        let before = fs::read(&receipt.path).unwrap();
        assert!(execute(&GentleEngine::new(), command, &mut |_| true).is_err());
        assert_eq!(fs::read(&receipt.path).unwrap(), before);
        fs::remove_file(input).unwrap();
        let inspected = execute(
            &GentleEngine::new(),
            Command::Inspect {
                report: receipt,
                view: view(),
            },
            &mut |_| true,
        )
        .unwrap();
        assert!(
            inspected["scene"]["rows"]
                .as_array()
                .unwrap()
                .iter()
                .any(|r| !r["runs"].as_array().unwrap().is_empty())
        );
    }
}

#[test]
fn read_end_service_tampered_original_bases_refuse_before_output() {
    let (root, mut request) = fixture();
    let path = root
        .path()
        .canonicalize()
        .unwrap()
        .join("changed bases.json");
    fs::write(&path, b"{}").unwrap();
    request.samples[0].original_bases =
        Some(ingress::fingerprint(&path, CAP, &mut || Ok(())).unwrap());
    fs::write(&path, b"{\"changed\":true}").unwrap();
    let out = request.output_dir.clone();
    assert!(
        execute(
            &GentleEngine::new(),
            Command::Build {
                request: Box::new(request)
            },
            &mut |_| true
        )
        .is_err()
    );
    assert!(!Path::new(&out).exists());
}

#[test]
fn read_end_service_build_inspect_export_is_saved_only_and_project_preserving() {
    let (root, request) = fixture();
    let engine = GentleEngine::new();
    let before = json_bytes(engine.state()).unwrap();
    let journal = engine.journal_len();
    let original = request.samples[0].model.path.clone();
    let report = run(request);
    fs::remove_file(original).unwrap();
    let inspected = execute(
        &engine,
        Command::Inspect {
            report: report.clone(),
            view: view(),
        },
        &mut |_| true,
    )
    .unwrap();
    assert!(
        inspected["page"]["coverage"]["samples"][0]["declared_units"]
            .as_u64()
            .unwrap()
            > 0
    );
    assert_eq!(json_bytes(engine.state()).unwrap(), before);
    assert_eq!(engine.journal_len(), journal);
    let output = root
        .path()
        .canonicalize()
        .unwrap()
        .join("fresh exports's α");
    let export = Command::Export {
        report,
        view: view(),
        formats: vec![Format::Json, Format::Tsv, Format::Svg],
        output_dir: output.to_string_lossy().into_owned(),
        max_output_bytes: 8 * 1024 * 1024,
        cancel_path: None,
    };
    execute(&engine, export.clone(), &mut |_| true).unwrap();
    let bytes = fs::read(output.join("artifacts/pileup.svg")).unwrap();
    assert!(bytes.starts_with(b"<svg"));
    assert!(execute(&engine, export, &mut |_| true).is_err());
    assert_eq!(
        fs::read(output.join("artifacts/pileup.svg")).unwrap(),
        bytes
    );
    fs::rename(&output, output.with_file_name("closed handles")).unwrap();
}

#[test]
fn read_end_service_raw_hashes_tampering_and_refused_model_geometry_are_not_weakened() {
    let (_root, request) = fixture();
    let model = request.samples[0].model.path.clone();
    fs::write(&model, b"{}").unwrap();
    assert!(
        execute(
            &GentleEngine::new(),
            Command::Build {
                request: Box::new(request.clone())
            },
            &mut |_| true
        )
        .is_err()
    );
    assert!(!Path::new(&request.output_dir).exists());
    let (_root, request) = fixture();
    let report = run(request);
    let artifact = Path::new(&report.path)
        .parent()
        .unwrap()
        .join("artifacts/pileup.json");
    fs::write(&artifact, b"{}").unwrap();
    assert!(
        execute(
            &GentleEngine::new(),
            Command::Inspect {
                report,
                view: view()
            },
            &mut |_| true
        )
        .is_err()
    );
}

#[test]
fn read_end_service_cancellation_is_request_bound_and_prevents_completion_publication() {
    let (_root, request) = fixture();
    let command = Command::Build {
        request: Box::new(request.clone()),
    };
    let result = execute(
        &GentleEngine::new(),
        Command::Cancel {
            command: Box::new(command.clone()),
        },
        &mut |_| true,
    )
    .unwrap();
    assert_eq!(result["termination_confirmed"], false);
    assert!(
        execute(&GentleEngine::new(), command, &mut |_| true)
            .unwrap_err()
            .message
            .contains("cancelled")
    );
    assert!(!Path::new(&request.output_dir).exists());
    let (_root, request) = fixture();
    assert!(
        execute(
            &GentleEngine::new(),
            Command::Build {
                request: Box::new(request.clone())
            },
            &mut |_| false
        )
        .is_err()
    );
    assert!(!Path::new(&request.output_dir).exists());
}

#[test]
fn read_end_service_shared_shell_parser_preserves_quote_paths_and_approval() {
    let (_root, request) = fixture();
    for command in [
        Command::Build {
            request: Box::new(request),
        },
        Command::Inspect {
            report: PrimerPayloadFile {
                path: "synthetic's report α.json".into(),
                bytes: 0,
                sha256: format!("sha256:{}", "a".repeat(64)),
            },
            view: view(),
        },
    ] {
        let text = serde_json::to_string(&command).unwrap();
        let shell = crate::engine_shell::parse_shell_line(&format!(
            "rna-reads read-ends {}",
            crate::engine_shell::quote_shell_arg(&text)
        ))
        .unwrap();
        assert_eq!(shell.is_state_mutating(), command.requires_confirmation());
        let crate::engine_shell::ShellCommand::Op { payload } = shell else {
            panic!("shared op");
        };
        let Operation::RnaReadEnds { command: decoded } = serde_json::from_str(&payload).unwrap()
        else {
            panic!("read-end op");
        };
        assert_eq!(
            json_bytes(&*decoded).unwrap(),
            json_bytes(&command).unwrap()
        );
    }
}

#[test]
fn read_end_service_lf_crlf_bind_raw_bytes_but_preserve_model_semantics() {
    let (_root, mut request) = fixture();
    let raw = serde_json::to_string_pretty(
        &serde_json::from_reader::<_, Value>(File::open(&request.samples[0].model.path).unwrap())
            .unwrap(),
    )
    .unwrap();
    let mut reports = vec![];
    for (i, bytes) in [raw.clone(), raw.replace('\n', "\r\n")]
        .into_iter()
        .enumerate()
    {
        fs::write(&request.samples[0].model.path, bytes).unwrap();
        request.samples[0].model =
            ingress::fingerprint(Path::new(&request.samples[0].model.path), CAP, &mut || {
                Ok(())
            })
            .unwrap();
        request.output_dir = Path::new(&request.output_dir)
            .with_file_name(format!("line-ending-{i}"))
            .to_string_lossy()
            .into_owned();
        let file = run(request.clone());
        let p = completed(&file, &mut || Ok(())).unwrap();
        reports.push(p);
    }
    assert_eq!(
        json_bytes(&reports[0].model).unwrap(),
        json_bytes(&reports[1].model).unwrap()
    );
    assert_ne!(
        reports[0].coverage.samples[0].bindings.model_sha256,
        reports[1].coverage.samples[0].bindings.model_sha256
    );
}

#[test]
fn read_end_service_raster_budget_checks_boundary_and_multiplication_overflow() {
    check_raster_pixel_budget(1000, 32_000).unwrap();
    for (width, height) in [(1000, 32_001), (8000, u64::MAX)] {
        assert!(
            check_raster_pixel_budget(width, height)
                .unwrap_err()
                .message
                .contains("raster_pixel_budget")
        );
    }
}

#[test]
fn read_end_service_raster_budget_uses_wrapped_headers_and_lanes_before_conversion() {
    use gentle_render::read_end_pileup::{self, Row, Scene};
    let long = "synthetic status ".repeat(300);
    for (long_headers, long_lanes) in [(true, false), (false, true), (true, true)] {
        let scene = Scene {
            contig: "synthetic_contig".into(),
            start: 100,
            end: 200,
            strand: 1,
            width: 8000,
            header: if long_headers {
                vec![long.clone(); 64]
            } else {
                vec!["synthetic header".into()]
            },
            model: vec![],
            top_tracks: vec![],
            bottom_tracks: vec![],
            axis_kind: None,
            rows: vec![
                Row {
                    label: if long_lanes {
                        "sample-unit-".repeat(16)
                    } else {
                        "synthetic unit".into()
                    },
                    details: "Hand-crafted layout budget control; not biological evidence".into(),
                    off_axis: if long_lanes {
                        long.clone()
                    } else {
                        "observed: 0 bp".into()
                    },
                    bands: vec![],
                    orientation: None,
                    sublanes: vec![],
                    runs: vec![],
                };
                if long_lanes { 32 } else { 1 }
            ],
        };
        let (width, height) = read_end_pileup::dimensions(&scene).unwrap();
        let old_height = 30 + 19 * scene.header.len() + 95 + 54 * scene.rows.len();
        assert!(u64::from(width) * old_height as u64 <= RASTER_PIXEL_LIMIT);
        assert!(u64::from(width) * height > RASTER_PIXEL_LIMIT);
        for format in [Format::Png, Format::Pdf] {
            assert!(
                export_svg(&scene, &[Format::Svg, format])
                    .unwrap_err()
                    .message
                    .contains("raster_pixel_budget")
            );
        }
        // Vector/text-only output remains allowed; no raster converter is called by this test.
        let svg = export_svg(&scene, &[Format::Json, Format::Tsv, Format::Svg]).unwrap();
        let root = svg.split('>').next().unwrap();
        assert!(root.contains(&format!(" height=\"{height}\"")));
    }
}

#[test]
fn read_end_service_png_pdf_derive_from_the_same_svg_with_explicit_dimensions() {
    let (root, mut request) = fixture();
    attach_original_bases(&root, &mut request, "\n");
    let report = run(request);
    let saved = completed(&report, &mut || Ok(())).unwrap();
    let page = pileup::page(&saved, &view()).unwrap();
    let scene = pileup::scene(&saved, &page).unwrap();
    let (width, height) = gentle_render::read_end_pileup::dimensions(&scene).unwrap();
    let output = root.path().canonicalize().unwrap().join("raster vector");
    execute(
        &GentleEngine::new(),
        Command::Export {
            report,
            view: view(),
            formats: vec![Format::Svg, Format::Png, Format::Pdf],
            output_dir: output.to_string_lossy().into_owned(),
            max_output_bytes: 16 * 1024 * 1024,
            cancel_path: None,
        },
        &mut |_| true,
    )
    .unwrap();
    let png = fs::read(output.join("artifacts/pileup.png")).unwrap();
    assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
    assert_eq!(u32::from_be_bytes(png[16..20].try_into().unwrap()), width);
    assert_eq!(
        u64::from(u32::from_be_bytes(png[20..24].try_into().unwrap())),
        height
    );
    let svg = fs::read_to_string(output.join("artifacts/pileup.svg")).unwrap();
    assert!(
        svg.split('>')
            .next()
            .unwrap()
            .contains(&format!(" height=\"{height}\""))
    );
    assert!(
        fs::read(output.join("artifacts/pileup.pdf"))
            .unwrap()
            .starts_with(b"%PDF-")
    );
}

#[test]
fn read_end_service_original_input_cannot_be_used_as_cancellation_marker() {
    let (root, mut request) = fixture();
    attach_original_bases(&root, &mut request, "\n");
    let file = request.samples[0].original_bases.as_ref().unwrap();
    let before = fs::read(&file.path).unwrap();
    request.cancel_path = Some(file.path.clone());
    let original_path = file.path.clone();
    let out = request.output_dir.clone();
    let command = Command::Cancel {
        command: Box::new(Command::Build {
            request: Box::new(request),
        }),
    };
    assert!(
        execute(&GentleEngine::new(), command, &mut |_| true)
            .unwrap_err()
            .message
            .contains("control_collision")
    );
    assert_eq!(fs::read(original_path).unwrap(), before);
    assert!(!Path::new(&out).exists());
}

#[test]
fn read_end_service_oriented_glyph_heights_participate_in_raster_guard() {
    use gentle_render::read_end_pileup::{self, Orientation, Row, RunGlyph, Scene};
    let scene = Scene {
        contig: "synthetic".into(),
        start: 100,
        end: 200,
        strand: 1,
        width: 8000,
        header: vec![],
        model: vec![],
        top_tracks: vec![],
        bottom_tracks: vec![],
        axis_kind: None,
        rows: vec![
            Row {
                label: "synthetic_read".into(),
                details: "synthetic original reads".into(),
                off_axis: "Not pore direction or physical tail length".into(),
                bands: vec![],
                orientation: Some(Orientation::Same),
                sublanes: vec![],
                runs: vec![
                    RunGlyph {
                        base: 'T',
                        length: 18,
                        reference_anchor: 120,
                        aligned_to_reference: false,
                        label: "18T original query [112,130), unaligned, not a tail length".into(),
                    };
                    2
                ],
            };
            64
        ],
    };
    assert!(8000 * (30 + 95 + 54 * 64) < RASTER_PIXEL_LIMIT);
    let (width, height) = read_end_pileup::dimensions(&scene).unwrap();
    assert!(u64::from(width) * height > RASTER_PIXEL_LIMIT);
    for format in [Format::Png, Format::Pdf] {
        assert!(
            export_svg(&scene, &[format])
                .unwrap_err()
                .message
                .contains("raster_pixel_budget")
        );
    }
    assert!(export_svg(&scene, &[Format::Svg]).unwrap().contains("18T"));
}
