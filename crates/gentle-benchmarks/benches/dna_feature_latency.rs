//! Synthetic feature-density audit. CPU preparation/paint only, never native latency.
//! See benches/README.md for provenance, prebuilt execution and auditor ownership.

use criterion::{BatchSize, BenchmarkId, Criterion};
use gentle::{
    about::{GENTLE_GIT_COMMIT, GENTLE_SOURCE_REVISION},
    dna_sequence::DNAsequence,
    main_area_dna::{
        MainAreaDna,
        latency_benchmark::{
            DNA_LATENCY_WORKLOAD, DnaLatencyInteraction, INTERACTIVE_BOUNDARY,
            feature_density_boundary_fixture, feature_density_fixture, feature_density_workloads,
        },
    },
};
use ring::digest::{SHA256, digest};
use serde_json::json;
use std::hint::black_box;

fn content_hash(dna: &DNAsequence, derived_layers: bool) -> String {
    let bytes = if derived_layers {
        serde_json::to_vec(&(
            dna.forward_bytes(),
            dna.features(),
            dna.restriction_enzymes(),
            dna.methylation_mode(),
        ))
        .unwrap()
    } else {
        // Preserve the existing nine sequence/annotation fixture identities.
        serde_json::to_vec(&(dna.forward_bytes(), dna.features())).unwrap()
    };
    digest(&SHA256, &bytes)
        .as_ref()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

fn area(dna: &DNAsequence, tree_loaded: bool, derived_layers: bool) -> MainAreaDna {
    let mut area = MainAreaDna::new(dna.clone(), None, None);
    area.prepare_latency_benchmark(tree_loaded);
    if derived_layers {
        area.enable_latency_benchmark_derived_layers();
    }
    area
}

fn frame(ctx: &egui::Context, area: &mut MainAreaDna, size: egui::Vec2, hover: bool) -> usize {
    let mut input = egui::RawInput {
        screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, size)),
        ..Default::default()
    };
    if hover {
        input.events.push(egui::Event::PointerMoved(
            area.latency_benchmark_hover_position().unwrap(),
        ));
    }
    let mut output = ctx.run_ui(input, |ui| area.render_inside(ui));
    let shapes = output.shapes.len();
    output.textures_delta.clear();
    assert!(shapes > 0);
    shapes
}

fn warm(dna: &DNAsequence, derived_layers: bool) -> (egui::Context, MainAreaDna) {
    let ctx = egui::Context::default();
    let mut area = area(dna, true, derived_layers);
    for _ in 0..2 {
        frame(&ctx, &mut area, egui::vec2(1200.0, 800.0), false);
    }
    (ctx, area)
}

fn benchmark(c: &mut Criterion) {
    for (length, count) in feature_density_workloads() {
        let derived_layers = (length, count) == INTERACTIVE_BOUNDARY;
        let dna = if derived_layers {
            feature_density_boundary_fixture()
        } else {
            feature_density_fixture(length, count)
        };
        let hash = content_hash(&dna, derived_layers);
        let id = format!("{length}bp_{count}features_{}", &hash[..12]);
        let mut group = c.benchmark_group("dna_feature_latency");
        group.bench_function(BenchmarkId::new(&id, "constructor"), |b| {
            b.iter_batched(
                || dna.clone(),
                |dna| black_box(MainAreaDna::new(dna, None, None)),
                BatchSize::LargeInput,
            );
        });
        group.bench_function(BenchmarkId::new(&id, "hydration"), |b| {
            b.iter_batched(
                || {
                    (
                        area(
                            &DNAsequence::from_sequence("").unwrap(),
                            false,
                            derived_layers,
                        ),
                        dna.clone(),
                    )
                },
                |(mut area, dna)| {
                    area.replace_loaded_sequence(dna);
                    black_box(area)
                },
                BatchSize::LargeInput,
            );
        });
        for tree_loaded in [false, true] {
            let name = if tree_loaded {
                "first_frame_tree_loaded"
            } else {
                "first_frame_tree_deferred"
            };
            group.bench_function(BenchmarkId::new(&id, name), |b| {
                b.iter_batched(
                    || {
                        (
                            egui::Context::default(),
                            area(&dna, tree_loaded, derived_layers),
                        )
                    },
                    |(ctx, mut area)| {
                        black_box(frame(&ctx, &mut area, egui::vec2(1200.0, 800.0), false))
                    },
                    BatchSize::LargeInput,
                );
            });
        }
        for (name, interaction, size, hover) in [
            (
                "steady",
                DnaLatencyInteraction::Steady,
                [1200.0, 800.0],
                false,
            ),
            (
                "pan_1bp",
                DnaLatencyInteraction::Pan,
                [1200.0, 800.0],
                false,
            ),
            ("zoom", DnaLatencyInteraction::Zoom, [1200.0, 800.0], false),
            (
                "toggle_mrna",
                DnaLatencyInteraction::ToggleMrna,
                [1200.0, 800.0],
                false,
            ),
            (
                "select",
                DnaLatencyInteraction::Select,
                [1200.0, 800.0],
                false,
            ),
            (
                "hover",
                DnaLatencyInteraction::Steady,
                [1200.0, 800.0],
                true,
            ),
            (
                "resize_compact",
                DnaLatencyInteraction::Steady,
                [820.0, 520.0],
                false,
            ),
            (
                "resize_desktop",
                DnaLatencyInteraction::Steady,
                [1600.0, 1000.0],
                false,
            ),
            (
                "resize_fullhd",
                DnaLatencyInteraction::Steady,
                [1920.0, 1080.0],
                false,
            ),
        ] {
            // Counter evidence and content checks are outside the timed sample.
            let (ctx, mut verification) = warm(&dna, derived_layers);
            let layers_before = verification.latency_benchmark_layer_inventory().unwrap();
            let before = verification.cache_diagnostics();
            verification.apply_latency_benchmark_interaction(interaction);
            frame(&ctx, &mut verification, size.into(), hover);
            let after = verification.cache_diagnostics();
            let layers_after = verification.latency_benchmark_layer_inventory().unwrap();
            assert_eq!(
                hash,
                content_hash(&verification.dna().read().unwrap(), derived_layers)
            );
            println!(
                "GENTLE_DNA_LATENCY_DIAGNOSTICS={}",
                json!({
                    "schema":"gentle.dna_feature_latency_observation.v1",
                    "revision":GENTLE_GIT_COMMIT, "fixture_sha256":hash,
                    "length_bp":length, "feature_count":count, "interaction":name,
                    "tree_loaded":true, "viewport_size":size, "before":before, "after":after,
                    "workload":DNA_LATENCY_WORKLOAD, "derived_layers":derived_layers,
                    "layers_before":layers_before, "layers_after":layers_after,
                })
            );
            group.bench_function(BenchmarkId::new(&id, name), |b| {
                b.iter_batched(
                    || warm(&dna, derived_layers),
                    |(ctx, mut area)| {
                        area.apply_latency_benchmark_interaction(interaction);
                        black_box(frame(&ctx, &mut area, size.into(), hover))
                    },
                    BatchSize::LargeInput,
                );
            });
        }
        group.finish();
    }
}

fn main() {
    if std::env::args().any(|arg| arg == "--gentle-identity") {
        println!(
            "{}",
            json!({"schema":"gentle.dna_feature_latency_binary.v1", "revision":GENTLE_GIT_COMMIT, "source_revision":GENTLE_SOURCE_REVISION, "workload":DNA_LATENCY_WORKLOAD})
        );
        return;
    }
    let mut criterion = Criterion::default().configure_from_args();
    benchmark(&mut criterion);
    criterion.final_summary();
}
