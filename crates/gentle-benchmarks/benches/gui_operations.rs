//! Headless Criterion proxies for GUI-critical DNA-window work.
//!
//! These benchmarks execute the real `WindowDna` constructor and embedded egui
//! render path. They measure CPU-side preparation and painting, not native
//! window-system, GPU, input-dispatch, or compositor latency.

use criterion::{BatchSize, BenchmarkId, Criterion, criterion_group, criterion_main};
use egui;
use gentle::{
    about::GENTLE_SOURCE_REVISION,
    dna_sequence::{DNAsequence, load_from_file},
    engine::{GentleEngine, ProjectState},
    main_area_dna::MainAreaDna,
    window_dna::WindowDna,
};
use ring::digest::{SHA256, digest};
use serde_json::Value;
use std::{env, fs, hint::black_box, path::PathBuf, sync::Arc, sync::RwLock};

const SYNTHETIC_CONTROL_BP: usize = 120_000;
const PATZ1_STATE_ENV: &str = "GENTLE_GUI_BENCH_PATZ1_STATE";
const PATZ1_REPORT_ENV: &str = "GENTLE_GUI_BENCH_PATZ1_REPORT";

#[derive(Clone, Copy)]
struct ScreenSize {
    id: &'static str,
    width: f32,
    height: f32,
}

const SCREEN_SIZES: [ScreenSize; 4] = [
    ScreenSize {
        id: "compact_820x520",
        width: 820.0,
        height: 520.0,
    },
    ScreenSize {
        id: "laptop_1200x800",
        width: 1_200.0,
        height: 800.0,
    },
    ScreenSize {
        id: "desktop_1600x1000",
        width: 1_600.0,
        height: 1_000.0,
    },
    ScreenSize {
        id: "fullhd_1920x1080",
        width: 1_920.0,
        height: 1_080.0,
    },
];

struct GuiFixture {
    id: String,
    seq_id: String,
    dna: DNAsequence,
    state: ProjectState,
    source_sha256: String,
    transcript_count: usize,
    comparison_record_count: usize,
    locus_report_path: Option<PathBuf>,
}

impl GuiFixture {
    fn benchmark_id(&self) -> String {
        let digest_token = self
            .source_sha256
            .strip_prefix("sha256:")
            .unwrap_or(&self.source_sha256);
        format!(
            "{}_{}bp_{}features_{}transcripts_{}source-records_{}",
            self.id,
            self.dna.len(),
            self.dna.features().len(),
            self.transcript_count,
            self.comparison_record_count,
            &digest_token[..12.min(digest_token.len())]
        )
    }
}

fn sha256_prefixed(bytes: &[u8]) -> String {
    let hex = digest(&SHA256, bytes)
        .as_ref()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<String>();
    format!("sha256:{hex}")
}

fn project_fixture(id: &str, seq_id: &str, dna: DNAsequence, source_sha256: String) -> GuiFixture {
    let mut state = ProjectState::default();
    state.sequences.insert(seq_id.to_string(), dna.clone());
    GuiFixture {
        id: id.to_string(),
        seq_id: seq_id.to_string(),
        dna,
        state,
        source_sha256,
        transcript_count: 0,
        comparison_record_count: 0,
        locus_report_path: None,
    }
}

fn count_patz1_transcripts(state_json: &Value) -> usize {
    state_json
        .pointer("/metadata/isoform_panels/records/0/resource/evidence")
        .and_then(Value::as_array)
        .map_or(0, Vec::len)
}

fn load_patz1_fixture(state_path: PathBuf, report_path: PathBuf) -> GuiFixture {
    let state_bytes = fs::read(&state_path).expect("read prepared public PATZ1 project");
    let state_json: Value =
        serde_json::from_slice(&state_bytes).expect("parse prepared public PATZ1 project JSON");
    let state = ProjectState::load_from_path(&state_path.to_string_lossy())
        .expect("load prepared public PATZ1 project");
    let seq_id = "patz1_ensembl_116";
    let dna = state
        .sequences
        .get(seq_id)
        .cloned()
        .expect("prepared public PATZ1 project contains expected sequence");
    let report_json: Value = serde_json::from_slice(
        &fs::read(&report_path).expect("read prepared public PATZ1 locus report"),
    )
    .expect("parse prepared public PATZ1 locus report JSON");
    let transcript_count = count_patz1_transcripts(&state_json);
    let comparison_record_count = report_json
        .pointer("/transcript_presentation/records")
        .and_then(Value::as_array)
        .map_or(0, Vec::len);
    assert_eq!(dna.len(), 20_802, "public PATZ1 locus length drifted");
    assert_eq!(
        transcript_count, 13,
        "public PATZ1 Ensembl transcript count drifted"
    );
    assert_eq!(
        comparison_record_count, 17,
        "public PATZ1 source-record count drifted"
    );
    assert!(
        dna.features().len() >= 75,
        "public PATZ1 feature set is incomplete"
    );
    let fixture_manifest_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../test_files/fixtures/transcript_assay_panel/patz1_reference/manifest.json");
    let fixture_manifest_bytes =
        fs::read(fixture_manifest_path).expect("read pinned public PATZ1 fixture manifest");
    GuiFixture {
        id: "patz1_full_annotation".to_string(),
        seq_id: seq_id.to_string(),
        dna,
        state,
        // The generated project contains an import timestamp. Benchmark IDs
        // must remain stable across equivalent preparations, so bind the
        // statistical series to the pinned source manifest instead. The audit
        // runner records the exact generated-project hash separately.
        source_sha256: sha256_prefixed(&fixture_manifest_bytes),
        transcript_count,
        comparison_record_count,
        locus_report_path: Some(report_path),
    }
}

fn benchmark_fixtures() -> Vec<GuiFixture> {
    let pattern = b"ACGTGCAATTCG";
    let mut synthetic_sequence = String::with_capacity(SYNTHETIC_CONTROL_BP);
    for index in 0..SYNTHETIC_CONTROL_BP {
        synthetic_sequence.push(pattern[index % pattern.len()] as char);
    }
    let synthetic = DNAsequence::from_sequence(&synthetic_sequence)
        .expect("construct feature-free GUI benchmark control");

    let tp73_path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../test_files/tp73.ncbi.gb");
    let tp73_bytes = fs::read(&tp73_path).expect("read public TP73 GUI benchmark fixture");
    let tp73 = load_from_file(&tp73_path.to_string_lossy())
        .expect("parse public TP73 GUI benchmark fixture");

    let mut fixtures = vec![
        project_fixture(
            "synthetic_control",
            "gui_benchmark_synthetic",
            synthetic,
            sha256_prefixed(synthetic_sequence.as_bytes()),
        ),
        project_fixture(
            "tp73_locus",
            "gui_benchmark_tp73",
            tp73,
            sha256_prefixed(&tp73_bytes),
        ),
    ];
    match (env::var_os(PATZ1_STATE_ENV), env::var_os(PATZ1_REPORT_ENV)) {
        (Some(state), Some(report)) => {
            fixtures.push(load_patz1_fixture(state.into(), report.into()));
        }
        (None, None) => eprintln!(
            "PATZ1 fixture omitted; set {PATZ1_STATE_ENV} and {PATZ1_REPORT_ENV} for the full public-data benchmark"
        ),
        _ => panic!("{PATZ1_STATE_ENV} and {PATZ1_REPORT_ENV} must be supplied together"),
    }
    fixtures
}

fn benchmark_engine(fixture: &GuiFixture) -> Arc<RwLock<GentleEngine>> {
    Arc::new(RwLock::new(GentleEngine::from_state(fixture.state.clone())))
}

fn benchmark_window(fixture: &GuiFixture) -> WindowDna {
    let mut window = WindowDna::new(
        fixture.dna.clone(),
        fixture.seq_id.clone(),
        benchmark_engine(fixture),
    );
    if let Some(report_path) = &fixture.locus_report_path {
        window
            .load_locus_report_for_benchmark(report_path)
            .expect("bind prepared public PATZ1 locus report to benchmark window");
    }
    window
}

fn deferred_main_area(fixture: &GuiFixture) -> MainAreaDna {
    let placeholder = DNAsequence::from_sequence("").expect("construct empty DNA placeholder");
    let mut main_area = MainAreaDna::new(
        placeholder,
        Some(fixture.seq_id.clone()),
        Some(benchmark_engine(fixture)),
    );
    main_area.defer_feature_tree_until_interaction();
    main_area
}

fn raw_input(size: ScreenSize) -> egui::RawInput {
    egui::RawInput {
        screen_rect: Some(egui::Rect::from_min_size(
            egui::Pos2::ZERO,
            egui::vec2(size.width, size.height),
        )),
        ..egui::RawInput::default()
    }
}

fn render_embedded_frame(
    context: &egui::Context,
    window: &mut WindowDna,
    size: ScreenSize,
) -> usize {
    let mut output = context.run_ui(raw_input(size), |ui| window.update_embedded(ui));
    let shape_count = output.shapes.len();
    // This headless benchmark has no renderer to apply texture uploads.
    output.textures_delta.clear();
    shape_count
}

fn assert_nonempty_frame(fixture: &GuiFixture, shape_count: usize) {
    assert!(
        shape_count > 0,
        "{} GUI fixture produced no paint shapes",
        fixture.id
    );
}

fn benchmark_gui_operations(c: &mut Criterion) {
    let fixtures = benchmark_fixtures();
    eprintln!("GENtle source revision: {GENTLE_SOURCE_REVISION}");

    let mut constructor_group = c.benchmark_group("dna_window_constructor");
    for fixture in &fixtures {
        constructor_group.bench_with_input(
            BenchmarkId::from_parameter(fixture.benchmark_id()),
            fixture,
            |b, fixture| {
                b.iter_batched(
                    || (fixture.dna.clone(), benchmark_engine(fixture)),
                    |(dna, engine)| black_box(WindowDna::new(dna, fixture.seq_id.clone(), engine)),
                    BatchSize::SmallInput,
                );
            },
        );
    }
    constructor_group.finish();

    let mut hydration_group = c.benchmark_group("dna_window_deferred_hydration");
    for fixture in &fixtures {
        hydration_group.bench_with_input(
            BenchmarkId::from_parameter(fixture.benchmark_id()),
            fixture,
            |b, fixture| {
                b.iter_batched(
                    || (deferred_main_area(fixture), fixture.dna.clone()),
                    |(mut main_area, dna)| {
                        main_area.replace_loaded_sequence(dna);
                        black_box(main_area)
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }
    hydration_group.finish();

    let mut first_frame_group = c.benchmark_group("dna_window_first_embedded_frame");
    for fixture in &fixtures {
        let mut verification_window = benchmark_window(fixture);
        let verification_context = egui::Context::default();
        for size in SCREEN_SIZES {
            let verification_shape_count =
                render_embedded_frame(&verification_context, &mut verification_window, size);
            assert_nonempty_frame(fixture, verification_shape_count);
            first_frame_group.bench_with_input(
                BenchmarkId::new(fixture.benchmark_id(), size.id),
                &(fixture, size),
                |b, (fixture, size)| {
                    b.iter_batched(
                        || (egui::Context::default(), benchmark_window(fixture)),
                        |(context, mut window)| {
                            black_box(render_embedded_frame(&context, &mut window, *size))
                        },
                        BatchSize::SmallInput,
                    );
                },
            );
        }
    }
    first_frame_group.finish();

    let mut steady_frame_group = c.benchmark_group("dna_window_steady_embedded_frame");
    for fixture in &fixtures {
        let context = egui::Context::default();
        let mut window = benchmark_window(fixture);
        for size in SCREEN_SIZES {
            let first_shape_count = render_embedded_frame(&context, &mut window, size);
            assert_nonempty_frame(fixture, first_shape_count);
            steady_frame_group.bench_with_input(
                BenchmarkId::new(fixture.benchmark_id(), size.id),
                &(fixture, size),
                |b, (_fixture, size)| {
                    b.iter(|| black_box(render_embedded_frame(&context, &mut window, *size)));
                },
            );
        }
    }
    steady_frame_group.finish();

    let mut resize_group = c.benchmark_group("dna_window_first_frame_after_resize");
    let source_size = SCREEN_SIZES[1];
    for fixture in &fixtures {
        for target_size in [SCREEN_SIZES[0], SCREEN_SIZES[2], SCREEN_SIZES[3]] {
            resize_group.bench_with_input(
                BenchmarkId::new(
                    fixture.benchmark_id(),
                    format!("{}_to_{}", source_size.id, target_size.id),
                ),
                &(fixture, target_size),
                |b, (fixture, target_size)| {
                    b.iter_batched(
                        || {
                            let context = egui::Context::default();
                            let mut window = benchmark_window(fixture);
                            let shape_count =
                                render_embedded_frame(&context, &mut window, source_size);
                            assert_nonempty_frame(fixture, shape_count);
                            (context, window)
                        },
                        |(context, mut window)| {
                            black_box(render_embedded_frame(&context, &mut window, *target_size))
                        },
                        BatchSize::SmallInput,
                    );
                },
            );
        }
    }
    resize_group.finish();
}

criterion_group!(benches, benchmark_gui_operations);
criterion_main!(benches);
