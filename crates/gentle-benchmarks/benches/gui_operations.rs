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
    engine::GentleEngine,
    main_area_dna::MainAreaDna,
    window_dna::WindowDna,
};
use ring::digest::{SHA256, digest};
use std::{fs, hint::black_box, path::PathBuf, sync::Arc, sync::RwLock};

const SCREEN_WIDTH: f32 = 1_600.0;
const SCREEN_HEIGHT: f32 = 1_000.0;
const SYNTHETIC_CONTROL_BP: usize = 120_000;

struct GuiFixture {
    id: String,
    seq_id: String,
    dna: DNAsequence,
    source_sha256: String,
}

impl GuiFixture {
    fn benchmark_id(&self) -> String {
        let digest_token = self
            .source_sha256
            .strip_prefix("sha256:")
            .unwrap_or(&self.source_sha256);
        format!(
            "{}_{}bp_{}features_{}",
            self.id,
            self.dna.len(),
            self.dna.features().len(),
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

    vec![
        GuiFixture {
            id: "synthetic_control".to_string(),
            seq_id: "gui_benchmark_synthetic".to_string(),
            dna: synthetic,
            source_sha256: sha256_prefixed(synthetic_sequence.as_bytes()),
        },
        GuiFixture {
            id: "tp73_locus".to_string(),
            seq_id: "gui_benchmark_tp73".to_string(),
            dna: tp73,
            source_sha256: sha256_prefixed(&tp73_bytes),
        },
    ]
}

fn benchmark_engine(fixture: &GuiFixture) -> Arc<RwLock<GentleEngine>> {
    let mut engine = GentleEngine::new();
    engine
        .state_mut()
        .sequences
        .insert(fixture.seq_id.clone(), fixture.dna.clone());
    Arc::new(RwLock::new(engine))
}

fn benchmark_window(fixture: &GuiFixture) -> WindowDna {
    WindowDna::new(
        fixture.dna.clone(),
        fixture.seq_id.clone(),
        benchmark_engine(fixture),
    )
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

fn raw_input() -> egui::RawInput {
    egui::RawInput {
        screen_rect: Some(egui::Rect::from_min_size(
            egui::Pos2::ZERO,
            egui::vec2(SCREEN_WIDTH, SCREEN_HEIGHT),
        )),
        ..egui::RawInput::default()
    }
}

fn render_embedded_frame(context: &egui::Context, window: &mut WindowDna) -> usize {
    let mut output = context.run_ui(raw_input(), |ui| window.update_embedded(ui));
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
        let verification_shape_count =
            render_embedded_frame(&verification_context, &mut verification_window);
        assert_nonempty_frame(fixture, verification_shape_count);

        first_frame_group.bench_with_input(
            BenchmarkId::from_parameter(fixture.benchmark_id()),
            fixture,
            |b, fixture| {
                b.iter_batched(
                    || (egui::Context::default(), benchmark_window(fixture)),
                    |(context, mut window)| black_box(render_embedded_frame(&context, &mut window)),
                    BatchSize::SmallInput,
                );
            },
        );
    }
    first_frame_group.finish();

    let mut steady_frame_group = c.benchmark_group("dna_window_steady_embedded_frame");
    for fixture in &fixtures {
        let context = egui::Context::default();
        let mut window = benchmark_window(fixture);
        let first_shape_count = render_embedded_frame(&context, &mut window);
        assert_nonempty_frame(fixture, first_shape_count);
        steady_frame_group.bench_with_input(
            BenchmarkId::from_parameter(fixture.benchmark_id()),
            fixture,
            |b, _fixture| {
                b.iter(|| black_box(render_embedded_frame(&context, &mut window)));
            },
        );
    }
    steady_frame_group.finish();
}

criterion_group!(benches, benchmark_gui_operations);
criterion_main!(benches);
