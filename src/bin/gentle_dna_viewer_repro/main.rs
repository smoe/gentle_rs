//! Standalone DNA-viewer repro harness.
//!
//! This binary deliberately starts only the DNA viewer, with selectable render
//! modes, so macOS embedded-window resize bugs can be bisected outside the full
//! GENtle project workspace.

use eframe::{NativeOptions, egui};
use gentle::dna_viewer_repro::{
    DnaViewerReproApp, load_repro_sequence, mode_help, parse_mode_list,
};
use std::path::PathBuf;

const DEFAULT_SEQUENCE_PATH: &str = "test_files/tp73.ncbi.gb";

#[derive(Debug)]
struct Args {
    sequence: Option<PathBuf>,
    synthetic_len_bp: Option<usize>,
    modes: Vec<String>,
    list_modes: bool,
}

impl Args {
    fn parse() -> Result<Self, String> {
        let mut args = std::env::args().skip(1);
        let mut parsed = Self {
            sequence: None,
            synthetic_len_bp: None,
            modes: Vec::new(),
            list_modes: false,
        };
        while let Some(arg) = args.next() {
            match arg.as_str() {
                "--help" | "-h" => return Err(help_text()),
                "--list-modes" => parsed.list_modes = true,
                "--mode" | "-m" => {
                    let value = args
                        .next()
                        .ok_or_else(|| "--mode requires a value".to_string())?;
                    parsed.modes.push(value);
                }
                "--sequence" | "--file" | "-s" => {
                    let value = args
                        .next()
                        .ok_or_else(|| "--sequence requires a path".to_string())?;
                    parsed.sequence = Some(PathBuf::from(value));
                }
                "--synthetic-len" => {
                    let value = args
                        .next()
                        .ok_or_else(|| "--synthetic-len requires a bp count".to_string())?;
                    let len = value
                        .parse::<usize>()
                        .map_err(|_| format!("invalid --synthetic-len value '{value}'"))?;
                    parsed.synthetic_len_bp = Some(len);
                }
                _ if arg.starts_with('-') => {
                    return Err(format!("unknown option '{arg}'\n\n{}", help_text()));
                }
                _ => {
                    if parsed.sequence.is_none() {
                        parsed.sequence = Some(PathBuf::from(arg));
                    } else {
                        return Err(format!("unexpected positional argument '{arg}'"));
                    }
                }
            }
        }
        Ok(parsed)
    }
}

fn main() -> eframe::Result<()> {
    let args = match Args::parse() {
        Ok(args) => args,
        Err(message) => {
            eprintln!("{message}");
            return Ok(());
        }
    };
    if args.list_modes {
        println!("{}", mode_help());
        return Ok(());
    }

    let modes = match parse_mode_list(&args.modes) {
        Ok(modes) => modes,
        Err(message) => {
            eprintln!("{message}");
            return Ok(());
        }
    };

    let default_path = PathBuf::from(DEFAULT_SEQUENCE_PATH);
    let sequence_path = if args.synthetic_len_bp.is_some() && args.sequence.is_none() {
        None
    } else {
        args.sequence
            .as_deref()
            .or_else(|| default_path.exists().then_some(default_path.as_path()))
    };
    let (dna, label) = match load_repro_sequence(sequence_path, args.synthetic_len_bp) {
        Ok(loaded) => loaded,
        Err(message) => {
            eprintln!("{message}");
            return Ok(());
        }
    };

    let options = NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1500.0, 980.0])
            .with_min_inner_size([640.0, 420.0]),
        ..Default::default()
    };
    eframe::run_native(
        "GENtle DNA viewer repro",
        options,
        Box::new(move |_cc| Ok(Box::new(DnaViewerReproApp::new(dna, label, modes)))),
    )
}

fn help_text() -> String {
    format!(
        "Usage: cargo run --bin gentle_dna_viewer_repro -- [OPTIONS] [SEQUENCE]\n\n\
         Options:\n\
           -s, --sequence PATH       Load a GenBank/FASTA/EMBL/SnapGene/XML sequence file\n\
           -m, --mode MODE           Enable a mode; repeat or comma-separate values\n\
               --synthetic-len BP    Use synthetic DNA if no sequence file is available\n\
               --list-modes          Print mode names and descriptions\n\
           -h, --help                Print this help\n\n\
         Default sequence: {DEFAULT_SEQUENCE_PATH} when present; otherwise synthetic DNA.\n\n\
         Modes:\n{}",
        mode_help()
    )
}
