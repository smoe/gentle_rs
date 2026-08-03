use std::{env, path::Path};

use gentle::gene_set_publication::{
    generate_gene_isoform_assay_publication, generate_gene_set_publication,
};
use gentle_protocol::GENE_ISOFORM_ASSAY_PUBLICATION_REQUEST_SCHEMA;

fn main() {
    let args = env::args().skip(1).collect::<Vec<_>>();
    if args.len() < 2 {
        eprintln!(
            "Usage: gentle_publication_report REQUEST.json OUTPUT_DIR [--profile ID] [--blocks ID,ID] [--pdf]"
        );
        std::process::exit(2);
    }
    let mut generate_pdf = false;
    let mut profile = None;
    let mut blocks = vec![];
    let mut idx = 2usize;
    while idx < args.len() {
        match args[idx].as_str() {
            "--pdf" => {
                generate_pdf = true;
                idx += 1;
            }
            "--profile" if idx + 1 < args.len() => {
                profile = Some(args[idx + 1].clone());
                idx += 2;
            }
            "--blocks" if idx + 1 < args.len() => {
                blocks = args[idx + 1]
                    .split(',')
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(ToString::to_string)
                    .collect();
                idx += 2;
            }
            other => {
                eprintln!("Unknown or incomplete option '{other}'");
                std::process::exit(2);
            }
        }
    }
    let request_path = Path::new(&args[0]);
    let schema = std::fs::read_to_string(request_path)
        .ok()
        .and_then(|text| serde_json::from_str::<serde_json::Value>(&text).ok())
        .and_then(|value| value.get("schema")?.as_str().map(ToString::to_string));
    let result = if schema.as_deref() == Some(GENE_ISOFORM_ASSAY_PUBLICATION_REQUEST_SCHEMA) {
        generate_gene_isoform_assay_publication(
            request_path,
            Path::new(&args[1]),
            profile.as_deref(),
            &blocks,
            generate_pdf,
        )
        .and_then(|report| serde_json::to_value(report).map_err(|error| error.to_string()))
    } else if profile.is_some() || !blocks.is_empty() {
        Err(
            "--profile and --blocks require gentle.gene_isoform_assay_publication_request.v1"
                .to_string(),
        )
    } else {
        generate_gene_set_publication(request_path, Path::new(&args[1]), generate_pdf)
            .and_then(|report| serde_json::to_value(report).map_err(|error| error.to_string()))
    };
    match result {
        Ok(report) => println!(
            "{}",
            serde_json::to_string_pretty(&report).expect("serialize generation report")
        ),
        Err(error) => {
            eprintln!("gentle_publication_report: {error}");
            std::process::exit(1);
        }
    }
}
