use std::{env, path::Path};

use gentle::gene_set_publication::generate_gene_set_publication;

fn main() {
    let args = env::args().skip(1).collect::<Vec<_>>();
    if args.len() < 2 || args.len() > 3 || (args.len() == 3 && args[2] != "--pdf") {
        eprintln!("Usage: gentle_publication_report REQUEST.json OUTPUT_DIR [--pdf]");
        std::process::exit(2);
    }
    match generate_gene_set_publication(Path::new(&args[0]), Path::new(&args[1]), args.len() == 3) {
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
