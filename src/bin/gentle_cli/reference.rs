//! Direct `genomes`/`helpers` CLI command handling kept out of `run()`.

use super::*;

pub(super) fn handle_reference_family(
    args: &[String],
    cmd_idx: usize,
    command: &str,
    state_path: &str,
    global: &GlobalCliArgs,
) -> Result<(), String> {
    let helper_mode = command == "helpers";
    let label = if helper_mode { "helpers" } else { "genomes" };
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err(format!("{label} requires a subcommand"));
    }
    match args[cmd_idx + 1].as_str() {
        "vocabulary" if helper_mode => {
            if args.len() <= cmd_idx + 2 {
                return Err(
                    "helpers vocabulary requires a subcommand (expected list or doctor)"
                        .to_string(),
                );
            }
            match args[cmd_idx + 2].as_str() {
                "list" => {
                    let mut vocabulary_path: Option<String> = None;
                    let mut filter: Option<String> = None;
                    let mut idx = cmd_idx + 3;
                    while idx < args.len() {
                        match args[idx].as_str() {
                            "--vocabulary" => {
                                vocabulary_path = Some(gentle_cli_args::required_value(
                                    args,
                                    idx,
                                    "--vocabulary",
                                    "PATH",
                                    "helpers vocabulary list",
                                )?);
                                idx += 2;
                            }
                            "--filter" => {
                                filter = Some(gentle_cli_args::required_value(
                                    args,
                                    idx,
                                    "--filter",
                                    "TEXT",
                                    "helpers vocabulary list",
                                )?);
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for helpers vocabulary list",
                                    other
                                ));
                            }
                        }
                    }
                    let terms = GentleEngine::list_helper_semantics_vocabulary_terms(
                        vocabulary_path.as_deref(),
                        filter.as_deref(),
                    )
                    .map_err(|e| e.to_string())?;
                    print_json(&json!({
                        "vocabulary_path": vocabulary_path
                            .clone()
                            .unwrap_or_else(|| default_helper_semantics_vocabulary_discovery_label().to_string()),
                        "filter": filter,
                        "term_count": terms.len(),
                        "terms": terms,
                    }))
                }
                "doctor" | "validate" => {
                    let mut vocabulary_path: Option<String> = None;
                    let mut routine_catalog_path: Option<String> = None;
                    let mut idx = cmd_idx + 3;
                    while idx < args.len() {
                        match args[idx].as_str() {
                            "--vocabulary" => {
                                vocabulary_path = Some(gentle_cli_args::required_value(
                                    args,
                                    idx,
                                    "--vocabulary",
                                    "PATH",
                                    "helpers vocabulary doctor",
                                )?);
                                idx += 2;
                            }
                            "--routine-catalog" => {
                                routine_catalog_path = Some(gentle_cli_args::required_value(
                                    args,
                                    idx,
                                    "--routine-catalog",
                                    "PATH",
                                    "helpers vocabulary doctor",
                                )?);
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for helpers vocabulary doctor",
                                    other
                                ));
                            }
                        }
                    }
                    let known_routine_families =
                        cloning_routine_families_from_catalog(routine_catalog_path.as_deref())?;
                    let report = GentleEngine::doctor_helper_semantics_vocabulary(
                        vocabulary_path.as_deref(),
                        &known_routine_families,
                    )
                    .map_err(|e| e.to_string())?;
                    print_json(&json!({
                        "vocabulary_path": vocabulary_path
                            .clone()
                            .unwrap_or_else(|| default_helper_semantics_vocabulary_discovery_label().to_string()),
                        "routine_catalog_path": routine_catalog_path
                            .clone()
                            .unwrap_or_else(|| DEFAULT_CLONING_ROUTINE_CATALOG_PATH.to_string()),
                        "report": report,
                    }))
                }
                other => Err(format!(
                    "Unknown helpers vocabulary subcommand '{}' (expected list or doctor)",
                    other
                )),
            }
        }
        "ensembl-available" => {
            let mut collection: Option<String> = None;
            let mut filter: Option<String> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--collection" => {
                        collection = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--collection",
                            "VALUE",
                            &format!("{label} ensembl-available"),
                        )?);
                        idx += 2;
                    }
                    "--filter" => {
                        filter = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--filter",
                            "TEXT",
                            &format!("{label} ensembl-available"),
                        )?);
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for {label} ensembl-available",
                            other
                        ));
                    }
                }
            }
            let report = GentleEngine::discover_ensembl_installable_genomes(
                collection.as_deref(),
                filter.as_deref(),
            )
            .map_err(|e| e.to_string())?;
            print_json(&json!({
                "scope": label,
                "report": report,
            }))
        }
        "list" => {
            let mut catalog_path: Option<String> = None;
            let mut filter: Option<String> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--catalog",
                            "PATH",
                            &format!("{label} list"),
                        )?);
                        idx += 2;
                    }
                    "--filter" => {
                        filter = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--filter",
                            "TEXT",
                            &format!("{label} list"),
                        )?);
                        idx += 2;
                    }
                    other => {
                        return Err(format!("Unknown option '{}' for {label} list", other));
                    }
                }
            }
            let resolved_catalog = explicit_catalog_arg(&catalog_path);
            let entries = if helper_mode {
                GentleEngine::list_helper_catalog_entries(resolved_catalog, filter.as_deref())
            } else {
                GentleEngine::list_reference_catalog_entries(resolved_catalog, filter.as_deref())
            }
            .map_err(|e| e.to_string())?;
            let genomes = entries
                .iter()
                .map(|entry| entry.genome_id.clone())
                .collect::<Vec<_>>();
            let effective_catalog = effective_catalog_label(&catalog_path, helper_mode);
            print_json(&json!({
                "catalog_path": effective_catalog,
                "filter": filter,
                "genome_count": genomes.len(),
                "genomes": genomes,
                "entries": entries,
            }))
        }
        "validate-catalog" => {
            let mut catalog_path: Option<String> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--catalog",
                            "PATH",
                            &format!("{label} validate-catalog"),
                        )?);
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for {label} validate-catalog",
                            other
                        ));
                    }
                }
            }
            let resolved_catalog = explicit_catalog_arg(&catalog_path);
            let genomes = if helper_mode {
                GentleEngine::list_helper_genomes(resolved_catalog)
            } else {
                GentleEngine::list_reference_genomes(resolved_catalog)
            }
            .map_err(|e| e.to_string())?;
            print_json(&json!({
                "catalog_path": effective_catalog_label(&catalog_path, helper_mode),
                "valid": true,
                "genome_count": genomes.len(),
            }))
        }
        "doctor-catalog" if helper_mode => {
            let mut catalog_path: Option<String> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--catalog",
                            "PATH",
                            &format!("{label} doctor-catalog"),
                        )?);
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for {label} doctor-catalog",
                            other
                        ));
                    }
                }
            }
            let resolved_catalog = explicit_catalog_arg(&catalog_path);
            let report = GentleEngine::doctor_helper_vector_catalog(resolved_catalog)
                .map_err(|e| e.to_string())?;
            print_json(&json!({
                "catalog_path": effective_catalog_label(&catalog_path, helper_mode),
                "report": report,
            }))
        }
        "show-card" if helper_mode => {
            let mut catalog_path: Option<String> = None;
            let mut filter: Option<String> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--catalog",
                            "PATH",
                            &format!("{label} show-card"),
                        )?);
                        idx += 2;
                    }
                    "--filter" | "--name" => {
                        let option = args[idx].clone();
                        filter = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            &option,
                            "TEXT",
                            &format!("{label} show-card"),
                        )?);
                        idx += 2;
                    }
                    other => {
                        return Err(format!("Unknown option '{}' for {label} show-card", other));
                    }
                }
            }
            let resolved_catalog = explicit_catalog_arg(&catalog_path);
            let report =
                GentleEngine::list_helper_vector_cards(resolved_catalog, filter.as_deref())
                    .map_err(|e| e.to_string())?;
            print_json(&json!({
                "catalog_path": effective_catalog_label(&catalog_path, helper_mode),
                "filter": filter,
                "report": report,
            }))
        }
        "update-ensembl-specs" => {
            let mut catalog_path: Option<String> = None;
            let mut output_catalog_path: Option<String> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--catalog",
                            "PATH",
                            &format!("{label} update-ensembl-specs"),
                        )?);
                        idx += 2;
                    }
                    "--output-catalog" => {
                        output_catalog_path = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--output-catalog",
                            "PATH",
                            &format!("{label} update-ensembl-specs"),
                        )?);
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for {label} update-ensembl-specs",
                            other
                        ));
                    }
                }
            }
            let resolved_catalog = explicit_catalog_arg(&catalog_path);
            let report = if helper_mode {
                GentleEngine::apply_helper_genome_ensembl_catalog_updates(
                    resolved_catalog,
                    output_catalog_path.as_deref(),
                )
            } else {
                GentleEngine::apply_reference_genome_ensembl_catalog_updates(
                    resolved_catalog,
                    output_catalog_path.as_deref(),
                )
            }
            .map_err(|e| e.to_string())?;
            print_json(&report)
        }
        "status" => {
            if args.len() <= cmd_idx + 2 {
                usage();
                return Err(format!(
                    "{label} status requires GENOME_ID [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = args[cmd_idx + 2].clone();
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = cmd_idx + 3;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--catalog",
                            "PATH",
                            &format!("{label} status"),
                        )?);
                        idx += 2;
                    }
                    "--cache-dir" => {
                        cache_dir = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--cache-dir",
                            "PATH",
                            &format!("{label} status"),
                        )?);
                        idx += 2;
                    }
                    other => {
                        return Err(format!("Unknown option '{}' for {label} status", other));
                    }
                }
            }
            let resolved_catalog = explicit_catalog_arg(&catalog_path);
            let prepared = if helper_mode {
                GentleEngine::is_helper_genome_prepared(
                    &genome_id,
                    resolved_catalog,
                    cache_dir.as_deref(),
                )
            } else {
                GentleEngine::is_reference_genome_prepared(
                    resolved_catalog,
                    &genome_id,
                    cache_dir.as_deref(),
                )
            }
            .map_err(|e| e.to_string())?;
            let source_plan = if helper_mode {
                GentleEngine::describe_helper_genome_sources(
                    &genome_id,
                    resolved_catalog,
                    cache_dir.as_deref(),
                )
            } else {
                GentleEngine::describe_reference_genome_sources(
                    resolved_catalog,
                    &genome_id,
                    cache_dir.as_deref(),
                )
            }
            .map_err(|e| e.to_string())?;
            let interpretation = if helper_mode {
                GentleEngine::interpret_helper_genome(&genome_id, resolved_catalog)
                    .map_err(|e| e.to_string())?
            } else {
                None
            };
            let effective_catalog = effective_catalog_label(&catalog_path, helper_mode);
            print_json(&json!({
                "genome_id": genome_id,
                "catalog_path": effective_catalog,
                "cache_dir": cache_dir,
                "prepared": prepared,
                "sequence_source_type": source_plan.sequence_source_type,
                "annotation_source_type": source_plan.annotation_source_type,
                "sequence_source": source_plan.sequence_source,
                "annotation_source": source_plan.annotation_source,
                "nucleotide_length_bp": source_plan.nucleotide_length_bp,
                "molecular_mass_da": source_plan.molecular_mass_da,
                "molecular_mass_source": source_plan.molecular_mass_source,
                "interpretation": interpretation,
            }))
        }
        "genes" => {
            if args.len() <= cmd_idx + 2 {
                usage();
                return Err(format!(
                    "{label} genes requires GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]"
                ));
            }
            let genome_id = args[cmd_idx + 2].clone();
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut filter = String::new();
            let mut biotype_filters: Vec<String> = vec![];
            let mut limit: usize = 200;
            let mut offset: usize = 0;
            let mut idx = cmd_idx + 3;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--catalog",
                            "PATH",
                            &format!("{label} genes"),
                        )?);
                        idx += 2;
                    }
                    "--cache-dir" => {
                        cache_dir = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--cache-dir",
                            "PATH",
                            &format!("{label} genes"),
                        )?);
                        idx += 2;
                    }
                    "--filter" => {
                        filter = gentle_cli_args::required_value(
                            args,
                            idx,
                            "--filter",
                            "TEXT",
                            &format!("{label} genes"),
                        )?;
                        idx += 2;
                    }
                    "--biotype" => {
                        biotype_filters.push(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--biotype",
                            "NAME",
                            &format!("{label} genes"),
                        )?);
                        idx += 2;
                    }
                    "--limit" => {
                        limit = gentle_cli_args::parse_required_value(
                            args,
                            idx,
                            "--limit",
                            "N",
                            &format!("{label} genes"),
                            None,
                        )?;
                        if limit == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                        idx += 2;
                    }
                    "--offset" => {
                        offset = gentle_cli_args::parse_required_value(
                            args,
                            idx,
                            "--offset",
                            "N",
                            &format!("{label} genes"),
                            None,
                        )?;
                        idx += 2;
                    }
                    other => {
                        return Err(format!("Unknown option '{}' for {label} genes", other));
                    }
                }
            }

            let resolved_catalog = explicit_catalog_arg(&catalog_path);
            let genes = if helper_mode {
                GentleEngine::list_helper_genome_features(
                    &genome_id,
                    resolved_catalog,
                    cache_dir.as_deref(),
                )
            } else {
                GentleEngine::list_reference_genome_genes(
                    resolved_catalog,
                    &genome_id,
                    cache_dir.as_deref(),
                )
            }
            .map_err(|e| e.to_string())?;
            let filter_regex = compile_gene_filter_regex(&filter)?;
            let available_biotypes = collect_biotypes(&genes);
            let allowed_biotypes_lower: Vec<String> = biotype_filters
                .iter()
                .map(|v| v.trim().to_ascii_lowercase())
                .filter(|v| !v.is_empty())
                .collect();
            let filtered: Vec<GenomeGeneRecord> = genes
                .into_iter()
                .filter(|g| {
                    genome_gene_matches_filter(g, filter_regex.as_ref(), &allowed_biotypes_lower)
                })
                .collect();
            let total = filtered.len();
            let offset = offset.min(total);
            let returned: Vec<GenomeGeneRecord> =
                filtered.into_iter().skip(offset).take(limit).collect();
            let effective_catalog = effective_catalog_label(&catalog_path, helper_mode);
            print_json(&json!({
                "genome_id": genome_id,
                "catalog_path": effective_catalog,
                "cache_dir": cache_dir,
                "filter": filter,
                "biotype_filter": biotype_filters,
                "available_biotypes": available_biotypes,
                "offset": offset,
                "limit": limit,
                "total_matches": total,
                "returned": returned.len(),
                "genes": returned,
            }))
        }
        "prepare" => {
            if args.len() <= cmd_idx + 2 {
                usage();
                return Err(format!(
                    "{label} prepare requires GENOME_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]"
                ));
            }
            let genome_id = args[cmd_idx + 2].clone();
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut timeout_seconds: Option<u64> = None;
            let mut idx = cmd_idx + 3;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--catalog",
                            "PATH",
                            &format!("{label} prepare"),
                        )?);
                        idx += 2;
                    }
                    "--cache-dir" => {
                        cache_dir = Some(gentle_cli_args::required_value(
                            args,
                            idx,
                            "--cache-dir",
                            "PATH",
                            &format!("{label} prepare"),
                        )?);
                        idx += 2;
                    }
                    "--timeout-secs" => {
                        let parsed = gentle_cli_args::parse_required_value(
                            args,
                            idx,
                            "--timeout-secs",
                            "N",
                            &format!("{label} prepare"),
                            Some(&format!("{label} prepare")),
                        )?;
                        if parsed == 0 {
                            timeout_seconds = None;
                        } else {
                            timeout_seconds = Some(parsed);
                        }
                        idx += 2;
                    }
                    other => {
                        return Err(format!("Unknown option '{}' for {label} prepare", other));
                    }
                }
            }
            let op_catalog_path = operation_catalog_arg(&catalog_path, helper_mode);
            let mut engine = GentleEngine::from_state(load_state(state_path)?);
            let op = Operation::PrepareGenome {
                genome_id,
                catalog_path: op_catalog_path,
                cache_dir,
                timeout_seconds,
            };
            let result = if let Some(sink) = global.progress_sink {
                let mut printer = ProgressPrinter::new(sink);
                engine
                    .apply_with_progress(op, |p| {
                        printer.on_progress(p);
                        true
                    })
                    .map_err(|e| e.to_string())?
            } else {
                engine.apply(op).map_err(|e| e.to_string())?
            };
            engine
                .state()
                .save_to_path(state_path)
                .map_err(|e| e.to_string())?;
            print_json(&result)
        }
        "remove-prepared" => {
            if args.len() <= cmd_idx + 2 {
                usage();
                return Err(format!(
                    "{label} remove-prepared requires GENOME_ID [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = args[cmd_idx + 2].clone();
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = cmd_idx + 3;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--catalog" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing PATH after --catalog for {label} remove-prepared"
                            ));
                        }
                        catalog_path = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--cache-dir" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing PATH after --cache-dir for {label} remove-prepared"
                            ));
                        }
                        cache_dir = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for {label} remove-prepared",
                            other
                        ));
                    }
                }
            }
            let resolved_catalog = explicit_catalog_arg(&catalog_path);
            let report = if helper_mode {
                GentleEngine::remove_prepared_helper_genome(
                    resolved_catalog,
                    &genome_id,
                    cache_dir.as_deref(),
                )
            } else {
                GentleEngine::remove_prepared_reference_genome(
                    resolved_catalog,
                    &genome_id,
                    cache_dir.as_deref(),
                )
            }
            .map_err(|e| e.to_string())?;
            print_json(&report)
        }
        "remove-catalog-entry" => {
            if args.len() <= cmd_idx + 2 {
                usage();
                return Err(format!(
                    "{label} remove-catalog-entry requires GENOME_ID [--catalog PATH] [--output-catalog PATH]"
                ));
            }
            let genome_id = args[cmd_idx + 2].clone();
            let mut catalog_path: Option<String> = None;
            let mut output_catalog_path: Option<String> = None;
            let mut idx = cmd_idx + 3;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--catalog" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing PATH after --catalog for {label} remove-catalog-entry"
                            ));
                        }
                        catalog_path = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--output-catalog" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing PATH after --output-catalog for {label} remove-catalog-entry"
                            ));
                        }
                        output_catalog_path = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for {label} remove-catalog-entry",
                            other
                        ));
                    }
                }
            }
            let resolved_catalog = explicit_catalog_arg(&catalog_path);
            let report = if helper_mode {
                GentleEngine::remove_helper_genome_catalog_entry(
                    resolved_catalog,
                    &genome_id,
                    output_catalog_path.as_deref(),
                )
            } else {
                GentleEngine::remove_reference_genome_catalog_entry(
                    resolved_catalog,
                    &genome_id,
                    output_catalog_path.as_deref(),
                )
            }
            .map_err(|e| e.to_string())?;
            print_json(&report)
        }
        "extract-region" => {
            if args.len() <= cmd_idx + 6 {
                usage();
                return Err(format!(
                    "{label} extract-region requires GENOME_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = args[cmd_idx + 2].clone();
            let chromosome = args[cmd_idx + 3].clone();
            let start_1based = args[cmd_idx + 4]
                .parse::<usize>()
                .map_err(|e| format!("Invalid START coordinate '{}': {}", args[cmd_idx + 4], e))?;
            let end_1based = args[cmd_idx + 5]
                .parse::<usize>()
                .map_err(|e| format!("Invalid END coordinate '{}': {}", args[cmd_idx + 5], e))?;
            let mut output_id: Option<String> = None;
            let mut annotation_scope: Option<GenomeAnnotationScope> = None;
            let mut max_annotation_features: Option<usize> = None;
            let mut include_genomic_annotation: Option<bool> = None;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = cmd_idx + 6;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--output-id" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing ID after --output-id for {label} extract-region"
                            ));
                        }
                        output_id = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--catalog" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing PATH after --catalog for {label} extract-region"
                            ));
                        }
                        catalog_path = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--cache-dir" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing PATH after --cache-dir for {label} extract-region"
                            ));
                        }
                        cache_dir = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--annotation-scope" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing VALUE after --annotation-scope for {label} extract-region"
                            ));
                        }
                        let value = args[idx + 1].trim().to_ascii_lowercase();
                        let parsed = match value.as_str() {
                            "none" => GenomeAnnotationScope::None,
                            "core" => GenomeAnnotationScope::Core,
                            "full" => GenomeAnnotationScope::Full,
                            other => {
                                return Err(format!(
                                    "Invalid --annotation-scope value '{}' for {label} extract-region (expected none|core|full)",
                                    other
                                ));
                            }
                        };
                        annotation_scope = Some(parsed);
                        idx += 2;
                    }
                    "--max-annotation-features" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing N after --max-annotation-features for {label} extract-region"
                            ));
                        }
                        let raw = args[idx + 1].trim().to_string();
                        let parsed = raw.parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --max-annotation-features value '{}' for {label} extract-region: {}",
                                raw, e
                            )
                        })?;
                        max_annotation_features = Some(parsed);
                        idx += 2;
                    }
                    "--include-genomic-annotation" => {
                        include_genomic_annotation = Some(true);
                        idx += 1;
                    }
                    "--no-include-genomic-annotation" => {
                        include_genomic_annotation = Some(false);
                        idx += 1;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for {label} extract-region",
                            other
                        ));
                    }
                }
            }
            if let Some(include) = include_genomic_annotation {
                let mapped_scope = if include {
                    GenomeAnnotationScope::Core
                } else {
                    GenomeAnnotationScope::None
                };
                if let Some(explicit_scope) = annotation_scope {
                    if explicit_scope != mapped_scope {
                        return Err(format!(
                            "Conflicting annotation options for {label} extract-region: --annotation-scope={} with legacy include/no-include flag",
                            match explicit_scope {
                                GenomeAnnotationScope::None => "none",
                                GenomeAnnotationScope::Core => "core",
                                GenomeAnnotationScope::Full => "full",
                            }
                        ));
                    }
                } else {
                    annotation_scope = Some(mapped_scope);
                }
            }
            let op_catalog_path = operation_catalog_arg(&catalog_path, helper_mode);
            let mut engine = GentleEngine::from_state(load_state(state_path)?);
            let result = engine
                .apply(Operation::ExtractGenomeRegion {
                    genome_id,
                    chromosome,
                    start_1based,
                    end_1based,
                    output_id,
                    annotation_scope,
                    max_annotation_features,
                    include_genomic_annotation,
                    catalog_path: op_catalog_path,
                    cache_dir,
                })
                .map_err(|e| e.to_string())?;
            engine
                .state()
                .save_to_path(state_path)
                .map_err(|e| e.to_string())?;
            print_json(&result)
        }
        "extract-gene" => {
            if args.len() <= cmd_idx + 3 {
                usage();
                return Err(format!(
                    "{label} extract-gene requires GENOME_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = args[cmd_idx + 2].clone();
            let gene_query = args[cmd_idx + 3].clone();
            let mut occurrence: Option<usize> = None;
            let mut output_id: Option<String> = None;
            let mut extract_mode: Option<GenomeGeneExtractMode> = None;
            let mut promoter_upstream_bp: Option<usize> = None;
            let mut annotation_scope: Option<GenomeAnnotationScope> = None;
            let mut max_annotation_features: Option<usize> = None;
            let mut include_genomic_annotation: Option<bool> = None;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = cmd_idx + 4;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--occurrence" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing N after --occurrence for {label} extract-gene"
                            ));
                        }
                        let occ = args[idx + 1].parse::<usize>().map_err(|e| {
                            format!("Invalid --occurrence value '{}': {}", args[idx + 1], e)
                        })?;
                        if occ == 0 {
                            return Err("--occurrence must be >= 1".to_string());
                        }
                        occurrence = Some(occ);
                        idx += 2;
                    }
                    "--output-id" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing ID after --output-id for {label} extract-gene"
                            ));
                        }
                        output_id = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--extract-mode" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing VALUE after --extract-mode for {label} extract-gene"
                            ));
                        }
                        let parsed = match args[idx + 1].trim().to_ascii_lowercase().as_str() {
                            "gene" => GenomeGeneExtractMode::Gene,
                            "coding_with_promoter" => GenomeGeneExtractMode::CodingWithPromoter,
                            other => {
                                return Err(format!(
                                    "Invalid --extract-mode value '{}' for {label} extract-gene (expected gene|coding_with_promoter)",
                                    other
                                ));
                            }
                        };
                        extract_mode = Some(parsed);
                        idx += 2;
                    }
                    "--promoter-upstream-bp" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing N after --promoter-upstream-bp for {label} extract-gene"
                            ));
                        }
                        let parsed = args[idx + 1].parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --promoter-upstream-bp value '{}' for {label} extract-gene: {}",
                                args[idx + 1],
                                e
                            )
                        })?;
                        promoter_upstream_bp = Some(parsed);
                        idx += 2;
                    }
                    "--annotation-scope" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing VALUE after --annotation-scope for {label} extract-gene"
                            ));
                        }
                        let parsed = match args[idx + 1].trim().to_ascii_lowercase().as_str() {
                            "none" => GenomeAnnotationScope::None,
                            "core" => GenomeAnnotationScope::Core,
                            "full" => GenomeAnnotationScope::Full,
                            other => {
                                return Err(format!(
                                    "Invalid --annotation-scope value '{}' for {label} extract-gene (expected none|core|full)",
                                    other
                                ));
                            }
                        };
                        annotation_scope = Some(parsed);
                        idx += 2;
                    }
                    "--max-annotation-features" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing N after --max-annotation-features for {label} extract-gene"
                            ));
                        }
                        let parsed = args[idx + 1].parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --max-annotation-features value '{}' for {label} extract-gene: {}",
                                args[idx + 1],
                                e
                            )
                        })?;
                        max_annotation_features = Some(parsed);
                        idx += 2;
                    }
                    "--include-genomic-annotation" => {
                        include_genomic_annotation = Some(true);
                        idx += 1;
                    }
                    "--no-include-genomic-annotation" => {
                        include_genomic_annotation = Some(false);
                        idx += 1;
                    }
                    "--catalog" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing PATH after --catalog for {label} extract-gene"
                            ));
                        }
                        catalog_path = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--cache-dir" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing PATH after --cache-dir for {label} extract-gene"
                            ));
                        }
                        cache_dir = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for {label} extract-gene",
                            other
                        ));
                    }
                }
            }
            if let Some(include) = include_genomic_annotation {
                let mapped_scope = if include {
                    GenomeAnnotationScope::Core
                } else {
                    GenomeAnnotationScope::None
                };
                if let Some(explicit_scope) = annotation_scope {
                    if explicit_scope != mapped_scope {
                        return Err(format!(
                            "Conflicting annotation options for {label} extract-gene: --annotation-scope={} with legacy include/no-include flag",
                            explicit_scope.as_str()
                        ));
                    }
                } else {
                    annotation_scope = Some(mapped_scope);
                }
            }
            if promoter_upstream_bp.is_some() && extract_mode.is_none() {
                extract_mode = Some(GenomeGeneExtractMode::CodingWithPromoter);
            }
            if matches!(extract_mode, Some(GenomeGeneExtractMode::Gene))
                && promoter_upstream_bp.unwrap_or(0) > 0
            {
                return Err(format!(
                    "--promoter-upstream-bp requires --extract-mode coding_with_promoter for {label} extract-gene"
                ));
            }
            let op_catalog_path = operation_catalog_arg(&catalog_path, helper_mode);
            let mut engine = GentleEngine::from_state(load_state(state_path)?);
            let result = engine
                .apply(Operation::ExtractGenomeGene {
                    genome_id,
                    gene_query,
                    occurrence,
                    output_id,
                    extract_mode,
                    promoter_upstream_bp,
                    annotation_scope,
                    max_annotation_features,
                    include_genomic_annotation,
                    catalog_path: op_catalog_path,
                    cache_dir,
                })
                .map_err(|e| e.to_string())?;
            engine
                .state()
                .save_to_path(state_path)
                .map_err(|e| e.to_string())?;
            print_json(&result)
        }
        "extract-promoter" => {
            if args.len() <= cmd_idx + 3 {
                usage();
                return Err(format!(
                    "{label} extract-promoter requires GENOME_ID QUERY [--occurrence N] [--transcript-id ID] [--output-id ID] [--upstream-bp N] [--downstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = args[cmd_idx + 2].clone();
            let gene_query = args[cmd_idx + 3].clone();
            let mut occurrence: Option<usize> = None;
            let mut transcript_id: Option<String> = None;
            let mut output_id: Option<String> = None;
            let mut upstream_bp: Option<usize> = None;
            let mut downstream_bp: Option<usize> = None;
            let mut annotation_scope: Option<GenomeAnnotationScope> = None;
            let mut max_annotation_features: Option<usize> = None;
            let mut include_genomic_annotation: Option<bool> = None;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = cmd_idx + 4;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--occurrence" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing N after --occurrence for {label} extract-promoter"
                            ));
                        }
                        let occ = args[idx + 1].parse::<usize>().map_err(|e| {
                            format!("Invalid --occurrence value '{}': {}", args[idx + 1], e)
                        })?;
                        if occ == 0 {
                            return Err("--occurrence must be >= 1".to_string());
                        }
                        occurrence = Some(occ);
                        idx += 2;
                    }
                    "--transcript-id" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing ID after --transcript-id for {label} extract-promoter"
                            ));
                        }
                        transcript_id = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--output-id" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing ID after --output-id for {label} extract-promoter"
                            ));
                        }
                        output_id = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--upstream-bp" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing N after --upstream-bp for {label} extract-promoter"
                            ));
                        }
                        let parsed = args[idx + 1].parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --upstream-bp value '{}' for {label} extract-promoter: {}",
                                args[idx + 1],
                                e
                            )
                        })?;
                        upstream_bp = Some(parsed);
                        idx += 2;
                    }
                    "--downstream-bp" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing N after --downstream-bp for {label} extract-promoter"
                            ));
                        }
                        let parsed = args[idx + 1].parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --downstream-bp value '{}' for {label} extract-promoter: {}",
                                args[idx + 1],
                                e
                            )
                        })?;
                        downstream_bp = Some(parsed);
                        idx += 2;
                    }
                    "--annotation-scope" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing VALUE after --annotation-scope for {label} extract-promoter"
                            ));
                        }
                        let parsed = match args[idx + 1].trim().to_ascii_lowercase().as_str() {
                            "none" => GenomeAnnotationScope::None,
                            "core" => GenomeAnnotationScope::Core,
                            "full" => GenomeAnnotationScope::Full,
                            other => {
                                return Err(format!(
                                    "Invalid --annotation-scope value '{}' for {label} extract-promoter (expected none|core|full)",
                                    other
                                ));
                            }
                        };
                        annotation_scope = Some(parsed);
                        idx += 2;
                    }
                    "--max-annotation-features" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing N after --max-annotation-features for {label} extract-promoter"
                            ));
                        }
                        let parsed = args[idx + 1].parse::<usize>().map_err(|e| {
                            format!(
                                "Invalid --max-annotation-features value '{}' for {label} extract-promoter: {}",
                                args[idx + 1],
                                e
                            )
                        })?;
                        max_annotation_features = Some(parsed);
                        idx += 2;
                    }
                    "--include-genomic-annotation" => {
                        include_genomic_annotation = Some(true);
                        idx += 1;
                    }
                    "--no-include-genomic-annotation" => {
                        include_genomic_annotation = Some(false);
                        idx += 1;
                    }
                    "--catalog" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing PATH after --catalog for {label} extract-promoter"
                            ));
                        }
                        catalog_path = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--cache-dir" => {
                        if idx + 1 >= args.len() {
                            return Err(format!(
                                "Missing PATH after --cache-dir for {label} extract-promoter"
                            ));
                        }
                        cache_dir = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for {label} extract-promoter",
                            other
                        ));
                    }
                }
            }
            if let Some(include) = include_genomic_annotation {
                let mapped_scope = if include {
                    GenomeAnnotationScope::Core
                } else {
                    GenomeAnnotationScope::None
                };
                if let Some(explicit_scope) = annotation_scope {
                    if explicit_scope != mapped_scope {
                        return Err(format!(
                            "Conflicting annotation options for {label} extract-promoter: --annotation-scope={} with legacy include/no-include flag",
                            explicit_scope.as_str()
                        ));
                    }
                } else {
                    annotation_scope = Some(mapped_scope);
                }
            }
            let op_catalog_path = operation_catalog_arg(&catalog_path, helper_mode);
            let mut engine = GentleEngine::from_state(load_state(state_path)?);
            let result = engine
                .apply(Operation::ExtractGenomePromoterSlice {
                    genome_id,
                    gene_query,
                    occurrence,
                    transcript_id,
                    output_id,
                    upstream_bp: upstream_bp.unwrap_or(DEFAULT_PROMOTER_EXTRACT_UPSTREAM_BP),
                    downstream_bp: downstream_bp.unwrap_or(DEFAULT_PROMOTER_EXTRACT_DOWNSTREAM_BP),
                    annotation_scope,
                    max_annotation_features,
                    include_genomic_annotation,
                    catalog_path: op_catalog_path,
                    cache_dir,
                })
                .map_err(|e| e.to_string())?;
            engine
                .state()
                .save_to_path(state_path)
                .map_err(|e| e.to_string())?;
            print_json(&result)
        }
        other => {
            let helper_only = if helper_mode {
                ", doctor-catalog, show-card"
            } else {
                ""
            };
            Err(format!(
                "Unknown {label} subcommand '{}' (expected ensembl-available, list, validate-catalog{helper_only}, update-ensembl-specs, status, genes, prepare, remove-prepared, remove-catalog-entry, extract-region, extract-gene, extract-promoter)",
                other
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn global_args() -> GlobalCliArgs {
        GlobalCliArgs {
            state_path: DEFAULT_STATE_PATH.to_string(),
            cmd_idx: 1,
            progress_sink: None,
            allow_screenshots: false,
        }
    }

    fn argv(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    #[test]
    fn direct_reference_family_requires_subcommand() {
        let args = argv(&["gentle_cli", "genomes"]);
        let err = handle_reference_family(&args, 1, "genomes", DEFAULT_STATE_PATH, &global_args())
            .expect_err("missing subcommand should fail");
        assert_eq!(err, "genomes requires a subcommand");
    }

    #[test]
    fn direct_helper_vocabulary_parser_reports_missing_filter_value() {
        let args = argv(&["gentle_cli", "helpers", "vocabulary", "list", "--filter"]);
        let err = handle_reference_family(&args, 1, "helpers", DEFAULT_STATE_PATH, &global_args())
            .expect_err("missing filter value should fail");
        assert_eq!(
            err,
            "Missing TEXT after --filter for helpers vocabulary list"
        );
    }

    #[test]
    fn direct_reference_genes_parser_reports_invalid_limit() {
        let args = argv(&["gentle_cli", "genomes", "genes", "demo", "--limit", "abc"]);
        let err = handle_reference_family(&args, 1, "genomes", DEFAULT_STATE_PATH, &global_args())
            .expect_err("invalid limit should fail");
        assert!(err.starts_with("Invalid --limit value 'abc': "));
    }
}
