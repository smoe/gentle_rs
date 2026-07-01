//! Direct pool, gel, and arrangement CLI command handling kept out of `run()`.
//!
//! These routes share pool/arrangement state and gel export behavior but are
//! still legacy direct CLI commands rather than one subcommand tree. Keeping
//! them together removes another sizeable command island from the binary entry
//! point while preserving the existing top-level command names.

use gentle::{
    dna_sequence::DNAsequence,
    engine::{GelBufferModel, GelRunConditions, GelTopologyForm, GentleEngine, Operation},
};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{collections::HashMap, fs};

use super::*;

const POOL_COMMANDS: &[&str] = &[
    "render-pool-gel-svg",
    "render-gel-svg",
    "arrange-serial",
    "arrange-set-ladders",
    "export-pool",
    "export-run-bundle",
    "import-pool",
];

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PoolEnd {
    end_type: String,
    forward_5: String,
    forward_3: String,
    reverse_5: String,
    reverse_3: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PoolMember {
    seq_id: String,
    #[serde(default)]
    human_id: String,
    name: Option<String>,
    sequence: String,
    length_bp: usize,
    topology: String,
    ends: PoolEnd,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PoolExport {
    schema: String,
    pool_id: String,
    #[serde(default)]
    human_id: String,
    member_count: usize,
    members: Vec<PoolMember>,
}

pub(super) fn is_pool_command(command: &str) -> bool {
    POOL_COMMANDS.contains(&command)
}

pub(super) fn handle_pool_family(
    command: &str,
    args: &[String],
    cmd_idx: usize,
    state_path: &str,
) -> Result<(), String> {
    match command {
        "render-pool-gel-svg" | "render-gel-svg" => {
            handle_render_pool_gel_svg(command, args, cmd_idx, state_path)
        }
        "arrange-serial" => handle_arrange_serial(args, cmd_idx, state_path),
        "arrange-set-ladders" => handle_arrange_set_ladders(args, cmd_idx, state_path),
        "export-pool" => handle_export_pool(args, cmd_idx, state_path),
        "export-run-bundle" => handle_export_run_bundle(args, cmd_idx, state_path),
        "import-pool" => handle_import_pool(args, cmd_idx, state_path),
        _ => Err(format!("Expected pool command, got '{command}'")),
    }
}

fn handle_render_pool_gel_svg(
    cmd_name: &str,
    args: &[String],
    cmd_idx: usize,
    state_path: &str,
) -> Result<(), String> {
    if args.len() <= cmd_idx + 2 {
        usage();
        return Err(format!(
            "{cmd_name} requires: IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID] [--agarose-pct FLOAT] [--buffer tae|tbe] [--topology-aware true|false]"
        ));
    }
    let ids = match args[cmd_idx + 1].trim() {
        "-" | "_" => vec![],
        raw => comma_items(raw),
    };
    let output = &args[cmd_idx + 2];
    let mut ladders: Option<Vec<String>> = None;
    let mut container_ids: Option<Vec<String>> = None;
    let mut arrangement_id: Option<String> = None;
    let mut conditions = GelRunConditions::default();
    let mut idx = cmd_idx + 3;
    while idx < args.len() {
        match args[idx].as_str() {
            "--ladders" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --ladders".to_string());
                }
                let parsed = comma_items(&args[idx + 1]);
                if !parsed.is_empty() {
                    ladders = Some(parsed);
                }
                idx += 2;
            }
            "--containers" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --containers".to_string());
                }
                let parsed = comma_items(&args[idx + 1]);
                if !parsed.is_empty() {
                    container_ids = Some(parsed);
                }
                idx += 2;
            }
            "--arrangement" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --arrangement".to_string());
                }
                let value = args[idx + 1].trim();
                if !value.is_empty() {
                    arrangement_id = Some(value.to_string());
                }
                idx += 2;
            }
            "--agarose-pct" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --agarose-pct".to_string());
                }
                conditions.agarose_percent = args[idx + 1]
                    .trim()
                    .parse::<f32>()
                    .map_err(|e| format!("Invalid agarose percent '{}': {e}", args[idx + 1]))?;
                idx += 2;
            }
            "--buffer" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --buffer".to_string());
                }
                conditions.buffer_model = match args[idx + 1].trim().to_ascii_lowercase().as_str() {
                    "tae" => GelBufferModel::Tae,
                    "tbe" => GelBufferModel::Tbe,
                    other => {
                        return Err(format!(
                            "Unknown buffer '{other}' for {cmd_name} (expected tae|tbe)"
                        ));
                    }
                };
                idx += 2;
            }
            "--topology-aware" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --topology-aware".to_string());
                }
                conditions.topology_aware = parse_bool_flag(&args[idx + 1]).ok_or_else(|| {
                    format!(
                        "Boolean value '{}' is invalid (expected true|false|1|0|yes|no|on|off)",
                        args[idx + 1]
                    )
                })?;
                idx += 2;
            }
            other => {
                return Err(format!(
                    "Unknown argument '{other}' for {cmd_name} (expected --ladders/--containers/--arrangement/--agarose-pct/--buffer/--topology-aware)"
                ));
            }
        }
    }
    if ids.is_empty()
        && container_ids.as_ref().is_none_or(|v| v.is_empty())
        && arrangement_id.is_none()
    {
        return Err(format!(
            "{cmd_name} requires sequence ids, --containers, or --arrangement"
        ));
    }
    let mut engine = GentleEngine::from_state(load_state(state_path)?);
    let result = engine
        .apply(Operation::RenderPoolGelSvg {
            inputs: ids,
            path: output.to_string(),
            ladders,
            container_ids,
            arrangement_id,
            conditions: Some(conditions.normalized()),
        })
        .map_err(|e| e.to_string())?;
    save_state_and_print_first_message(&engine, state_path, &result.messages)
}

fn handle_arrange_serial(args: &[String], cmd_idx: usize, state_path: &str) -> Result<(), String> {
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err(
            "arrange-serial requires: CONTAINER_IDS [--id ARR_ID] [--name TEXT] [--ladders NAME[,NAME]]"
                .to_string(),
        );
    }
    let container_ids = comma_items(&args[cmd_idx + 1]);
    if container_ids.is_empty() {
        return Err("arrange-serial requires at least one container id".to_string());
    }
    let mut arrangement_id: Option<String> = None;
    let mut name: Option<String> = None;
    let mut ladders: Option<Vec<String>> = None;
    let mut idx = cmd_idx + 2;
    while idx < args.len() {
        match args[idx].as_str() {
            "--id" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --id".to_string());
                }
                let value = args[idx + 1].trim();
                if !value.is_empty() {
                    arrangement_id = Some(value.to_string());
                }
                idx += 2;
            }
            "--name" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --name".to_string());
                }
                let value = args[idx + 1].trim();
                if !value.is_empty() {
                    name = Some(value.to_string());
                }
                idx += 2;
            }
            "--ladders" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --ladders".to_string());
                }
                let parsed = comma_items(&args[idx + 1]);
                if !parsed.is_empty() {
                    ladders = Some(parsed);
                }
                idx += 2;
            }
            other => {
                return Err(format!(
                    "Unknown argument '{other}' for arrange-serial (expected --id/--name/--ladders)"
                ));
            }
        }
    }
    let mut engine = GentleEngine::from_state(load_state(state_path)?);
    let result = engine
        .apply(Operation::CreateArrangementSerial {
            container_ids,
            arrangement_id,
            name,
            ladders,
        })
        .map_err(|e| e.to_string())?;
    save_state_and_print_first_message(&engine, state_path, &result.messages)
}

fn handle_arrange_set_ladders(
    args: &[String],
    cmd_idx: usize,
    state_path: &str,
) -> Result<(), String> {
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err("arrange-set-ladders requires: ARR_ID [--ladders NAME[,NAME]]".to_string());
    }
    let arrangement_id = args[cmd_idx + 1].trim().to_string();
    if arrangement_id.is_empty() {
        return Err("arrange-set-ladders requires a non-empty ARR_ID".to_string());
    }
    let mut ladders: Option<Vec<String>> = None;
    let mut idx = cmd_idx + 2;
    while idx < args.len() {
        match args[idx].as_str() {
            "--ladders" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --ladders".to_string());
                }
                ladders = Some(comma_items(&args[idx + 1]));
                idx += 2;
            }
            other => {
                return Err(format!(
                    "Unknown argument '{other}' for arrange-set-ladders (expected --ladders)"
                ));
            }
        }
    }
    let mut engine = GentleEngine::from_state(load_state(state_path)?);
    let result = engine
        .apply(Operation::SetArrangementLadders {
            arrangement_id,
            ladders: ladders.filter(|values| !values.is_empty()),
        })
        .map_err(|e| e.to_string())?;
    save_state_and_print_first_message(&engine, state_path, &result.messages)
}

fn handle_export_pool(args: &[String], cmd_idx: usize, state_path: &str) -> Result<(), String> {
    if args.len() <= cmd_idx + 2 {
        usage();
        return Err("export-pool requires: IDS OUTPUT.pool.gentle.json [HUMAN_ID]".to_string());
    }
    let ids = comma_items(&args[cmd_idx + 1]);
    if ids.is_empty() {
        return Err("export-pool requires at least one sequence id".to_string());
    }
    let output = &args[cmd_idx + 2];
    let human_id = args.get(cmd_idx + 3).cloned();
    let mut engine = GentleEngine::from_state(load_state(state_path)?);
    let result = engine
        .apply(Operation::ExportPool {
            inputs: ids,
            path: output.to_string(),
            pool_id: Some("pool_export".to_string()),
            human_id,
        })
        .map_err(|e| e.to_string())?;
    engine
        .state()
        .save_to_path(state_path)
        .map_err(|e| e.to_string())?;
    print_json(&result)
}

fn handle_export_run_bundle(
    args: &[String],
    cmd_idx: usize,
    state_path: &str,
) -> Result<(), String> {
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err(
            "export-run-bundle requires: OUTPUT.run_bundle.json [--run-id RUN_ID]".to_string(),
        );
    }
    let output = args[cmd_idx + 1].clone();
    let mut run_id: Option<String> = None;
    let mut idx = cmd_idx + 2;
    while idx < args.len() {
        match args[idx].as_str() {
            "--run-id" => {
                if idx + 1 >= args.len() {
                    return Err("Missing value after --run-id".to_string());
                }
                let value = args[idx + 1].trim();
                if !value.is_empty() {
                    run_id = Some(value.to_string());
                }
                idx += 2;
            }
            other => {
                return Err(format!(
                    "Unknown argument '{other}' for export-run-bundle (expected --run-id)"
                ));
            }
        }
    }
    let mut engine = GentleEngine::from_state(load_state(state_path)?);
    let result = engine
        .apply(Operation::ExportProcessRunBundle {
            path: output,
            run_id,
        })
        .map_err(|e| e.to_string())?;
    engine
        .state()
        .save_to_path(state_path)
        .map_err(|e| e.to_string())?;
    print_json(&result)
}

fn handle_import_pool(args: &[String], cmd_idx: usize, state_path: &str) -> Result<(), String> {
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err("import-pool requires: INPUT.pool.gentle.json [PREFIX]".to_string());
    }
    let input = &args[cmd_idx + 1];
    let prefix = args
        .get(cmd_idx + 2)
        .cloned()
        .unwrap_or_else(|| "pool".to_string());
    let text = fs::read_to_string(input)
        .map_err(|e| format!("Could not read pool file '{input}': {e}"))?;
    let pool: PoolExport =
        serde_json::from_str(&text).map_err(|e| format!("Invalid pool JSON '{input}': {e}"))?;
    if pool.schema != "gentle.pool.v1" {
        return Err(format!(
            "Unsupported pool schema '{}', expected 'gentle.pool.v1'",
            pool.schema
        ));
    }

    let mut state = load_state(state_path)?;
    for (idx, member) in pool.members.iter().enumerate() {
        let mut dna = DNAsequence::from_sequence(&member.sequence)
            .map_err(|e| format!("Invalid DNA in pool member '{}': {e}", member.seq_id))?;
        if let Some(name) = &member.name {
            let mut value = serde_json::to_value(&dna)
                .map_err(|e| format!("Could not serialize sequence: {e}"))?;
            if let Some(obj) = value.as_object_mut()
                && let Some(seq_obj) = obj.get_mut("seq").and_then(|v| v.as_object_mut())
            {
                seq_obj.insert("name".to_string(), json!(name));
            }
            dna = serde_json::from_value(value)
                .map_err(|e| format!("Could not set sequence name: {e}"))?;
        }
        apply_pool_member_topology_hint(&member.topology, &mut dna)?;
        apply_member_overhang(member, &mut dna)?;
        dna.update_computed_features();
        let base = format!("{prefix}_{}", idx + 1);
        let id = unique_id(&state.sequences, &base);
        state.sequences.insert(id, dna);
    }
    state.save_to_path(state_path).map_err(|e| e.to_string())?;
    println!(
        "Imported pool '{}' ({} members) into '{}'",
        pool.pool_id, pool.member_count, state_path
    );
    Ok(())
}

fn comma_items(raw: &str) -> Vec<String> {
    raw.split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(str::to_string)
        .collect()
}

fn unique_id(existing: &HashMap<String, DNAsequence>, base: &str) -> String {
    if !existing.contains_key(base) {
        return base.to_string();
    }
    let mut i = 2usize;
    loop {
        let candidate = format!("{base}_{i}");
        if !existing.contains_key(&candidate) {
            return candidate;
        }
        i += 1;
    }
}

fn apply_pool_member_topology_hint(
    member_topology: &str,
    dna: &mut DNAsequence,
) -> Result<(), String> {
    let Some(form) = GelTopologyForm::from_hint(member_topology) else {
        return Ok(());
    };
    if form.is_circular() {
        dna.set_circular(true);
    }
    if !matches!(form, GelTopologyForm::Linear | GelTopologyForm::Circular) {
        let mut value = serde_json::to_value(&*dna)
            .map_err(|e| format!("Could not serialize sequence for topology hint: {e}"))?;
        if let Some(obj) = value.as_object_mut()
            && let Some(seq_obj) = obj.get_mut("seq").and_then(|v| v.as_object_mut())
        {
            let comments = seq_obj
                .entry("comments".to_string())
                .or_insert_with(|| json!([]));
            if let Some(arr) = comments.as_array_mut() {
                arr.push(json!(format!("gel_topology={}", form.as_str())));
            }
        }
        *dna = serde_json::from_value(value)
            .map_err(|e| format!("Could not apply topology hint to sequence: {e}"))?;
    }
    Ok(())
}

fn apply_member_overhang(member: &PoolMember, dna: &mut DNAsequence) -> Result<(), String> {
    let mut value =
        serde_json::to_value(&*dna).map_err(|e| format!("Could not serialize sequence: {e}"))?;
    let obj = value
        .as_object_mut()
        .ok_or_else(|| "Internal serialization shape error".to_string())?;
    obj.insert(
        "overhang".to_string(),
        json!({
            "forward_3": member.ends.forward_3.as_bytes(),
            "forward_5": member.ends.forward_5.as_bytes(),
            "reverse_3": member.ends.reverse_3.as_bytes(),
            "reverse_5": member.ends.reverse_5.as_bytes(),
        }),
    );
    let patched: DNAsequence =
        serde_json::from_value(value).map_err(|e| format!("Could not restore overhang: {e}"))?;
    *dna = patched;
    Ok(())
}

fn save_state_and_print_first_message(
    engine: &GentleEngine,
    state_path: &str,
    messages: &[String],
) -> Result<(), String> {
    engine
        .state()
        .save_to_path(state_path)
        .map_err(|e| e.to_string())?;
    if let Some(msg) = messages.first() {
        println!("{msg}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use gentle::engine::ProjectState;
    use tempfile::tempdir;

    fn argv(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    #[test]
    fn pool_command_filter_recognizes_moved_routes() {
        assert!(is_pool_command("render-gel-svg"));
        assert!(is_pool_command("arrange-serial"));
        assert!(is_pool_command("export-run-bundle"));
        assert!(!is_pool_command("protocol-cartoon"));
    }

    #[test]
    fn render_pool_gel_missing_args_reports_legacy_message() {
        let args = argv(&["gentle_cli", "render-pool-gel-svg"]);
        let err = handle_pool_family("render-pool-gel-svg", &args, 1, ".gentle_state.json")
            .expect_err("missing args should fail");
        assert_eq!(
            err,
            "render-pool-gel-svg requires: IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID] [--agarose-pct FLOAT] [--buffer tae|tbe] [--topology-aware true|false]"
        );
    }

    #[test]
    fn render_pool_gel_rejects_empty_inputs_without_alternates() {
        let args = argv(&["gentle_cli", "render-pool-gel-svg", "-", "out.svg"]);
        let err = handle_pool_family("render-pool-gel-svg", &args, 1, ".gentle_state.json")
            .expect_err("empty inputs should fail");
        assert_eq!(
            err,
            "render-pool-gel-svg requires sequence ids, --containers, or --arrangement"
        );
    }

    #[test]
    fn arrange_serial_rejects_empty_container_list() {
        let args = argv(&["gentle_cli", "arrange-serial", ","]);
        let err = handle_pool_family("arrange-serial", &args, 1, ".gentle_state.json")
            .expect_err("empty container list should fail");
        assert_eq!(err, "arrange-serial requires at least one container id");
    }

    #[test]
    fn export_run_bundle_rejects_unknown_option() {
        let args = argv(&[
            "gentle_cli",
            "export-run-bundle",
            "bundle.json",
            "--unknown",
        ]);
        let err = handle_pool_family("export-run-bundle", &args, 1, ".gentle_state.json")
            .expect_err("unknown option should fail");
        assert_eq!(
            err,
            "Unknown argument '--unknown' for export-run-bundle (expected --run-id)"
        );
    }

    #[test]
    fn import_pool_rejects_unsupported_schema() {
        let td = tempdir().expect("tempdir");
        let input = td.path().join("bad.pool.gentle.json");
        fs::write(
            &input,
            r#"{"schema":"other","pool_id":"pool","member_count":0,"members":[]}"#,
        )
        .expect("write pool");
        let input_arg = input.to_string_lossy().to_string();
        let args = argv(&["gentle_cli", "import-pool", &input_arg]);
        let err = handle_pool_family("import-pool", &args, 1, ".gentle_state.json")
            .expect_err("unsupported schema should fail");
        assert_eq!(
            err,
            "Unsupported pool schema 'other', expected 'gentle.pool.v1'"
        );
    }

    #[test]
    fn import_pool_writes_member_to_state() {
        let td = tempdir().expect("tempdir");
        let input = td.path().join("demo.pool.gentle.json");
        let state_path = td.path().join("state.json");
        fs::write(
            &input,
            r#"{
  "schema": "gentle.pool.v1",
  "pool_id": "demo_pool",
  "member_count": 1,
  "members": [
    {
      "seq_id": "member_a",
      "human_id": "tube A",
      "name": "Member A",
      "sequence": "ATGC",
      "length_bp": 4,
      "topology": "linear",
      "ends": {
        "end_type": "blunt",
        "forward_5": "",
        "forward_3": "",
        "reverse_5": "",
        "reverse_3": ""
      }
    }
  ]
}"#,
        )
        .expect("write pool");
        let input_arg = input.to_string_lossy().to_string();
        let args = argv(&["gentle_cli", "import-pool", &input_arg, "imported"]);
        handle_pool_family("import-pool", &args, 1, &state_path.to_string_lossy())
            .expect("pool import should execute");
        let state =
            ProjectState::load_from_path(&state_path.to_string_lossy()).expect("load state");
        assert!(
            state.sequences.contains_key("imported_1"),
            "pool import should create prefixed sequence"
        );
    }
}
