//! Validated synthetic bench allowances, separate from computational macros.

use super::CloningRoutineCatalog;
use crate::engine::{
    PlanningBenchStepSource, PlanningEstimate, PlanningTimeBasis, ScheduledStep, TimeRange,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;

pub(super) const SCHEMA: &str = "gentle.cloning_wetlab_steps.v1";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Catalog {
    pub schema: String,
    pub routines: Vec<Routine>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Routine {
    pub routine_id: String,
    pub steps: Vec<Step>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Step {
    pub step_id: String,
    pub label: String,
    pub hands_on_hours: TimeRange,
    pub unattended_hours: TimeRange,
    pub depends_on: Vec<String>,
    pub source: PlanningBenchStepSource,
}

/// A load/admission failure, serialized intact through the shared shell error.
#[derive(Debug, Clone, Serialize)]
pub(super) struct CatalogError {
    pub code: &'static str,
    pub routine_id: Option<String>,
    pub step_id: Option<String>,
    pub detail: String,
}

impl std::fmt::Display for CatalogError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            serde_json::to_string(self).map_err(|_| std::fmt::Error)?
        )
    }
}

fn error(code: &'static str, routine: &str, step: &str, detail: &str) -> CatalogError {
    CatalogError {
        code,
        routine_id: (!routine.is_empty()).then(|| routine.into()),
        step_id: (!step.is_empty()).then(|| step.into()),
        detail: detail.into(),
    }
}

fn valid_range(range: TimeRange) -> bool {
    [range.min, range.typical, range.max]
        .iter()
        .all(|value| value.is_finite() && *value >= 0.0)
        && range.min <= range.typical
        && range.typical <= range.max
}

fn add(left: TimeRange, right: TimeRange) -> TimeRange {
    TimeRange {
        min: left.min + right.min,
        typical: left.typical + right.typical,
        max: left.max + right.max,
    }
}

fn components(range: TimeRange) -> [f64; 3] {
    [range.min, range.typical, range.max]
}

pub(super) struct Schedule {
    pub hands_on: TimeRange,
    pub unattended: TimeRange,
    pub elapsed: TimeRange,
    pub rows: Vec<ScheduledStep>,
    pub critical_paths: [Vec<String>; 3],
}

impl Catalog {
    pub(super) fn parse(raw: &str, known_ids: &BTreeSet<String>) -> Result<Self, CatalogError> {
        let catalog: Self = serde_json::from_str(raw)
            .map_err(|err| error("invalid_wetlab_catalog", "", "", &err.to_string()))?;
        if catalog.schema != SCHEMA {
            return Err(error("unsupported_wetlab_schema", "", "", &catalog.schema));
        }
        let mut ids = BTreeSet::new();
        for routine in &catalog.routines {
            if !known_ids.contains(&routine.routine_id) {
                return Err(error(
                    "unknown_routine_id",
                    &routine.routine_id,
                    "",
                    "No exact bundled routine ID",
                ));
            }
            if !ids.insert(&routine.routine_id) {
                return Err(error(
                    "duplicate_wetlab_routine",
                    &routine.routine_id,
                    "",
                    "Duplicate routine model",
                ));
            }
            routine.schedule(0.0)?;
        }
        Ok(catalog)
    }

    pub(super) fn model(&self, id: &str) -> Option<&Routine> {
        self.routines.iter().find(|row| row.routine_id == id)
    }
}

impl Routine {
    pub(super) fn schedule(&self, delay_hours: f64) -> Result<Schedule, CatalogError> {
        let fail = |code, step: &str, detail: &str| error(code, &self.routine_id, step, detail);
        if !delay_hours.is_finite() || delay_hours < 0.0 {
            return Err(fail(
                "invalid_planning_delay",
                "",
                "Delay must be finite and nonnegative",
            ));
        }
        if self.steps.is_empty() {
            return Err(fail(
                "empty_wetlab_model",
                "",
                "An empty model cannot assert zero bench work",
            ));
        }
        let mut pending = BTreeMap::new();
        for step in &self.steps {
            let mut bytes = step.step_id.bytes();
            let valid_id = bytes.next().is_some_and(|b| b.is_ascii_lowercase())
                && bytes.all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'_');
            if !valid_id || step.label.trim().is_empty() {
                return Err(fail(
                    "invalid_wetlab_step",
                    &step.step_id,
                    "A snake_case ID and label are required",
                ));
            }
            if pending.insert(step.step_id.clone(), step).is_some() {
                return Err(fail(
                    "duplicate_wetlab_step",
                    &step.step_id,
                    "Duplicate step ID",
                ));
            }
            if !valid_range(step.hands_on_hours) || !valid_range(step.unattended_hours) {
                return Err(fail(
                    "invalid_wetlab_duration",
                    &step.step_id,
                    "Require finite 0 <= min <= typical <= max",
                ));
            }
            if step.source.kind != "synthetic_planning_default"
                || step.source.note.trim().is_empty()
            {
                return Err(fail(
                    "invalid_wetlab_source",
                    &step.step_id,
                    "This catalog admits documented synthetic planning defaults only",
                ));
            }
        }
        for step in &self.steps {
            let mut seen = BTreeSet::new();
            for dependency in &step.depends_on {
                if !pending.contains_key(dependency) {
                    return Err(fail("unknown_wetlab_dependency", &step.step_id, dependency));
                }
                if !seen.insert(dependency) {
                    return Err(fail(
                        "duplicate_wetlab_dependency",
                        &step.step_id,
                        dependency,
                    ));
                }
            }
        }
        let mut result = Schedule {
            hands_on: TimeRange::default(),
            unattended: TimeRange::default(),
            elapsed: TimeRange {
                min: delay_hours,
                typical: delay_hours,
                max: delay_hours,
            },
            rows: vec![],
            critical_paths: Default::default(),
        };
        let mut completed: BTreeMap<String, (TimeRange, [Vec<String>; 3])> = BTreeMap::new();
        while !pending.is_empty() {
            let Some(id) = pending
                .iter()
                .find(|(_, step)| step.depends_on.iter().all(|id| completed.contains_key(id)))
                .map(|(id, _)| id.clone())
            else {
                return Err(fail(
                    "wetlab_dependency_cycle",
                    "",
                    "Dependency graph contains a cycle",
                ));
            };
            let step = pending.remove(&id).expect("Selected pending step exists");
            let mut start = [delay_hours; 3];
            let mut paths: [Vec<String>; 3] = Default::default();
            // Each scenario has its own longest predecessor path; equal paths
            // use lexical ID order, independent of JSON row/dependency order.
            let mut dependencies = step.depends_on.clone();
            dependencies.sort();
            for dependency in &dependencies {
                let (finish, previous_paths) = &completed[dependency];
                for (index, value) in components(*finish).into_iter().enumerate() {
                    if value > start[index]
                        || (value == start[index]
                            && (paths[index].is_empty() || previous_paths[index] < paths[index]))
                    {
                        start[index] = value;
                        paths[index] = previous_paths[index].clone();
                    }
                }
            }
            let start_range = TimeRange {
                min: start[0],
                typical: start[1],
                max: start[2],
            };
            let finish = add(start_range, add(step.hands_on_hours, step.unattended_hours));
            result.hands_on = add(result.hands_on, step.hands_on_hours);
            result.unattended = add(result.unattended, step.unattended_hours);
            if !valid_range(finish)
                || !valid_range(result.hands_on)
                || !valid_range(result.unattended)
            {
                return Err(fail(
                    "wetlab_duration_overflow",
                    &id,
                    "Accumulated durations are not finite",
                ));
            }
            for path in &mut paths {
                path.push(id.clone());
            }
            let elapsed = components(result.elapsed);
            let mut maxima = elapsed;
            for (index, value) in components(finish).into_iter().enumerate() {
                if value > elapsed[index]
                    || (value == elapsed[index]
                        && (result.critical_paths[index].is_empty()
                            || paths[index] < result.critical_paths[index]))
                {
                    maxima[index] = value;
                    result.critical_paths[index] = paths[index].clone();
                }
            }
            result.elapsed = TimeRange {
                min: maxima[0],
                typical: maxima[1],
                max: maxima[2],
            };
            result.rows.push(ScheduledStep {
                step_id: id.clone(),
                label: step.label.clone(),
                hands_on_hours: step.hands_on_hours,
                unattended_hours: step.unattended_hours,
                earliest_start_hours: start_range,
                earliest_finish_hours: finish,
                depends_on: dependencies,
                source: step.source.clone(),
            });
            completed.insert(id, (finish, paths));
        }
        Ok(result)
    }
}

pub(super) fn bundled() -> Result<&'static Catalog, String> {
    static CATALOG: OnceLock<Result<Catalog, CatalogError>> = OnceLock::new();
    CATALOG
        .get_or_init(|| {
            let routines: CloningRoutineCatalog =
                serde_json::from_str(include_str!("../../assets/cloning_routines.json"))
                    .map_err(|err| error("invalid_bundled_routines", "", "", &err.to_string()))?;
            let known = routines
                .routines
                .into_iter()
                .map(|row| row.routine_id)
                .collect();
            Catalog::parse(
                include_str!("../../assets/cloning_wetlab_steps.json"),
                &known,
            )
        })
        .as_ref()
        .map_err(ToString::to_string)
}

fn range_text(range: Option<TimeRange>) -> String {
    range
        .map(|range| {
            format!(
                "{:.2} h (min {:.2}, max {:.2})",
                range.typical, range.min, range.max
            )
        })
        .unwrap_or_else(|| "unknown".into())
}

fn basis_key(estimate: &PlanningEstimate) -> &'static str {
    match estimate.estimated_time_hours_basis {
        PlanningTimeBasis::StepModelCriticalPath => "step_model_critical_path",
        PlanningTimeBasis::LegacyFamilyScalar => "legacy_family_scalar",
    }
}

fn policy_time_hours(estimate: &PlanningEstimate) -> Option<f64> {
    let policy = estimate.explanation.get("time_policy")?;
    if policy.get("basis")?.as_str()? != "legacy_family_scalar" {
        return None;
    }
    policy
        .get("hours")?
        .as_f64()
        .filter(|hours| hours.is_finite() && *hours >= 0.0)
}

/// Thin adapters display these shared labels without deriving another estimate.
pub(crate) fn estimate_display_lines(estimate: &PlanningEstimate) -> Vec<String> {
    let basis = match estimate.estimated_time_hours_basis {
        PlanningTimeBasis::StepModelCriticalPath => "synthetic step-model critical path",
        PlanningTimeBasis::LegacyFamilyScalar => "legacy family scalar; breakdown unknown",
    };
    let mut lines = vec![
        format!(
            "Elapsed: {:.2} h ({basis}); range {}",
            estimate.estimated_time_hours,
            range_text(estimate.elapsed_hours)
        ),
        format!(
            "Hands-on: {} | Unattended work total (may overlap): {}",
            range_text(estimate.hands_on_hours),
            range_text(estimate.unattended_hours)
        ),
    ];
    lines.push(match policy_time_hours(estimate) {
        Some(hours) => format!("Ranking/time-guardrail proxy: {hours:.2} h (legacy family/base scalar + delays; not measured laboratory time)"),
        None => "Ranking/time-guardrail proxy: unknown (not recorded)".into(),
    });
    if estimate.schema != crate::engine::PLANNING_ESTIMATE_SCHEMA {
        lines.push(format!(
            "Legacy total cost: {:.2} (currency and assumptions as recorded)",
            estimate.estimated_cost
        ));
        lines.push(
            "Cost components: unknown (legacy estimate); no v2 cost breakdown was recorded".into(),
        );
        return lines;
    }
    lines.push(format!(
        "Stock-drawn material: {:.2} | Missing-material purchase proxy: {:.2}",
        estimate.stock_drawn_material_cost, estimate.missing_material_purchase_proxy_cost
    ));
    lines.push(format!(
        "Machine run: {:.2} | Routine base: {:.2} | Known-unit total cost: {:.2} ({})",
        estimate.machine_run_cost,
        estimate.routine_base_cost,
        estimate.estimated_cost,
        estimate
            .explanation
            .get("currency")
            .and_then(|value| value.as_str())
            .unwrap_or("currency unspecified")
    ));
    lines.push("Material subtotals are incomplete estimates, not package prices or orderable quotes; machine occupancy uses the legacy duration basis.".into());
    lines
}

pub(super) fn comparison_axes(
    left: &PlanningEstimate,
    right: &PlanningEstimate,
) -> Vec<serde_json::Value> {
    use serde_json::json;
    let mut rows = vec![
        json!({"axis":"estimated_time_hours_basis", "left":basis_key(left), "right":basis_key(right),
        "same":left.estimated_time_hours_basis == right.estimated_time_hours_basis}),
    ];
    let (a, b) = (policy_time_hours(left), policy_time_hours(right));
    rows.push(json!({
        "axis":"time_policy_proxy_hours", "basis":"legacy_family_scalar",
        "left":a.map(|hours| format!("{hours:.2} h (legacy policy proxy)")).unwrap_or_else(|| "unknown".into()),
        "right":b.map(|hours| format!("{hours:.2} h (legacy policy proxy)")).unwrap_or_else(|| "unknown".into()),
        "left_hours":a, "right_hours":b, "comparable":a.is_some() && b.is_some(),
        "same":matches!((a,b), (Some(a),Some(b)) if a == b)
    }));
    for (axis, a, b) in [
        ("hands_on_hours", left.hands_on_hours, right.hands_on_hours),
        (
            "unattended_hours",
            left.unattended_hours,
            right.unattended_hours,
        ),
        ("elapsed_hours", left.elapsed_hours, right.elapsed_hours),
    ] {
        rows.push(
            json!({"axis":axis, "left":range_text(a), "right":range_text(b),
            "left_range":a, "right_range":b, "comparable":a.is_some() && b.is_some(),
            "same": matches!((a,b), (Some(a),Some(b)) if a == b)}),
        );
    }
    for (axis, a, b) in [
        (
            "stock_drawn_material_cost",
            left.stock_drawn_material_cost,
            right.stock_drawn_material_cost,
        ),
        (
            "missing_material_purchase_proxy_cost",
            left.missing_material_purchase_proxy_cost,
            right.missing_material_purchase_proxy_cost,
        ),
        (
            "machine_run_cost",
            left.machine_run_cost,
            right.machine_run_cost,
        ),
        (
            "routine_base_cost",
            left.routine_base_cost,
            right.routine_base_cost,
        ),
    ] {
        rows.push(json!({"axis":axis, "left":format!("{a:.2}"), "right":format!("{b:.2}"), "same":(a-b).abs() < 1e-9}));
    }
    rows
}
