//! Native reader for reduced promoter-cofactor packages, not a scan provider.
//! Source SQL, databases and scripts are validated as files, never executed.

use crate::digest_utils::{sha256_file_hex, sha256_hex_bytes, short_sha256_id};
use crate::genomic_motif_evidence::query_local_parquet_json;
use gentle_protocol::promoter_cofactors::*;
use serde::{Deserialize, de::DeserializeOwned};
use serde_json::Value;
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::Read,
    path::{Component, Path, PathBuf},
    time::{Duration, Instant},
};

const BANDS: [&str; 6] = [
    "overlap",
    "adjacent_0_5",
    "gap_6_20",
    "gap_21_50",
    "gap_51_100",
    "gap_101_150",
];
const MAX_METADATA: u64 = 4 * 1024 * 1024;

#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
struct FileEntry {
    path: String,
    bytes: u64,
    sha256: String,
}
#[derive(Deserialize)]
struct PanelEntry {
    motif_id: String,
}
#[derive(Deserialize)]
struct Manifest {
    kind: String,
    schema_version: u32,
    state: String,
    assembly: String,
    taxon_id: u64,
    coordinate_mode: String,
    complete_genome_scan: bool,
    chromosomes: Vec<String>,
    distance_bands: Vec<String>,
    panel: Vec<PanelEntry>,
    source_score_floor: f64,
    positive_threshold: f64,
    score_configuration: BTreeMap<String, Value>,
    retention: String,
    h3k4me3_model_effects: String,
    requested_genes: Vec<CofactorCandidate>,
    identity: BTreeMap<String, Value>,
    files: Vec<FileEntry>,
}
#[derive(Deserialize)]
struct Complete {
    identity: BTreeMap<String, Value>,
    files: Vec<FileEntry>,
}

struct Package {
    root: PathBuf,
    manifest: Manifest,
    files: BTreeMap<String, FileEntry>,
}

fn read_metadata(path: &Path) -> Result<Vec<u8>, String> {
    let mut bytes = Vec::new();
    fs::File::open(path)
        .map_err(|e| format!("{}: {e}", path.display()))?
        .take(MAX_METADATA + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
    if bytes.len() as u64 > MAX_METADATA {
        return Err("Metadata exceeds 4 MiB".into());
    }
    Ok(bytes)
}

fn checked_path(root: &Path, name: &str) -> Result<PathBuf, String> {
    if name.is_empty()
        || Path::new(name)
            .components()
            .any(|c| !matches!(c, Component::Normal(_)))
    {
        return Err(format!("Unsafe package path: {name}"));
    }
    let path = root
        .join(name)
        .canonicalize()
        .map_err(|e| format!("Missing package file {name}: {e}"))?;
    if !path.starts_with(root) || !path.is_file() {
        return Err(format!(
            "Package path escapes root or is not a file: {name}"
        ));
    }
    Ok(path)
}

fn validate_package(
    request: &PromoterCofactorRequest,
    report: &mut PromoterCofactorReport,
    deadline: Instant,
) -> Result<Package, String> {
    let root = Path::new(&request.package_path)
        .canonicalize()
        .map_err(|e| e.to_string())?;
    let manifest_bytes = read_metadata(&checked_path(&root, "manifest.json")?)?;
    let complete_bytes = read_metadata(&checked_path(&root, "complete.json")?)?;
    let manifest: Manifest = serde_json::from_slice(&manifest_bytes).map_err(|e| e.to_string())?;
    let complete: Complete = serde_json::from_slice(&complete_bytes).map_err(|e| e.to_string())?;
    if manifest.kind != "tp73_promoter_collaboration"
        || manifest.schema_version != 1
        || manifest.state != "complete"
        || manifest.complete_genome_scan
        || manifest.coordinate_mode != "bed_0based_half_open"
        || manifest.distance_bands != BANDS
        || manifest.panel.is_empty()
        || manifest.source_score_floor != -1.0
        || manifest.positive_threshold != 0.0
        || manifest.identity != complete.identity
        || !manifest.identity.contains_key("plan_sha256")
    {
        return Err("Unsupported collaborator package kind/version/coverage or inconsistent completion identity".into());
    }
    let mut files = BTreeMap::new();
    let mut total = 0u64;
    if complete.files.len() > 4096 {
        return Err("Package inventory exceeds 4096 files".into());
    }
    for entry in complete.files {
        total = total
            .checked_add(entry.bytes)
            .ok_or("Package size overflow")?;
        if total > 2_000_000_000
            || entry.sha256.len() != 64
            || !entry.sha256.bytes().all(|b| b.is_ascii_hexdigit())
        {
            return Err("Invalid file digest or package exceeds 2 GB contract".into());
        }
        if Instant::now() >= deadline {
            return Err("Package verification timed out".into());
        }
        let path = checked_path(&root, &entry.path)?;
        if fs::metadata(&path).map_err(|e| e.to_string())?.len() != entry.bytes
            || sha256_file_hex(&path).map_err(|e| e.to_string())? != entry.sha256
        {
            return Err(format!("Package integrity mismatch: {}", entry.path));
        }
        if files.insert(entry.path.clone(), entry).is_some() {
            return Err("Duplicate inventory path".into());
        }
    }
    if files.get("manifest.json").map(|f| f.sha256.as_str())
        != Some(sha256_hex_bytes(&manifest_bytes).as_str())
        || manifest.files.iter().any(|f| files.get(&f.path) != Some(f))
    {
        return Err("Completion inventory does not bind the manifest and its files".into());
    }
    for name in [
        "anchors.parquet",
        "cofactor_distance_isoform_comparison.parquet",
        "anchor_promoter.parquet",
        "promoter.parquet",
        "promoter_gene.parquet",
    ] {
        if !files.contains_key(name) {
            return Err(format!("Missing required inventory entry {name}"));
        }
    }
    let mut motifs = BTreeSet::new();
    for entry in &manifest.panel {
        if !entry
            .motif_id
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'.')
            || !motifs.insert(entry.motif_id.clone())
            || !files.contains_key(&format!("feature_{}.parquet", entry.motif_id))
        {
            return Err("Invalid, duplicate or missing detailed motif inventory".into());
        }
    }
    report.package_manifest_sha256 = Some(sha256_hex_bytes(&manifest_bytes));
    report.completion_sha256 = Some(sha256_hex_bytes(&complete_bytes));
    report.verified_file_sha256 = files
        .iter()
        .map(|(k, v)| (k.clone(), v.sha256.clone()))
        .collect();
    report.coverage = Some(CofactorPackageCoverage {
        assembly: manifest.assembly.clone(),
        taxon_id: manifest.taxon_id,
        chromosomes: manifest.chromosomes.clone(),
        detailed_motif_ids: motifs.into_iter().collect(),
        distance_bands: manifest.distance_bands.clone(),
        source_score_floor: manifest.source_score_floor,
        positive_threshold: manifest.positive_threshold,
        score_configuration: manifest.score_configuration.clone(),
        retention: manifest.retention.clone(),
        complete_genome_scan: false,
        h3k4me3_model_effects: manifest.h3k4me3_model_effects.clone(),
        requested_candidates: manifest.requested_genes.clone(),
    });
    Ok(Package {
        root,
        manifest,
        files,
    })
}

fn literal(s: &str) -> String {
    format!("'{}'", s.replace('\'', "''"))
}
impl Package {
    fn table(&self, name: &str) -> Result<String, String> {
        if !self.files.contains_key(name) {
            return Err(format!("File not inventoried: {name}"));
        }
        Ok(format!(
            "read_parquet({}, hive_partitioning=false)",
            literal(&checked_path(&self.root, name)?.to_string_lossy())
        ))
    }
}

pub fn validate_request(r: &PromoterCofactorRequest) -> Result<(), String> {
    if !matches!(
        r.query,
        CofactorQuery::Rankings | CofactorQuery::AnchorDetail
    ) && (r.motif.is_some() || r.distance_band.is_some())
    {
        return Err("Motif/band filters apply to rankings and anchor_detail; other queries enumerate the anchor cohort".into());
    }
    if matches!(r.query, CofactorQuery::Inspect | CofactorQuery::Candidates)
        && r.presence_threshold != 0.0
    {
        return Err("Inspection does not evaluate presence thresholds".into());
    }
    if !(1..=2000).contains(&r.max_rows)
        || !(1..=120).contains(&r.timeout_seconds)
        || r.assembly.trim().is_empty()
        || !r.presence_threshold.is_finite()
        || r.presence_threshold < -1.0
        || r.max_q_value
            .is_some_and(|v| !v.is_finite() || !(0.0..=1.0).contains(&v))
    {
        return Err("Require assembly, max_rows 1..2000, timeout_seconds 1..120, presence_threshold >= -1 and q in [0,1]".into());
    }
    if r.region.as_ref().is_some_and(|v| {
        v.chromosome.is_empty()
            || v.start_0based >= v.end_0based_exclusive
            || v.end_0based_exclusive > i64::MAX as u64
    }) {
        return Err("Region requires chromosome and BED 0 <= start < end <= i64::MAX".into());
    }
    if r.distance_band
        .as_deref()
        .is_some_and(|b| !BANDS.contains(&b))
    {
        return Err("Unknown exclusive distance band".into());
    }
    if r.query == CofactorQuery::AnchorDetail
        && (r.anchor_id.is_none() || r.motif.as_deref().is_none_or(str::is_empty))
    {
        return Err(
            "anchor_detail requires package-local anchor_id and exact motif accession".into(),
        );
    }
    if matches!(r.query, CofactorQuery::Anchors | CofactorQuery::Promoters)
        && r.region.is_none()
        && r.gene_id.is_none()
        && r.anchor_id.is_none()
    {
        return Err("Anchors/promoters require a region, gene_id or anchor_id".into());
    }
    if r.query != CofactorQuery::Rankings && (r.source_species.is_some() || r.max_q_value.is_some())
    {
        return Err("Species/q filters apply only to cohort rankings".into());
    }
    if matches!(
        r.query,
        CofactorQuery::Rankings | CofactorQuery::Inspect | CofactorQuery::Candidates
    ) && (r.region.is_some() || r.gene_id.is_some() || r.anchor_id.is_some())
    {
        return Err(
            "Cohort statistics cannot be reinterpreted as region/gene-local statistics".into(),
        );
    }
    Ok(())
}

fn decode<T: DeserializeOwned>(rows: Vec<Value>) -> Result<Vec<T>, String> {
    rows.into_iter()
        .map(|v| serde_json::from_value(v).map_err(|e| format!("Unsupported source row: {e}")))
        .collect()
}

fn run_queries(
    p: &Package,
    r: &PromoterCofactorRequest,
    report: &mut PromoterCofactorReport,
    executable: &str,
    deadline: Instant,
) -> Result<(), String> {
    let run = |sql: &str| -> Result<Vec<Value>, String> {
        let remaining = deadline
            .checked_duration_since(Instant::now())
            .ok_or("Query timed out")?;
        query_local_parquet_json(executable, Some(sql), remaining)
    };
    let limit = r.max_rows + 1;
    if r.query == CofactorQuery::Rankings {
        let (sort, q, status) = match r.ranking {
            CofactorRanking::TaEnriched => (
                "ta_adjusted_odds_ratio DESC",
                "ta_q_value_bh_tax_group",
                "ta_evaluation_status",
            ),
            CofactorRanking::DnEnriched => (
                "dn_adjusted_odds_ratio DESC",
                "dn_q_value_bh_tax_group",
                "dn_evaluation_status",
            ),
            CofactorRanking::TaDepleted => (
                "ta_adjusted_odds_ratio ASC",
                "ta_q_value_bh_tax_group",
                "ta_evaluation_status",
            ),
            CofactorRanking::DnDepleted => (
                "dn_adjusted_odds_ratio ASC",
                "dn_q_value_bh_tax_group",
                "dn_evaluation_status",
            ),
            CofactorRanking::IsoformDifference => (
                "abs(ta_vs_dn_log_odds_difference) DESC",
                "q_value_bh_tax_group",
                "evaluation_status",
            ),
        };
        let mut filters = vec!["true".to_string()];
        if let Some(m) = &r.motif {
            filters.push(format!(
                "(motif_id={} OR contains(lower(motif_name),lower({})))",
                literal(m),
                literal(m)
            ));
        }
        if let Some(b) = &r.distance_band {
            filters.push(format!("distance_band={}", literal(b)));
        }
        if let Some(s) = &r.source_species {
            filters.push(format!(
                "contains(lower(source_species),lower({}))",
                literal(s)
            ));
        }
        if let Some(v) = r.max_q_value {
            filters.push(format!("{q}<={v} AND {status}='ok'"));
        }
        // Retain non-estimable rows, sorted after estimable ones unless explicitly filtered.
        report.rankings = decode(run(&format!(
            "SELECT * FROM {} WHERE {} ORDER BY ({status}='ok') DESC, {sort} NULLS LAST, motif_id,distance_band_order LIMIT {limit}",
            p.table("cofactor_distance_isoform_comparison.parquet")?,
            filters.join(" AND ")
        ))?)?;
        return Ok(());
    }
    let anchors = p.table("anchors.parquet")?;
    let ap = p.table("anchor_promoter.parquet")?;
    let pg = p.table("promoter_gene.parquet")?;
    let mut predicates = vec!["true".to_string()];
    if let Some(v) = &r.region {
        predicates.push(format!(
            "a.chrom={} AND a.anchor_start<{} AND a.anchor_end>{}",
            literal(&v.chromosome),
            v.end_0based_exclusive,
            v.start_0based
        ));
    }
    if let Some(id) = r.anchor_id {
        predicates.push(format!("a.anchor_id={id}"));
    }
    if let Some(g) = &r.gene_id {
        predicates.push(format!("EXISTS(SELECT 1 FROM {ap} ap JOIN {pg} g USING(regulatory_feature_id) WHERE ap.anchor_id=a.anchor_id AND g.gene_id={})",literal(g)));
    }
    let selected = format!(
        "SELECT a.* FROM {anchors} a WHERE {}",
        predicates.join(" AND ")
    );
    report.anchors = decode(run(&format!(
        "{selected} ORDER BY chrom,anchor_start,anchor_end,anchor_id LIMIT {limit}"
    ))?)?;
    if report.anchors.len() > r.max_rows {
        return Ok(());
    }
    let mut anchor_ids = BTreeSet::new();
    for a in &report.anchors {
        if a.anchor_start >= a.anchor_end
            || a.anchor_end > i64::MAX as u64
            || !a.anchor_score.is_finite()
            || !p.manifest.chromosomes.contains(&a.chrom)
            || !anchor_ids.insert(a.anchor_id)
        {
            return Err("Invalid or duplicate anchor geometry/identity".into());
        }
    }
    if r.query == CofactorQuery::AnchorDetail {
        let motif = r.motif.as_deref().ok_or("Missing motif")?;
        let file = p.table(&format!("feature_{motif}.parquet"))?;
        let bands = BANDS
            .iter()
            .enumerate()
            .filter(|(_, b)| r.distance_band.as_deref().is_none_or(|f| f == **b))
            .map(|(i, b)| format!("({}, {i})", literal(b)))
            .collect::<Vec<_>>()
            .join(",");
        report.details=decode(run(&format!("WITH selected AS ({selected}) SELECT a.anchor_id,{} AS motif_id,b.distance_band,
            f.hit_start,f.hit_end,f.best_score,f.plus_score,f.minus_score,f.best_strand,f.interval_distance_bp,f.genomic_side,
            coalesce(f.n_source_loci,0)::BIGINT AS n_source_loci,coalesce(f.n_score_zero_loci,0)::BIGINT AS n_score_zero_loci,
            coalesce(f.best_score>={},false) AS present_at_requested_threshold
            FROM selected a CROSS JOIN (VALUES {bands}) b(distance_band,ordinal)
            LEFT JOIN {file} f ON f.anchor_id=a.anchor_id AND f.distance_band=b.distance_band AND f.motif_id={}
            ORDER BY a.anchor_id,b.ordinal LIMIT {limit}", literal(motif),r.presence_threshold,literal(motif)))?)?;
        for d in &report.details {
            let a = report
                .anchors
                .iter()
                .find(|a| a.anchor_id == d.anchor_id)
                .ok_or("Unknown detail anchor")?;
            validate_detail(a, d)?;
        }
        let mut keys = BTreeSet::new();
        if report
            .details
            .iter()
            .any(|d| !keys.insert((d.anchor_id, d.distance_band.clone())))
        {
            return Err("Duplicate strongest-per-band rows".into());
        }
        let band_filter = r
            .distance_band
            .as_ref()
            .map(|b| format!(" AND distance_band={}", literal(b)))
            .unwrap_or_default();
        report.rankings = decode(run(&format!(
            "SELECT * FROM {} WHERE motif_id={} {band_filter} ORDER BY distance_band_order LIMIT {limit}",
            p.table("cofactor_distance_isoform_comparison.parquet")?,
            literal(motif)
        ))?)?;
    }
    report.memberships = decode(run(&format!(
        "WITH selected AS ({selected}) SELECT ap.* FROM {ap} ap JOIN selected a USING(anchor_id) ORDER BY ap.anchor_id,ap.regulatory_feature_id LIMIT {limit}"
    ))?)?;
    let promoters = p.table("promoter.parquet")?;
    report.promoters=decode(run(&format!("WITH selected AS ({selected}) SELECT p.*,
        coalesce((SELECT list(struct_pack(gene_id:=g.gene_id,link_source:=g.link_source,annotation_release:=g.annotation_release) ORDER BY g.gene_id,g.link_source,g.annotation_release) FROM {pg} g WHERE g.regulatory_feature_id=p.regulatory_feature_id),[]) AS gene_links
        FROM {promoters} p WHERE EXISTS(SELECT 1 FROM {ap} ap JOIN selected a USING(anchor_id) WHERE ap.regulatory_feature_id=p.regulatory_feature_id)
        ORDER BY p.chrom,p.extended_start,p.extended_end,p.regulatory_feature_id LIMIT {limit}"))?)?;
    if report.promoters.iter().any(|v| {
        v.extended_start >= v.extended_end
            || v.extended_end > i64::MAX as u64
            || !p.manifest.chromosomes.contains(&v.chrom)
            || v.source_fields
                .get("assembly")
                .is_some_and(|a| a.as_str() != Some(p.manifest.assembly.as_str()))
    }) {
        return Err("Invalid promoter geometry/reference".into());
    }
    Ok(())
}

fn validate_detail(a: &CofactorAnchor, d: &CofactorDetail) -> Result<(), String> {
    if d.n_score_zero_loci > d.n_source_loci {
        return Err("Invalid occurrence counts".into());
    }
    if d.n_source_loci == 0 {
        if d.best_score.is_some()
            || d.hit_start.is_some()
            || d.hit_end.is_some()
            || d.plus_score.is_some()
            || d.minus_score.is_some()
            || d.best_strand.is_some()
            || d.interval_distance_bp.is_some()
            || d.present_at_requested_threshold
        {
            return Err("Zero row has a fabricated hit".into());
        }
        return Ok(());
    }
    let (start, end, score) = match (d.hit_start, d.hit_end, d.best_score) {
        (Some(s), Some(e), Some(v))
            if s < e && e <= i64::MAX as u64 && v.is_finite() && v >= -1.0 =>
        {
            (s, e, v)
        }
        _ => return Err("Invalid strongest-hit geometry or score".into()),
    };
    let gap = a.anchor_start.max(start) as i64 - a.anchor_end.min(end) as i64;
    if (score >= 0.0) != (d.n_score_zero_loci > 0) {
        return Err("Positive count disagrees with maximum".into());
    }
    let band = match gap {
        i64::MIN..=-1 => BANDS[0],
        0..=5 => BANDS[1],
        6..=20 => BANDS[2],
        21..=50 => BANDS[3],
        51..=100 => BANDS[4],
        101..=150 => BANDS[5],
        _ => return Err("Hit exceeds covered distance".into()),
    };
    if d.interval_distance_bp != Some(gap) || d.distance_band != band {
        return Err("Distance/band disagrees with hit coordinates".into());
    }
    let strand = match (d.plus_score, d.minus_score) {
        (Some(p), Some(m)) if p == m => ".",
        (Some(p), Some(m)) if p < m => "-",
        (Some(_), _) => "+",
        (_, Some(_)) => "-",
        _ => return Err("Missing strand scores".into()),
    };
    if d.best_strand.as_deref() != Some(strand)
        || d.plus_score
            .into_iter()
            .chain(d.minus_score)
            .any(|s| !s.is_finite() || s < -1.0 || s > score)
        || !d
            .plus_score
            .into_iter()
            .chain(d.minus_score)
            .any(|s| s == score)
    {
        return Err("Strand scores disagree with strongest hit".into());
    }
    Ok(())
}

/// Query immutable package files. No project state, downloads, or source writes.
pub fn query(r: &PromoterCofactorRequest) -> Result<PromoterCofactorReport, String> {
    validate_request(r)?;
    let mut report=PromoterCofactorReport {
        schema:PROMOTER_COFACTOR_SCHEMA.into(),report_id:String::new(),request:r.clone(),
        availability:CofactorAvailability::Available,diagnostic:None,package_manifest_sha256:None,
        completion_sha256:None,verified_file_sha256:BTreeMap::new(),duckdb_version:None,coverage:None,
        rankings:vec![],more_rankings_available:false,anchors:vec![],details:vec![],promoters:vec![],memberships:vec![],
        non_claims:vec![
            "Adjusted odds ratios and BH q-values describe cohort associations, not site significance or confirmed cofactor occupancy.".into(),
            "TA/DN is a ratio of adjusted odds ratios, not a binding-probability ratio.".into(),
            "Strongest-per-band retention is not a complete scan, all occurrences or pair architecture; missing outside covered scope is unavailable.".into(),
            "CUT&RUN depths are conditional motif-span maxima, not raw coverage tracks. SAOS2 and SK-MEL-29_2 are distinct series; SK-MEL-29_1 is excluded.".into(),
            "Genomic left/right is not transcriptional upstream/downstream. Matrix source species is provenance, not species exclusivity.".into(),
        ],
    };
    let deadline = Instant::now() + Duration::from_secs(r.timeout_seconds);
    let outcome = (|| -> Result<(), (CofactorAvailability, String)> {
        if !Path::new(&r.package_path).is_dir() {
            return Err((
                CofactorAvailability::PackageMissing,
                "Select the complete local package directory".into(),
            ));
        }
        let p = validate_package(r, &mut report, deadline)
            .map_err(|e| (CofactorAvailability::InvalidPackage, e))?;
        if p.manifest.assembly != r.assembly {
            return Err((
                CofactorAvailability::AssemblyMismatch,
                format!("Requested {}, package {}", r.assembly, p.manifest.assembly),
            ));
        }
        if r.region
            .as_ref()
            .is_some_and(|v| !p.manifest.chromosomes.contains(&v.chromosome))
        {
            return Err((
                CofactorAvailability::UnsupportedCoverage,
                "Chromosome is outside the included anchor cohort".into(),
            ));
        }
        if r.query == CofactorQuery::AnchorDetail
            && !p
                .manifest
                .panel
                .iter()
                .any(|m| Some(&m.motif_id) == r.motif.as_ref())
        {
            return Err((
                CofactorAvailability::UnsupportedCoverage,
                "Motif has no positional detail in this package; overview availability is separate"
                    .into(),
            ));
        }
        if matches!(r.query, CofactorQuery::Inspect | CofactorQuery::Candidates) {
            return Ok(());
        }
        let executable = r
            .duckdb_executable
            .clone()
            .or_else(|| std::env::var("GENTLE_DUCKDB_BIN").ok())
            .unwrap_or_else(|| "duckdb".into());
        report.request.duckdb_executable = Some(executable.clone());
        let remaining = deadline.checked_duration_since(Instant::now()).ok_or((
            CofactorAvailability::QueryFailed,
            "Verification exceeded timeout".into(),
        ))?;
        let version = query_local_parquet_json(&executable, None, remaining)
            .map_err(|e| (CofactorAvailability::RuntimeUnavailable, e))?;
        report.duckdb_version = version.first().and_then(Value::as_str).map(str::to_string);
        run_queries(&p, r, &mut report, &executable, deadline)
            .map_err(|e| (CofactorAvailability::QueryFailed, e))?;
        if r.query == CofactorQuery::Rankings && report.rankings.len() > r.max_rows {
            report.more_rankings_available = true;
            report.rankings.truncate(r.max_rows);
        }
        if report.rankings.len()
            + report.anchors.len()
            + report.details.len()
            + report.promoters.len()
            + report.memberships.len()
            > r.max_rows
        {
            return Err((
                CofactorAvailability::RowLimitExceeded,
                "Combined result exceeds max_rows; narrow the query or explicitly raise its limit"
                    .into(),
            ));
        }
        if r.query != CofactorQuery::Rankings && report.anchors.is_empty() {
            return Err((
                CofactorAvailability::UnsupportedCoverage,
                "No included anchor matches the requested ID, region or gene; cofactor absence is not established".into(),
            ));
        }
        // Detect replacement while querying, including replacements at identical paths.
        for (name, entry) in &p.files {
            if Instant::now() >= deadline {
                return Err((
                    CofactorAvailability::QueryFailed,
                    "Integrity recheck exceeded timeout".into(),
                ));
            }
            let path = checked_path(&p.root, name)
                .map_err(|e| (CofactorAvailability::InvalidPackage, e))?;
            if sha256_file_hex(&path)
                .map_err(|e| (CofactorAvailability::InvalidPackage, e.to_string()))?
                != entry.sha256
            {
                return Err((
                    CofactorAvailability::InvalidPackage,
                    format!("Package changed during query: {name}"),
                ));
            }
        }
        let completion = read_metadata(
            &checked_path(&p.root, "complete.json")
                .map_err(|e| (CofactorAvailability::InvalidPackage, e))?,
        )
        .map_err(|e| (CofactorAvailability::InvalidPackage, e))?;
        if report.completion_sha256.as_deref() != Some(sha256_hex_bytes(&completion).as_str()) {
            return Err((
                CofactorAvailability::InvalidPackage,
                "Completion marker changed during query".into(),
            ));
        }
        Ok(())
    })();
    if let Err((status, message)) = outcome {
        report.availability = status;
        report.diagnostic = Some(message);
        report.rankings.clear();
        report.anchors.clear();
        report.details.clear();
        report.promoters.clear();
        report.memberships.clear();
        report.more_rankings_available = false;
    }
    // Location is transport metadata, not content identity.
    let mut identity_request = r.clone();
    identity_request.package_path.clear();
    identity_request.duckdb_executable = None;
    report.report_id = short_sha256_id(
        "promoter_cofactors",
        &serde_json::to_string(&(
            identity_request,
            &report.package_manifest_sha256,
            &report.completion_sha256,
            &report.duckdb_version,
            report.availability,
            &report.rankings,
            &report.anchors,
            &report.details,
            &report.promoters,
            &report.memberships,
        ))
        .map_err(|e| e.to_string())?,
    );
    Ok(report)
}

#[cfg(test)]
#[path = "promoter_cofactors_tests.rs"]
mod tests;
