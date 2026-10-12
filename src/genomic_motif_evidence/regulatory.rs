//! Read-only regulatory/TSS-subset adapter for the shared genomic motif query.
//! Package inspection never opens hit payloads. Regional queries verify only
//! selected payloads and annotation dimensions, under one wall-clock budget.

use super::*;
use gentle_protocol::genomic_motif_evidence::*;

mod package;
use package::{Package, checked_file, verify};

const MAX_OWNERS: usize = 20_000;
const MAX_WINDOWS: usize = 256;
const MAX_SELECTED_BYTES: u64 = 1024 * 1024 * 1024;

pub(super) fn validate_target(target: &GenomicMotifEvidenceTarget) -> Result<(), String> {
    match target {
        GenomicMotifEvidenceTarget::PackageCatalog {
            search,
            offset,
            limit,
        } => {
            if search.len() > 256 || *offset > 250_000 || *limit == 0 || *limit > 1000 {
                return Err(
                    "Catalog requires search <=256 bytes, offset <=250000, limit 1..1000".into(),
                );
            }
        }
        GenomicMotifEvidenceTarget::PackageTssWindows { gene_query, tss_id } => {
            if usize::from(gene_query.is_some()) + usize::from(tss_id.is_some()) != 1
                || gene_query
                    .iter()
                    .chain(tss_id)
                    .any(|s| s.trim().is_empty() || s.len() > 256)
            {
                return Err(
                    "Specify exactly one nonempty package gene or TSS ID (<=256 bytes)".into(),
                );
            }
        }
        _ => {}
    }
    Ok(())
}

struct Runtime {
    executable: String,
    start: Instant,
    timeout: Duration,
}

impl Runtime {
    fn remaining(&self) -> Result<Duration, String> {
        self.timeout
            .checked_sub(self.start.elapsed())
            .filter(|d| !d.is_zero())
            .ok_or_else(|| {
                "query_failed:Regulatory motif query exceeded its total wall-clock budget".into()
            })
    }

    fn sql(&self, database: &Path, sql: &str) -> Result<Vec<Value>, String> {
        let mut command = Command::new(&self.executable);
        command.args(["-no-init", "-batch", "-bail", "-json", ":memory:", "-c"]);
        command.arg(format!(
            "SET threads=2; SET memory_limit='512MB'; SET max_temp_directory_size='0B'; \
             SET autoinstall_known_extensions=false; SET autoload_known_extensions=false; \
             ATTACH {} AS evidence (READ_ONLY); {sql}",
            sql_string(&database.to_string_lossy())
        ));
        let output =
            bounded_command_output(&mut command, self.remaining()?, MAX_DUCKDB_OUTPUT_BYTES)
                .map_err(|e| format!("query_failed:{e}"))?;
        if output.stdout_truncated || output.stderr_truncated {
            return Err("query_failed:DuckDB output exceeded 64 MiB; narrow the query".into());
        }
        if !output.output.status.success() {
            return Err(format!(
                "query_failed:DuckDB: {}",
                String::from_utf8_lossy(&output.output.stderr)
            ));
        }
        if output.output.stdout.iter().all(u8::is_ascii_whitespace) {
            return Ok(vec![]);
        }
        serde_json::from_slice(&output.output.stdout)
            .map_err(|e| format!("query_failed:Invalid DuckDB JSON: {e}"))
    }
}

fn number(row: &Value, key: &str) -> Result<u64, String> {
    optional_u64(row, key).ok_or_else(|| format!("Missing nonnegative integer {key}"))
}

fn finite(row: &Value, key: &str) -> Result<f64, String> {
    optional_f64(row, key)
        .filter(|v| v.is_finite())
        .ok_or_else(|| format!("Missing finite {key}"))
}

fn values(strings: impl IntoIterator<Item = impl AsRef<str>>) -> String {
    strings
        .into_iter()
        .map(|s| sql_string(s.as_ref()))
        .collect::<Vec<_>>()
        .join(",")
}

fn parquet(path: &Path) -> String {
    format!(
        "read_parquet({}, hive_partitioning=false)",
        sql_string(&path.to_string_lossy())
    )
}

fn region(id: String, chrom: String, start: u64, end: u64) -> GenomicMotifQueryRegion {
    GenomicMotifQueryRegion {
        interval_id: id,
        chromosome: chrom,
        start_0based: start,
        end_0based_exclusive: end,
        source_reference: None,
        label: None,
        source_seq_id: None,
        source_start_0based: None,
        source_end_0based_exclusive: None,
        source_sequence_length_bp: None,
        source_anchor_start_1based: None,
        source_anchor_end_1based: None,
        source_anchor_reverse: false,
    }
}

pub(super) fn query(
    request: &GenomicMotifEvidenceRequest,
    regions: &[GenomicMotifQueryRegion],
    paths: &PackagePaths,
) -> GenomicMotifEvidenceReport {
    match query_inner(request, regions, paths) {
        Ok(report) => report,
        Err(error) => {
            let (availability, message) = if let Some(e) = error.strip_prefix("duckdb_unavailable:")
            {
                (GenomicMotifEvidenceAvailability::DuckdbUnavailable, e)
            } else if let Some(e) = error.strip_prefix("incompatible_package:") {
                (GenomicMotifEvidenceAvailability::IncompatiblePackage, e)
            } else if let Some(e) = error.strip_prefix("query_failed:") {
                (GenomicMotifEvidenceAvailability::QueryFailed, e)
            } else {
                (
                    GenomicMotifEvidenceAvailability::InvalidPackage,
                    error.as_str(),
                )
            };
            unavailable_report(
                request,
                regions,
                availability,
                format!("{message}; local GENtle TFBS scoring remains available"),
            )
        }
    }
}

fn query_inner(
    request: &GenomicMotifEvidenceRequest,
    input_regions: &[GenomicMotifQueryRegion],
    paths: &PackagePaths,
) -> Result<GenomicMotifEvidenceReport, String> {
    let runtime = Runtime {
        executable: request
            .duckdb_executable
            .clone()
            .filter(|s| !s.is_empty())
            .or_else(|| env::var(DUCKDB_BIN_ENV).ok())
            .unwrap_or_else(|| "duckdb".into()),
        start: Instant::now(),
        timeout: Duration::from_secs(request.timeout_seconds),
    };
    let mut package = Package::open(paths, &runtime)?;
    let m = &package.manifest;
    if request
        .expected_genome_id
        .as_ref()
        .is_some_and(|g| g != &m.genome_id)
    {
        return Err("incompatible_package:Exact requested genome ID differs from package".into());
    }
    let version = run_duckdb(&runtime.executable, None, None, runtime.remaining()?)
        .map_err(|e| format!("duckdb_unavailable:{e}"))?;
    if !version.output.status.success() || version.stdout_truncated {
        return Err("duckdb_unavailable:DuckDB version probe failed".into());
    }
    let database = &paths.database_path;
    let tables = runtime.sql(database, "SELECT table_name FROM duckdb_tables() WHERE database_name='evidence' AND schema_name='main';")?;
    for name in [
        "genome",
        "motif_metadata",
        "sequence_region",
        "scan_file_inventory",
        "file_inventory",
    ] {
        if tables.iter().filter(|r| r["table_name"] == name).count() != 1 {
            return Err(format!(
                "Expected materialized table {name}; views/macros are not accepted"
            ));
        }
    }
    let genomes = runtime.sql(database, "SELECT * FROM evidence.main.genome LIMIT 2;")?;
    if genomes.len() != 1
        || genomes[0]["genome_id"] != m.genome_id
        || genomes[0]["assembly_name"] != m.assembly
    {
        return Err("Genome catalog and manifest disagree".into());
    }
    let genome = &genomes[0];
    let counts = runtime.sql(database, "SELECT count(*) AS n, count(DISTINCT motif_id) AS unique_n FROM evidence.main.motif_metadata;")?;
    if counts.len() != 1
        || number(&counts[0], "n")? != m.motif_count as u64
        || number(&counts[0], "unique_n")? != m.motif_count as u64
    {
        return Err("Motif catalog is incomplete or duplicated".into());
    }
    let ids = runtime
        .sql(
            database,
            "SELECT motif_id FROM evidence.main.motif_metadata ORDER BY motif_id;",
        )?
        .iter()
        .map(|r| required_string(r, "motif_id"))
        .collect::<Result<BTreeSet<_>, _>>()?;
    if ids.iter().collect::<BTreeSet<_>>()
        != package.entries.keys().map(|(_, motif)| motif).collect()
    {
        return Err("Motif catalog and delivered inventory contain different accessions".into());
    }
    let mut subset = RegulatoryMotifSubset {
        scope: m.scope.clone(),
        score_selection: m.score_selection.clone(),
        complete_genome_scan: false,
        requires_tp73: false,
        annotation_release: m.annotation_release.clone(),
        regulatory_release: m.regulatory_release.clone(),
        promoter_definition_id: m.promoter_definition_id.clone(),
        upstream_bp: m.upstream_bp,
        downstream_bp: m.downstream_bp,
        production_plan_sha256: m.production_plan_sha256.clone(),
        source_commit: m.source_commit.clone(),
        catalog_total_motifs: m.motif_count,
        ..Default::default()
    };
    let mut provider = GenomicMotifEvidenceProviderProvenance {
        provider_kind: REGULATORY_MOTIF_PROVIDER.into(),
        package_root: paths.root.display().to_string(),
        manifest_path: paths.manifest_path.display().to_string(),
        manifest_sha256: sha256_prefixed_bytes(&paths.manifest_bytes),
        declared_content_fingerprint_sha256: sha256_prefixed_bytes(&paths.manifest_bytes),
        manifest_schema_version: 1,
        database_path: database.display().to_string(),
        duckdb_executable: runtime.executable.clone(),
        duckdb_version: Some(
            String::from_utf8_lossy(&version.output.stdout)
                .trim()
                .into(),
        ),
        genome_id: m.genome_id.clone(),
        assembly_name: Some(m.assembly.clone()),
        assembly_accession: optional_string(genome, "assembly_accession"),
        ensembl_release: optional_string(genome, "ensembl_release"),
        coordinate_mode: m.coordinate_mode.clone(),
        ..Default::default()
    };
    let mut report = GenomicMotifEvidenceReport {
        schema: GENOMIC_MOTIF_EVIDENCE_SCHEMA.into(), request: request.clone(), availability: GenomicMotifEvidenceAvailability::Available,
        warnings: vec![subset.summary(), "Manifest/hash agreement is package consistency, not independent reference authentication. All source taxa are retained; matrix names are not TSS-owning genes.".into()],
        ..Default::default()
    };
    if let GenomicMotifEvidenceTarget::PackageCatalog {
        search,
        offset,
        limit,
    } = &request.target
    {
        let filter = format!(
            "contains(lower(motif_id),lower({0})) OR contains(lower(coalesce(motif_name,'')),lower({0}))",
            sql_string(search)
        );
        let count = runtime.sql(
            database,
            &format!("SELECT count(*) AS n FROM evidence.main.motif_metadata WHERE {filter};"),
        )?;
        subset.catalog_matched_motifs =
            number(count.first().ok_or("Missing catalog count")?, "n")? as usize;
        subset.catalog_offset = *offset;
        subset.catalog_has_more = offset.saturating_add(*limit) < subset.catalog_matched_motifs;
        subset.motif_metadata = runtime.sql(database, &format!("SELECT * FROM evidence.main.motif_metadata WHERE {filter} ORDER BY motif_id LIMIT {limit} OFFSET {offset};"))?;
        report.query_complete = true;
        report.warnings.push("Catalog inspection only; no hit payloads were opened. Score policy and genomic absence were not assessed.".into());
    } else {
        let mut regions = input_regions.to_vec();
        let exact_selection = match &request.target {
            GenomicMotifEvidenceTarget::PackageTssWindows { gene_query, tss_id } => {
                let p = package.annotation_file(paths, "promoter.parquet", &runtime)?;
                let filter = if let Some(gene) = gene_query {
                    let owners =
                        package.annotation_file(paths, "transcript_tss.parquet", &runtime)?;
                    format!(
                        "EXISTS(SELECT 1 FROM {} o WHERE o.tss_id=p.tss_id AND (o.gene_id={g} OR o.gene_name={g}))",
                        parquet(&owners),
                        g = sql_string(gene)
                    )
                } else {
                    format!(
                        "p.tss_id={}",
                        sql_string(tss_id.as_deref().ok_or("Missing TSS")?)
                    )
                };
                let windows = runtime.sql(database, &format!("SELECT p.* FROM {} p WHERE {filter} ORDER BY chrom,promoter_start,promoter_id LIMIT {};", parquet(&p), MAX_WINDOWS+1))?;
                if windows.is_empty() {
                    return Err("incompatible_package:No exact annotated gene/TSS match; not evidence of absent motifs".into());
                }
                subset.tss_windows = parse_windows(&windows, &package)?;
                regions = subset
                    .tss_windows
                    .iter()
                    .map(|w| {
                        region(
                            w.promoter_id.clone(),
                            w.chromosome.clone(),
                            w.start_0based,
                            w.end_0based_exclusive,
                        )
                    })
                    .collect();
                true
            }
            _ => false,
        };
        let contigs = runtime.sql(
            database,
            "SELECT * FROM evidence.main.sequence_region WHERE included_in_scan ORDER BY chrom;",
        )?;
        let mut compatible = vec![];
        let mut all_compatible = true;
        for r in &regions {
            let c = package
                .annotation
                .chromosomes
                .iter()
                .find(|c| c.chrom == r.chromosome);
            let status = if let Some(reference) = &r.source_reference {
                check_assembly_reference(
                    reference,
                    provider.assembly_name.as_deref(),
                    provider.assembly_accession.as_deref(),
                )
                .map(|()| GenomicMotifEvidenceCompatibilityStatus::AssemblyAndContigGeometryMatched)
                .unwrap_or_else(|e| e)
            } else {
                GenomicMotifEvidenceCompatibilityStatus::ContigGeometryMatchedOnly
            };
            if matches!(
                status,
                GenomicMotifEvidenceCompatibilityStatus::AssemblyMismatch
                    | GenomicMotifEvidenceCompatibilityStatus::AssemblyNotVerified
            ) {
                return Err("incompatible_package:Saved-region assembly disagrees with the subset; no payloads read".into());
            }
            let valid = c.is_some_and(|c| {
                r.start_0based < r.end_0based_exclusive && r.end_0based_exclusive <= c.length
            });
            if let Some(c) = c {
                if contigs
                    .iter()
                    .filter(|row| {
                        row["chrom"] == c.chrom && row["length"].as_u64() == Some(c.length)
                    })
                    .count()
                    != 1
                {
                    return Err("Sequence-region and annotation lengths disagree".into());
                }
            }
            report.regions.push(GenomicMotifEvidenceResolvedRegion {
                interval_id: r.interval_id.clone(),
                source_reference: r.source_reference.clone(),
                label: r.label.clone(),
                requested_chromosome: r.chromosome.clone(),
                resolved_chromosome: c.map(|c| c.chrom.clone()),
                start_0based: r.start_0based,
                end_0based_exclusive: r.end_0based_exclusive,
                source_seq_id: r.source_seq_id.clone(),
                source_start_0based: r.source_start_0based,
                source_end_0based_exclusive: r.source_end_0based_exclusive,
                package_contig_length_bp: c.map(|c| c.length),
                compatibility_status: if valid {
                    status
                } else {
                    GenomicMotifEvidenceCompatibilityStatus::ContigGeometryMismatch
                },
            });
            if valid {
                compatible.push(r.clone());
            } else {
                all_compatible = false;
            }
        }
        if !exact_selection {
            subset.tss_windows = region_windows(&mut package, paths, &runtime, &compatible)?;
        }
        attach_owners(&mut subset, &mut package, paths, &runtime)?;
        attach_regulatory_context(&mut subset, &mut package, paths, &runtime, &compatible)?;
        let motifs = values(&request.motif_ids);
        subset.motif_metadata = runtime.sql(database, &format!("SELECT * FROM evidence.main.motif_metadata WHERE motif_id IN ({motifs}) ORDER BY motif_id;"))?;
        subset.catalog_matched_motifs = subset.motif_metadata.len();
        let chromosomes = regions
            .iter()
            .map(|r| r.chromosome.clone())
            .collect::<BTreeSet<_>>();
        let mut selected = vec![];
        for chrom in &chromosomes {
            for motif in &request.motif_ids {
                let entry = package.entries.get(&(chrom.clone(), motif.clone()));
                let state = match entry {
                    Some(e) if e.state == "known_empty_intersection" => {
                        RegulatoryMotifCoverageState::KnownEmptyIntersection
                    }
                    Some(e) => {
                        if compatible.iter().any(|r| &r.chromosome == chrom) {
                            selected.push(e.clone());
                        }
                        RegulatoryMotifCoverageState::Available
                    }
                    None if !package.manifest.chromosomes.contains(chrom) => {
                        RegulatoryMotifCoverageState::ChromosomeNotInPackage
                    }
                    None => RegulatoryMotifCoverageState::MotifNotInPackage,
                };
                if state != RegulatoryMotifCoverageState::Available {
                    report.warnings.push(format!(
                        "{chrom}/{motif}: {state:?}; not evidence of no genomic motif matches"
                    ));
                }
                subset.coverage.push(RegulatoryMotifCoverage {
                    chromosome: chrom.clone(),
                    motif_id: motif.clone(),
                    state,
                });
            }
        }
        if selected.len() > request.max_payload_files {
            return Err("query_failed:Selected payload count exceeds max_payload_files".into());
        }
        let chrom_sql = values(&chromosomes);
        let original = runtime.sql(database, &format!("SELECT * FROM evidence.main.scan_file_inventory WHERE motif_id IN ({motifs}) AND chrom IN ({chrom_sql}) ORDER BY chrom,motif_id,strand LIMIT {};", MAX_GENOMIC_MOTIF_EVIDENCE_QUERY_MOTIFS * MAX_REGIONS * 2 + 1))?;
        let thresholds = source_policy(&original, &mut subset, &mut provider)?;
        if !chromosomes.is_empty() {
            let inventory = runtime.sql(database, &format!("SELECT * FROM evidence.main.file_inventory WHERE motif_id IN ({motifs}) AND chrom IN ({chrom_sql}) ORDER BY chrom,motif_id;"))?;
            for e in package.entries.values().filter(|e| {
                chromosomes.contains(&e.chrom) && request.motif_ids.contains(&e.motif_id)
            }) {
                let rows = inventory
                    .iter()
                    .filter(|r| r["chrom"] == e.chrom && r["motif_id"] == e.motif_id)
                    .collect::<Vec<_>>();
                if rows.len() != 1
                    || rows[0]["path"].as_str() != e.path.as_deref()
                    || rows[0]["sha256"].as_str() != e.sha256.as_deref()
                    || rows[0]["rows"].as_u64() != Some(e.rows)
                    || rows[0]["bytes"].as_u64() != Some(e.bytes)
                    || rows[0]["state"] != e.state
                {
                    return Err("Materialized inventory disagrees with file_inventory.json".into());
                }
            }
        }
        if !selected.is_empty() {
            let mut payloads = vec![];
            let mut bytes = 0_u64;
            for e in &selected {
                bytes = bytes.checked_add(e.bytes).ok_or("Payload byte overflow")?;
                if bytes > MAX_SELECTED_BYTES {
                    return Err(
                        "query_failed:Selected payloads exceed 1 GiB; narrow the query".into(),
                    );
                }
                let record = RegulatoryMotifFileBinding {
                    path: e.path.clone().ok_or("Missing payload")?,
                    bytes: e.bytes,
                    sha256: e.sha256.clone().ok_or("Missing payload hash")?,
                };
                payloads.push(verify(&paths.root, &record, &runtime, MAX_SELECTED_BYTES)?);
                package.verified.push(record);
                provider
                    .selected_payloads
                    .push(GenomicMotifEvidencePayloadProvenance {
                        task_id: String::new(),
                        output_relative_path: e.path.clone().unwrap(),
                        declared_sha256: e.sha256.clone().unwrap(),
                        emitted_hits: e.rows,
                    });
            }
            query_hits(
                &runtime,
                paths,
                request,
                &compatible,
                &payloads,
                &thresholds,
                &mut provider,
                &mut report,
                &mut subset,
            )?;
        }
        report.selected_payload_file_count = selected.len();
        report.selected_payload_emitted_hit_count = selected.iter().map(|e| e.rows).sum();
        report.motif_coverage =
            motif_coverage_rows(request, &thresholds, None, &report.hits, report.truncated);
        if compatible.is_empty() {
            for coverage in &mut report.motif_coverage {
                if let Some(m) = subset
                    .motif_metadata
                    .iter()
                    .find(|m| m["motif_id"] == coverage.motif_id)
                {
                    coverage.status = GenomicMotifEvidenceCoverageStatus::NotAssessed;
                    coverage.motif_name = optional_string(m, "motif_name");
                }
            }
            report.warnings.push(
                "No compatible genomic interval was queried; catalog presence is not hit coverage."
                    .into(),
            );
        }
        report.query_complete =
            all_compatible && motif_coverage_is_query_complete(&report.motif_coverage);
        report
            .warnings
            .extend(motif_coverage_warnings(&report.motif_coverage));
        if report.truncated {
            report.warnings.push(
                "TRUNCATED at max_rows; narrow the query before interpreting completeness".into(),
            );
        }
    }
    runtime.remaining()?;
    if subset.score_pseudocount.is_none() {
        report.warnings.push("Scoring pseudocount not available from the inspected rows; regulatory_subset.score_pseudocount=null is unknown, not a zero-pseudocount claim.".into());
    }
    // No reusable verification cache: replacement at the same path is rechecked.
    for record in &package.verified {
        let path = checked_file(&paths.root, &record.path)?;
        if fs::metadata(path).map_err(|e| e.to_string())?.len() != record.bytes {
            return Err("Package changed during query".into());
        }
    }
    subset.verified_files = package.verified;
    report.returned_hit_count = report.hits.len();
    report.matched_hit_count = report.hits.len();
    report.report_id = report_id(request, input_regions, Some(&provider.manifest_sha256));
    report.provider = Some(provider);
    report.regulatory_subset = Some(subset);
    validate_regulatory_subset(&report)?;
    Ok(report)
}

fn parse_windows(
    rows: &[Value],
    package: &Package,
) -> Result<Vec<RegulatoryMotifTssWindow>, String> {
    if rows.len() > MAX_WINDOWS {
        return Err(
            "query_failed:More than 256 physical TSS windows; select a TSS ID or narrower region"
                .into(),
        );
    }
    let mut ids = BTreeSet::new();
    rows.iter()
        .map(|r| {
            let w = RegulatoryMotifTssWindow {
                promoter_id: required_string(r, "promoter_id")?,
                tss_id: required_string(r, "tss_id")?,
                chromosome: required_string(r, "chrom")?,
                start_0based: number(r, "promoter_start")?,
                end_0based_exclusive: number(r, "promoter_end")?,
                tss_0based: number(r, "tss_start")?,
                strand: required_string(r, "strand")?,
            };
            if w.start_0based > w.tss_0based
                || w.tss_0based >= w.end_0based_exclusive
                || !matches!(w.strand.as_str(), "+" | "-")
                || !ids.insert(w.promoter_id.clone())
                || r["genome_id"] != package.manifest.genome_id
                || r["annotation_release"] != package.manifest.annotation_release
                || r["promoter_definition_id"] != package.manifest.promoter_definition_id
                || package
                    .annotation
                    .chromosomes
                    .iter()
                    .find(|c| c.chrom == w.chromosome)
                    .is_none_or(|c| w.end_0based_exclusive > c.length)
            {
                return Err("Invalid physical TSS window or source identity".into());
            }
            Ok(w)
        })
        .collect()
}

fn region_windows(
    package: &mut Package,
    paths: &PackagePaths,
    runtime: &Runtime,
    regions: &[GenomicMotifQueryRegion],
) -> Result<Vec<RegulatoryMotifTssWindow>, String> {
    let mut rows = vec![];
    let chroms = regions
        .iter()
        .map(|r| r.chromosome.clone())
        .collect::<BTreeSet<_>>();
    for chrom in chroms {
        let c = package
            .annotation
            .chromosomes
            .iter()
            .find(|c| c.chrom == chrom)
            .ok_or("Unknown chromosome")?;
        let name = c.promoters.path.clone();
        let p = package.annotation_file(paths, &name, runtime)?;
        let filter = regions
            .iter()
            .filter(|r| r.chromosome == chrom)
            .map(|r| {
                format!(
                    "(promoter_start<{} AND promoter_end>{})",
                    r.end_0based_exclusive, r.start_0based
                )
            })
            .collect::<Vec<_>>()
            .join(" OR ");
        rows.extend(runtime.sql(
            &paths.database_path,
            &format!(
                "SELECT * FROM {} WHERE {filter} ORDER BY promoter_start,promoter_id LIMIT {};",
                parquet(&p),
                MAX_WINDOWS + 1
            ),
        )?);
    }
    parse_windows(&rows, package)
}

fn attach_owners(
    subset: &mut RegulatoryMotifSubset,
    package: &mut Package,
    paths: &PackagePaths,
    runtime: &Runtime,
) -> Result<(), String> {
    if subset.tss_windows.is_empty() {
        return Ok(());
    }
    let owner_file = package.annotation_file(paths, "transcript_tss.parquet", runtime)?;
    let ids = values(subset.tss_windows.iter().map(|w| &w.tss_id));
    let rows=runtime.sql(&paths.database_path,&format!("SELECT * FROM {} WHERE tss_id IN ({ids}) ORDER BY tss_id,gene_id,transcript_id LIMIT {};",parquet(&owner_file),MAX_OWNERS+1))?;
    if rows.len() > MAX_OWNERS {
        return Err(
            "query_failed:Transcript ownership exceeds 20000 rows; narrow the query".into(),
        );
    }
    let mut seen = BTreeSet::new();
    for r in rows {
        let owner = RegulatoryMotifTranscriptOwner {
            tss_id: required_string(&r, "tss_id")?,
            transcript_id: required_string(&r, "transcript_id")?,
            gene_id: required_string(&r, "gene_id")?,
            gene_name: optional_string(&r, "gene_name"),
        };
        if r["genome_id"] != package.manifest.genome_id
            || r["annotation_release"] != package.manifest.annotation_release
            || !seen.insert((
                owner.tss_id.clone(),
                owner.transcript_id.clone(),
                owner.gene_id.clone(),
            ))
        {
            return Err("Invalid/duplicate transcript ownership".into());
        }
        subset.transcript_owners.push(owner);
    }
    if subset.tss_windows.iter().any(|w| {
        !subset
            .transcript_owners
            .iter()
            .any(|o| o.tss_id == w.tss_id)
    }) {
        return Err("Physical TSS lacks transcript ownership".into());
    }
    Ok(())
}

fn attach_regulatory_context(
    subset: &mut RegulatoryMotifSubset,
    package: &mut Package,
    paths: &PackagePaths,
    runtime: &Runtime,
    regions: &[GenomicMotifQueryRegion],
) -> Result<(), String> {
    for chrom in regions
        .iter()
        .map(|r| &r.chromosome)
        .collect::<BTreeSet<_>>()
    {
        let c = package
            .annotation
            .chromosomes
            .iter()
            .find(|c| &c.chrom == chrom)
            .ok_or("Missing annotation chromosome")?;
        let name = c.features.path.clone();
        let file = package.annotation_file(paths, &name, runtime)?;
        let filter = regions
            .iter()
            .filter(|r| &r.chromosome == chrom)
            .map(|r| {
                format!(
                    "(coalesce(extended_start,start)<{} AND coalesce(extended_end,\"end\")>{})",
                    r.end_0based_exclusive, r.start_0based
                )
            })
            .collect::<Vec<_>>()
            .join(" OR ");
        let rows = runtime.sql(&paths.database_path, &format!("SELECT * FROM {} WHERE {filter} ORDER BY start,\"end\",regulatory_feature_id LIMIT {};",parquet(&file),MAX_OWNERS+1))?;
        if subset.regulatory_features.len() + rows.len() > MAX_OWNERS {
            return Err(
                "query_failed:Regulatory annotations exceed 20000 rows; narrow the query".into(),
            );
        }
        for row in &rows {
            required_string(row, "regulatory_feature_id")?;
            if row["assembly"] != package.manifest.assembly
                || row["chrom"] != *chrom
                || number(row, "start")? >= number(row, "end")?
            {
                return Err("Invalid regulatory annotation identity/geometry".into());
            }
        }
        subset.regulatory_features.extend(rows);
    }
    if !subset.regulatory_features.is_empty() {
        let ids = values(
            subset
                .regulatory_features
                .iter()
                .map(|r| r["regulatory_feature_id"].as_str().unwrap()),
        );
        let file = package.annotation_file(paths, "regulatory_feature_gene.parquet", runtime)?;
        subset.regulatory_gene_links = runtime.sql(&paths.database_path,&format!("SELECT * FROM {} WHERE regulatory_feature_id IN ({ids}) ORDER BY regulatory_feature_id,gene_id LIMIT {};",parquet(&file),MAX_OWNERS+1))?;
        if subset.regulatory_gene_links.len() > MAX_OWNERS {
            return Err("query_failed:Regulatory gene links exceed 20000 rows".into());
        }
    }
    Ok(())
}

fn source_policy(
    rows: &[Value],
    subset: &mut RegulatoryMotifSubset,
    provider: &mut GenomicMotifEvidenceProviderProvenance,
) -> Result<BTreeMap<String, MotifThresholdRow>, String> {
    let mut thresholds = BTreeMap::new();
    for coverage in &subset.coverage {
        if !matches!(
            coverage.state,
            RegulatoryMotifCoverageState::Available
                | RegulatoryMotifCoverageState::KnownEmptyIntersection
        ) {
            continue;
        }
        let selected = rows
            .iter()
            .filter(|r| r["chrom"] == coverage.chromosome && r["motif_id"] == coverage.motif_id)
            .collect::<Vec<_>>();
        if selected.len() != 2
            || !["+", "-"]
                .iter()
                .all(|s| selected.iter().filter(|r| r["strand"] == *s).count() == 1)
        {
            return Err("Incomplete original scan policy for motif/chromosome/strands".into());
        }
        let floor = finite(selected[0], "minimum_score")?;
        for r in selected {
            if r["state"] != "complete"
                || r["coordinate_mode"] != "bed"
                || finite(r, "minimum_score")? != floor
                || !r["minimum_pwm_relative_score"].is_null()
                || !r["maximum_pwm_relative_score"].is_null()
            {
                return Err("Inconsistent/censored source scan policy".into());
            }
            merge_policy(provider, &mut subset.score_pseudocount, r)?;
        }
        if let Some(old) = thresholds.get(&coverage.motif_id) {
            let old: &MotifThresholdRow = old;
            if old.final_minimum_score != Some(floor) {
                return Err("Different chromosome source floors require separate queries".into());
            }
        }
        thresholds.insert(
            coverage.motif_id.clone(),
            MotifThresholdRow {
                motif_id: coverage.motif_id.clone(),
                motif_name: subset
                    .motif_metadata
                    .iter()
                    .find(|r| r["motif_id"] == coverage.motif_id)
                    .and_then(|r| optional_string(r, "motif_name")),
                threshold_set_id: None,
                informative_threshold: None,
                final_minimum_score: Some(floor),
                density_limited: None,
            },
        );
    }
    Ok(thresholds)
}

fn merge_policy(
    provider: &mut GenomicMotifEvidenceProviderProvenance,
    pseudocount: &mut Option<f64>,
    row: &Value,
) -> Result<(), String> {
    for (key, target) in [
        ("run_id", &mut provider.run_id),
        ("motif_set_id", &mut provider.motif_set_id),
        ("score_mode", &mut provider.score_mode),
        ("pseudocount_scheme", &mut provider.pseudocount_scheme),
        ("background_model_id", &mut provider.background_model_id),
        ("n_policy", &mut provider.n_policy),
    ] {
        if let Some(value) = optional_string(row, key) {
            if !target.is_empty() && *target != value {
                return Err(format!("Mixed {key} within one evidence query"));
            }
            *target = value;
        }
    }
    if let Some(p) = optional_f64(row, "pseudocount") {
        if !p.is_finite() || p < 0.0 || pseudocount.is_some_and(|old| old != p) {
            return Err("Mixed/invalid pseudocount".into());
        }
        provider.pseudocount = p;
        *pseudocount = Some(p);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn query_hits(
    runtime: &Runtime,
    paths: &PackagePaths,
    request: &GenomicMotifEvidenceRequest,
    regions: &[GenomicMotifQueryRegion],
    payloads: &[PathBuf],
    thresholds: &BTreeMap<String, MotifThresholdRow>,
    provider: &mut GenomicMotifEvidenceProviderProvenance,
    report: &mut GenomicMotifEvidenceReport,
    subset: &mut RegulatoryMotifSubset,
) -> Result<(), String> {
    let filter = regions
        .iter()
        .map(|r| {
            format!(
                "(h.chrom={} AND h.start<{} AND h.\"end\">{})",
                sql_string(&r.chromosome),
                r.end_0based_exclusive,
                r.start_0based
            )
        })
        .collect::<Vec<_>>()
        .join(" OR ");
    let mut predicate = format!("({filter})");
    if let Some(n) = request.minimum_score {
        predicate.push_str(&format!(" AND h.score>={n:.17}"));
    }
    if let Some(n) = request.minimum_pwm_relative_score {
        predicate.push_str(&format!(" AND h.pwm_relative_score>={n:.17}"));
    }
    let path_strings = payloads
        .iter()
        .map(|p| p.to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    let mut rows=runtime.sql(&paths.database_path,&format!("SELECT h.* FROM read_parquet({},hive_partitioning=false) h WHERE {predicate} ORDER BY chrom,start,\"end\",motif_id,strand LIMIT {};",sql_string_list(path_strings.iter().map(String::as_str)),request.max_rows+1))?;
    report.truncated = rows.len() > request.max_rows;
    rows.truncate(request.max_rows);
    let mut seen = BTreeSet::new();
    for row in rows {
        let chrom = required_string(&row, "chrom")?;
        let motif = required_string(&row, "motif_id")?;
        let start = number(&row, "start")?;
        let end = number(&row, "end")?;
        let strand = required_string(&row, "strand")?;
        let score = finite(&row, "score")?;
        let floor = finite(&row, "source_minimum_score")?;
        let tags = number(&row, "regulation_tags")?;
        let width = subset
            .motif_metadata
            .iter()
            .find(|m| m["motif_id"] == motif)
            .and_then(|m| m["motif_length"].as_u64())
            .ok_or("Missing exact motif width")?;
        let relative = finite(&row, "pwm_relative_score")?;
        finite(&row, "pseudocount")?;
        for name in [
            "motif_set_id",
            "score_mode",
            "background_model_id",
            "pseudocount_scheme",
            "n_policy",
        ] {
            required_string(&row, name)?;
        }
        if start >= end
            || end > i64::MAX as u64
            || end - start != width
            || !(0.0..=1.0).contains(&relative)
            || report
                .regions
                .iter()
                .find(|r| r.resolved_chromosome.as_deref() == Some(&chrom))
                .and_then(|r| r.package_contig_length_bp)
                .is_none_or(|length| end > length)
            || !matches!(strand.as_str(), "+" | "-")
            || !request.motif_ids.contains(&motif)
            || row["genome_id"] != provider.genome_id
            || row["overlaps_regulatory_tss_intersection"] != true
            || tags > 255
            || tags & 127 == 0
            || tags & 128 == 0
            || score < floor
            || thresholds
                .get(&motif)
                .is_none_or(|t| t.final_minimum_score != Some(floor))
            || !seen.insert((chrom.clone(), start, end, motif.clone(), strand.clone()))
        {
            return Err("Invalid/duplicate regulatory intersection hit or source floor".into());
        }
        for (bit, name) in [
            "overlaps_promoter_core",
            "overlaps_promoter_extended",
            "overlaps_enhancer",
            "overlaps_open_chromatin",
            "overlaps_ctcf",
            "overlaps_emar",
            "overlaps_other_regulatory",
            "overlaps_tss_window",
        ]
        .iter()
        .enumerate()
        {
            if row[*name].as_bool() != Some(tags & (1 << bit) != 0) {
                return Err("Regulatory bitmask disagrees with its source overlap flags".into());
            }
        }
        merge_policy(provider, &mut subset.score_pseudocount, &row)?;
        let matched = regions
            .iter()
            .filter(|r| {
                r.chromosome == chrom && r.start_0based < end && r.end_0based_exclusive > start
            })
            .collect::<Vec<_>>();
        let r = matched.first().ok_or("Hit outside requested intervals")?;
        let (local_start, local_end, local_forward) = source_coordinates(r, start, end, &strand);
        subset.hit_annotations.push(RegulatoryMotifHitAnnotation {
            hit_index: report.hits.len(),
            query_interval_ids: matched.iter().map(|r| r.interval_id.clone()).collect(),
            regulation_tags: tags as u16,
            overlaps_regulatory_tss_intersection: true,
            promoter_ids: subset
                .tss_windows
                .iter()
                .filter(|w| {
                    w.chromosome == chrom && w.start_0based < end && w.end_0based_exclusive > start
                })
                .map(|w| w.promoter_id.clone())
                .collect(),
            regulatory_feature_ids: subset
                .regulatory_features
                .iter()
                .filter(|r| {
                    r["chrom"] == chrom
                        && r["extended_start"]
                            .as_u64()
                            .or_else(|| r["start"].as_u64())
                            .is_some_and(|s| s < end)
                        && r["extended_end"]
                            .as_u64()
                            .or_else(|| r["end"].as_u64())
                            .is_some_and(|e| e > start)
                })
                .filter_map(|r| r["regulatory_feature_id"].as_str().map(str::to_string))
                .collect(),
        });
        report.hits.push(GenomicMotifEvidenceHit {
            interval_id: r.interval_id.clone(),
            chromosome: chrom,
            start_0based: start,
            end_0based_exclusive: end,
            motif_id: motif,
            motif_name: optional_string(&row, "motif_name"),
            strand,
            score,
            pwm_relative_score: optional_f64(&row, "pwm_relative_score"),
            score_mode: required_string(&row, "score_mode")?,
            minimum_score: Some(floor),
            source_seq_id: r.source_seq_id.clone(),
            source_start_0based: local_start,
            source_end_0based_exclusive: local_end,
            source_forward_strand: local_forward,
            ..Default::default()
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests;
