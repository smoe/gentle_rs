//! Bound saved artifacts, original annotation catalogs and streaming reference windows.

use super::*;

pub(super) struct Catalog {
    pub source_id: String,
    pub gene_id: Option<String>,
    pub reason: String,
    pub models: Vec<PrimerPayloadTranscript>,
}
pub(super) struct Loaded {
    pub report: ObservedTranscriptStructuresReport,
    pub data: ObservedTranscriptStructuresData,
    pub catalogs: Vec<Catalog>,
    pub contexts: BTreeMap<s1::ReadRecord, Vec<ObservedSampleContext>>,
}

pub(super) fn load(
    r: &AssessObservedSpliceJunctionsRequest,
    m: &mut Meter<'_>,
    ctl: &mut Control<'_>,
) -> Result<Loaded, EngineError> {
    let ObservedJunctionSource::S1Report {
        result,
        report: file,
    } = &r.source;
    let raw: ObservedTranscriptStructuresReport = read_bound(file, 128_000_000, m, ctl)?;
    m.bytes(file.bytes)?;
    for a in &raw.artifacts {
        m.bytes(a.bytes)?;
    }
    ctl.check()?;
    let inspected = GentleEngine::inspect_observed_transcript_structures(
        &file.path,
        Some(result),
        "read_structures.tsv",
        0,
        true,
    )?;
    if encode(&raw)? != encode(&inspected.report)? {
        return Err(bad("junction_assessment_source_changed"));
    }
    ctl.check()?;
    let report = inspected.report;
    let root = Path::new(&file.path)
        .parent()
        .unwrap_or(Path::new("."))
        .join("artifacts");
    let a = report
        .artifacts
        .iter()
        .find(|a| a.path == "structures.json")
        .ok_or_else(|| bad("junction_assessment_source_data"))?;
    let data: ObservedTranscriptStructuresData =
        read_bound(&artifact_file(&root, a)?, 128_000_000, m, ctl)?;
    if data.schema != "gentle.observed_transcript_structures_data.v1"
        || data.request_sha256 != report.request_sha256
    {
        return Err(bad("junction_assessment_source_data_binding"));
    }
    let mut loaded = Loaded {
        report,
        data,
        catalogs: Vec::new(),
        contexts: BTreeMap::new(),
    };
    if loaded.report.status != "completed" {
        return Ok(loaded);
    }
    let summary = loaded
        .data
        .analysis
        .as_ref()
        .ok_or_else(|| bad("junction_assessment_source_missing_analysis"))?;
    let counts = loaded
        .report
        .counts
        .as_ref()
        .ok_or_else(|| bad("junction_assessment_source_missing_counts"))?;
    if counts.records != summary.records.len()
        || counts.observations != summary.observations.len()
        || counts.chain_families != summary.families.len()
        || counts.unavailable_views != loaded.data.unavailable.len()
        || summary.reference_gene.source_id != loaded.report.request.reference_gene.source_id
        || summary.reference_gene.gene_id != loaded.report.request.reference_gene.gene_id
    {
        return Err(bad("junction_assessment_source_population"));
    }
    let mut record_map = BTreeMap::new();
    let mut identities = BTreeMap::new();
    for record in &summary.records {
        if record_map.insert(record.record.clone(), record).is_some()
            || identities
                .insert(
                    (&record.record.source_sha256, record.record.record_index),
                    &record.record.read_id,
                )
                .is_some_and(|id| id != &record.record.read_id)
        {
            return Err(bad("junction_assessment_source_record_identity"));
        }
    }
    let artifact = loaded
        .report
        .artifacts
        .iter()
        .find(|a| a.path == "read_structures.tsv")
        .ok_or_else(|| bad("junction_assessment_source_read_rows"))?;
    let file = artifact_file(&root, artifact)?;
    m.bytes(file.bytes)?;
    let mut csv = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(VerifiedReader::open(&file, 128_000_000)?);
    if csv
        .headers()
        .map_err(|e| bad(e.to_string()))?
        .iter()
        .collect::<Vec<_>>()
        != ["row_json"]
    {
        return Err(bad("junction_assessment_source_columns"));
    }
    let mut seen = BTreeSet::new();
    let mut file_contexts = BTreeMap::<String, ObservedSampleContext>::new();
    let mut file_sources = BTreeMap::<String, String>::new();
    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    struct Row {
        support: s1::RecordSupport,
        queries: Vec<PrimerPayloadQuery>,
        observations: Vec<s1::Observation>,
    }
    for row in csv.records() {
        ctl.check()?;
        let row = row.map_err(|e| bad(e.to_string()))?;
        if row.len() != 1 || row[0].len() > loaded.report.request.budgets.max_line_bytes {
            return Err(bad("junction_assessment_source_row_budget"));
        }
        let row: Row = serde_json::from_str(&row[0]).map_err(|e| bad(e.to_string()))?;
        let key = &row.support.record;
        if !seen.insert(key.clone())
            || record_map
                .get(key)
                .is_none_or(|s| encode(*s).ok() != encode(&row.support).ok())
        {
            return Err(bad("junction_assessment_source_read_binding"));
        }
        let expected = summary
            .observations
            .iter()
            .filter(|o| o.record == *key)
            .collect::<Vec<_>>();
        if encode(&expected)? != encode(&row.observations)? {
            return Err(bad("junction_assessment_source_observation_binding"));
        }
        let mut contexts = BTreeSet::new();
        for q in row.queries {
            if q.record_index != key.record_index || q.read_id != key.read_id {
                return Err(bad("junction_assessment_source_context_binding"));
            }
            let c = ObservedSampleContext {
                sample_id: q.sample_id,
                barcode: q.barcode,
                priming: q.priming,
                file_id: q.file_id,
            };
            if file_contexts
                .insert(c.file_id.clone(), c.clone())
                .is_some_and(|old| old != c)
                || file_sources
                    .insert(c.file_id.clone(), key.source_sha256.clone())
                    .is_some_and(|old| old != key.source_sha256)
            {
                return Err(bad("junction_assessment_source_file_context"));
            }
            contexts.insert(c);
        }
        loaded
            .contexts
            .insert(key.clone(), contexts.into_iter().collect());
    }
    csv.into_inner().finish()?;
    if seen.len() != summary.records.len() {
        return Err(bad("junction_assessment_source_read_coverage"));
    }
    // A file-wide context can be reused only when retained query evidence established it.
    for record in &summary.records {
        let contexts = loaded.contexts.entry(record.record.clone()).or_default();
        if contexts.is_empty() {
            for u in &loaded.data.unavailable {
                if u.record_index == record.record.record_index
                    && u.read_id == record.record.read_id
                    && file_sources.get(&u.file_id) == Some(&record.record.source_sha256)
                {
                    if let Some(c) = file_contexts.get(&u.file_id) {
                        contexts.push(c.clone());
                    }
                }
            }
            contexts.sort();
            contexts.dedup();
        }
    }
    for resolution in &summary.gene_resolutions {
        ctl.check()?;
        let mut catalog = Catalog {
            source_id: resolution.source_id.clone(),
            gene_id: resolution.gene_id.clone(),
            reason: match resolution.resolution {
                s1::GeneResolution::Resolved => "resolved",
                s1::GeneResolution::MissingGene => "missing_gene",
                s1::GeneResolution::AmbiguousLocus => "ambiguous_locus",
                s1::GeneResolution::UnmappedSource => "unmapped_source",
            }
            .into(),
            models: Vec::new(),
        };
        let Some(a) = loaded
            .report
            .request
            .annotation_sources
            .iter()
            .find(|a| a.source_id == resolution.source_id)
        else {
            loaded.catalogs.push(catalog);
            continue;
        };
        if a.reference.assembly != summary.reference.assembly
            || a.reference.fasta.sha256 != summary.reference.fasta_sha256
        {
            return Err(bad("junction_assessment_annotation_reference"));
        }
        let binding = summary
            .annotation_bindings
            .iter()
            .find(|b| b.source_id == a.source_id)
            .ok_or_else(|| bad("junction_assessment_annotation_source_binding"))?;
        if binding.annotation_sha256 != a.reference.annotation.sha256
            || binding.assembly != a.reference.assembly
            || binding.reference_id != a.reference.reference_id
            || binding.release != a.reference.release
            || binding.reference_fasta_sha256 != a.reference.fasta.sha256
        {
            return Err(bad("junction_assessment_annotation_snapshot_binding"));
        }
        match fs::symlink_metadata(&a.reference.annotation.path) {
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
                catalog.reason = "catalog_unavailable".into();
                loaded.catalogs.push(catalog);
                continue;
            }
            Err(e) => return Err(bad(e.to_string())),
            Ok(_) => {}
        }
        m.bytes(a.reference.annotation.bytes)?;
        let original = reference::catalog(&a.reference, &loaded.report.request.budgets)?;
        catalog.models = original
            .annotation
            .transcripts
            .into_iter()
            .filter(|t| {
                catalog
                    .gene_id
                    .as_ref()
                    .is_some_and(|g| t.gene_ids.contains(g))
            })
            .collect();
        if catalog
            .models
            .iter()
            .any(|t| t.target_length > i64::MAX as usize)
        {
            return Err(bad("junction_assessment_coordinate_range"));
        }
        if resolution.resolution == s1::GeneResolution::Resolved
            && (catalog.models.is_empty()
                || catalog
                    .models
                    .iter()
                    .map(|m| (&m.target_name, m.genomic_strand))
                    .collect::<BTreeSet<_>>()
                    .len()
                    != 1)
        {
            return Err(bad("junction_assessment_resolution_binding"));
        }
        if resolution.resolution == s1::GeneResolution::MissingGene && !catalog.models.is_empty() {
            return Err(bad("junction_assessment_resolution_binding"));
        }
        catalog
            .models
            .sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));
        loaded.catalogs.push(catalog);
    }
    let mut lengths = BTreeMap::new();
    for model in loaded.catalogs.iter().flat_map(|c| &c.models) {
        if lengths
            .insert(&model.target_name, model.target_length)
            .is_some_and(|old| old != model.target_length)
        {
            return Err(bad("junction_assessment_annotation_contig_length"));
        }
    }
    loaded
        .catalogs
        .sort_by(|a, b| a.source_id.cmp(&b.source_id));
    Ok(loaded)
}

pub(super) fn compare_source(
    contig: &str,
    skip: s1::Interval,
    c: &Catalog,
    window: usize,
) -> ObservedJunctionSourceAssessment {
    let mut out = ObservedJunctionSourceAssessment {
        source_id: c.source_id.clone(),
        gene_id: c.gene_id.clone(),
        relation: "not_assessable".into(),
        reason: c.reason.clone(),
        model_strand: None,
        shared_end: None,
        near_hits: Vec::new(),
        events: Vec::new(),
        model_events: Vec::new(),
    };
    if c.reason != "resolved" {
        return out;
    }
    if c.models.iter().any(|m| m.target_name != contig) {
        out.reason = "outside_gene_locus".into();
        return out;
    }
    out.model_strand = c.models.first().and_then(|m| m.genomic_strand);
    out.near_hits = near_hits_for_models(skip, &c.models, window);
    out.relation = near_relation(&out.near_hits).into();
    out.reason = "bound_gene_models".into();
    out.shared_end = match out.relation.as_str() {
        "shared_start_only" => Some(
            if out.model_strand == Some('+') {
                "donor"
            } else {
                "acceptor"
            }
            .into(),
        ),
        "shared_end_only" => Some(
            if out.model_strand == Some('+') {
                "acceptor"
            } else {
                "donor"
            }
            .into(),
        ),
        _ => None,
    };
    let pairs = c
        .models
        .iter()
        .flat_map(|m| m.genomic_exons.windows(2).map(|w| (w[0].1, w[1].0)))
        .collect::<BTreeSet<_>>();
    let start = pairs.iter().any(|p| p.0 == skip.0);
    let end = pairs.iter().any(|p| p.1 == skip.1);
    let novel = start && end && !pairs.contains(&skip);
    let mut labels = BTreeSet::new();
    for model in &c.models {
        let mut events = geometry_events(skip, model);
        if novel {
            events.push("novel_combination".into());
            events.sort();
        }
        labels.extend(events.iter().cloned());
        out.model_events.push(ObservedModelJunctionEvents {
            transcript_id: model.transcript_id.clone(),
            labels: events,
        });
    }
    if start || end {
        labels.remove("unannotated_both_ends");
    }
    out.events = labels.into_iter().collect();
    out.model_events
        .sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));
    out
}

pub(super) type Manifest = BTreeMap<ObservedSampleContext, ObservedSampleBinding>;
pub(super) fn manifest(
    r: &AssessObservedSpliceJunctionsRequest,
    m: &mut Meter<'_>,
    ctl: &mut Control<'_>,
) -> Result<Manifest, EngineError> {
    let Some(file) = &r.sample_manifest else {
        return Ok(BTreeMap::new());
    };
    let data: ObservedSampleManifest = read_bound(file, r.budgets.max_manifest_bytes, m, ctl)?;
    if data.schema != OBSERVED_SAMPLE_MANIFEST_SCHEMA || data.rows.len() > 100_000 {
        return Err(bad("junction_assessment_manifest_schema"));
    }
    let (mut rows, mut libraries, mut preps) = (BTreeMap::new(), BTreeMap::new(), BTreeMap::new());
    for row in data.rows {
        ctl.check()?;
        m.charge(1)?;
        if [
            &row.library_id,
            &row.cdna_prep_id,
            &row.biological_sample_id,
        ]
        .iter()
        .any(|s| s.trim().is_empty() || s.len() > 512 || s.chars().any(char::is_control))
            || [
                &row.context.sample_id,
                &row.context.barcode,
                &row.context.priming,
                &row.context.file_id,
            ]
            .into_iter()
            .chain(row.cell_line.iter())
            .any(|s| s.len() > 4096 || s.chars().any(char::is_control))
        {
            return Err(bad("junction_assessment_manifest_identifier"));
        }
        if libraries
            .insert(row.library_id.clone(), row.cdna_prep_id.clone())
            .is_some_and(|v| v != row.cdna_prep_id)
            || preps
                .insert(row.cdna_prep_id.clone(), row.biological_sample_id.clone())
                .is_some_and(|v| v != row.biological_sample_id)
        {
            return Err(bad("junction_assessment_manifest_hierarchy"));
        }
        if rows.insert(row.context.clone(), row).is_some() {
            return Err(bad("junction_assessment_manifest_duplicate_context"));
        }
    }
    Ok(rows)
}

pub(super) struct Bases {
    pub authority: String,
    pub values: BTreeMap<(String, usize), Vec<u8>>,
}
pub(super) fn reference_bases(
    r: &AssessObservedSpliceJunctionsRequest,
    summary: &s1::Summary,
    physical: &BTreeSet<(String, s1::Interval)>,
    catalogs: &[Catalog],
    m: &mut Meter<'_>,
    ctl: &mut Control<'_>,
) -> Result<Bases, EngineError> {
    let mut out = Bases {
        authority: "not_assessable".into(),
        values: BTreeMap::new(),
    };
    let Some(reference) = &r.reference_sequence else {
        return Ok(out);
    };
    let (declared, base_offset) = match &reference.authority {
        ObservedJunctionReferenceAuthority::FullAssembly => {
            if reference.file.sha256 != summary.reference.fasta_sha256 {
                return Err(bad("junction_assessment_reference_hash_mismatch"));
            }
            (None, 0)
        }
        ObservedJunctionReferenceAuthority::DeclaredExcerpt {
            contig,
            offset_0based,
            length,
        } => (Some((contig.as_str(), *length)), *offset_0based),
    };
    let mut ranges = BTreeMap::<String, BTreeSet<usize>>::new();
    for (contig, (a, z)) in physical {
        if z - a >= 2 {
            ranges
                .entry(contig.clone())
                .or_default()
                .extend([*a, z - 2]);
        }
    }
    let ranges = ranges
        .into_iter()
        .map(|(c, v)| (c, v.into_iter().collect::<Vec<_>>()))
        .collect::<BTreeMap<_, _>>();
    let mut reader = BufReader::new(VerifiedReader::open(
        &reference.file,
        r.budgets.max_reference_bytes,
    )?);
    m.bytes(reference.file.bytes)?;
    let (mut current, mut offset, mut seen, mut name_bytes) =
        (String::new(), 0usize, BTreeSet::new(), 0usize);
    let lengths = catalogs
        .iter()
        .flat_map(|c| &c.models)
        .map(|m| (m.target_name.as_str(), m.target_length))
        .collect::<BTreeMap<_, _>>();
    let finish = |name: &str, n: usize| -> Result<(), EngineError> {
        if let Some((contig, len)) = declared {
            if name != contig || n != len {
                return Err(bad("junction_assessment_excerpt_extent"));
            }
        } else if lengths.get(name).is_some_and(|length| *length != n) {
            return Err(bad("junction_assessment_reference_contig_length"));
        }
        Ok(())
    };
    while let Some(text) = line(&mut reader, 1_048_576)? {
        ctl.check()?;
        if let Some(header) = text.strip_prefix('>') {
            if !current.is_empty() {
                finish(&current, offset)?;
            }
            current = header.split_ascii_whitespace().next().unwrap_or("").into();
            name_bytes = name_bytes
                .checked_add(current.len())
                .filter(|n| *n <= 8_000_000)
                .ok_or_else(|| bad("junction_assessment_reference_header_budget"))?;
            if current.is_empty()
                || current.len() > 4096
                || !seen.insert(current.clone())
                || seen.len() > 100_000
                || declared.is_some_and(|(c, _)| c != current)
            {
                return Err(bad("junction_assessment_reference_header"));
            }
            offset = 0;
        } else {
            if current.is_empty()
                || text.is_empty()
                || !text
                    .bytes()
                    .all(|b| b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(&b))
            {
                return Err(bad("junction_assessment_reference_sequence"));
            }
            let end = offset
                .checked_add(text.len())
                .ok_or_else(|| bad("junction_assessment_reference_length"))?;
            let absolute = base_offset
                .checked_add(offset)
                .ok_or_else(|| bad("junction_assessment_reference_length"))?;
            let absolute_end = base_offset
                .checked_add(end)
                .ok_or_else(|| bad("junction_assessment_reference_length"))?;
            if let Some(wanted) = ranges.get(&current) {
                let lo = wanted.partition_point(|p| p.saturating_add(2) <= absolute);
                let hi = wanted.partition_point(|p| *p < absolute_end);
                for &start in &wanted[lo..hi] {
                    let a = start.max(absolute);
                    let z = (start + 2).min(absolute_end);
                    out.values
                        .entry((current.clone(), start))
                        .or_default()
                        .extend_from_slice(&text.as_bytes()[a - absolute..z - absolute]);
                }
            }
            offset = end;
        }
    }
    if current.is_empty() {
        return Err(bad("junction_assessment_reference_empty"));
    }
    finish(&current, offset)?;
    reader.into_inner().finish()?;
    out.values.retain(|_, v| v.len() == 2);
    out.authority = if declared.is_some() {
        "declared"
    } else {
        "verified"
    }
    .into();
    Ok(out)
}
