//! Pure validation of accession-pinned panels and transcript-oriented TSS inputs.
//!
//! These checks establish structural consistency only. Registry resolution owns
//! exact, case-sensitive factor-name/PFM verification; filesystem readers own
//! containment, byte limits, hashes, and sequence/reference verification.

use gentle_protocol::EngineError;
use gentle_protocol::tss_profiles::{
    BUNDLE_SCHEMA, JasparTargetPanel, PANEL_SCHEMA, SELECTION_SCHEMA, TssBundleManifest,
    TssCalibrationState, TssGeometry, TssRecord, TssReference, TssScaleMode, TssSelection,
    TssStrand,
};
use std::collections::{BTreeMap, BTreeSet};

/// Maximum number of distinct windows admitted by the v1 bundle reader.
pub const MAX_BUNDLE_RECORDS: usize = 10_000;
/// Maximum number of FASTA members admitted by the v1 bundle reader.
pub const MAX_BUNDLE_FASTA_FILES: usize = 1_024;

fn text_field(value: &str, field: &str) -> Result<(), EngineError> {
    if value.is_empty()
        || value.trim() != value
        || value.len() > 4_096
        || value.chars().any(char::is_control)
    {
        return Err(EngineError::invalid_input(format!(
            "{field} must be nonempty, trimmed text without control characters (at most 4096 bytes)"
        )));
    }
    Ok(())
}

fn header_token(value: &str, field: &str) -> Result<(), EngineError> {
    text_field(value, field)?;
    if value
        .chars()
        .any(|c| c.is_whitespace() || matches!(c, '|' | '=' | ','))
    {
        return Err(EngineError::invalid_input(format!(
            "{field} cannot contain whitespace or FASTA header separators"
        )));
    }
    Ok(())
}

/// Require an unprefixed, lowercase SHA-256, matching the bundle byte convention.
pub fn validate_sha256(value: &str) -> Result<(), EngineError> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(EngineError::invalid_input(
            "SHA-256 must contain exactly 64 lowercase hexadecimal digits",
        ));
    }
    Ok(())
}

/// Reference identifiers are independent, exact identities, not aliases.
pub fn validate_reference(reference: &TssReference) -> Result<(), EngineError> {
    text_field(&reference.genome_id, "reference.genome_id")?;
    text_field(&reference.assembly, "reference.assembly")?;
    if let Some(release) = &reference.annotation_release {
        text_field(release, "reference.annotation_release")?;
    }
    Ok(())
}

/// Validate an unclipped inclusive genomic interval against its oriented extents.
///
/// Minus-strand upstream bases lie above the TSS, but stored intervals remain
/// ascending. Arithmetic must also fit the protocol's signed relative axis.
pub fn validate_geometry(geometry: &TssGeometry) -> Result<(), EngineError> {
    header_token(&geometry.chromosome, "geometry.chromosome")?;
    let invalid = || {
        EngineError::invalid_input(
            "TSS geometry must be an unclipped, positive 1-based interval matching strand and upstream/downstream extents without overflow",
        )
    };
    let length = geometry.length().ok_or_else(invalid)?;
    i64::try_from(length).map_err(|_| invalid())?;
    let upstream = u64::try_from(geometry.upstream_bp).map_err(|_| invalid())?;
    let downstream = u64::try_from(geometry.downstream_bp).map_err(|_| invalid())?;
    let (left, right) = match geometry.strand {
        TssStrand::Plus => (upstream, downstream),
        TssStrand::Minus => (downstream, upstream),
    };
    let start = geometry.tss_1based.checked_sub(left).ok_or_else(invalid)?;
    let end = geometry.tss_1based.checked_add(right).ok_or_else(invalid)?;
    if start == 0 || geometry.start_1based != start || geometry.end_1based != end {
        return Err(invalid());
    }
    Ok(())
}

/// Validate record identities, geometry, transcript uniqueness and digest syntax.
pub fn validate_record(record: &TssRecord) -> Result<(), EngineError> {
    header_token(&record.promoter_id, "promoter_id")?;
    header_token(&record.gene_id, "gene_id")?;
    header_token(&record.gene_symbol, "gene_symbol")?;
    validate_geometry(&record.geometry)?;
    validate_sha256(&record.sequence_sha256)?;
    if record.transcripts.is_empty() || record.transcripts.len() > 4_096 {
        return Err(EngineError::invalid_input(
            "Each TSS must have 1..=4096 transcript memberships",
        ));
    }
    let mut transcripts = BTreeSet::new();
    for transcript in &record.transcripts {
        header_token(transcript, "transcript ID")?;
        if !transcripts.insert(transcript) {
            return Err(EngineError::invalid_input(format!(
                "Duplicate transcript ID {transcript} for promoter {}",
                record.promoter_id
            )));
        }
    }
    Ok(())
}

/// Validate a homogeneous manifest without filesystem access or reordering rows.
///
/// One physical TSS is `(reference, chromosome, strand, tss_1based)`, not one
/// transcript or sequence digest. Repeated sequences at distinct loci are valid.
/// The v1 single-gene record cannot express multiple genes sharing a physical
/// TSS: such duplicates fail rather than dropping or merging gene memberships.
pub fn validate_manifest(manifest: &TssBundleManifest) -> Result<(), EngineError> {
    if manifest.schema != BUNDLE_SCHEMA {
        return Err(EngineError::invalid_input(format!(
            "Expected explicit {BUNDLE_SCHEMA} manifest; legacy compatibility is not inferred"
        )));
    }
    validate_reference(&manifest.reference)?;
    if manifest.fasta_files.is_empty() || manifest.fasta_files.len() > MAX_BUNDLE_FASTA_FILES {
        return Err(EngineError::invalid_input(format!(
            "Bundle requires 1..={MAX_BUNDLE_FASTA_FILES} FASTA files"
        )));
    }
    if manifest.records.is_empty() || manifest.records.len() > MAX_BUNDLE_RECORDS {
        return Err(EngineError::invalid_input(format!(
            "Bundle requires 1..={MAX_BUNDLE_RECORDS} TSS records"
        )));
    }
    for digest in manifest.fasta_files.values() {
        validate_sha256(digest)?;
    }
    let mut ids = BTreeSet::new();
    let mut loci = BTreeSet::new();
    for record in &manifest.records {
        validate_record(record)?;
        if !ids.insert(&record.promoter_id) {
            return Err(EngineError::invalid_input(format!(
                "Duplicate promoter ID {}",
                record.promoter_id
            )));
        }
        let geometry = &record.geometry;
        if !loci.insert((
            &geometry.chromosome,
            geometry.strand.as_str(),
            geometry.tss_1based,
        )) {
            return Err(EngineError::invalid_input(format!(
                "Duplicate physical TSS at {}:{}:{}; v1 single-gene records do not support multiple rows or multi-gene memberships at one TSS; no memberships were merged or discarded",
                geometry.chromosome,
                geometry.strand.as_str(),
                geometry.tss_1based
            )));
        }
    }
    Ok(())
}

/// Validate exact selection references; no selection is inferred from gene names.
pub fn validate_selection(
    selection: &TssSelection,
    manifest: &TssBundleManifest,
) -> Result<(), EngineError> {
    if selection.schema != SELECTION_SCHEMA {
        return Err(EngineError::invalid_input(format!(
            "Expected {SELECTION_SCHEMA} selection"
        )));
    }
    validate_reference(&selection.reference)?;
    if selection.reference != manifest.reference {
        return Err(EngineError::invalid_input(
            "Selection reference does not exactly match the bundle reference",
        ));
    }
    let records: BTreeMap<_, _> = manifest
        .records
        .iter()
        .map(|record| (record.promoter_id.as_str(), record.gene_id.as_str()))
        .collect();
    let mut ids = BTreeSet::new();
    for selected in &selection.selected {
        if !ids.insert(&selected.promoter_id) {
            return Err(EngineError::invalid_input(format!(
                "Duplicate selected promoter ID {}",
                selected.promoter_id
            )));
        }
        if records.get(selected.promoter_id.as_str()).copied() != Some(selected.gene_id.as_str()) {
            return Err(EngineError::invalid_input(format!(
                "Selection promoter/gene reference does not match the manifest: {}/{}",
                selected.promoter_id, selected.gene_id
            )));
        }
    }
    Ok(())
}

fn canonical_accession(value: &str) -> bool {
    let Some((base, version)) = value.split_once('.') else {
        return false;
    };
    base.len() == 6
        && base.starts_with("MA")
        && base.as_bytes()[2..].iter().all(u8::is_ascii_digit)
        && !version.starts_with('0')
        && !version.is_empty()
        && version.bytes().all(|b| b.is_ascii_digit())
        && version.parse::<u32>().is_ok()
}

/// Validate v1 panel policies without resolving names or mutating declared order.
///
/// Accessions use `MA` plus four digits and a canonical positive version. Array
/// order must agree with strictly increasing display orders (zero/one origins
/// and gaps are allowed). Expected factor names retain their case for the parent
/// registry's exact comparison; this function does not establish matrix/taxon identity.
pub fn validate_panel(panel: &JasparTargetPanel) -> Result<(), EngineError> {
    if panel.schema != PANEL_SCHEMA {
        return Err(EngineError::invalid_input(format!(
            "Expected explicit {PANEL_SCHEMA} panel"
        )));
    }
    text_field(&panel.panel_id, "panel_id")?;
    text_field(&panel.label, "panel label")?;
    text_field(&panel.calibration_statement, "calibration_statement")?;
    if !matches!(
        panel.score_kind.as_str(),
        "llr_bits"
            | "llr_quantile"
            | "llr_background_quantile"
            | "llr_background_tail_log10"
            | "true_log_odds_bits"
            | "true_log_odds_quantile"
            | "true_log_odds_background_quantile"
            | "true_log_odds_background_tail_log10"
    ) {
        return Err(EngineError::invalid_input(format!(
            "Unsupported TSS panel score kind {}",
            panel.score_kind
        )));
    }
    if panel.factors.is_empty() || panel.factors.len() > 1_024 {
        return Err(EngineError::invalid_input(
            "TSS panels require 1..=1024 exact matrices",
        ));
    }
    if panel.top_hit_count == 0 || panel.top_hit_count > 10_000 {
        return Err(EngineError::invalid_input(
            "TSS panel top_hit_count must be in 1..=10000",
        ));
    }
    match (&panel.calibration_id, &panel.calibration_sha256) {
        (Some(id), Some(digest)) => {
            text_field(id, "calibration_id")?;
            validate_sha256(digest)?;
        }
        (None, None) if panel.calibration_state == TssCalibrationState::MatrixSpecific => {}
        _ => {
            return Err(EngineError::invalid_input(
                "Calibration bindings require both an ID and SHA-256; cross-source calibration requires a binding",
            ));
        }
    }
    if panel.scale_mode == TssScaleMode::Shared
        && panel.calibration_state != TssCalibrationState::CrossSourceCalibrated
    {
        return Err(EngineError::invalid_input(
            "Shared TSS scales require typed cross-source calibration, not a prose claim",
        ));
    }
    let mut accessions = BTreeSet::new();
    let mut track_ids = BTreeSet::new();
    let mut previous_order = None;
    for track in &panel.factors {
        if !canonical_accession(&track.source_id) {
            return Err(EngineError::invalid_input(format!(
                "TSS panel requires a canonical versioned MA accession, not an alias: {}",
                track.source_id
            )));
        }
        if !accessions.insert(&track.source_id) {
            return Err(EngineError::invalid_input(format!(
                "Duplicate TSS panel accession {}",
                track.source_id
            )));
        }
        text_field(&track.factor_id, "expected factor_id")?;
        text_field(&track.label, "track label")?;
        if let Some(id) = &track.track_id {
            text_field(id, "track_id")?;
            if !track_ids.insert(id) {
                return Err(EngineError::invalid_input(format!(
                    "Duplicate track_id {id}"
                )));
            }
        }
        if let Some(label) = &track.factor_label {
            text_field(label, "factor_label")?;
        }
        if track
            .provider_kind
            .as_deref()
            .is_some_and(|kind| kind != "jaspar_pwm")
        {
            return Err(EngineError::invalid_input(format!(
                "Matrix {} requires provider_kind jaspar_pwm",
                track.source_id
            )));
        }
        if previous_order.is_some_and(|previous| track.display_order <= previous) {
            return Err(EngineError::invalid_input(format!(
                "Panel array order conflicts with display_order at {}",
                track.source_id
            )));
        }
        previous_order = Some(track.display_order);
        if track
            .score_kind
            .as_ref()
            .is_some_and(|kind| kind != &panel.score_kind)
        {
            return Err(EngineError::invalid_input(format!(
                "Conflicting score kind for {}: v1 requires the panel's single score kind {}",
                track.source_id, panel.score_kind
            )));
        }
        if let Some(color) = &track.color_hint {
            if color.len() != 7
                || !color.starts_with('#')
                || !color.as_bytes()[1..].iter().all(u8::is_ascii_hexdigit)
            {
                return Err(EngineError::invalid_input(format!(
                    "Unsupported color_hint for {}: expected #RRGGBB",
                    track.source_id
                )));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use gentle_protocol::tss_profiles::{JasparPanelTrack, TssSelectedRecord, TssStrandPolicy};

    fn geometry(strand: TssStrand, upstream: usize, downstream: usize) -> TssGeometry {
        let (left, right) = match strand {
            TssStrand::Plus => (upstream, downstream),
            TssStrand::Minus => (downstream, upstream),
        };
        TssGeometry {
            chromosome: "synthetic_chr".into(),
            strand,
            tss_1based: 1_000,
            start_1based: 1_000 - left as u64,
            end_1based: 1_000 + right as u64,
            upstream_bp: upstream,
            downstream_bp: downstream,
        }
    }

    fn panel() -> JasparTargetPanel {
        JasparTargetPanel {
            schema: PANEL_SCHEMA.into(),
            panel_id: "synthetic-panel".into(),
            label: "Synthetic validator input, not registry evidence".into(),
            score_kind: "llr_bits".into(),
            clip_negative: false,
            scale_mode: TssScaleMode::Independent,
            strand_policy: TssStrandPolicy::Both,
            calibration_state: TssCalibrationState::MatrixSpecific,
            calibration_statement: "Matrix-specific scales are not comparable".into(),
            calibration_id: None,
            calibration_sha256: None,
            top_hit_count: 5,
            factors: ["MA0001.1", "MA0002.2"]
                .into_iter()
                .enumerate()
                .map(|(index, source_id)| JasparPanelTrack {
                    source_id: source_id.into(),
                    factor_id: "CaseSensitive".into(),
                    label: format!("Synthetic matrix {index}"),
                    display_order: index,
                    color_hint: Some("#aBc123".into()),
                    score_kind: None,
                    track_id: None,
                    provider_kind: None,
                    factor_label: None,
                })
                .collect(),
        }
    }

    fn manifest() -> TssBundleManifest {
        // Inline synthetic identities and digest syntax, not sequence evidence.
        TssBundleManifest {
            schema: BUNDLE_SCHEMA.into(),
            reference: TssReference {
                genome_id: "synthetic-genome".into(),
                assembly: "synthetic-assembly".into(),
                annotation_release: Some("synthetic-release".into()),
            },
            fasta_files: BTreeMap::from([("synthetic.fa".into(), "a".repeat(64))]),
            records: vec![TssRecord {
                promoter_id: "synthetic-promoter".into(),
                gene_id: "synthetic-gene".into(),
                gene_symbol: "SyntheticGene".into(),
                geometry: geometry(TssStrand::Plus, 3, 7),
                transcripts: vec!["synthetic.tx2".into(), "synthetic.tx1".into()],
                sequence_sha256: "b".repeat(64),
            }],
        }
    }

    #[test]
    fn manifest_rejects_duplicate_ids_and_loci_but_not_repeated_digests() {
        let mut m = manifest();
        validate_manifest(&m).unwrap();
        m.records.push(m.records[0].clone());
        assert!(
            validate_manifest(&m)
                .unwrap_err()
                .message
                .contains("Duplicate promoter")
        );
        m.records[1].promoter_id = "another-promoter".into();
        m.records[1].gene_id = "another-gene".into();
        let error = validate_manifest(&m).unwrap_err();
        assert!(error.message.contains("Duplicate physical TSS"));
        assert!(error.message.contains("multi-gene memberships"));
        m.records[1].geometry.tss_1based += 100;
        m.records[1].geometry.start_1based += 100;
        m.records[1].geometry.end_1based += 100;
        validate_manifest(&m).unwrap();
        assert_eq!(m.records[0].sequence_sha256, m.records[1].sequence_sha256);
        m.records[1].geometry = geometry(TssStrand::Minus, 3, 7);
        validate_manifest(&m).unwrap();
        m.records[1].transcripts.push("synthetic.tx1".into());
        assert!(
            validate_manifest(&m)
                .unwrap_err()
                .message
                .contains("Duplicate transcript")
        );
    }

    #[test]
    fn selection_references_require_matching_gene_reference_and_unique_promoters() {
        let m = manifest();
        let mut selection = TssSelection {
            schema: SELECTION_SCHEMA.into(),
            reference: m.reference.clone(),
            selected: vec![],
        };
        validate_selection(&selection, &m).unwrap();
        let selected = TssSelectedRecord {
            promoter_id: m.records[0].promoter_id.clone(),
            gene_id: m.records[0].gene_id.clone(),
        };
        selection.selected.push(selected.clone());
        validate_selection(&selection, &m).unwrap();
        selection.selected[0].gene_id = "another-gene".into();
        assert!(validate_selection(&selection, &m).is_err());
        selection.selected[0] = selected.clone();
        selection.reference.annotation_release = Some("another-release".into());
        assert!(validate_selection(&selection, &m).is_err());
        selection.reference = m.reference.clone();
        selection.selected.push(selected);
        assert!(
            validate_selection(&selection, &m)
                .unwrap_err()
                .message
                .contains("Duplicate selected")
        );
    }

    #[test]
    fn digest_and_reference_strings_are_not_silently_normalized() {
        validate_sha256(&"ab0123456789cdef".repeat(4)).unwrap();
        for digest in [
            "a".repeat(63),
            "A".repeat(64),
            format!("sha256:{}", "a".repeat(64)),
        ] {
            assert!(validate_sha256(&digest).is_err());
        }
        let mut m = manifest();
        m.reference.annotation_release = Some(String::new());
        assert!(validate_manifest(&m).is_err());
        m = manifest();
        m.schema = "legacy-unknown-bundle".into();
        assert!(
            validate_manifest(&m)
                .unwrap_err()
                .message
                .contains("legacy compatibility is not inferred")
        );
    }

    #[test]
    fn geometry_maps_plus_minus_and_nondefault_windows() {
        for strand in [TssStrand::Plus, TssStrand::Minus] {
            for (upstream, downstream) in [(500, 200), (3, 7), (0, 0)] {
                let g = geometry(strand, upstream, downstream);
                validate_geometry(&g).unwrap();
                let length = g.length().unwrap();
                assert_eq!(g.genomic_at(upstream), Some(1_000));
                assert_eq!(g.relative_at(upstream), Some(0));
                assert_eq!(g.relative_at(0), Some(-(upstream as i64)));
                assert_eq!(g.relative_at(length - 1), Some(downstream as i64));
                let endpoints = match strand {
                    TssStrand::Plus => (g.start_1based, g.end_1based),
                    TssStrand::Minus => (g.end_1based, g.start_1based),
                };
                assert_eq!(g.genomic_at(0), Some(endpoints.0));
                assert_eq!(g.genomic_at(length - 1), Some(endpoints.1));
                assert_eq!(g.genomic_at(length), None);
            }
        }
    }

    #[test]
    fn geometry_rejects_clipping_zero_and_overflow() {
        let valid = geometry(TssStrand::Minus, 3, 7);
        let mut g = valid.clone();
        g.end_1based -= 1;
        assert!(validate_geometry(&g).is_err());
        g = valid.clone();
        g.tss_1based = 0;
        assert!(validate_geometry(&g).is_err());
        g = valid.clone();
        g.tss_1based = u64::MAX;
        assert!(validate_geometry(&g).is_err());
        g = valid.clone();
        g.upstream_bp = usize::MAX;
        assert!(validate_geometry(&g).is_err());
        g = valid;
        g.tss_1based = 7;
        g.start_1based = 0;
        g.end_1based = 10;
        assert!(validate_geometry(&g).is_err());
    }

    #[test]
    fn panel_preserves_case_order_and_settings_for_registry_validation() {
        let mut p = panel();
        validate_panel(&p).unwrap();
        p.factors[0].factor_id = "casesensitive".into();
        p.factors[0].display_order = 10;
        p.factors[1].display_order = 20;
        validate_panel(&p).unwrap();
        assert_eq!(p.factors[0].factor_id, "casesensitive");
        assert_eq!(p.factors[0].source_id, "MA0001.1");
        assert_eq!(p.top_hit_count, 5);
    }

    #[test]
    fn panel_rejects_noncanonical_accessions_duplicates_and_order_conflicts() {
        for accession in [
            "TP53",
            "MA0001",
            "ma0001.1",
            "MA001.1",
            "MA00001.1",
            "MA0001.0",
            "MA0001.01",
            "MA0001.+1",
            "MA0001.1.2",
            " MA0001.1",
            "MA0001.4294967296",
        ] {
            let mut p = panel();
            p.factors[0].source_id = accession.into();
            assert!(validate_panel(&p).is_err(), "{accession}");
        }
        let mut p = panel();
        p.factors[1].source_id = p.factors[0].source_id.clone();
        assert!(validate_panel(&p).is_err());
        p = panel();
        p.factors[1].display_order = 0;
        assert!(validate_panel(&p).is_err());
        p = panel();
        p.factors.swap(0, 1);
        assert!(validate_panel(&p).is_err());
    }

    #[test]
    fn panel_rejects_mixed_scores_unsupported_settings_and_unbound_shared_scales() {
        let mut p = panel();
        p.factors[1].score_kind = Some("true_log_odds_bits".into());
        let error = validate_panel(&p).unwrap_err();
        assert!(error.message.contains("MA0002.2"));
        p = panel();
        p.score_kind = "approximate_affinity".into();
        assert!(validate_panel(&p).is_err());
        p = panel();
        p.factors[0].color_hint = Some("red".into());
        assert!(validate_panel(&p).is_err());
        p = panel();
        p.top_hit_count = 0;
        assert!(validate_panel(&p).is_err());
        p = panel();
        p.scale_mode = TssScaleMode::Shared;
        assert!(validate_panel(&p).is_err());
        p.calibration_state = TssCalibrationState::CrossSourceCalibrated;
        assert!(validate_panel(&p).is_err());
        p.calibration_id = Some("synthetic-calibration".into());
        p.calibration_sha256 = Some("a".repeat(64));
        validate_panel(&p).unwrap();
        p.calibration_sha256 = Some("A".repeat(64));
        assert!(validate_panel(&p).is_err());
    }
}
