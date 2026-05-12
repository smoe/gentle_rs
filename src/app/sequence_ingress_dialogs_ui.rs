//! Sequence and protein ingress dialog helpers.
//!
//! This module is a move-only extraction from `app.rs`: it keeps the GenBank,
//! dbSNP, UniProt, Ensembl protein, reverse-translation, protease digest, and
//! protein-to-DNA handoff dialog helpers close to `GENtleApp` while reducing
//! the top-level app monolith.

use super::*;

impl GENtleApp {
    pub(super) fn uniprot_optional_trimmed(value: &str) -> Option<String> {
        let trimmed = value.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_string())
        }
    }

    pub(super) fn parse_protein_feature_key_list(value: &str) -> Vec<String> {
        let mut seen = HashSet::new();
        let mut keys = Vec::new();
        for part in value.split(|c: char| c == ',' || c == ';' || c.is_whitespace()) {
            let trimmed = part.trim();
            if trimmed.is_empty() {
                continue;
            }
            if seen.insert(trimmed.to_string()) {
                keys.push(trimmed.to_string());
            }
        }
        keys
    }

    pub(super) fn protein_feature_filter_from_dialog(&self) -> gentle_protocol::ProteinFeatureFilter {
        gentle_protocol::ProteinFeatureFilter {
            include_feature_keys: Self::parse_protein_feature_key_list(
                &self.protein_feature_key_include,
            ),
            exclude_feature_keys: Self::parse_protein_feature_key_list(
                &self.protein_feature_key_exclude,
            ),
        }
    }

    pub(super) fn fetch_genbank_accession_from_dialog(&mut self) {
        let accession = self.genbank_accession.trim().to_string();
        if accession.is_empty() {
            self.genbank_status = "GenBank accession cannot be empty".to_string();
            return;
        }
        let op = Operation::FetchGenBankAccession {
            accession: accession.clone(),
            as_id: Self::uniprot_optional_trimmed(&self.genbank_as_id),
        };
        let result = { self.engine.write().unwrap().apply(op) };
        match result {
            Ok(result) => {
                for seq_id in &result.created_seq_ids {
                    self.open_sequence_window(seq_id);
                }
                if self.genbank_as_id.trim().is_empty() {
                    self.genbank_as_id = accession;
                }
                self.genbank_status = Self::format_op_result_status(
                    "GenBank fetch: ok",
                    &result.created_seq_ids,
                    &result.warnings,
                    &result.messages,
                );
            }
            Err(err) => {
                self.genbank_status = format!("GenBank fetch failed: {}", err.message);
            }
        }
    }

    pub(super) fn fetch_dbsnp_region_from_dialog(&mut self) {
        if self.dbsnp_fetch_task.is_some() {
            self.dbsnp_status = "A dbSNP fetch is already running".to_string();
            return;
        }
        let rs_id = self.dbsnp_rs_id.trim().to_string();
        if rs_id.is_empty() {
            self.dbsnp_status = "dbSNP rs_id cannot be empty".to_string();
            return;
        }
        let genome_id = self.dbsnp_genome_id.trim().to_string();
        if genome_id.is_empty() {
            self.dbsnp_status = "Prepared genome id cannot be empty".to_string();
            return;
        }
        let flank_bp = if self.dbsnp_flank_bp.trim().is_empty() {
            3000
        } else {
            match self.dbsnp_flank_bp.trim().parse::<usize>() {
                Ok(value) => value,
                Err(_) => {
                    self.dbsnp_status =
                        "flank bp must be a non-negative integer (empty = 3000)".to_string();
                    return;
                }
            }
        };
        let op = Operation::FetchDbSnpRegion {
            rs_id: rs_id.clone(),
            genome_id: genome_id.clone(),
            flank_bp: Some(flank_bp),
            output_id: Self::uniprot_optional_trimmed(&self.dbsnp_output_id),
            annotation_scope: Some(GenomeAnnotationScope::Full),
            max_annotation_features: Some(0),
            catalog_path: self.dbsnp_catalog_path_opt(),
            cache_dir: self.dbsnp_cache_dir_opt(),
        };
        let (tx, rx) = mpsc::channel::<DbSnpFetchTaskMessage>();
        let engine = self.engine.clone();
        self.dbsnp_status = format!(
            "Starting dbSNP fetch for '{}' against '{}'",
            rs_id, genome_id
        );
        self.dbsnp_fetch_task = Some(DbSnpFetchTask {
            started: Instant::now(),
            receiver: rx,
        });
        std::thread::spawn(move || {
            let tx_progress = tx.clone();
            let result = {
                let mut guard = engine.write().expect("Engine lock poisoned");
                guard.apply_with_progress(op, move |progress| {
                    if let OperationProgress::DbSnpFetch(progress) = progress {
                        let _ = tx_progress.send(DbSnpFetchTaskMessage::Progress(progress));
                    }
                    true
                })
            };
            let _ = tx.send(DbSnpFetchTaskMessage::Done(result));
        });
    }

    pub(super) fn sort_uniprot_entries_recent(mut rows: Vec<UniprotEntrySummary>) -> Vec<UniprotEntrySummary> {
        rows.sort_by(|left, right| {
            right
                .imported_at_unix_ms
                .cmp(&left.imported_at_unix_ms)
                .then_with(|| left.entry_id.cmp(&right.entry_id))
        });
        rows
    }

    pub(super) fn recent_uniprot_entries_for_dialog(&self, limit: usize) -> Vec<UniprotEntrySummary> {
        let rows = self.engine.read().unwrap().list_uniprot_entries();
        let mut rows = Self::sort_uniprot_entries_recent(rows);
        let cap = limit.max(1);
        if rows.len() > cap {
            rows.truncate(cap);
        }
        rows
    }

    pub(super) fn sort_ensembl_protein_entries_recent(
        mut rows: Vec<EnsemblProteinEntrySummary>,
    ) -> Vec<EnsemblProteinEntrySummary> {
        rows.sort_by(|left, right| {
            right
                .imported_at_unix_ms
                .cmp(&left.imported_at_unix_ms)
                .then_with(|| left.entry_id.cmp(&right.entry_id))
        });
        rows
    }

    pub(super) fn recent_ensembl_protein_entries_for_dialog(
        &self,
        limit: usize,
    ) -> Vec<EnsemblProteinEntrySummary> {
        let rows = self.engine.read().unwrap().list_ensembl_protein_entries();
        let mut rows = Self::sort_ensembl_protein_entries_recent(rows);
        let cap = limit.max(1);
        if rows.len() > cap {
            rows.truncate(cap);
        }
        rows
    }

    pub(super) fn selected_ensembl_protein_entry_for_dialog(&self) -> Option<EnsemblProteinEntry> {
        let entry_id = self.ensembl_protein_entry_id.trim();
        if entry_id.is_empty() {
            return None;
        }
        self.engine
            .read()
            .unwrap()
            .get_ensembl_protein_entry(entry_id)
            .ok()
    }

    pub(super) fn summarize_ensembl_protein_feature_keys(
        entry: &EnsemblProteinEntry,
        limit: usize,
    ) -> Vec<String> {
        let mut counts: HashMap<String, usize> = HashMap::new();
        for feature in &entry.features {
            let label = if !feature.feature_key.trim().is_empty() {
                feature.feature_key.trim().to_string()
            } else if !feature.feature_type.trim().is_empty() {
                feature.feature_type.trim().to_string()
            } else {
                "feature".to_string()
            };
            *counts.entry(label).or_insert(0) += 1;
        }
        let mut rows = counts.into_iter().collect::<Vec<_>>();
        rows.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
        rows.into_iter()
            .take(limit.max(1))
            .map(|(label, count)| format!("{label} x{count}"))
            .collect()
    }

    pub(super) fn sort_uniprot_projections_recent(
        mut rows: Vec<UniprotGenomeProjectionSummary>,
    ) -> Vec<UniprotGenomeProjectionSummary> {
        rows.sort_by(|left, right| {
            right
                .created_at_unix_ms
                .cmp(&left.created_at_unix_ms)
                .then_with(|| left.projection_id.cmp(&right.projection_id))
        });
        rows
    }

    pub(super) fn resolve_uniprot_projection_id_from_dialog_fields(&self) -> Option<String> {
        let seq_id = self.uniprot_map_seq_id.trim();
        let entry_id = self.uniprot_entry_id.trim();
        if seq_id.is_empty() || entry_id.is_empty() {
            return None;
        }
        Some(
            Self::uniprot_optional_trimmed(&self.uniprot_map_projection_id)
                .unwrap_or_else(|| format!("{entry_id}@{seq_id}")),
        )
    }

    pub(super) fn recent_uniprot_projections_for_dialog(
        &self,
        limit: usize,
    ) -> Vec<UniprotGenomeProjectionSummary> {
        let seq_filter = Self::uniprot_optional_trimmed(&self.uniprot_map_seq_id);
        let entry_filter = Self::uniprot_optional_trimmed(&self.uniprot_entry_id);
        let mut rows = self
            .engine
            .read()
            .unwrap()
            .list_uniprot_genome_projections(seq_filter.as_deref());
        if let Some(entry_id) = entry_filter.as_deref() {
            rows.retain(|row| row.entry_id == entry_id);
        }
        let mut rows = Self::sort_uniprot_projections_recent(rows);
        let cap = limit.max(1);
        if rows.len() > cap {
            rows.truncate(cap);
        }
        rows
    }

    pub(super) fn recent_uniprot_audit_reports_for_dialog(
        &self,
        limit: usize,
    ) -> Vec<crate::engine::UniprotProjectionAuditReportSummary> {
        let seq_filter = Self::uniprot_optional_trimmed(&self.uniprot_map_seq_id);
        let entry_filter = Self::uniprot_optional_trimmed(&self.uniprot_entry_id);
        let mut rows = self
            .engine
            .read()
            .unwrap()
            .list_uniprot_projection_audit_reports(seq_filter.as_deref());
        if let Some(entry_id) = entry_filter.as_deref() {
            rows.retain(|row| row.entry_id == entry_id);
        }
        rows.sort_by(|left, right| {
            right
                .generated_at_unix_ms
                .cmp(&left.generated_at_unix_ms)
                .then_with(|| left.report_id.cmp(&right.report_id))
        });
        rows.truncate(limit.max(1));
        rows
    }

    pub(super) fn recent_uniprot_audit_parity_reports_for_dialog(
        &self,
        limit: usize,
    ) -> Vec<crate::engine::UniprotProjectionAuditParityReportSummary> {
        let seq_filter = Self::uniprot_optional_trimmed(&self.uniprot_map_seq_id);
        let entry_filter = Self::uniprot_optional_trimmed(&self.uniprot_entry_id);
        let mut rows = self
            .engine
            .read()
            .unwrap()
            .list_uniprot_projection_audit_parity_reports(seq_filter.as_deref());
        if let Some(entry_id) = entry_filter.as_deref() {
            rows.retain(|row| row.entry_id == entry_id);
        }
        rows.sort_by(|left, right| {
            right
                .generated_at_unix_ms
                .cmp(&left.generated_at_unix_ms)
                .then_with(|| left.report_id.cmp(&right.report_id))
        });
        rows.truncate(limit.max(1));
        rows
    }

    pub(super) fn open_sequence_window_for_uniprot_projection_expert(
        &mut self,
        seq_id: &str,
        projection_id: &str,
        protein_feature_filter: gentle_protocol::ProteinFeatureFilter,
    ) {
        if projection_id.trim().is_empty() {
            self.open_sequence_window(seq_id);
            return;
        }
        if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
            if let Some(window) = self.windows.get(&viewport_id) {
                if let Ok(mut window) = window.write() {
                    window.focus_uniprot_projection_expert(
                        projection_id,
                        protein_feature_filter.clone(),
                    );
                }
            }
            self.queue_focus_viewport(viewport_id);
            return;
        }
        if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
            window.focus_uniprot_projection_expert(projection_id, protein_feature_filter.clone());
            return;
        }
        let exists = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key(seq_id);
        if exists {
            let mut window = Window::new_dna_lazy(seq_id.to_string(), self.engine.clone());
            window.focus_uniprot_projection_expert(projection_id, protein_feature_filter);
            self.new_windows.push(window);
        }
    }

    pub(super) fn open_sequence_window_for_transcript_protein_expert(
        &mut self,
        seq_id: &str,
        transcript_id_filter: Option<&str>,
    ) {
        if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
            if let Some(window) = self.windows.get(&viewport_id) {
                if let Ok(mut window) = window.write() {
                    window.focus_transcript_protein_expert(transcript_id_filter);
                }
            }
            self.queue_focus_viewport(viewport_id);
            return;
        }
        if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
            window.focus_transcript_protein_expert(transcript_id_filter);
            return;
        }
        let exists = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key(seq_id);
        if exists {
            let mut window = Window::new_dna_lazy(seq_id.to_string(), self.engine.clone());
            window.focus_transcript_protein_expert(transcript_id_filter);
            self.new_windows.push(window);
        }
    }

    pub(super) fn open_sequence_window_for_ensembl_protein_expert(
        &mut self,
        seq_id: &str,
        transcript_id_filter: Option<&str>,
        entry_id: &str,
        protein_feature_filter: gentle_protocol::ProteinFeatureFilter,
    ) {
        if let Some(viewport_id) = self.find_open_sequence_viewport_id(seq_id) {
            if let Some(window) = self.windows.get(&viewport_id) {
                if let Ok(mut window) = window.write() {
                    window.focus_ensembl_entry_protein_expert(
                        transcript_id_filter,
                        entry_id,
                        protein_feature_filter.clone(),
                    );
                }
            }
            self.queue_focus_viewport(viewport_id);
            return;
        }
        if let Some(window) = self.find_pending_sequence_window_mut(seq_id) {
            window.focus_ensembl_entry_protein_expert(
                transcript_id_filter,
                entry_id,
                protein_feature_filter.clone(),
            );
            return;
        }
        let exists = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key(seq_id);
        if exists {
            let mut window = Window::new_dna_lazy(seq_id.to_string(), self.engine.clone());
            window.focus_ensembl_entry_protein_expert(
                transcript_id_filter,
                entry_id,
                protein_feature_filter,
            );
            self.new_windows.push(window);
        }
    }

    pub(super) fn open_uniprot_projection_expert_from_dialog(&mut self, seq_id: &str, projection_id: &str) {
        let seq_id = seq_id.trim();
        let projection_id = projection_id.trim();
        if seq_id.is_empty() {
            self.uniprot_status = "seq_id is required for UniProt protein expert".to_string();
            return;
        }
        if projection_id.is_empty() {
            self.uniprot_status =
                "projection_id is required for UniProt protein expert".to_string();
            return;
        }
        match self
            .engine
            .read()
            .unwrap()
            .get_uniprot_genome_projection(projection_id)
        {
            Ok(projection) => {
                if projection.seq_id != seq_id {
                    self.uniprot_status = format!(
                        "Projection '{}' belongs to '{}' rather than '{}'",
                        projection_id, projection.seq_id, seq_id
                    );
                    return;
                }
            }
            Err(err) => {
                self.uniprot_status =
                    format!("Could not open UniProt protein expert: {}", err.message);
                return;
            }
        }
        let protein_feature_filter = self.protein_feature_filter_from_dialog();
        self.open_sequence_window_for_uniprot_projection_expert(
            seq_id,
            projection_id,
            protein_feature_filter,
        );
        self.uniprot_map_seq_id = seq_id.to_string();
        self.uniprot_map_projection_id = projection_id.to_string();
        self.uniprot_status = format!(
            "Opening Protein Expert for UniProt projection '{}' on '{}'",
            projection_id, seq_id
        );
    }

    pub(super) fn open_transcript_protein_expert_from_dialog(
        &mut self,
        seq_id: &str,
        transcript_id_filter: Option<&str>,
    ) {
        let seq_id = seq_id.trim();
        if seq_id.is_empty() {
            self.uniprot_status =
                "seq_id is required for transcript-native Protein Expert".to_string();
            return;
        }
        self.open_sequence_window_for_transcript_protein_expert(seq_id, transcript_id_filter);
        self.uniprot_map_seq_id = seq_id.to_string();
        self.uniprot_status = if let Some(transcript_id_filter) = transcript_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            format!(
                "Opening transcript-native Protein Expert for '{}' on '{}'",
                transcript_id_filter, seq_id
            )
        } else {
            format!("Opening transcript-native Protein Expert on '{}'", seq_id)
        };
    }

    pub(super) fn open_ensembl_protein_expert_from_dialog(
        &mut self,
        seq_id: &str,
        transcript_id_filter: Option<&str>,
        entry_id: &str,
    ) {
        let seq_id = seq_id.trim();
        let entry_id = entry_id.trim();
        if seq_id.is_empty() {
            self.uniprot_status = "seq_id is required for Ensembl Protein Expert".to_string();
            return;
        }
        if entry_id.is_empty() {
            self.uniprot_status =
                "entry_id is required before opening an Ensembl Protein Expert".to_string();
            return;
        }
        match self
            .engine
            .read()
            .unwrap()
            .get_ensembl_protein_entry(entry_id)
        {
            Ok(entry) => {
                if transcript_id_filter
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .is_none()
                {
                    self.uniprot_map_transcript_id = entry.transcript_id.clone();
                }
            }
            Err(err) => {
                self.uniprot_status =
                    format!("Could not open Ensembl Protein Expert: {}", err.message);
                return;
            }
        }
        let protein_feature_filter = self.protein_feature_filter_from_dialog();
        self.open_sequence_window_for_ensembl_protein_expert(
            seq_id,
            transcript_id_filter,
            entry_id,
            protein_feature_filter,
        );
        self.uniprot_map_seq_id = seq_id.to_string();
        self.ensembl_protein_entry_id = entry_id.to_string();
        self.uniprot_status = if let Some(transcript_id_filter) = transcript_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            format!(
                "Opening Ensembl Protein Expert for '{}' / '{}' on '{}'",
                entry_id, transcript_id_filter, seq_id
            )
        } else {
            format!(
                "Opening Ensembl Protein Expert for '{}' on '{}'",
                entry_id, seq_id
            )
        };
    }

    pub(super) fn fetch_uniprot_entry_from_dialog(&mut self) {
        let query = self.uniprot_query.trim().to_string();
        if query.is_empty() {
            self.uniprot_status = "UniProt query cannot be empty".to_string();
            return;
        }
        let op = Operation::FetchUniprotSwissProt {
            query: query.clone(),
            entry_id: Self::uniprot_optional_trimmed(&self.uniprot_entry_id),
        };
        match self.engine.write().unwrap().apply(op) {
            Ok(result) => {
                if self.uniprot_entry_id.trim().is_empty() {
                    self.uniprot_entry_id = query;
                }
                self.uniprot_status = Self::format_op_result_status(
                    "UniProt fetch: ok",
                    &result.created_seq_ids,
                    &result.warnings,
                    &result.messages,
                );
            }
            Err(err) => {
                self.uniprot_status = format!("UniProt fetch failed: {}", err.message);
            }
        }
    }

    pub(super) fn import_uniprot_swiss_prot_from_dialog(&mut self) {
        let path = self.uniprot_swiss_path.trim().to_string();
        if path.is_empty() {
            self.uniprot_status = "SWISS-PROT file path cannot be empty".to_string();
            return;
        }
        let op = Operation::ImportUniprotSwissProt {
            path,
            entry_id: Self::uniprot_optional_trimmed(&self.uniprot_entry_id),
        };
        match self.engine.write().unwrap().apply(op) {
            Ok(result) => {
                self.uniprot_status = Self::format_op_result_status(
                    "UniProt import: ok",
                    &result.created_seq_ids,
                    &result.warnings,
                    &result.messages,
                );
            }
            Err(err) => {
                self.uniprot_status = format!("UniProt import failed: {}", err.message);
            }
        }
    }

    pub(super) fn fetch_ensembl_protein_from_dialog(&mut self) {
        let query = self.ensembl_protein_query.trim().to_string();
        if query.is_empty() {
            self.uniprot_status = "Ensembl protein query cannot be empty".to_string();
            return;
        }
        let entry_id_override = Self::uniprot_optional_trimmed(&self.ensembl_protein_entry_id);
        let op = Operation::FetchEnsemblProtein {
            query,
            entry_id: entry_id_override.clone(),
        };
        match self.engine.write().unwrap().apply(op) {
            Ok(result) => {
                if let Some(row) = self.recent_ensembl_protein_entries_for_dialog(1).first() {
                    self.ensembl_protein_entry_id = row.entry_id.clone();
                    if self.ensembl_protein_query.trim().is_empty() {
                        self.ensembl_protein_query = row.protein_id.clone();
                    }
                    if self.uniprot_map_transcript_id.trim().is_empty() {
                        self.uniprot_map_transcript_id = row.transcript_id.clone();
                    }
                } else if let Some(entry_id_override) = entry_id_override {
                    self.ensembl_protein_entry_id = entry_id_override;
                }
                self.uniprot_status = Self::format_op_result_status(
                    "Ensembl protein fetch: ok",
                    &result.created_seq_ids,
                    &result.warnings,
                    &result.messages,
                );
            }
            Err(err) => {
                self.uniprot_status = format!("Ensembl protein fetch failed: {}", err.message);
            }
        }
    }

    pub(super) fn import_ensembl_protein_sequence_from_dialog(&mut self) {
        let entry_id = self.ensembl_protein_entry_id.trim().to_string();
        if entry_id.is_empty() {
            self.uniprot_status =
                "entry_id is required before importing an Ensembl protein sequence".to_string();
            return;
        }
        let op = Operation::ImportEnsemblProteinSequence {
            entry_id,
            output_id: Self::uniprot_optional_trimmed(&self.ensembl_protein_output_id),
        };
        let result = { self.engine.write().unwrap().apply(op) };
        match result {
            Ok(result) => {
                for seq_id in &result.created_seq_ids {
                    self.open_sequence_window(seq_id);
                }
                self.uniprot_status = Self::format_op_result_status(
                    "Ensembl protein sequence import: ok",
                    &result.created_seq_ids,
                    &result.warnings,
                    &result.messages,
                );
            }
            Err(err) => {
                self.uniprot_status =
                    format!("Ensembl protein sequence import failed: {}", err.message);
            }
        }
    }

    pub(super) fn fetch_uniprot_linked_genbank_from_dialog(&mut self) {
        let entry_id = self.uniprot_entry_id.trim().to_string();
        if entry_id.is_empty() {
            self.uniprot_status = "entry_id is required for linked GenBank fetch".to_string();
            return;
        }
        let op = Operation::FetchUniprotLinkedGenBank {
            entry_id,
            accession: Self::uniprot_optional_trimmed(&self.uniprot_linked_accession),
            as_id: Self::uniprot_optional_trimmed(&self.uniprot_linked_as_id),
        };
        let result = { self.engine.write().unwrap().apply(op) };
        match result {
            Ok(result) => {
                for seq_id in &result.created_seq_ids {
                    self.open_sequence_window(seq_id);
                }
                if let Some(seq_id) = result.created_seq_ids.first() {
                    self.uniprot_map_seq_id = seq_id.clone();
                }
                self.uniprot_status = Self::format_op_result_status(
                    "UniProt linked GenBank fetch: ok",
                    &result.created_seq_ids,
                    &result.warnings,
                    &result.messages,
                );
            }
            Err(err) => {
                self.uniprot_status =
                    format!("UniProt linked GenBank fetch failed: {}", err.message);
            }
        }
    }

    pub(super) fn project_uniprot_entry_from_dialog(&mut self) {
        let entry_id = self.uniprot_entry_id.trim().to_string();
        if entry_id.is_empty() {
            self.uniprot_status = "entry_id is required for projection".to_string();
            return;
        }
        let seq_ids = self.project_sequence_ids_for_blast();
        if seq_ids.is_empty() {
            self.uniprot_status =
                "No sequence is loaded. Load/import sequence data before mapping.".to_string();
            return;
        }
        let seq_id = if self.uniprot_map_seq_id.trim().is_empty() {
            seq_ids[0].clone()
        } else {
            self.uniprot_map_seq_id.trim().to_string()
        };
        if !seq_ids.iter().any(|candidate| candidate == &seq_id) {
            self.uniprot_status = format!(
                "Sequence '{}' is not available anymore; pick an existing sequence",
                seq_id
            );
            return;
        }
        self.uniprot_map_seq_id = seq_id.clone();
        let resolved_projection_id =
            Self::uniprot_optional_trimmed(&self.uniprot_map_projection_id)
                .unwrap_or_else(|| format!("{entry_id}@{seq_id}"));
        let op = Operation::ProjectUniprotToGenome {
            seq_id,
            entry_id,
            projection_id: Some(resolved_projection_id.clone()),
            transcript_id: Self::uniprot_optional_trimmed(&self.uniprot_map_transcript_id),
        };
        match self.engine.write().unwrap().apply(op) {
            Ok(result) => {
                self.uniprot_map_projection_id = resolved_projection_id;
                self.uniprot_status = Self::format_op_result_status(
                    "UniProt projection: ok",
                    &result.created_seq_ids,
                    &result.warnings,
                    &result.messages,
                );
            }
            Err(err) => {
                self.uniprot_status = format!("UniProt projection failed: {}", err.message);
            }
        }
    }

    pub(super) fn default_uniprot_projection_svg_file_name(seq_id: &str, projection_id: &str) -> String {
        let stem = Self::sanitize_file_stem(
            &format!("{}_{}", seq_id.trim(), projection_id.trim()),
            "uniprot_projection",
        );
        format!("{stem}.protein_mapping.svg")
    }

    pub(super) fn default_transcript_protein_svg_file_name(
        seq_id: &str,
        transcript_id_filter: Option<&str>,
    ) -> String {
        let descriptor = transcript_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| format!("{}_{}", seq_id.trim(), value))
            .unwrap_or_else(|| format!("{}_derived_proteins", seq_id.trim()));
        let stem = Self::sanitize_file_stem(&descriptor, "protein_compare");
        format!("{stem}.protein_compare.svg")
    }

    pub(super) fn default_ensembl_protein_svg_file_name(
        seq_id: &str,
        transcript_id_filter: Option<&str>,
        entry_id: &str,
    ) -> String {
        let descriptor = transcript_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| format!("{}_{}_{}", seq_id.trim(), value, entry_id.trim()))
            .unwrap_or_else(|| format!("{}_{}", seq_id.trim(), entry_id.trim()));
        let stem = Self::sanitize_file_stem(&descriptor, "ensembl_protein_compare");
        format!("{stem}.protein_compare.svg")
    }

    pub(super) fn export_uniprot_projection_svg_from_dialog(&mut self, seq_id: &str, projection_id: &str) {
        let seq_id = seq_id.trim();
        let projection_id = projection_id.trim();
        if seq_id.is_empty() || projection_id.is_empty() {
            self.uniprot_status =
                "Project an entry first or provide seq_id + projection_id before exporting SVG"
                    .to_string();
            return;
        }
        let default_file_name =
            Self::default_uniprot_projection_svg_file_name(seq_id, projection_id);
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("SVG", &["svg"])
            .save_file()
        else {
            self.uniprot_status = "UniProt protein mapping SVG export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let protein_feature_filter = self.protein_feature_filter_from_dialog();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::RenderFeatureExpertSvg {
                seq_id: seq_id.to_string(),
                target: FeatureExpertTarget::UniprotProjection {
                    projection_id: projection_id.to_string(),
                    protein_feature_filter,
                },
                path: path_text.clone(),
            });
        match result {
            Ok(op_result) => {
                self.uniprot_status = op_result.messages.first().cloned().unwrap_or_else(|| {
                    format!("Wrote UniProt protein mapping SVG to '{path_text}'")
                });
            }
            Err(err) => {
                self.uniprot_status = format!(
                    "Could not export UniProt protein mapping SVG: {}",
                    err.message
                );
            }
        }
    }

    pub(super) fn export_transcript_protein_svg_from_dialog(
        &mut self,
        seq_id: &str,
        transcript_id_filter: Option<&str>,
    ) {
        let seq_id = seq_id.trim();
        if seq_id.is_empty() {
            self.uniprot_status =
                "seq_id is required before exporting transcript-native Protein Expert SVG"
                    .to_string();
            return;
        }
        let normalized_filter = transcript_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let default_file_name =
            Self::default_transcript_protein_svg_file_name(seq_id, normalized_filter.as_deref());
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("SVG", &["svg"])
            .save_file()
        else {
            self.uniprot_status =
                "Transcript-native Protein Expert SVG export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::RenderFeatureExpertSvg {
                seq_id: seq_id.to_string(),
                target: FeatureExpertTarget::ProteinComparison {
                    transcript_id_filter: normalized_filter.clone(),
                    protein_feature_filter: Default::default(),
                    external_source: None,
                    external_entry_id: None,
                },
                path: path_text.clone(),
            });
        match result {
            Ok(op_result) => {
                self.uniprot_status = op_result.messages.first().cloned().unwrap_or_else(|| {
                    if let Some(transcript_id_filter) = normalized_filter.as_deref() {
                        format!(
                            "Wrote transcript-native Protein Expert SVG for '{}' to '{}'",
                            transcript_id_filter, path_text
                        )
                    } else {
                        format!(
                            "Wrote transcript-native Protein Expert SVG to '{}'",
                            path_text
                        )
                    }
                });
            }
            Err(err) => {
                self.uniprot_status = format!(
                    "Could not export transcript-native Protein Expert SVG: {}",
                    err.message
                );
            }
        }
    }

    pub(super) fn export_ensembl_protein_svg_from_dialog(
        &mut self,
        seq_id: &str,
        transcript_id_filter: Option<&str>,
        entry_id: &str,
    ) {
        let seq_id = seq_id.trim();
        let entry_id = entry_id.trim();
        if seq_id.is_empty() {
            self.uniprot_status =
                "seq_id is required before exporting Ensembl Protein Expert SVG".to_string();
            return;
        }
        if entry_id.is_empty() {
            self.uniprot_status =
                "entry_id is required before exporting Ensembl Protein Expert SVG".to_string();
            return;
        }
        let normalized_filter = transcript_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let default_file_name = Self::default_ensembl_protein_svg_file_name(
            seq_id,
            normalized_filter.as_deref(),
            entry_id,
        );
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_file_name)
            .add_filter("SVG", &["svg"])
            .save_file()
        else {
            self.uniprot_status = "Ensembl Protein Expert SVG export canceled".to_string();
            return;
        };
        let path_text = path.display().to_string();
        let protein_feature_filter = self.protein_feature_filter_from_dialog();
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::RenderFeatureExpertSvg {
                seq_id: seq_id.to_string(),
                target: FeatureExpertTarget::ProteinComparison {
                    transcript_id_filter: normalized_filter.clone(),
                    protein_feature_filter,
                    external_source: Some(gentle_protocol::ProteinExternalOpinionSource::Ensembl),
                    external_entry_id: Some(entry_id.to_string()),
                },
                path: path_text.clone(),
            });
        match result {
            Ok(op_result) => {
                self.uniprot_status = op_result.messages.first().cloned().unwrap_or_else(|| {
                    if let Some(transcript_id_filter) = normalized_filter.as_deref() {
                        format!(
                            "Wrote Ensembl Protein Expert SVG for '{}' / '{}' to '{}'",
                            entry_id, transcript_id_filter, path_text
                        )
                    } else {
                        format!(
                            "Wrote Ensembl Protein Expert SVG for '{}' to '{}'",
                            entry_id, path_text
                        )
                    }
                });
            }
            Err(err) => {
                self.uniprot_status = format!(
                    "Could not export Ensembl Protein Expert SVG: {}",
                    err.message
                );
            }
        }
    }

    pub(super) fn uniprot_feature_speed_profile_label(
        profile: Option<TranslationSpeedProfile>,
    ) -> &'static str {
        match profile {
            None => "Auto",
            Some(TranslationSpeedProfile::Human) => "Human",
            Some(TranslationSpeedProfile::Mouse) => "Mouse",
            Some(TranslationSpeedProfile::Yeast) => "Yeast",
            Some(TranslationSpeedProfile::Ecoli) => "E. coli",
        }
    }

    pub(super) fn reverse_translation_speed_mark_label(mark: Option<TranslationSpeedMark>) -> &'static str {
        match mark {
            None => "Auto",
            Some(TranslationSpeedMark::Fast) => "Fast",
            Some(TranslationSpeedMark::Slow) => "Slow",
        }
    }

    pub(super) fn format_reverse_translation_speed_resolution_summary(
        report: &ReverseTranslationReport,
    ) -> String {
        let requested = report
            .requested_speed_profile
            .map(|profile| profile.as_str())
            .unwrap_or("auto");
        let resolved = report
            .resolved_speed_profile
            .map(|profile| profile.as_str())
            .unwrap_or("auto-unresolved");
        let source = report
            .resolved_speed_profile_source
            .map(|source| source.as_str())
            .unwrap_or("unspecified");
        let reference = report
            .translation_speed_reference_species
            .as_deref()
            .unwrap_or("-");
        format!(
            "requested={} | resolved={} | source={} | ref={}",
            requested, resolved, source, reference
        )
    }

    pub(super) fn format_reverse_translation_table_resolution_summary(
        report: &ReverseTranslationReport,
    ) -> String {
        let organism = report
            .translation_context_organism
            .as_deref()
            .unwrap_or("-");
        let organelle = report
            .translation_context_organelle
            .as_deref()
            .unwrap_or("-");
        format!(
            "{} ({}) | source={} | organism={} | organelle={}",
            report.translation_table,
            report.translation_table_label,
            report.translation_table_source,
            organism,
            organelle
        )
    }

    pub(super) fn protein_sequence_ids_for_dialog(&self) -> Vec<String> {
        let engine = self.engine.read().unwrap();
        let mut ids = engine
            .state()
            .sequences
            .iter()
            .filter_map(|(seq_id, dna)| dna.is_protein_sequence().then_some(seq_id.clone()))
            .collect::<Vec<_>>();
        ids.sort();
        ids
    }

    pub(super) fn format_uniprot_feature_speed_resolution_summary(
        report: &UniprotFeatureCodingDnaQueryReport,
    ) -> String {
        let requested = report
            .requested_translation_speed_profile
            .map(|profile| profile.as_str())
            .unwrap_or("auto");
        let resolved = report
            .resolved_translation_speed_profile
            .map(|profile| profile.as_str())
            .unwrap_or("auto-unresolved");
        let source = report
            .resolved_translation_speed_profile_source
            .map(|source| source.as_str())
            .unwrap_or("unspecified");
        let reference = report
            .resolved_translation_speed_reference_species
            .as_deref()
            .unwrap_or("-");
        format!(
            "requested={} | resolved={} | source={} | ref={}",
            requested, resolved, source, reference
        )
    }

    pub(super) fn reverse_translate_protein_from_dialog(&mut self) {
        let seq_id = self.reverse_translate_protein_seq_id.trim().to_string();
        if seq_id.is_empty() {
            self.uniprot_status =
                "Choose a protein sequence before reverse translation".to_string();
            return;
        }
        let translation_table =
            Self::uniprot_optional_trimmed(&self.reverse_translate_translation_table)
                .map(|raw| {
                    raw.parse::<usize>().map_err(|_| {
                        format!(
                            "translation table '{}' is not a valid positive integer",
                            raw
                        )
                    })
                })
                .transpose();
        let translation_table = match translation_table {
            Ok(value) => value.filter(|value| *value > 0),
            Err(message) => {
                self.uniprot_status = message;
                return;
            }
        };
        let target_anneal_tm_c =
            Self::uniprot_optional_trimmed(&self.reverse_translate_target_anneal_tm_c)
                .map(|raw| {
                    raw.parse::<f64>()
                        .map_err(|_| format!("target anneal Tm '{}' is not a valid number", raw))
                })
                .transpose();
        let target_anneal_tm_c = match target_anneal_tm_c {
            Ok(value) => value,
            Err(message) => {
                self.uniprot_status = message;
                return;
            }
        };
        let anneal_window_bp =
            Self::uniprot_optional_trimmed(&self.reverse_translate_anneal_window_bp)
                .map(|raw| {
                    raw.parse::<usize>().map_err(|_| {
                        format!("anneal window '{}' is not a valid positive integer", raw)
                    })
                })
                .transpose();
        let anneal_window_bp = match anneal_window_bp {
            Ok(value) => value.filter(|value| *value > 0),
            Err(message) => {
                self.uniprot_status = message;
                return;
            }
        };
        let result =
            self.engine
                .write()
                .unwrap()
                .apply(Operation::ReverseTranslateProteinSequence {
                    seq_id,
                    output_id: Self::uniprot_optional_trimmed(&self.reverse_translate_output_id),
                    speed_profile: self.reverse_translate_speed_profile,
                    speed_mark: self.reverse_translate_speed_mark,
                    translation_table,
                    target_anneal_tm_c,
                    anneal_window_bp,
                });
        match result {
            Ok(result) => {
                for seq_id in &result.created_seq_ids {
                    self.open_sequence_window(seq_id);
                }
                let status = if let Some(report) = result.reverse_translation_report.as_ref() {
                    self.reverse_translation_report = Some(report.clone());
                    format!(
                        "Reverse translation: ok\ncreated: {}\ntranslation table: {}\nspeed profile: {}\nspeed mark: {}\nanneal heuristic: {}\nwarnings: {}\nmessages: {}",
                        report.coding_seq_id,
                        Self::format_reverse_translation_table_resolution_summary(report),
                        Self::format_reverse_translation_speed_resolution_summary(report),
                        report
                            .speed_mark
                            .map(|mark| mark.as_str())
                            .unwrap_or("auto"),
                        report
                            .target_anneal_tm_c
                            .map(|tm| format!("{tm:.1} °C / {} bp", report.anneal_window_bp))
                            .unwrap_or_else(|| "-".to_string()),
                        if result.warnings.is_empty() {
                            "-".to_string()
                        } else {
                            result.warnings.join(" | ")
                        },
                        if result.messages.is_empty() {
                            "-".to_string()
                        } else {
                            result.messages.join(" | ")
                        }
                    )
                } else {
                    self.reverse_translation_report = None;
                    Self::format_op_result_status(
                        "Reverse translation: ok",
                        &result.created_seq_ids,
                        &result.warnings,
                        &result.messages,
                    )
                };
                self.uniprot_status = status;
            }
            Err(err) => {
                self.reverse_translation_report = None;
                self.uniprot_status = format!("Reverse translation failed: {}", err.message);
            }
        }
    }

    pub(super) fn protease_digest_names_from_dialog(&self) -> Vec<String> {
        self.protease_digest_names
            .split(',')
            .map(str::trim)
            .filter(|name| !name.is_empty())
            .map(ToOwned::to_owned)
            .collect()
    }

    pub(super) fn digest_selected_protein_from_dialog(&mut self) {
        let seq_id = self.reverse_translate_protein_seq_id.trim().to_string();
        if seq_id.is_empty() {
            self.uniprot_status = "Choose a protein sequence before protease digest".to_string();
            return;
        }
        let proteases = self.protease_digest_names_from_dialog();
        if proteases.is_empty() {
            self.uniprot_status =
                "Enter at least one protease name, for example Trypsin".to_string();
            return;
        }
        let min_length_aa = Self::uniprot_optional_trimmed(&self.protease_digest_min_length_aa)
            .map(|raw| {
                raw.parse::<usize>()
                    .map_err(|_| format!("minimum peptide length '{}' is not a valid integer", raw))
            })
            .transpose();
        let min_length_aa = match min_length_aa {
            Ok(value) => value.map(|value| value.max(1)),
            Err(message) => {
                self.uniprot_status = message;
                return;
            }
        };
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::ProteaseDigestProteinSequence {
                seq_id,
                proteases,
                output_prefix: Self::uniprot_optional_trimmed(&self.protease_digest_output_prefix),
                min_length_aa,
                materialize: self.protease_digest_materialize,
            });
        match result {
            Ok(result) => {
                for seq_id in &result.created_seq_ids {
                    self.open_sequence_window(seq_id);
                }
                if let Some(report) = result.protease_digest_report.as_ref() {
                    self.protease_digest_report = Some(report.clone());
                    self.uniprot_status = format!(
                        "Protease digest: ok\nsource: {}\nproteases: {}\ncleavage sites: {}\npeptides: {}\ncreated: {}\nwarnings: {}\nmessages: {}",
                        report.source_seq_id,
                        report
                            .resolved_proteases
                            .iter()
                            .map(|protease| protease.name.as_str())
                            .collect::<Vec<_>>()
                            .join(", "),
                        report.cleavage_site_count,
                        report.peptide_count,
                        if report.created_seq_ids.is_empty() {
                            "-".to_string()
                        } else {
                            report.created_seq_ids.join(", ")
                        },
                        if result.warnings.is_empty() {
                            "-".to_string()
                        } else {
                            result.warnings.join(" | ")
                        },
                        if result.messages.is_empty() {
                            "-".to_string()
                        } else {
                            result.messages.join(" | ")
                        }
                    );
                } else {
                    self.protease_digest_report = None;
                    self.uniprot_status = Self::format_op_result_status(
                        "Protease digest: ok",
                        &result.created_seq_ids,
                        &result.warnings,
                        &result.messages,
                    );
                }
            }
            Err(err) => {
                self.protease_digest_report = None;
                self.uniprot_status = format!("Protease digest failed: {}", err.message);
            }
        }
    }

    pub(super) fn protein_handoff_transcript_filter_from_dialog(&self) -> Option<String> {
        Self::uniprot_optional_trimmed(&self.uniprot_feature_transcript_id)
            .or_else(|| Self::uniprot_optional_trimmed(&self.uniprot_map_transcript_id))
    }

    pub(super) fn sync_protein_handoff_candidate_selection(&mut self) {
        let Some(graph) = self.protein_handoff_graph.as_ref() else {
            self.protein_handoff_selected_candidate_id.clear();
            return;
        };
        let current = self.protein_handoff_selected_candidate_id.trim();
        if current.is_empty()
            || !graph
                .candidates
                .iter()
                .any(|candidate| candidate.candidate_id == current)
        {
            self.protein_handoff_selected_candidate_id = graph
                .candidates
                .first()
                .map(|candidate| candidate.candidate_id.clone())
                .unwrap_or_default();
        }
    }

    pub(super) fn build_protein_to_dna_handoff_reasoning_from_dialog(&mut self) {
        let seq_id = self.uniprot_map_seq_id.trim().to_string();
        if seq_id.is_empty() {
            self.protein_handoff_status =
                "Choose a target sequence before building protein-to-DNA handoff reasoning"
                    .to_string();
            return;
        }
        let protein_seq_id = self.reverse_translate_protein_seq_id.trim().to_string();
        if protein_seq_id.is_empty() {
            self.protein_handoff_status =
                "Choose a protein sequence before building protein-to-DNA handoff reasoning"
                    .to_string();
            return;
        }
        let translation_table =
            Self::uniprot_optional_trimmed(&self.reverse_translate_translation_table)
                .map(|raw| {
                    raw.parse::<usize>().map_err(|_| {
                        format!(
                            "translation table '{}' is not a valid positive integer",
                            raw
                        )
                    })
                })
                .transpose();
        let translation_table = match translation_table {
            Ok(value) => value.filter(|value| *value > 0),
            Err(message) => {
                self.protein_handoff_status = message;
                return;
            }
        };
        let target_anneal_tm_c =
            Self::uniprot_optional_trimmed(&self.reverse_translate_target_anneal_tm_c)
                .map(|raw| {
                    raw.parse::<f64>()
                        .map_err(|_| format!("target anneal Tm '{}' is not a valid number", raw))
                })
                .transpose();
        let target_anneal_tm_c = match target_anneal_tm_c {
            Ok(value) => value,
            Err(message) => {
                self.protein_handoff_status = message;
                return;
            }
        };
        let anneal_window_bp =
            Self::uniprot_optional_trimmed(&self.reverse_translate_anneal_window_bp)
                .map(|raw| {
                    raw.parse::<usize>().map_err(|_| {
                        format!("anneal window '{}' is not a valid positive integer", raw)
                    })
                })
                .transpose();
        let anneal_window_bp = match anneal_window_bp {
            Ok(value) => value.filter(|value| *value > 0),
            Err(message) => {
                self.protein_handoff_status = message;
                return;
            }
        };
        let result =
            self.engine
                .write()
                .unwrap()
                .apply(Operation::BuildProteinToDnaHandoffReasoning {
                    seq_id: seq_id.clone(),
                    protein_seq_id: protein_seq_id.clone(),
                    transcript_filter: self.protein_handoff_transcript_filter_from_dialog(),
                    projection_id: self.resolve_uniprot_projection_id_from_dialog_fields(),
                    ensembl_entry_id: Self::uniprot_optional_trimmed(
                        &self.ensembl_protein_entry_id,
                    ),
                    feature_query: Self::uniprot_optional_trimmed(&self.uniprot_feature_query),
                    ranking_goal: self.protein_handoff_ranking_goal,
                    speed_profile: self.reverse_translate_speed_profile,
                    speed_mark: self.reverse_translate_speed_mark,
                    translation_table,
                    target_anneal_tm_c,
                    anneal_window_bp,
                    objective_id: None,
                    graph_id: self
                        .protein_handoff_graph
                        .as_ref()
                        .map(|graph| graph.graph_id.clone()),
                });
        match result {
            Ok(result) => {
                if let Some(graph) = result.construct_reasoning_graph.as_ref() {
                    self.protein_handoff_graph = Some((**graph).clone());
                    self.sync_protein_handoff_candidate_selection();
                    let handoff_candidate_count = graph
                        .candidates
                        .iter()
                        .filter(|candidate| candidate.protein_to_dna_handoff.is_some())
                        .count();
                    self.protein_handoff_status = format!(
                        "Protein-to-DNA handoff reasoning: ok (graph='{}', candidates={}, warnings={}, messages={})",
                        graph.graph_id,
                        handoff_candidate_count,
                        result.warnings.len(),
                        result.messages.len()
                    );
                } else {
                    self.protein_handoff_graph = None;
                    self.protein_handoff_selected_candidate_id.clear();
                    self.protein_handoff_status = Self::format_op_result_status(
                        "Protein-to-DNA handoff reasoning: ok",
                        &result.created_seq_ids,
                        &result.warnings,
                        &result.messages,
                    );
                }
            }
            Err(err) => {
                self.protein_handoff_graph = None;
                self.protein_handoff_selected_candidate_id.clear();
                self.protein_handoff_status =
                    format!("Protein-to-DNA handoff reasoning failed: {}", err.message);
            }
        }
    }

    pub(super) fn render_protein_handoff_graph(&mut self, ui: &mut Ui) {
        self.sync_protein_handoff_candidate_selection();
        let Some(graph) = self.protein_handoff_graph.clone() else {
            return;
        };
        ui.separator();
        ui.label("DNA handoff candidates");
        ui.small(format!(
            "graph={} | evidence={} | decisions={} | candidates={}",
            graph.graph_id,
            graph.evidence.len(),
            graph.decisions.len(),
            graph.candidates.len()
        ));
        egui::Grid::new("protein_handoff_candidate_grid")
            .striped(true)
            .show(ui, |ui| {
                ui.strong("rank");
                ui.strong("candidate");
                ui.strong("strategy");
                ui.strong("coverage");
                ui.strong("provenance");
                ui.strong("warnings");
                ui.strong("action");
                ui.end_row();
                for (idx, candidate) in graph
                    .candidates
                    .iter()
                    .filter(|candidate| candidate.protein_to_dna_handoff.is_some())
                    .enumerate()
                {
                    let detail = candidate
                        .protein_to_dna_handoff
                        .as_ref()
                        .expect("handoff detail");
                    let selected =
                        self.protein_handoff_selected_candidate_id == candidate.candidate_id;
                    ui.label((idx + 1).to_string());
                    ui.label(&candidate.title);
                    ui.monospace(detail.strategy.as_str());
                    ui.monospace(format!(
                        "{} / {} aa",
                        detail.coverage.covered_amino_acids, detail.coverage.requested_amino_acids
                    ));
                    ui.monospace(
                        detail
                            .provenance_score
                            .map(|value| format!("{value:.2}"))
                            .unwrap_or_else(|| "-".to_string()),
                    );
                    ui.label(candidate.warnings.len().to_string());
                    if ui
                        .selectable_label(selected, if selected { "Selected" } else { "Select" })
                        .clicked()
                    {
                        self.protein_handoff_selected_candidate_id = candidate.candidate_id.clone();
                    }
                    ui.end_row();
                }
            });
        if let Some(candidate) = graph
            .candidates
            .iter()
            .find(|candidate| candidate.candidate_id == self.protein_handoff_selected_candidate_id)
        {
            let detail = candidate
                .protein_to_dna_handoff
                .as_ref()
                .expect("handoff detail");
            ui.separator();
            ui.label("Selected handoff candidate");
            egui::Grid::new("protein_handoff_candidate_detail")
                .num_columns(2)
                .spacing([12.0, 4.0])
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("title");
                    ui.label(&candidate.title);
                    ui.end_row();
                    ui.strong("strategy");
                    ui.monospace(detail.strategy.as_str());
                    ui.end_row();
                    ui.strong("protein");
                    ui.monospace(&detail.source_protein_seq_id);
                    ui.end_row();
                    ui.strong("transcript");
                    ui.monospace(detail.transcript_id.as_deref().unwrap_or("-"));
                    ui.end_row();
                    ui.strong("projection");
                    ui.monospace(detail.projection_id.as_deref().unwrap_or("-"));
                    ui.end_row();
                    ui.strong("ensembl entry");
                    ui.monospace(detail.ensembl_entry_id.as_deref().unwrap_or("-"));
                    ui.end_row();
                    ui.strong("feature query");
                    ui.label(detail.feature_query.as_deref().unwrap_or("-"));
                    ui.end_row();
                    ui.strong("feature match");
                    ui.label(detail.matched_feature_label.as_deref().unwrap_or("-"));
                    ui.end_row();
                    ui.strong("coverage");
                    ui.monospace(format!(
                        "{} / {} aa",
                        detail.coverage.covered_amino_acids, detail.coverage.requested_amino_acids
                    ));
                    ui.end_row();
                    ui.strong("translation table");
                    ui.monospace(
                        detail
                            .translation_table
                            .map(|value| value.to_string())
                            .unwrap_or_else(|| "-".to_string()),
                    );
                    ui.end_row();
                    ui.strong("speed profile");
                    ui.monospace(
                        detail
                            .speed_profile
                            .map(TranslationSpeedProfile::as_str)
                            .unwrap_or("-"),
                    );
                    ui.end_row();
                    ui.strong("speed mark");
                    ui.monospace(
                        detail
                            .speed_mark
                            .map(TranslationSpeedMark::as_str)
                            .unwrap_or("-"),
                    );
                    ui.end_row();
                    ui.strong("codon policy");
                    ui.label(&detail.codon_policy_summary);
                    ui.end_row();
                });
            if !detail.preserved_constraints.is_empty() {
                ui.small(format!(
                    "preserved: {}",
                    detail.preserved_constraints.join(" | ")
                ));
            }
            if !detail.relaxed_constraints.is_empty() {
                ui.small(format!(
                    "relaxed: {}",
                    detail.relaxed_constraints.join(" | ")
                ));
            }
            if !detail.next_step_recommendations.is_empty() {
                ui.small(format!(
                    "next: {}",
                    detail.next_step_recommendations.join(" | ")
                ));
            }
            if !candidate.warnings.is_empty() {
                ui.small(format!("warnings: {}", candidate.warnings.join(" | ")));
            }
            if !candidate.notes.is_empty() {
                ui.small(format!("notes: {}", candidate.notes.join(" | ")));
            }
        }
    }

    pub(super) fn query_uniprot_feature_coding_dna_from_dialog(&mut self) {
        let Some(projection_id) = self.resolve_uniprot_projection_id_from_dialog_fields() else {
            self.uniprot_status = "Project an entry first or provide seq_id + entry_id to resolve a projection before querying coding DNA.".to_string();
            return;
        };
        let feature_query = self.uniprot_feature_query.trim().to_string();
        if feature_query.is_empty() {
            self.uniprot_status = "Protein feature query cannot be empty".to_string();
            return;
        }
        let transcript_filter = Self::uniprot_optional_trimmed(&self.uniprot_feature_transcript_id);
        let result = self
            .engine
            .read()
            .unwrap()
            .query_uniprot_feature_coding_dna(
                &projection_id,
                &feature_query,
                transcript_filter.as_deref(),
                self.uniprot_feature_query_mode,
                self.uniprot_feature_speed_profile,
            );
        match result {
            Ok(report) => {
                let speed_summary = Self::format_uniprot_feature_speed_resolution_summary(&report);
                self.uniprot_feature_report = Some(report);
                self.uniprot_status = format!(
                    "UniProt feature coding DNA query: ok (projection='{}', matches={}, mode={}, {})",
                    projection_id,
                    self.uniprot_feature_report
                        .as_ref()
                        .map(|report| report.match_count)
                        .unwrap_or(0),
                    self.uniprot_feature_query_mode.as_str(),
                    speed_summary,
                );
            }
            Err(err) => {
                self.uniprot_status =
                    format!("UniProt feature coding DNA query failed: {}", err.message);
            }
        }
    }

    pub(super) fn run_uniprot_projection_audit_from_dialog(&mut self) {
        let Some(projection_id) = self.resolve_uniprot_projection_id_from_dialog_fields() else {
            self.uniprot_status = "Project an entry first or provide seq_id + entry_id to resolve a projection before running the UniProt audit.".to_string();
            return;
        };
        let transcript_id = Self::uniprot_optional_trimmed(&self.uniprot_audit_transcript_id)
            .or_else(|| Self::uniprot_optional_trimmed(&self.uniprot_map_transcript_id));
        let report_id = Self::uniprot_optional_trimmed(&self.uniprot_audit_report_id);
        let ensembl_entry_id = Self::uniprot_optional_trimmed(&self.ensembl_protein_entry_id);
        let result =
            self.engine
                .write()
                .unwrap()
                .apply(Operation::AuditUniprotProjectionConsistency {
                    projection_id,
                    transcript_id,
                    report_id,
                    ensembl_entry_id,
                });
        match result {
            Ok(result) => {
                self.uniprot_audit_report = result.uniprot_projection_audit.as_deref().cloned();
                if let Some(report) = self.uniprot_audit_report.as_ref() {
                    self.uniprot_audit_report_id = report.report_id.clone();
                    self.uniprot_audit_transcript_id =
                        report.transcript_id_filter.clone().unwrap_or_default();
                }
                self.uniprot_status = Self::format_op_result_status(
                    "UniProt audit: ok",
                    &result.created_seq_ids,
                    &result.warnings,
                    &result.messages,
                );
            }
            Err(err) => {
                self.uniprot_status = format!("UniProt audit failed: {}", err.message);
            }
        }
    }

    pub(super) fn run_uniprot_projection_audit_parity_from_dialog(&mut self) {
        let Some(projection_id) = self.resolve_uniprot_projection_id_from_dialog_fields() else {
            self.uniprot_status = "Project an entry first or provide seq_id + entry_id to resolve a projection before running the UniProt audit parity report.".to_string();
            return;
        };
        let transcript_id = Self::uniprot_optional_trimmed(&self.uniprot_audit_transcript_id)
            .or_else(|| Self::uniprot_optional_trimmed(&self.uniprot_map_transcript_id));
        let report_id = Self::uniprot_optional_trimmed(&self.uniprot_audit_parity_report_id);
        let ensembl_entry_id = Self::uniprot_optional_trimmed(&self.ensembl_protein_entry_id);
        let result = self
            .engine
            .write()
            .unwrap()
            .apply(Operation::AuditUniprotProjectionParity {
                projection_id,
                transcript_id,
                report_id,
                ensembl_entry_id,
            });
        match result {
            Ok(result) => {
                self.uniprot_audit_parity_report =
                    result.uniprot_projection_audit_parity.as_deref().cloned();
                if let Some(report) = self.uniprot_audit_parity_report.as_ref() {
                    self.uniprot_audit_parity_report_id = report.report_id.clone();
                }
                self.uniprot_status = Self::format_op_result_status(
                    "UniProt audit parity: ok",
                    &result.created_seq_ids,
                    &result.warnings,
                    &result.messages,
                );
            }
            Err(err) => {
                self.uniprot_status = format!("UniProt audit parity failed: {}", err.message);
            }
        }
    }

    pub(super) fn load_uniprot_audit_report_from_dialog(&mut self, report_id: &str) {
        match self
            .engine
            .read()
            .unwrap()
            .get_uniprot_projection_audit_report(report_id)
        {
            Ok(report) => {
                self.uniprot_audit_report_id = report.report_id.clone();
                self.uniprot_map_seq_id = report.seq_id.clone();
                self.uniprot_entry_id = report.entry_id.clone();
                self.uniprot_audit_transcript_id =
                    report.transcript_id_filter.clone().unwrap_or_default();
                self.uniprot_audit_report = Some(report);
                self.uniprot_status = format!("Loaded stored UniProt audit report '{}'", report_id);
            }
            Err(err) => {
                self.uniprot_status = format!(
                    "Could not load UniProt audit report '{}': {}",
                    report_id, err.message
                );
            }
        }
    }

    pub(super) fn load_uniprot_audit_parity_report_from_dialog(&mut self, report_id: &str) {
        match self
            .engine
            .read()
            .unwrap()
            .get_uniprot_projection_audit_parity_report(report_id)
        {
            Ok(report) => {
                self.uniprot_audit_parity_report_id = report.report_id.clone();
                self.uniprot_map_seq_id = report.seq_id.clone();
                self.uniprot_entry_id = report.entry_id.clone();
                self.uniprot_audit_transcript_id =
                    report.transcript_id_filter.clone().unwrap_or_default();
                self.uniprot_audit_parity_report = Some(report);
                self.uniprot_status =
                    format!("Loaded stored UniProt audit parity report '{}'", report_id);
            }
            Err(err) => {
                self.uniprot_status = format!(
                    "Could not load UniProt audit parity report '{}': {}",
                    report_id, err.message
                );
            }
        }
    }

    pub(super) fn render_uniprot_feature_coding_report(&self, ui: &mut Ui) {
        let Some(report) = self.uniprot_feature_report.as_ref() else {
            return;
        };
        ui.separator();
        ui.label(format!(
            "Feature coding DNA matches ({})",
            report.match_count
        ));
        egui::Grid::new("uniprot_feature_coding_report_summary")
            .num_columns(2)
            .spacing([12.0, 4.0])
            .striped(true)
            .show(ui, |ui| {
                ui.strong("projection");
                ui.monospace(&report.projection_id);
                ui.end_row();
                ui.strong("mode");
                ui.monospace(report.query_mode.as_str());
                ui.end_row();
                ui.strong("speed profile");
                ui.monospace(Self::format_uniprot_feature_speed_resolution_summary(
                    report,
                ));
                ui.end_row();
                ui.strong("optimized DNA");
                ui.label(
                    if report
                        .matches
                        .iter()
                        .any(|row| row.translation_speed_optimized_dna.is_some())
                    {
                        "available for at least one match"
                    } else {
                        "not available"
                    },
                );
                ui.end_row();
            });
        if !report.warnings.is_empty() {
            ui.small(report.warnings.join(" | "));
        }
        egui::ScrollArea::vertical()
            .max_height(240.0)
            .show(ui, |ui| {
                scroll_input_policy::apply_scrollarea_keyboard_navigation(
                    ui,
                    scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                );
                for (idx, row) in report.matches.iter().enumerate() {
                    let header = format!(
                        "{}. {} {} aa {}..{}",
                        idx + 1,
                        row.transcript_id,
                        row.feature_key,
                        row.aa_start,
                        row.aa_end
                    );
                    ui.collapsing(header, |ui| {
                        if let Some(note) = row.feature_note.as_deref() {
                            ui.label(format!("note: {note}"));
                        }
                        ui.small(format!(
                            "matched in {} | strand={} | feature fields={}",
                            row.transcript_id,
                            row.strand,
                            if row.matched_feature_fields.is_empty() {
                                "-".to_string()
                            } else {
                                row.matched_feature_fields.join(", ")
                            }
                        ));
                        ui.small(format!("AA sequence: {}", row.amino_acid_sequence));
                        let exon_summary = if let Some(exon) = row.primary_exon_ordinal {
                            format!("exon {exon}")
                        } else if let Some(pair) = row.primary_exon_pair.as_ref() {
                            format!(
                                "exon pair {} -> {}",
                                pair.from_exon_ordinal, pair.to_exon_ordinal
                            )
                        } else if row.exon_spans.is_empty() {
                            "exon attribution unavailable".to_string()
                        } else {
                            format!(
                                "exons {}",
                                row.exon_spans
                                    .iter()
                                    .map(|span| span.exon_ordinal.to_string())
                                    .collect::<Vec<_>>()
                                    .join(", ")
                            )
                        };
                        ui.small(exon_summary);
                        if !row.exon_spans.is_empty() {
                            ui.small(format!(
                                "coding spans: {}",
                                row.exon_spans
                                    .iter()
                                    .map(|span| format!(
                                        "exon {} [{}..{}] coding {}..{}",
                                        span.exon_ordinal,
                                        span.exon_start_1based,
                                        span.exon_end_1based,
                                        span.coding_start_1based,
                                        span.coding_end_1based
                                    ))
                                    .collect::<Vec<_>>()
                                    .join(" | ")
                            ));
                        }
                        if !row.genomic_segments.is_empty() {
                            ui.small(format!(
                                "genomic segments: {}",
                                row.genomic_segments
                                    .iter()
                                    .map(|segment| format!(
                                        "{}..{} ({}) aa {}..{}",
                                        segment.genomic_start_1based,
                                        segment.genomic_end_1based,
                                        segment.strand,
                                        segment.aa_start,
                                        segment.aa_end
                                    ))
                                    .collect::<Vec<_>>()
                                    .join(" | ")
                            ));
                        }
                        ui.label("Genomic coding DNA");
                        let mut genomic_text = row.genomic_coding_dna.clone();
                        ui.add(
                            egui::TextEdit::multiline(&mut genomic_text)
                                .desired_width(f32::INFINITY)
                                .font(egui::TextStyle::Monospace)
                                .interactive(false),
                        );
                        if let Some(optimized) = row.translation_speed_optimized_dna.as_ref() {
                            ui.label("Translation-speed optimized DNA");
                            let mut optimized_text = optimized.clone();
                            ui.add(
                                egui::TextEdit::multiline(&mut optimized_text)
                                    .desired_width(f32::INFINITY)
                                    .font(egui::TextStyle::Monospace)
                                    .interactive(false),
                            );
                        }
                        if !row.warnings.is_empty() {
                            ui.small(format!("warnings: {}", row.warnings.join(" | ")));
                        }
                    });
                }
            });
    }

    pub(super) fn render_reverse_translation_report(&self, ui: &mut Ui) {
        let Some(report) = self.reverse_translation_report.as_ref() else {
            return;
        };
        ui.separator();
        ui.label("Reverse-translated coding DNA");
        egui::Grid::new("reverse_translation_report_summary")
            .num_columns(2)
            .spacing([12.0, 4.0])
            .striped(true)
            .show(ui, |ui| {
                ui.strong("protein");
                ui.monospace(&report.protein_seq_id);
                ui.end_row();
                ui.strong("output");
                ui.monospace(&report.coding_seq_id);
                ui.end_row();
                ui.strong("length");
                ui.monospace(format!(
                    "{} aa -> {} bp",
                    report.protein_length_aa, report.coding_length_bp
                ));
                ui.end_row();
                ui.strong("translation table");
                ui.monospace(Self::format_reverse_translation_table_resolution_summary(
                    report,
                ));
                ui.end_row();
                ui.strong("speed profile");
                ui.monospace(Self::format_reverse_translation_speed_resolution_summary(
                    report,
                ));
                ui.end_row();
                ui.strong("speed mark");
                ui.monospace(
                    report
                        .speed_mark
                        .map(|mark| mark.as_str())
                        .unwrap_or("auto"),
                );
                ui.end_row();
                ui.strong("anneal heuristic");
                ui.monospace(
                    report
                        .target_anneal_tm_c
                        .map(|tm| format!("{tm:.1} °C / {} bp", report.anneal_window_bp))
                        .unwrap_or_else(|| "-".to_string()),
                );
                ui.end_row();
            });
        if let Some(sequence) = self
            .engine
            .read()
            .unwrap()
            .state()
            .sequences
            .get(&report.coding_seq_id)
        {
            ui.label("Coding DNA");
            let mut dna_text = sequence.get_forward_string();
            ui.add(
                egui::TextEdit::multiline(&mut dna_text)
                    .desired_width(f32::INFINITY)
                    .font(egui::TextStyle::Monospace)
                    .interactive(false),
            );
        }
        if !report.warnings.is_empty() {
            ui.small(format!("warnings: {}", report.warnings.join(" | ")));
        }
    }

    pub(super) fn render_protease_digest_report(&self, ui: &mut Ui) {
        let Some(report) = self.protease_digest_report.as_ref() else {
            return;
        };
        ui.separator();
        ui.label("Protease digest report");
        egui::Grid::new("protease_digest_report_summary")
            .num_columns(2)
            .spacing([12.0, 4.0])
            .striped(true)
            .show(ui, |ui| {
                ui.strong("source");
                ui.monospace(format!(
                    "{} ({} aa)",
                    report.source_seq_id, report.source_length_aa
                ));
                ui.end_row();
                ui.strong("proteases");
                ui.monospace(
                    report
                        .resolved_proteases
                        .iter()
                        .map(|protease| protease.name.as_str())
                        .collect::<Vec<_>>()
                        .join(", "),
                );
                ui.end_row();
                ui.strong("sites / peptides");
                ui.monospace(format!(
                    "{} cleavage site(s), {} peptide(s)",
                    report.cleavage_site_count, report.peptide_count
                ));
                ui.end_row();
                ui.strong("materialized");
                ui.monospace(if report.materialized {
                    format!("yes ({})", report.created_seq_ids.len())
                } else {
                    "prediction only".to_string()
                });
                ui.end_row();
                if let Some(transcript_id) = report.source_transcript_id.as_deref() {
                    ui.strong("transcript");
                    ui.monospace(transcript_id);
                    ui.end_row();
                }
            });
        if !report.sites.is_empty() {
            let site_preview = report
                .sites
                .iter()
                .take(16)
                .map(|site| {
                    format!(
                        "{}@{} ({})",
                        site.protease_name, site.cleavage_after_aa_1based, site.context
                    )
                })
                .collect::<Vec<_>>()
                .join(" | ");
            ui.small(format!(
                "cleavage preview: {}{}",
                site_preview,
                if report.sites.len() > 16 {
                    " | ..."
                } else {
                    ""
                }
            ));
        }
        if !report.peptides.is_empty() {
            ui.label("Peptide products");
            egui::Grid::new("protease_digest_peptide_preview")
                .num_columns(5)
                .spacing([10.0, 3.0])
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("#");
                    ui.strong("aa span");
                    ui.strong("len");
                    ui.strong("seq_id");
                    ui.strong("sequence");
                    ui.end_row();
                    for peptide in report.peptides.iter().take(12) {
                        ui.monospace(peptide.peptide_index.to_string());
                        ui.monospace(format!(
                            "{}..{}",
                            peptide.source_start_aa_1based, peptide.source_end_aa_1based
                        ));
                        ui.monospace(peptide.length_aa.to_string());
                        ui.monospace(peptide.created_seq_id.as_deref().unwrap_or("-"));
                        let mut peptide_sequence = peptide.sequence.clone();
                        if peptide_sequence.len() > 42 {
                            peptide_sequence.truncate(42);
                            peptide_sequence.push_str("...");
                        }
                        ui.monospace(peptide_sequence);
                        ui.end_row();
                    }
                });
            if report.peptides.len() > 12 {
                ui.small(format!(
                    "showing first 12 of {} peptide products",
                    report.peptides.len()
                ));
            }
        }
    }

    pub(super) fn render_genbank_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_genbank_dialog {
            return;
        }
        const WINDOW_TITLE: &str = "GenBank / dbSNP Fetch";
        let viewport_id = Self::genbank_viewport_id();
        let mut open = self.show_genbank_dialog;
        let mut close_requested = false;
        let spec = self.hosted_window_spec_for_viewport(
            WINDOW_TITLE,
            egui::Id::new(("hosted_genbank_window", viewport_id)),
            viewport_id,
            Vec2::new(820.0, 420.0),
            Vec2::new(680.0, 320.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            let close_hover = Self::specialist_window_close_hover_text(WINDOW_TITLE);
            if self
                .render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str())))
            {
                close_requested = true;
            }
            ui.label(
                    "Fetch one GenBank accession as a project sequence, or resolve a dbSNP rsID into an annotated genomic region from a prepared reference genome.",
                );
            ui.separator();
            ui.label("GenBank accession");
            ui.horizontal(|ui| {
                ui.label("accession");
                ui.text_edit_singleline(&mut self.genbank_accession)
                    .on_hover_text("GenBank accession (for example AY738222 or NC_000001)");
                if ui
                    .button("Fetch")
                    .on_hover_text("Fetch accession and import sequence into the project")
                    .clicked()
                {
                    self.fetch_genbank_accession_from_dialog();
                }
            });
            ui.horizontal(|ui| {
                ui.label("as_id");
                ui.text_edit_singleline(&mut self.genbank_as_id)
                    .on_hover_text(
                        "Optional project sequence id override (auto-uses accession when empty)",
                    );
            });
            ui.small("Examples: AY738222, NC_000001");
            ui.small(
                    "If network fetch is unavailable, use File -> Open Sequence... with a local GenBank file.",
                );
            if !self.genbank_status.trim().is_empty() {
                ui.separator();
                ui.monospace(self.genbank_status.clone());
            }
            ui.separator();
            ui.label("dbSNP locus extraction");
            ui.small(
                    "Resolve one rsID through NCBI Variation, then extract +/- flank bp with full feature annotation from the selected prepared reference genome.",
                );
            let dbsnp_fetch_running = self.dbsnp_fetch_task.is_some();
            ui.add_enabled_ui(!dbsnp_fetch_running, |ui| {
                ui.horizontal(|ui| {
                    ui.label("rs_id");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.dbsnp_rs_id)
                            .hint_text(DEFAULT_DBSNP_TUTORIAL_RS_ID),
                    )
                    .on_hover_text("dbSNP identifier such as rs9923231");
                    ui.label("genome");
                    ui.text_edit_singleline(&mut self.dbsnp_genome_id)
                        .on_hover_text(
                            "Prepared reference genome ID from the reference genome catalog",
                        );
                    if ui
                        .button("Fetch Region")
                        .on_hover_text(
                            "Resolve the rsID and extract the annotated genomic interval",
                        )
                        .clicked()
                    {
                        self.fetch_dbsnp_region_from_dialog();
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("flank_bp");
                    ui.text_edit_singleline(&mut self.dbsnp_flank_bp)
                        .on_hover_text("Bases to include on each side of the SNP (default 3000)");
                    ui.label("output_id");
                    ui.text_edit_singleline(&mut self.dbsnp_output_id)
                        .on_hover_text(
                            "Optional project sequence id override (auto-generated when empty)",
                        );
                });
            });
            ui.small(format!(
                "Quick try: {} (VKORC1 warfarin-sensitivity locus)",
                DEFAULT_DBSNP_TUTORIAL_RS_ID
            ));
            let catalog_label = self
                .dbsnp_catalog_path_opt()
                .unwrap_or_else(|| DEFAULT_GENOME_CATALOG_PATH.to_string());
            let cache_label = self
                .dbsnp_cache_dir_opt()
                .unwrap_or_else(configured_reference_genome_cache_dir);
            ui.small(format!("Examples: {}, rs334", DEFAULT_DBSNP_TUTORIAL_RS_ID));
            ui.small(format!(
                "Uses reference genome catalog '{}' and cache '{}'.",
                catalog_label, cache_label
            ));
            if dbsnp_fetch_running {
                ui.horizontal(|ui| {
                    ui.add(egui::Spinner::new());
                    ui.label("dbSNP fetch in progress...");
                });
            }
            if !self.dbsnp_status.trim().is_empty() {
                ui.separator();
                ui.monospace(self.dbsnp_status.clone());
            }
        });
        if close_requested {
            open = false;
        }
        if Self::viewport_close_requested_or_shortcut(ctx) {
            open = false;
        }
        self.show_genbank_dialog = open;
    }

    pub(super) fn render_uniprot_dialog_contents(&mut self, ui: &mut Ui) -> bool {
        let mut close_requested = false;
        let seq_ids = self.project_sequence_ids_for_blast();
        let protein_seq_ids = self.protein_sequence_ids_for_dialog();
        if !seq_ids.is_empty() {
            let selected = self.uniprot_map_seq_id.trim();
            if selected.is_empty() || !seq_ids.iter().any(|candidate| candidate == selected) {
                self.uniprot_map_seq_id = seq_ids[0].clone();
            }
        }
        if !protein_seq_ids.is_empty() {
            let selected = self.reverse_translate_protein_seq_id.trim();
            if selected.is_empty()
                || !protein_seq_ids
                    .iter()
                    .any(|candidate| candidate == selected)
            {
                self.reverse_translate_protein_seq_id = protein_seq_ids[0].clone();
            }
        }

        let close_hover = Self::specialist_window_close_hover_text("Protein Evidence");
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            close_requested = true;
        }
        ui.label(
            "Fetch/import UniProt or Ensembl protein evidence and compare it against transcript-native protein derivation.",
        );
        ui.horizontal(|ui| {
            ui.label("entry_id");
            ui.text_edit_singleline(&mut self.uniprot_entry_id)
                .on_hover_text("Stable entry identifier in this project (used by projection)");
        });
        ui.small(
            "Optional override for fetch/import. Required for projection; if empty after fetch, query is used.",
        );
        ui.small(
            "UniProt evidence in this window stays compare/projection-oriented; the Ensembl section below can also import one standalone protein sequence into the project.",
        );
        ui.separator();
        ui.label("Online fetch");
        ui.horizontal(|ui| {
            ui.label("query");
            ui.text_edit_singleline(&mut self.uniprot_query)
                .on_hover_text(
                    "UniProt accession or entry ID to fetch (for example P04637 or P53_HUMAN)",
                );
            if ui
                .button("Fetch")
                .on_hover_text("Fetch UniProt SWISS-PROT text by accession or UniProt entry ID")
                .clicked()
            {
                self.fetch_uniprot_entry_from_dialog();
            }
        });
        ui.small("Examples: P04637, P53_HUMAN");
        ui.separator();
        ui.label("Offline import (SWISS-PROT text)");
        ui.horizontal(|ui| {
            ui.label("path");
            ui.text_edit_singleline(&mut self.uniprot_swiss_path)
                .on_hover_text("Local path to a SWISS-PROT text file (.txt/.dat)");
            if ui
                .button("Browse...")
                .on_hover_text("Pick a local SWISS-PROT text file")
                .clicked()
            {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("SWISS-PROT text", &["txt", "dat"])
                    .add_filter("Text", &["txt"])
                    .pick_file()
                {
                    self.uniprot_swiss_path = path.display().to_string();
                }
            }
            if ui
                .button("Import")
                .on_hover_text("Import local SWISS-PROT text into project metadata")
                .clicked()
            {
                self.import_uniprot_swiss_prot_from_dialog();
            }
        });
        ui.separator();
        ui.label("UniProt projection to sequence");
        if seq_ids.is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                "No project sequence available. Load/import one first.",
            );
            ui.horizontal(|ui| {
                ui.label("seq_id");
                ui.text_edit_singleline(&mut self.uniprot_map_seq_id)
                    .on_hover_text("Target sequence ID for projection (must exist in project)");
            });
        } else {
            egui::ComboBox::from_label("seq_id")
                .selected_text(self.uniprot_map_seq_id.clone())
                .show_ui(ui, |ui| {
                    for seq_id in &seq_ids {
                        ui.selectable_value(&mut self.uniprot_map_seq_id, seq_id.clone(), seq_id);
                    }
                })
                .response
                .on_hover_text("Select the target sequence for UniProt projection");
        }
        ui.horizontal(|ui| {
            ui.label("projection_id");
            ui.text_edit_singleline(&mut self.uniprot_map_projection_id)
                .on_hover_text("Optional projection record ID (auto-generated when empty)");
            ui.label("transcript");
            ui.text_edit_singleline(&mut self.uniprot_map_transcript_id)
                .on_hover_text("Optional transcript filter (for example ENST... )");
        });
        ui.horizontal(|ui| {
            if ui
                .button("Project To Sequence")
                .on_hover_text(
                    "Map UniProt feature intervals through transcript/CDS onto genomic coordinates",
                )
                .clicked()
            {
                self.project_uniprot_entry_from_dialog();
            }
            if ui
                .button("Open Protein Expert")
                .on_hover_text(
                    "Open the Protein Expert with transcript-native translation as primary and the stored UniProt projection as optional external evidence",
                )
                .clicked()
            {
                let seq_id = self.uniprot_map_seq_id.trim().to_string();
                if let Some(projection_id) = self.resolve_uniprot_projection_id_from_dialog_fields()
                {
                    self.open_uniprot_projection_expert_from_dialog(&seq_id, &projection_id);
                } else {
                    self.uniprot_status =
                        "Project an entry first or provide seq_id + entry_id to resolve a projection"
                            .to_string();
                }
            }
            if ui
                .button("Open Derived Protein Expert")
                .on_hover_text(
                    "Open the same Protein Expert without requiring any external UniProt projection; optional transcript filter uses the current transcript field",
                )
                .clicked()
            {
                let seq_id = self.uniprot_map_seq_id.trim().to_string();
                self.open_transcript_protein_expert_from_dialog(
                    &seq_id,
                    Self::uniprot_optional_trimmed(&self.uniprot_map_transcript_id).as_deref(),
                );
            }
            if ui
                .button("Render Derived Protein SVG...")
                .on_hover_text(
                    "Export the transcript-native Protein Expert directly as SVG without requiring any stored UniProt projection; optional transcript filter uses the current transcript field",
                )
                .clicked()
            {
                let seq_id = self.uniprot_map_seq_id.trim().to_string();
                let transcript_filter =
                    Self::uniprot_optional_trimmed(&self.uniprot_map_transcript_id);
                self.export_transcript_protein_svg_from_dialog(
                    &seq_id,
                    transcript_filter.as_deref(),
                );
            }
            if ui
                .button("Render Protein Mapping SVG...")
                .on_hover_text(
                    "Export the current stored UniProt projection directly as an SVG from the shared protein-mapping expert route",
                )
                .clicked()
            {
                let seq_id = self.uniprot_map_seq_id.trim().to_string();
                if let Some(projection_id) = self.resolve_uniprot_projection_id_from_dialog_fields()
                {
                    self.export_uniprot_projection_svg_from_dialog(&seq_id, &projection_id);
                } else {
                    self.uniprot_status =
                        "Project an entry first or provide seq_id + entry_id to resolve a projection"
                            .to_string();
                }
            }
        });
        ui.horizontal(|ui| {
            ui.label("feature include");
            ui.text_edit_singleline(&mut self.protein_feature_key_include)
                .on_hover_text(
                    "Optional external protein feature keys to include (comma or whitespace separated, for example DOMAIN DNA_BIND PF02196)",
                );
            ui.label("exclude");
            ui.text_edit_singleline(&mut self.protein_feature_key_exclude)
                .on_hover_text(
                    "Optional external protein feature keys to exclude (comma or whitespace separated, for example CONFLICT REGION)",
                );
        });
        ui.small(
            "Open Protein Expert uses transcript-native translation as the primary product model and treats UniProt as optional external evidence; Open Derived Protein Expert omits the external opinion entirely. Render Derived Protein SVG... exports the transcript-native view directly, while Render Protein Mapping SVG... exports the shared Protein Expert view for the stored UniProt projection. The Ensembl section below reuses that same compare window with fetched Ensembl protein entries as optional external evidence.",
        );
        ui.small(
            "The feature include/exclude fields above apply to external protein annotation overlays used by Open Protein Expert, Render Protein Mapping SVG..., Open Ensembl Protein Expert, and Render Ensembl Protein SVG.... Derived-only transcript Protein Expert views ignore them. CONFLICT stays hidden by default unless explicitly re-included.",
        );
        ui.separator();
        ui.label("Feature coding DNA query");
        ui.small(
            "Query one mapped UniProt feature and show the exact coding DNA as encoded in the genome, plus an optional preferred-codon translational-speed version.",
        );
        ui.horizontal(|ui| {
            ui.label("feature");
            ui.text_edit_singleline(&mut self.uniprot_feature_query).on_hover_text(
                "Case-insensitive substring against mapped feature key/note (for example DOMAIN, DNA-binding, activation)",
            );
            if ui
                .button("Query Coding DNA")
                .on_hover_text(
                    "Look up coding DNA and exon attribution for the selected mapped UniProt feature",
                )
                .clicked()
            {
                self.query_uniprot_feature_coding_dna_from_dialog();
            }
        });
        ui.horizontal(|ui| {
            ui.label("feature transcript");
            ui.text_edit_singleline(&mut self.uniprot_feature_transcript_id)
                .on_hover_text(
                    "Optional stored projection transcript filter for feature DNA lookup",
                );
            egui::ComboBox::from_label("mode")
                .selected_text(self.uniprot_feature_query_mode.as_str())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.uniprot_feature_query_mode,
                        UniprotFeatureCodingDnaQueryMode::GenomicAsEncoded,
                        "genomic_as_encoded",
                    );
                    ui.selectable_value(
                        &mut self.uniprot_feature_query_mode,
                        UniprotFeatureCodingDnaQueryMode::TranslationSpeedOptimized,
                        "translation_speed_optimized",
                    );
                    ui.selectable_value(
                        &mut self.uniprot_feature_query_mode,
                        UniprotFeatureCodingDnaQueryMode::Both,
                        "both",
                    );
                });
            egui::ComboBox::from_label("speed profile")
                .selected_text(Self::uniprot_feature_speed_profile_label(
                    self.uniprot_feature_speed_profile,
                ))
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.uniprot_feature_speed_profile, None, "Auto");
                    ui.selectable_value(
                        &mut self.uniprot_feature_speed_profile,
                        Some(TranslationSpeedProfile::Human),
                        "Human",
                    );
                    ui.selectable_value(
                        &mut self.uniprot_feature_speed_profile,
                        Some(TranslationSpeedProfile::Mouse),
                        "Mouse",
                    );
                    ui.selectable_value(
                        &mut self.uniprot_feature_speed_profile,
                        Some(TranslationSpeedProfile::Yeast),
                        "Yeast",
                    );
                    ui.selectable_value(
                        &mut self.uniprot_feature_speed_profile,
                        Some(TranslationSpeedProfile::Ecoli),
                        "E. coli",
                    );
                });
        });
        ui.small(
            "Exon numbering follows transcript order. On reverse-strand transcripts, exon 1 is the transcript 5' exon rather than the lowest genomic coordinate.",
        );
        self.render_uniprot_feature_coding_report(ui);
        ui.separator();
        ui.label("UniProt projection audit");
        ui.small(
            "Run the shared audit that compares projected transcript accounting against stored UniProt and imported Ensembl evidence. The audit stores a reusable report plus a local unsent maintainer-email draft.",
        );
        ui.horizontal(|ui| {
            ui.label("audit transcript");
            ui.text_edit_singleline(&mut self.uniprot_audit_transcript_id)
                .on_hover_text(
                    "Optional transcript filter for the audit (empty = all projected transcripts)",
                );
            ui.label("audit report_id");
            ui.text_edit_singleline(&mut self.uniprot_audit_report_id)
                .on_hover_text("Optional stored audit report id override");
            if ui
                .button("Audit Mapped Isoforms")
                .on_hover_text("Run the high-level UniProt projection audit and store its report")
                .clicked()
            {
                self.run_uniprot_projection_audit_from_dialog();
            }
        });
        ui.horizontal(|ui| {
            ui.label("parity report_id");
            ui.text_edit_singleline(&mut self.uniprot_audit_parity_report_id)
                .on_hover_text("Optional stored audit parity report id override");
            if ui
                .button("Audit Parity")
                .on_hover_text(
                    "Compare the direct integrated audit against the same result reconstructed from the reusable primitives",
                )
                .clicked()
            {
                self.run_uniprot_projection_audit_parity_from_dialog();
            }
        });
        let recent_audits = self.recent_uniprot_audit_reports_for_dialog(8);
        if !recent_audits.is_empty() {
            ui.small(format!("Recent audits ({})", recent_audits.len()));
            egui::Grid::new("uniprot_audit_reports_grid")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("report_id");
                    ui.strong("projection");
                    ui.strong("generated_ms");
                    ui.strong("failing");
                    ui.strong("action");
                    ui.end_row();
                    for row in &recent_audits {
                        ui.monospace(&row.report_id);
                        ui.monospace(&row.projection_id);
                        ui.label(row.generated_at_unix_ms.to_string());
                        ui.label(row.failing_transcript_count.to_string());
                        if ui
                            .small_button("Show")
                            .on_hover_text("Load this stored UniProt audit report")
                            .clicked()
                        {
                            self.load_uniprot_audit_report_from_dialog(&row.report_id);
                        }
                        ui.end_row();
                    }
                });
        }
        if let Some(report) = self.uniprot_audit_report.as_ref() {
            ui.separator();
            ui.label(format!("Loaded audit: {}", report.report_id));
            ui.small(format!(
                "projection={} | rows={} | warnings={}",
                report.projection_id,
                report.rows.len(),
                report.warnings.len()
            ));
            for row in &report.rows {
                ui.collapsing(
                    format!(
                        "{} ({}) [{}]",
                        row.transcript_id,
                        row.transcript_label,
                        row.status.as_str()
                    ),
                    |ui| {
                        ui.small(format!(
                            "exon_nt_sum={} | untranslated_5'={} | untranslated_3'={} | translated_nt={} | divisible_by_3={} | expected_aa={} | uniprot_aa={}",
                            row.accounting.contributing_exon_nt_sum,
                            row.accounting.untranslated_5prime_nt,
                            row.accounting.untranslated_3prime_nt,
                            row.accounting.translated_nt,
                            row.accounting.translated_nt_divisible_by_3,
                            row.accounting.expected_aa_count,
                            row.accounting.uniprot_aa_count
                        ));
                        if !row.mismatch_reasons.is_empty() {
                            ui.small(format!("reasons: {}", row.mismatch_reasons.join(" | ")));
                        }
                    },
                );
            }
            if let Some(draft) = report.maintainer_email_draft.as_ref() {
                ui.label("Unsent maintainer email draft");
                let mut draft_body = draft.body.clone();
                ui.add(
                    egui::TextEdit::multiline(&mut draft_body)
                        .desired_width(f32::INFINITY)
                        .font(egui::TextStyle::Monospace)
                        .desired_rows(8)
                        .interactive(false),
                );
            }
        }
        let recent_parity = self.recent_uniprot_audit_parity_reports_for_dialog(6);
        if !recent_parity.is_empty() {
            ui.separator();
            ui.small(format!(
                "Recent audit parity reports ({})",
                recent_parity.len()
            ));
            egui::Grid::new("uniprot_audit_parity_reports_grid")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("report_id");
                    ui.strong("projection");
                    ui.strong("generated_ms");
                    ui.strong("rows");
                    ui.strong("action");
                    ui.end_row();
                    for row in &recent_parity {
                        ui.monospace(&row.report_id);
                        ui.monospace(&row.projection_id);
                        ui.label(row.generated_at_unix_ms.to_string());
                        ui.label(format!(
                            "{} / {}",
                            row.divergent_transcript_count, row.transcript_count
                        ));
                        if ui
                            .small_button("Show")
                            .on_hover_text("Load this stored UniProt audit parity report")
                            .clicked()
                        {
                            self.load_uniprot_audit_parity_report_from_dialog(&row.report_id);
                        }
                        ui.end_row();
                    }
                });
        }
        if let Some(report) = self.uniprot_audit_parity_report.as_ref() {
            ui.separator();
            ui.label(format!("Loaded audit parity: {}", report.report_id));
            ui.small(format!(
                "rows={} | email_draft_transcripts_match={}",
                report.rows.len(),
                report.email_draft_transcripts_match
            ));
        }
        ui.separator();
        ui.label("Linked nucleotide retrieval (EMBL/GenBank crossref)");
        ui.horizontal(|ui| {
            ui.label("accession");
            ui.text_edit_singleline(&mut self.uniprot_linked_accession)
                .on_hover_text(
                    "Optional accession override (empty = auto from UniProt DR EMBL/GenBank)",
                );
            ui.label("as_id");
            ui.text_edit_singleline(&mut self.uniprot_linked_as_id)
                .on_hover_text("Optional project sequence id override");
            if ui
                .button("Fetch Linked GenBank")
                .on_hover_text("Fetch nucleotide entry linked from UniProt and import as sequence")
                .clicked()
            {
                self.fetch_uniprot_linked_genbank_from_dialog();
            }
        });
        ui.separator();
        let recent_entries = self.recent_uniprot_entries_for_dialog(12);
        ui.label(format!(
            "Recent imported entries ({})",
            recent_entries.len()
        ));
        if recent_entries.is_empty() {
            ui.small("No UniProt entries imported in this project yet.");
        } else {
            egui::ScrollArea::vertical()
                .id_salt("protein_evidence_recent_uniprot_entries_scroll")
                .max_height(180.0)
                .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    egui::Grid::new("uniprot_recent_entries_grid")
                        .striped(true)
                        .show(ui, |ui| {
                            ui.strong("entry_id");
                            ui.strong("accession");
                            ui.strong("primary_id");
                            ui.strong("AA");
                            ui.strong("imported_ms");
                            ui.strong("source");
                            ui.strong("action");
                            ui.end_row();
                            for row in recent_entries.iter() {
                                ui.monospace(&row.entry_id);
                                ui.monospace(&row.accession);
                                ui.label(&row.primary_id);
                                ui.label(row.sequence_length.to_string());
                                ui.label(row.imported_at_unix_ms.to_string());
                                ui.label(&row.source);
                                if ui
                                    .small_button("Select")
                                    .on_hover_text(
                                        "Use this row as the active entry_id for projection",
                                    )
                                    .clicked()
                                {
                                    self.uniprot_entry_id = row.entry_id.clone();
                                    if self.uniprot_query.trim().is_empty()
                                        && !row.accession.trim().is_empty()
                                    {
                                        self.uniprot_query = row.accession.clone();
                                    }
                                }
                                ui.end_row();
                            }
                        });
                });
        }
        ui.separator();
        let recent_projections = self.recent_uniprot_projections_for_dialog(12);
        ui.label(format!(
            "Recent projections for current context ({})",
            recent_projections.len()
        ));
        if recent_projections.is_empty() {
            ui.small(
                "No stored projections match the current seq_id/entry_id filter yet. Project an entry to unlock the protein expert view.",
            );
        } else {
            egui::ScrollArea::vertical()
                .id_salt("protein_evidence_recent_projections_scroll")
                .max_height(160.0)
                .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    egui::Grid::new("uniprot_recent_projection_grid")
                        .striped(true)
                        .show(ui, |ui| {
                            ui.strong("projection_id");
                            ui.strong("entry_id");
                            ui.strong("seq_id");
                            ui.strong("transcripts");
                            ui.strong("created_ms");
                            ui.strong("actions");
                            ui.end_row();
                            for row in &recent_projections {
                                ui.monospace(&row.projection_id);
                                ui.monospace(&row.entry_id);
                                ui.monospace(&row.seq_id);
                                ui.label(row.transcript_projection_count.to_string());
                                ui.label(row.created_at_unix_ms.to_string());
                                ui.horizontal(|ui| {
                                    if ui
                                        .small_button("Use")
                                        .on_hover_text(
                                            "Use this stored projection as the active projection_id and seq_id",
                                        )
                                        .clicked()
                                    {
                                        self.uniprot_entry_id = row.entry_id.clone();
                                        self.uniprot_map_seq_id = row.seq_id.clone();
                                        self.uniprot_map_projection_id = row.projection_id.clone();
                                    }
                                    if ui
                                        .small_button("Open Protein Expert")
                                        .on_hover_text(
                                            "Open the dedicated UniProt protein mapping expert for this stored projection",
                                        )
                                        .clicked()
                                    {
                                        self.open_uniprot_projection_expert_from_dialog(
                                            &row.seq_id,
                                            &row.projection_id,
                                        );
                                    }
                                    if ui
                                        .small_button("Render SVG...")
                                        .on_hover_text(
                                            "Export this stored UniProt projection directly as an SVG expert-view artifact",
                                        )
                                        .clicked()
                                    {
                                        self.export_uniprot_projection_svg_from_dialog(
                                            &row.seq_id,
                                            &row.projection_id,
                                        );
                                    }
                                });
                                ui.end_row();
                            }
                        });
                });
        }
        ui.separator();
        ui.label("Ensembl protein evidence");
        ui.small(
            "Fetch one Ensembl protein or transcript through the shared REST-backed engine path, reuse it as optional external evidence in the Protein Expert, or import the protein as a standalone sequence.",
        );
        ui.horizontal(|ui| {
            ui.label("query");
            ui.text_edit_singleline(&mut self.ensembl_protein_query)
                .on_hover_text(
                    "Ensembl protein or transcript identifier to fetch (for example ENSP... or ENST...)",
                );
            if ui
                .button("Fetch Ensembl")
                .on_hover_text(
                    "Fetch one Ensembl protein entry into project metadata for later compare/import reuse",
                )
                .clicked()
            {
                self.fetch_ensembl_protein_from_dialog();
            }
        });
        ui.horizontal(|ui| {
            ui.label("entry_id");
            ui.text_edit_singleline(&mut self.ensembl_protein_entry_id).on_hover_text(
                "Stable project entry identifier for the fetched Ensembl protein evidence",
            );
            ui.label("output_id");
            ui.text_edit_singleline(&mut self.ensembl_protein_output_id).on_hover_text(
                "Optional sequence id override when importing the Ensembl protein as a standalone sequence",
            );
            if ui
                .button("Import Sequence")
                .on_hover_text(
                    "Import the selected Ensembl protein entry as a first-class protein sequence in the project",
                )
                .clicked()
            {
                self.import_ensembl_protein_sequence_from_dialog();
            }
        });
        ui.horizontal(|ui| {
            if ui
                .button("Open Ensembl Protein Expert")
                .on_hover_text(
                    "Open the shared Protein Expert with the selected Ensembl protein entry as optional external evidence",
                )
                .clicked()
            {
                let seq_id = self.uniprot_map_seq_id.trim().to_string();
                let transcript_filter =
                    Self::uniprot_optional_trimmed(&self.uniprot_map_transcript_id);
                let entry_id = self.ensembl_protein_entry_id.trim().to_string();
                self.open_ensembl_protein_expert_from_dialog(
                    &seq_id,
                    transcript_filter.as_deref(),
                    &entry_id,
                );
            }
            if ui
                .button("Render Ensembl Protein SVG...")
                .on_hover_text(
                    "Export the Ensembl-backed Protein Expert directly as an SVG through the shared feature-expert route",
                )
                .clicked()
            {
                let seq_id = self.uniprot_map_seq_id.trim().to_string();
                let transcript_filter =
                    Self::uniprot_optional_trimmed(&self.uniprot_map_transcript_id);
                let entry_id = self.ensembl_protein_entry_id.trim().to_string();
                self.export_ensembl_protein_svg_from_dialog(
                    &seq_id,
                    transcript_filter.as_deref(),
                    &entry_id,
                );
            }
        });
        ui.small(
            "Ensembl compare/export reuses the same seq_id / transcript filter controls above. If transcript is empty, the fetched entry's transcript will be used when available.",
        );
        let recent_ensembl_entries = self.recent_ensembl_protein_entries_for_dialog(12);
        ui.label(format!(
            "Recent Ensembl protein entries ({})",
            recent_ensembl_entries.len()
        ));
        if recent_ensembl_entries.is_empty() {
            ui.small("No Ensembl protein entries fetched in this project yet.");
        } else {
            egui::ScrollArea::vertical()
                .id_salt("protein_evidence_recent_ensembl_entries_scroll")
                .max_height(180.0)
                .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    egui::Grid::new("ensembl_recent_entries_grid")
                        .striped(true)
                        .show(ui, |ui| {
                            ui.strong("entry_id");
                            ui.strong("protein_id");
                            ui.strong("transcript_id");
                            ui.strong("gene");
                            ui.strong("AA");
                            ui.strong("imported_ms");
                            ui.strong("source_query");
                            ui.strong("actions");
                            ui.end_row();
                            for row in &recent_ensembl_entries {
                                ui.monospace(&row.entry_id);
                                ui.monospace(&row.protein_id);
                                ui.monospace(&row.transcript_id);
                                ui.label(row.gene_symbol.as_deref().unwrap_or("-"));
                                ui.label(row.sequence_length.to_string());
                                ui.label(row.imported_at_unix_ms.to_string());
                                ui.label(row.source_query.as_deref().unwrap_or("-"));
                                ui.horizontal(|ui| {
                                    if ui
                                        .small_button("Select")
                                        .on_hover_text(
                                            "Use this row as the active Ensembl protein evidence entry",
                                        )
                                        .clicked()
                                    {
                                        self.ensembl_protein_entry_id = row.entry_id.clone();
                                        self.ensembl_protein_query = row.protein_id.clone();
                                        if self.uniprot_map_transcript_id.trim().is_empty() {
                                            self.uniprot_map_transcript_id =
                                                row.transcript_id.clone();
                                        }
                                    }
                                    if ui
                                        .small_button("Open Protein Expert")
                                        .on_hover_text(
                                            "Open the shared Protein Expert with this Ensembl entry as optional external evidence",
                                        )
                                        .clicked()
                                    {
                                        self.ensembl_protein_entry_id = row.entry_id.clone();
                                        if self.uniprot_map_transcript_id.trim().is_empty() {
                                            self.uniprot_map_transcript_id =
                                                row.transcript_id.clone();
                                        }
                                        let seq_id = self.uniprot_map_seq_id.trim().to_string();
                                        let transcript_filter =
                                            Self::uniprot_optional_trimmed(
                                                &self.uniprot_map_transcript_id,
                                            );
                                        self.open_ensembl_protein_expert_from_dialog(
                                            &seq_id,
                                            transcript_filter.as_deref(),
                                            &row.entry_id,
                                        );
                                    }
                                    if ui
                                        .small_button("Import Sequence")
                                        .on_hover_text(
                                            "Import this Ensembl protein entry as a standalone sequence",
                                        )
                                        .clicked()
                                    {
                                        self.ensembl_protein_entry_id = row.entry_id.clone();
                                        self.import_ensembl_protein_sequence_from_dialog();
                                    }
                                    if ui
                                        .small_button("Render SVG...")
                                        .on_hover_text(
                                            "Export this Ensembl-backed Protein Expert directly as SVG",
                                        )
                                        .clicked()
                                    {
                                        self.ensembl_protein_entry_id = row.entry_id.clone();
                                        if self.uniprot_map_transcript_id.trim().is_empty() {
                                            self.uniprot_map_transcript_id =
                                                row.transcript_id.clone();
                                        }
                                        let seq_id = self.uniprot_map_seq_id.trim().to_string();
                                        let transcript_filter =
                                            Self::uniprot_optional_trimmed(
                                                &self.uniprot_map_transcript_id,
                                            );
                                        self.export_ensembl_protein_svg_from_dialog(
                                            &seq_id,
                                            transcript_filter.as_deref(),
                                            &row.entry_id,
                                        );
                                    }
                                });
                                ui.end_row();
                            }
                        });
                });
        }
        if let Some(entry) = self.selected_ensembl_protein_entry_for_dialog() {
            ui.separator();
            ui.label("Selected Ensembl evidence");
            ui.small(
                "Stored Ensembl protein evidence stays optional here: transcript-native translation remains the authoritative product model, while this record acts as one external comparison opinion.",
            );
            let current_transcript_filter =
                Self::uniprot_optional_trimmed(&self.uniprot_map_transcript_id);
            if current_transcript_filter
                .as_deref()
                .map(|value| value != entry.transcript_id)
                .unwrap_or(false)
            {
                ui.horizontal_wrapped(|ui| {
                    ui.colored_label(
                        egui::Color32::from_rgb(180, 83, 9),
                        format!(
                            "Current transcript filter '{}' differs from Ensembl transcript '{}'.",
                            current_transcript_filter.as_deref().unwrap_or("-"),
                            entry.transcript_id
                        ),
                    );
                    if ui
                        .small_button("Use entry transcript")
                        .on_hover_text(
                            "Replace the current transcript filter with the transcript recorded in this Ensembl protein entry",
                        )
                        .clicked()
                    {
                        self.uniprot_map_transcript_id = entry.transcript_id.clone();
                    }
                });
            }
            egui::Grid::new("ensembl_selected_entry_grid")
                .num_columns(2)
                .spacing([12.0, 4.0])
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("entry_id");
                    ui.monospace(&entry.entry_id);
                    ui.end_row();
                    ui.strong("protein_id");
                    ui.monospace(&entry.protein_id);
                    ui.end_row();
                    ui.strong("transcript_id");
                    ui.monospace(&entry.transcript_id);
                    ui.end_row();
                    ui.strong("gene");
                    ui.label(
                        entry
                            .gene_symbol
                            .as_deref()
                            .or(entry.gene_id.as_deref())
                            .unwrap_or("-"),
                    );
                    ui.end_row();
                    ui.strong("transcript label");
                    ui.label(entry.transcript_display_name.as_deref().unwrap_or("-"));
                    ui.end_row();
                    ui.strong("species");
                    ui.label(entry.species.as_deref().unwrap_or("-"));
                    ui.end_row();
                    ui.strong("derived aa");
                    ui.monospace(entry.sequence_length.to_string());
                    ui.end_row();
                    ui.strong("feature count");
                    ui.monospace(entry.features.len().to_string());
                    ui.end_row();
                    ui.strong("source query");
                    ui.label(entry.source_query.as_deref().unwrap_or("-"));
                    ui.end_row();
                });
            let feature_key_summary = Self::summarize_ensembl_protein_feature_keys(&entry, 8);
            if !feature_key_summary.is_empty() {
                ui.horizontal_wrapped(|ui| {
                    ui.strong("feature keys");
                    ui.monospace(feature_key_summary.join(", "));
                });
            }
            if !entry.aliases.is_empty() {
                let mut aliases = entry.aliases.clone();
                aliases.sort();
                if aliases.len() > 8 {
                    aliases.truncate(8);
                    aliases.push("...".to_string());
                }
                ui.horizontal_wrapped(|ui| {
                    ui.strong("aliases");
                    ui.monospace(aliases.join(", "));
                });
            }
        } else if !self.ensembl_protein_entry_id.trim().is_empty() {
            ui.separator();
            ui.colored_label(
                egui::Color32::from_rgb(180, 83, 9),
                format!(
                    "Selected Ensembl entry '{}' is not currently available in project metadata.",
                    self.ensembl_protein_entry_id.trim()
                ),
            );
        }
        ui.separator();
        ui.label("Reverse translate protein");
        ui.small(
            "Generate one synthetic coding DNA sequence from any first-class protein in the project and inspect the resolved translation table, codon-speed profile, and optional annealing heuristic.",
        );
        if protein_seq_ids.is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(180, 83, 9),
                "No first-class protein sequences are available yet. Import one from Ensembl or derive one first.",
            );
        } else {
            egui::ComboBox::from_label("protein seq_id")
                .selected_text(self.reverse_translate_protein_seq_id.clone())
                .show_ui(ui, |ui| {
                    for seq_id in &protein_seq_ids {
                        ui.selectable_value(
                            &mut self.reverse_translate_protein_seq_id,
                            seq_id.clone(),
                            seq_id,
                        );
                    }
                });
            ui.horizontal(|ui| {
                ui.label("output_id");
                ui.text_edit_singleline(&mut self.reverse_translate_output_id)
                    .on_hover_text(
                        "Optional coding-sequence id override (defaults to PROTEIN_ID__coding)",
                    );
                ui.label("translation table");
                ui.text_edit_singleline(&mut self.reverse_translate_translation_table)
                    .on_hover_text(
                        "Optional explicit NCBI translation table override (empty = resolve from protein/context)",
                    );
            });
            ui.horizontal(|ui| {
                egui::ComboBox::from_label("speed profile")
                    .selected_text(Self::uniprot_feature_speed_profile_label(
                        self.reverse_translate_speed_profile,
                    ))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.reverse_translate_speed_profile, None, "Auto");
                        ui.selectable_value(
                            &mut self.reverse_translate_speed_profile,
                            Some(TranslationSpeedProfile::Human),
                            "Human",
                        );
                        ui.selectable_value(
                            &mut self.reverse_translate_speed_profile,
                            Some(TranslationSpeedProfile::Mouse),
                            "Mouse",
                        );
                        ui.selectable_value(
                            &mut self.reverse_translate_speed_profile,
                            Some(TranslationSpeedProfile::Yeast),
                            "Yeast",
                        );
                        ui.selectable_value(
                            &mut self.reverse_translate_speed_profile,
                            Some(TranslationSpeedProfile::Ecoli),
                            "E. coli",
                        );
                    });
                egui::ComboBox::from_label("speed mark")
                    .selected_text(Self::reverse_translation_speed_mark_label(
                        self.reverse_translate_speed_mark,
                    ))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.reverse_translate_speed_mark, None, "Auto");
                        ui.selectable_value(
                            &mut self.reverse_translate_speed_mark,
                            Some(TranslationSpeedMark::Fast),
                            "Fast",
                        );
                        ui.selectable_value(
                            &mut self.reverse_translate_speed_mark,
                            Some(TranslationSpeedMark::Slow),
                            "Slow",
                        );
                    });
                ui.label("target anneal Tm");
                ui.text_edit_singleline(&mut self.reverse_translate_target_anneal_tm_c)
                    .on_hover_text(
                        "Optional lightweight local Tm heuristic target in °C (empty = disabled)",
                    );
                ui.label("window bp");
                ui.text_edit_singleline(&mut self.reverse_translate_anneal_window_bp)
                    .on_hover_text(
                        "Window size for the optional local annealing heuristic (default 20)",
                    );
                if ui
                    .button("Reverse Translate")
                    .on_hover_text(
                        "Create one synthetic coding DNA sequence from the selected protein and store the resolved translation provenance",
                    )
                    .clicked()
                {
                    self.reverse_translate_protein_from_dialog();
                }
            });
            self.render_reverse_translation_report(ui);
        }
        ui.separator();
        ui.label("Protease digest protein");
        ui.small(
            "Predict protease cleavage sites on the selected first-class protein and optionally materialize peptide products while keeping transcript-derived provenance.",
        );
        if protein_seq_ids.is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(180, 83, 9),
                "No protein sequence is available for protease digestion.",
            );
        } else {
            ui.horizontal(|ui| {
                ui.label("proteases");
                ui.text_edit_singleline(&mut self.protease_digest_names)
                    .on_hover_text("Comma-separated protease names or aliases, for example Trypsin or Trypsin,Lys-C");
                ui.label("min aa");
                ui.text_edit_singleline(&mut self.protease_digest_min_length_aa)
                    .on_hover_text("Minimum peptide product length to keep (default 1)");
                ui.checkbox(&mut self.protease_digest_materialize, "materialize")
                    .on_hover_text("Create peptide sequence entries for the predicted digest products");
            });
            ui.horizontal(|ui| {
                ui.label("output prefix");
                ui.text_edit_singleline(&mut self.protease_digest_output_prefix)
                    .on_hover_text(
                        "Optional peptide seq_id prefix (defaults to PROTEIN_ID_protease_digest)",
                    );
                if ui
                    .button("Digest Protein")
                    .on_hover_text(
                        "Run the shared ProteaseDigestProteinSequence operation on the selected protein sequence",
                    )
                    .clicked()
                {
                    self.digest_selected_protein_from_dialog();
                }
            });
            ui.small(
                "Uses the same protein seq_id selector as reverse translation above; uncheck materialize for report-only prediction.",
            );
            self.render_protease_digest_report(ui);
        }
        ui.separator();
        ui.label("DNA handoff reasoning");
        ui.small(
            "Rank transcript-native CDS reuse, mapped feature-coding DNA, and synthetic reverse-translation fallback without materializing DNA sequences yet.",
        );
        ui.horizontal(|ui| {
            egui::ComboBox::from_label("ranking goal")
                .selected_text(self.protein_handoff_ranking_goal.as_str())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.protein_handoff_ranking_goal,
                        ProteinToDnaHandoffRankingGoal::BalancedProvenance,
                        "Balanced provenance",
                    );
                    ui.selectable_value(
                        &mut self.protein_handoff_ranking_goal,
                        ProteinToDnaHandoffRankingGoal::NativeFidelity,
                        "Native fidelity",
                    );
                    ui.selectable_value(
                        &mut self.protein_handoff_ranking_goal,
                        ProteinToDnaHandoffRankingGoal::ExpressionOptimized,
                        "Expression optimized",
                    );
                });
            if ui
                .button("Build DNA Handoff")
                .on_hover_text(
                    "Build a construct-reasoning graph that ranks transcript-native, mapped-feature, and synthetic protein-to-DNA handoff candidates",
                )
                .clicked()
            {
                self.build_protein_to_dna_handoff_reasoning_from_dialog();
            }
        });
        ui.small(
            "Reuses the active protein sequence, transcript/projection filters, feature query, and reverse-translation controls above.",
        );
        self.render_protein_handoff_graph(ui);
        if !self.protein_handoff_status.trim().is_empty() {
            ui.small(self.protein_handoff_status.clone());
        }
        if !self.uniprot_status.trim().is_empty() {
            ui.separator();
            ui.monospace(self.uniprot_status.clone());
        }
        close_requested
    }

    pub(super) fn render_uniprot_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_uniprot_dialog {
            return;
        }

        let viewport_id = Self::uniprot_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "Protein Evidence",
            egui::Id::new(("hosted_uniprot_window", viewport_id)),
            viewport_id,
            Vec2::new(1040.0, 680.0),
            Vec2::new(820.0, 480.0),
        );
        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut open = self.show_uniprot_dialog;
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                    egui::ScrollArea::vertical()
                        .id_salt("protein_evidence_embedded_scroll")
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            scroll_input_policy::apply_scrollarea_keyboard_navigation(
                                ui,
                                scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                            );
                            close_requested = self.render_uniprot_dialog_contents(ui);
                        });
                });
                if close_requested {
                    open = false;
                }
                self.show_uniprot_dialog = open;
                return;
            }

            let mut close_requested = false;
            crate::egui_compat::show_central_panel(ctx, egui::CentralPanel::default(), |ui| {
                egui::ScrollArea::vertical()
                    .id_salt("protein_evidence_viewport_scroll")
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        scroll_input_policy::apply_scrollarea_keyboard_navigation(
                            ui,
                            scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                        );
                        close_requested = self.render_uniprot_dialog_contents(ui);
                    });
            });

            if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                self.show_uniprot_dialog = false;
            }
        });
    }
}
