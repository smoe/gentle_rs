//! Bind RNA reference panels to existing gene-set operands and persist opt-in graph snapshots.

use super::*;

const KEY: &str = "rna_read_panels.v1";

pub(super) fn validate_gene_set(panel: &RnaReadPanel) -> Result<(), EngineError> {
    let Some(report) = &panel.gene_set_resolution else {
        return Ok(());
    };
    if report.schema != "gentle.gene_set_resolution.v1"
        || report.unresolved_member_count != 0
        || !report.unresolved_members.is_empty()
        || report.resolved_member_count != report.resolved_members.len()
        || report.resolved_members.len() != panel.loci.len()
    {
        return Err(bad(
            "panel_gene_set_incomplete: every resolved member needs one explicit reference locus",
        ));
    }
    let mut used = BTreeSet::new();
    for member in &report.resolved_members {
        let id = member
            .gene_id
            .as_deref()
            .filter(|id| !id.trim().is_empty())
            .ok_or_else(|| {
                bad(format!(
                    "panel_gene_set_missing_stable_id: {}",
                    member.symbol
                ))
            })?;
        let matches = panel
            .loci
            .iter()
            .filter(|l| l.source.gene_id == id)
            .collect::<Vec<_>>();
        let [locus] = matches.as_slice() else {
            return Err(bad(format!(
                "panel_gene_set_ambiguous_or_missing_source: {id}"
            )));
        };
        if !used.insert(locus.source.locus_id.clone()) {
            return Err(bad("panel_gene_set_duplicate_binding"));
        }
        if let Some(context) = report
            .effective_biological_context(member)
            .map_err(|e| bad(e.to_string()))?
        {
            if context.genome_id.as_deref().is_some_and(|genome| {
                locus
                    .genome_anchor
                    .as_ref()
                    .is_none_or(|a| a.genome_id != genome)
            }) || context
                .assembly_accession
                .as_deref()
                .is_some_and(|a| a != locus.source.assembly)
                || context
                    .annotation_release
                    .as_deref()
                    .is_some_and(|r| r != locus.source.annotation_release)
            {
                return Err(bad(format!(
                    "panel_gene_set_reference_mismatch_or_unassessed: {id}"
                )));
            }
        }
    }
    Ok(())
}

impl GentleEngine {
    pub(crate) fn rna_read_panel_gene_set_binding(
        report: &GeneSetResolutionReport,
    ) -> Result<RnaReadPanelGeneSetBinding, EngineError> {
        Ok(RnaReadPanelGeneSetBinding {
            resolution_id: Self::gene_set_resolution_artifact_id(report),
            resolution_sha256: sha256_prefixed_bytes(&encode(report)?),
        })
    }

    pub fn rna_read_panel_gene_sets(&self) -> Result<Vec<RnaReadPanelGeneSetChoice>, EngineError> {
        Self::gene_set_resolution_artifacts_from_state(self.state())
            .into_iter()
            .map(|resolution| {
                Ok(RnaReadPanelGeneSetChoice {
                    binding: Self::rna_read_panel_gene_set_binding(&resolution)?,
                    resolution,
                })
            })
            .collect()
    }

    pub(crate) fn rna_read_panel_resolve_gene_set(
        &self,
        binding: &RnaReadPanelGeneSetBinding,
    ) -> Result<GeneSetResolutionReport, EngineError> {
        let report = self.get_gene_set_resolution_artifact(&binding.resolution_id)?;
        if Self::rna_read_panel_gene_set_binding(&report)? != *binding {
            return Err(bad("panel_gene_set_changed"));
        }
        Ok(report)
    }

    /// Exact stable-ID mapping only. Missing or competing source views require explicit selection.
    pub fn rna_read_panel_gene_set_sources(
        &self,
        id: &str,
    ) -> Result<(RnaReadPanelGeneSetBinding, Vec<RnaReadPanelSourceChoice>), EngineError> {
        let report = self.get_gene_set_resolution_artifact(id)?;
        if report.unresolved_member_count != 0 || !report.unresolved_members.is_empty() {
            return Err(bad("panel_gene_set_unresolved_members"));
        }
        let choices = self.rna_read_panel_source_choices()?;
        let mut selected = Vec::new();
        for member in &report.resolved_members {
            let gene = member.gene_id.as_deref().ok_or_else(|| {
                bad(format!(
                    "panel_gene_set_missing_stable_id: {}",
                    member.symbol
                ))
            })?;
            let matching = choices
                .iter()
                .filter(|c| c.diagnostic.is_none() && c.source.gene_id == gene)
                .collect::<Vec<_>>();
            let [source] = matching.as_slice() else {
                return Err(bad(format!(
                    "panel_gene_set_ambiguous_or_missing_source: {gene}"
                )));
            };
            selected.push((*source).clone());
        }
        Ok((Self::rna_read_panel_gene_set_binding(&report)?, selected))
    }

    pub(crate) fn register_rna_read_panel(
        &mut self,
        panel: &RnaReadPanel,
        output_path: &str,
        op_id: &str,
    ) -> Result<(), EngineError> {
        validate_panel(panel)?;
        let mut entries = Self::rna_read_panel_project_entries(self.state())?
            .into_iter()
            .map(|entry| (entry.panel.content_sha256.clone(), entry))
            .collect::<BTreeMap<_, _>>();
        entries
            .entry(panel.content_sha256.clone())
            .or_insert_with(|| RnaReadPanelProjectEntry {
                panel: panel.clone(),
                output_path: output_path.into(),
                op_id: op_id.into(),
                created_at_unix_ms: Self::now_unix_ms(),
            });
        self.state.metadata.insert(
            KEY.into(),
            serde_json::to_value(entries).map_err(|e| bad(e.to_string()))?,
        );
        Ok(())
    }

    pub fn rna_read_panel_project_entries(
        state: &ProjectState,
    ) -> Result<Vec<RnaReadPanelProjectEntry>, EngineError> {
        match state.metadata.get(KEY) {
            None => Ok(vec![]),
            Some(value) => {
                let entries: BTreeMap<String, RnaReadPanelProjectEntry> =
                    serde_json::from_value(value.clone())
                        .map_err(|e| bad(format!("panel_project_store: {e}")))?;
                for (id, entry) in &entries {
                    validate_panel(&entry.panel)?;
                    if id != &entry.panel.content_sha256 {
                        return Err(bad("panel_project_store_identity"));
                    }
                }
                Ok(entries.into_values().collect())
            }
        }
    }
}
