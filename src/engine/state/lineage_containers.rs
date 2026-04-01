//! Lineage graph and container/arrangement helper routines.
//!
//! Container-first project semantics and lineage bookkeeping live together here
//! because most mutations need to keep both models in sync.
//!
//! Look here for:
//! - lineage-node creation and lookup helpers
//! - container/arrangement bookkeeping tied to sequence mutations
//! - macro-instance/container-state synchronization utilities

use super::*;

impl GentleEngine {
    pub(super) fn add_lineage_node(
        &mut self,
        seq_id: &str,
        origin: SequenceOrigin,
        created_by_op: Option<&str>,
    ) -> NodeId {
        self.state.lineage.next_node_counter += 1;
        let node_id = format!("n-{}", self.state.lineage.next_node_counter);
        let node = LineageNode {
            node_id: node_id.clone(),
            seq_id: seq_id.to_string(),
            created_by_op: created_by_op.map(|s| s.to_string()),
            origin,
            created_at_unix_ms: Self::now_unix_ms(),
        };
        self.state
            .lineage
            .seq_to_node
            .insert(seq_id.to_string(), node_id.clone());
        self.state.lineage.nodes.insert(node_id.clone(), node);
        node_id
    }

    pub(super) fn ensure_lineage_node(&mut self, seq_id: &str) -> NodeId {
        if let Some(node_id) = self.state.lineage.seq_to_node.get(seq_id) {
            return node_id.clone();
        }
        self.add_lineage_node(seq_id, SequenceOrigin::ImportedUnknown, None)
    }

    pub(super) fn next_macro_instance_id(&mut self) -> String {
        loop {
            self.state.lineage.next_macro_instance_counter += 1;
            let candidate = format!("macro-{}", self.state.lineage.next_macro_instance_counter);
            if !self
                .state
                .lineage
                .macro_instances
                .iter()
                .any(|instance| instance.macro_instance_id == candidate)
            {
                return candidate;
            }
        }
    }

    /// Append one macro-instance audit record, allocating ids/defaults as
    /// needed.
    pub fn record_lineage_macro_instance(&mut self, mut instance: LineageMacroInstance) -> String {
        if instance.macro_instance_id.trim().is_empty() {
            instance.macro_instance_id = self.next_macro_instance_id();
        } else if self
            .state
            .lineage
            .macro_instances
            .iter()
            .any(|existing| existing.macro_instance_id == instance.macro_instance_id)
        {
            instance.macro_instance_id = self.next_macro_instance_id();
        }
        if instance.created_at_unix_ms == 0 {
            instance.created_at_unix_ms = Self::now_unix_ms();
        }
        instance.run_id = instance.run_id.trim().to_string();
        if instance.run_id.is_empty() {
            instance.run_id = "macro".to_string();
        }
        self.state.lineage.macro_instances.push(instance.clone());
        instance.macro_instance_id
    }

    /// Borrow the persistent macro-instance audit trail in insertion order.
    pub fn lineage_macro_instances(&self) -> &[LineageMacroInstance] {
        &self.state.lineage.macro_instances
    }

    pub(super) fn reconcile_lineage_nodes(&mut self) {
        let seq_ids: Vec<String> = self.state.sequences.keys().cloned().collect();
        for seq_id in seq_ids {
            let _ = self.ensure_lineage_node(&seq_id);
        }
    }

    pub(super) fn next_container_id(&mut self) -> ContainerId {
        loop {
            self.state.container_state.next_container_counter += 1;
            let id = format!(
                "container-{}",
                self.state.container_state.next_container_counter
            );
            if !self.state.container_state.containers.contains_key(&id) {
                return id;
            }
        }
    }

    pub(super) fn next_arrangement_id(&mut self) -> String {
        loop {
            self.state.container_state.next_arrangement_counter += 1;
            let id = format!(
                "arrangement-{}",
                self.state.container_state.next_arrangement_counter
            );
            if !self.state.container_state.arrangements.contains_key(&id) {
                return id;
            }
        }
    }

    pub(super) fn add_serial_arrangement(
        &mut self,
        container_ids: &[ContainerId],
        arrangement_id: Option<String>,
        name: Option<String>,
        ladders: Option<Vec<String>>,
        created_by_op: Option<&str>,
    ) -> Result<String, EngineError> {
        if container_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "CreateArrangementSerial requires at least one container id".to_string(),
            });
        }
        let mut lane_container_ids: Vec<ContainerId> = vec![];
        let mut seen: HashSet<String> = HashSet::new();
        for container_id in container_ids {
            let trimmed = container_id.trim();
            if trimmed.is_empty() {
                continue;
            }
            if !self.state.container_state.containers.contains_key(trimmed) {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Container '{trimmed}' not found"),
                });
            }
            if seen.insert(trimmed.to_string()) {
                lane_container_ids.push(trimmed.to_string());
            }
        }
        if lane_container_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "CreateArrangementSerial requires at least one non-empty container id"
                    .to_string(),
            });
        }
        let arrangement_id = arrangement_id
            .map(|v| v.trim().to_string())
            .filter(|v| !v.is_empty())
            .unwrap_or_else(|| self.next_arrangement_id());
        if self
            .state
            .container_state
            .arrangements
            .contains_key(&arrangement_id)
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Arrangement '{arrangement_id}' already exists"),
            });
        }
        let arrangement = Arrangement {
            arrangement_id: arrangement_id.clone(),
            mode: ArrangementMode::Serial,
            name: name.map(|v| v.trim().to_string()).filter(|v| !v.is_empty()),
            lane_container_ids,
            ladders: ladders
                .unwrap_or_default()
                .into_iter()
                .map(|v| v.trim().to_string())
                .filter(|v| !v.is_empty())
                .collect::<Vec<_>>(),
            created_by_op: created_by_op.map(ToString::to_string),
            created_at_unix_ms: Self::now_unix_ms(),
        };
        self.state
            .container_state
            .arrangements
            .insert(arrangement_id.clone(), arrangement);
        Ok(arrangement_id)
    }

    pub(super) fn add_container(
        &mut self,
        members: &[SeqId],
        kind: ContainerKind,
        name: Option<String>,
        created_by_op: Option<&str>,
    ) -> Option<ContainerId> {
        if members.is_empty() {
            return None;
        }
        let container_id = self.next_container_id();
        let container = Container {
            container_id: container_id.clone(),
            kind,
            name,
            members: members.to_vec(),
            created_by_op: created_by_op.map(ToString::to_string),
            created_at_unix_ms: Self::now_unix_ms(),
        };
        self.state
            .container_state
            .containers
            .insert(container_id.clone(), container);
        for seq_id in members {
            self.state
                .container_state
                .seq_to_latest_container
                .insert(seq_id.clone(), container_id.clone());
        }
        Some(container_id)
    }

    pub(super) fn exact_singleton_container_for_seq(&self, seq_id: &str) -> Option<ContainerId> {
        if let Some(container_id) = self
            .state
            .container_state
            .seq_to_latest_container
            .get(seq_id)
            .filter(|container_id| {
                self.state
                    .container_state
                    .containers
                    .get(*container_id)
                    .is_some_and(|container| {
                        matches!(container.kind, ContainerKind::Singleton)
                            && container.members == vec![seq_id.to_string()]
                    })
            })
        {
            return Some(container_id.clone());
        }

        let mut exact_matches = self
            .state
            .container_state
            .containers
            .values()
            .filter(|container| {
                matches!(container.kind, ContainerKind::Singleton)
                    && container.members == vec![seq_id.to_string()]
            })
            .map(|container| container.container_id.clone())
            .collect::<Vec<_>>();
        exact_matches.sort();
        exact_matches.into_iter().next()
    }

    pub(super) fn ensure_exact_singleton_container_for_seq(
        &mut self,
        seq_id: &str,
        created_by_op: Option<&str>,
    ) -> Result<ContainerId, EngineError> {
        if let Some(container_id) = self.exact_singleton_container_for_seq(seq_id) {
            return Ok(container_id);
        }
        if !self.state.sequences.contains_key(seq_id) {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found for singleton container creation"),
            });
        }
        let name = self
            .state
            .sequences
            .get(seq_id)
            .and_then(|dna| dna.name().clone())
            .filter(|name| !name.trim().is_empty())
            .unwrap_or_else(|| seq_id.to_string());
        self.add_container(
            &[seq_id.to_string()],
            ContainerKind::Singleton,
            Some(name),
            created_by_op,
        )
        .ok_or_else(|| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not create singleton container for sequence '{seq_id}'"),
        })
    }

    pub(super) fn reconcile_containers(&mut self) {
        let seq_ids: Vec<SeqId> = self.state.sequences.keys().cloned().collect();
        for seq_id in &seq_ids {
            if !self
                .state
                .container_state
                .seq_to_latest_container
                .contains_key(seq_id)
            {
                let _ = self.add_container(
                    std::slice::from_ref(seq_id),
                    ContainerKind::Singleton,
                    Some(format!("Imported sequence {seq_id}")),
                    None,
                );
            }
        }
    }

    pub(super) fn add_container_from_result(&mut self, op: &Operation, result: &OpResult) {
        if result.created_seq_ids.is_empty() {
            return;
        }
        // Primer-pair design operations materialize one explicit container per
        // accepted pair inside the operation handler; skip generic aggregate
        // container creation so seq_to_latest_container keeps pointing to the
        // pair-scoped container.
        if matches!(
            op,
            Operation::DesignPrimerPairs { .. }
                | Operation::DesignInsertionPrimerPairs { .. }
                | Operation::PcrOverlapExtensionMutagenesis { .. }
        ) {
            return;
        }
        // Gibson apply materializes explicit singleton containers and its
        // arrangement directly inside the operation handler, so the generic
        // fallback must stay out of the way here.
        if matches!(op, Operation::ApplyGibsonAssemblyPlan { .. }) {
            return;
        }
        let kind = if matches!(op, Operation::SelectCandidate { .. }) {
            ContainerKind::Selection
        } else if result.created_seq_ids.len() > 1 {
            ContainerKind::Pool
        } else {
            ContainerKind::Singleton
        };
        let name = match op {
            Operation::LoadFile { .. } => Some("Imported sequence".to_string()),
            Operation::ImportUniprotEntrySequence { .. } => {
                Some("Imported UniProt sequence".to_string())
            }
            Operation::Digest { .. } => Some("Digest products".to_string()),
            Operation::DigestContainer { .. } => Some("Digest products".to_string()),
            Operation::MergeContainers { .. } => Some("Merged container".to_string()),
            Operation::MergeContainersById { .. } => Some("Merged container".to_string()),
            Operation::Ligation { .. } => Some("Ligation products".to_string()),
            Operation::LigationContainer { .. } => Some("Ligation products".to_string()),
            Operation::Pcr { .. }
            | Operation::PcrAdvanced { .. }
            | Operation::PcrMutagenesis { .. } => Some("PCR products".to_string()),
            Operation::ExtractRegion { .. } => Some("Extracted region".to_string()),
            Operation::ExtractAnchoredRegion { .. } => Some("Extracted region".to_string()),
            Operation::DeriveSplicingReferences { .. } => {
                Some("Derived splicing references".to_string())
            }
            Operation::ExtractGenomeRegion { .. } => Some("Extracted genome region".to_string()),
            Operation::ExtractGenomeGene { .. } => Some("Extracted genome gene".to_string()),
            Operation::FetchGenBankAccession { .. } => Some("Fetched GenBank sequence".to_string()),
            Operation::FetchDbSnpRegion { .. } => Some("Fetched dbSNP region".to_string()),
            Operation::FetchUniprotLinkedGenBank { .. } => {
                Some("Fetched UniProt-linked GenBank sequence".to_string())
            }
            Operation::ExtendGenomeAnchor { .. } => Some("Extended genome anchor".to_string()),
            Operation::ImportGenomeBedTrack { .. } => Some("Imported BED track".to_string()),
            Operation::ImportGenomeBigWigTrack { .. } => Some("Imported BigWig track".to_string()),
            Operation::ImportGenomeVcfTrack { .. } => Some("Imported VCF track".to_string()),
            Operation::ImportBlastHitsTrack { .. } => Some("Imported BLAST hit track".to_string()),
            Operation::SelectCandidate { .. } => Some("Selected candidate".to_string()),
            Operation::FilterByMolecularWeight { .. } => {
                Some("Molecular-weight filtered".to_string())
            }
            Operation::FilterByDesignConstraints { .. } => {
                Some("Design constraints filtered".to_string())
            }
            Operation::FilterContainerByMolecularWeight { .. } => {
                Some("Molecular-weight filtered".to_string())
            }
            Operation::DeriveTranscriptSequences { .. } => {
                Some("Derived transcript sequence".to_string())
            }
            Operation::Reverse { .. }
            | Operation::Complement { .. }
            | Operation::ReverseComplement { .. }
            | Operation::Branch { .. } => Some("Derived sequence".to_string()),
            _ => None,
        };
        let _ = self.add_container(&result.created_seq_ids, kind, name, Some(&result.op_id));
    }

    pub(super) fn add_lineage_edges(
        &mut self,
        parent_seq_ids: &[SeqId],
        created_seq_ids: &[SeqId],
        op_id: &str,
        run_id: &str,
    ) {
        if parent_seq_ids.is_empty() || created_seq_ids.is_empty() {
            return;
        }
        let parent_nodes: Vec<NodeId> = parent_seq_ids
            .iter()
            .map(|seq_id| self.ensure_lineage_node(seq_id))
            .collect();
        let child_nodes: Vec<NodeId> = created_seq_ids
            .iter()
            .map(|seq_id| self.ensure_lineage_node(seq_id))
            .collect();
        for from_node_id in &parent_nodes {
            for to_node_id in &child_nodes {
                self.state.lineage.edges.push(LineageEdge {
                    from_node_id: from_node_id.clone(),
                    to_node_id: to_node_id.clone(),
                    op_id: op_id.to_string(),
                    run_id: run_id.to_string(),
                });
            }
        }
    }

    pub(super) fn unique_seq_id(&self, base: &str) -> SeqId {
        if !self.state.sequences.contains_key(base) {
            return base.to_string();
        }
        let mut i = 2usize;
        loop {
            let candidate = format!("{base}_{i}");
            if !self.state.sequences.contains_key(&candidate) {
                return candidate;
            }
            i += 1;
        }
    }

    pub(super) fn prepare_sequence(dna: &mut DNAsequence) {
        Self::prepare_sequence_light(dna);
        dna.update_computed_features();
    }

    pub(super) fn prepare_sequence_light(dna: &mut DNAsequence) {
        *dna.restriction_enzymes_mut() = active_restriction_enzymes();
        dna.set_max_restriction_enzyme_sites(None);
        dna.set_methylation_mode(MethylationMode::both());
    }
}
