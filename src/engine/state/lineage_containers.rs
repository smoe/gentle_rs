//! Lineage graph and container/arrangement helper routines.
//!
//! Container-first project semantics and lineage bookkeeping live together here
//! because most mutations need to keep both models in sync.
//!
//! Look here for:
//! - lineage-node creation and lookup helpers
//! - container/arrangement/rack bookkeeping tied to sequence mutations
//! - macro-instance/container-state synchronization utilities

use super::*;
use std::fs;

impl GentleEngine {
    pub(super) fn normalize_serial_gel_ladders_from_iter<I>(iter: I) -> Vec<String>
    where
        I: IntoIterator<Item = String>,
    {
        let mut ladders = Vec::new();
        let mut seen = HashSet::new();
        for value in iter {
            let trimmed = value.trim();
            if trimmed.is_empty() {
                continue;
            }
            let normalized = trimmed.to_string();
            if seen.insert(normalized.clone()) {
                ladders.push(normalized);
            }
        }
        ladders
    }

    pub(super) fn normalize_serial_gel_ladders_owned(ladders: Option<Vec<String>>) -> Vec<String> {
        Self::normalize_serial_gel_ladders_from_iter(ladders.unwrap_or_default())
    }

    pub(super) fn normalize_serial_gel_ladders_slice(ladders: &[String]) -> Vec<String> {
        Self::normalize_serial_gel_ladders_from_iter(ladders.iter().cloned())
    }

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

    pub(super) fn next_rack_id(&mut self) -> String {
        loop {
            self.state.container_state.next_rack_counter += 1;
            let id = format!("rack-{}", self.state.container_state.next_rack_counter);
            if !self.state.container_state.racks.contains_key(&id) {
                return id;
            }
        }
    }

    fn normalize_arrangement_lane_role_labels(
        lane_container_ids: &[ContainerId],
        lane_role_labels: Option<Vec<String>>,
    ) -> Vec<String> {
        let mut out = lane_container_ids
            .iter()
            .enumerate()
            .map(|(idx, _)| format!("lane_{}", idx + 1))
            .collect::<Vec<_>>();
        if let Some(input) = lane_role_labels {
            for (idx, value) in input.into_iter().enumerate().take(out.len()) {
                let trimmed = value.trim();
                if !trimmed.is_empty() {
                    out[idx] = trimmed.to_string();
                }
            }
        }
        out
    }

    fn arrangement_ladder_pair(arrangement: &Arrangement) -> Option<(String, String)> {
        match arrangement.ladders.as_slice() {
            [] => None,
            [name] => Some((name.clone(), name.clone())),
            [left, right, ..] => Some((left.clone(), right.clone())),
        }
    }

    fn choose_smallest_rack_profile_for_slots(slot_count: usize) -> Result<RackProfileKind, EngineError> {
        for kind in [
            RackProfileKind::SmallTube4x6,
            RackProfileKind::Plate96,
            RackProfileKind::Plate384,
        ] {
            if slot_count <= kind.capacity() {
                return Ok(kind);
            }
        }
        Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "No built-in rack profile can fit {slot_count} occupied positions (max supported: {})",
                RackProfileKind::Plate384.capacity()
            ),
        })
    }

    pub(crate) fn rack_coordinate_from_index(
        profile: &RackProfileSnapshot,
        index: usize,
    ) -> Result<String, EngineError> {
        if index >= profile.capacity() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack slot index {index} is outside profile capacity {}",
                    profile.capacity()
                ),
            });
        }
        let row = index / profile.columns;
        let column = (index % profile.columns) + 1;
        let row_letter = ((b'A' + row as u8) as char).to_string();
        Ok(format!("{row_letter}{column}"))
    }

    pub(crate) fn rack_index_from_coordinate(
        profile: &RackProfileSnapshot,
        coordinate: &str,
    ) -> Result<usize, EngineError> {
        let trimmed = coordinate.trim().to_ascii_uppercase();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Rack coordinate cannot be empty".to_string(),
            });
        }
        let mut letters = String::new();
        let mut digits = String::new();
        for ch in trimmed.chars() {
            if ch.is_ascii_alphabetic() && digits.is_empty() {
                letters.push(ch);
            } else if ch.is_ascii_digit() {
                digits.push(ch);
            } else {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Invalid rack coordinate '{coordinate}'"),
                });
            }
        }
        if letters.len() != 1 || digits.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Invalid rack coordinate '{coordinate}'"),
            });
        }
        let row = (letters.as_bytes()[0] - b'A') as usize;
        let column = digits.parse::<usize>().map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Invalid rack coordinate '{coordinate}': {e}"),
        })?;
        if row >= profile.rows || column == 0 || column > profile.columns {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack coordinate '{}' is outside profile {}x{}",
                    coordinate, profile.rows, profile.columns
                ),
            });
        }
        Ok(row * profile.columns + (column - 1))
    }

    fn sorted_rack_placements(
        &self,
        rack: &Rack,
    ) -> Result<Vec<(usize, RackPlacementEntry)>, EngineError> {
        let mut rows = Vec::with_capacity(rack.placements.len());
        for entry in &rack.placements {
            let index = Self::rack_index_from_coordinate(&rack.profile, &entry.coordinate)?;
            rows.push((index, entry.clone()));
        }
        rows.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.order_index.cmp(&b.1.order_index)));
        Ok(rows)
    }

    fn reflow_rack_placements(
        profile: &RackProfileSnapshot,
        entries: &mut [RackPlacementEntry],
    ) -> Result<(), EngineError> {
        for (idx, entry) in entries.iter_mut().enumerate() {
            entry.coordinate = Self::rack_coordinate_from_index(profile, idx)?;
        }
        Ok(())
    }

    fn arrangement_rack_payload(
        &self,
        arrangement_id: &str,
    ) -> Result<Vec<(RackOccupant, String)>, EngineError> {
        let arrangement = self
            .state
            .container_state
            .arrangements
            .get(arrangement_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Arrangement '{arrangement_id}' not found"),
            })?;
        let mut payload = vec![];
        if let Some((left, right)) = Self::arrangement_ladder_pair(arrangement) {
            payload.push((
                RackOccupant::LadderReference { ladder_name: left },
                "ladder_left".to_string(),
            ));
            for (idx, container_id) in arrangement.lane_container_ids.iter().enumerate() {
                payload.push((
                    RackOccupant::Container {
                        container_id: container_id.clone(),
                    },
                    arrangement
                        .lane_role_labels
                        .get(idx)
                        .cloned()
                        .unwrap_or_else(|| format!("lane_{}", idx + 1)),
                ));
            }
            payload.push((
                RackOccupant::LadderReference { ladder_name: right },
                "ladder_right".to_string(),
            ));
        } else {
            for (idx, container_id) in arrangement.lane_container_ids.iter().enumerate() {
                payload.push((
                    RackOccupant::Container {
                        container_id: container_id.clone(),
                    },
                    arrangement
                        .lane_role_labels
                        .get(idx)
                        .cloned()
                        .unwrap_or_else(|| format!("lane_{}", idx + 1)),
                ));
            }
        }
        Ok(payload)
    }

    pub(super) fn add_serial_arrangement(
        &mut self,
        container_ids: &[ContainerId],
        arrangement_id: Option<String>,
        name: Option<String>,
        ladders: Option<Vec<String>>,
        lane_role_labels: Option<Vec<String>>,
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
            ladders: Self::normalize_serial_gel_ladders_owned(ladders),
            lane_role_labels: Self::normalize_arrangement_lane_role_labels(
                &container_ids.iter().map(|id| id.trim().to_string()).filter(|id| !id.is_empty()).collect::<Vec<_>>(),
                lane_role_labels,
            ),
            default_rack_id: None,
            created_by_op: created_by_op.map(ToString::to_string),
            created_at_unix_ms: Self::now_unix_ms(),
        };
        self.state
            .container_state
            .arrangements
            .insert(arrangement_id.clone(), arrangement);
        if let Err(err) =
            self.create_rack_from_arrangement(&arrangement_id, None, None, None, created_by_op, true)
        {
            self.state
                .container_state
                .arrangements
                .remove(&arrangement_id);
            return Err(err);
        }
        Ok(arrangement_id)
    }

    pub(super) fn create_rack_from_arrangement(
        &mut self,
        arrangement_id: &str,
        rack_id: Option<String>,
        name: Option<String>,
        profile: Option<RackProfileKind>,
        created_by_op: Option<&str>,
        set_as_default: bool,
    ) -> Result<String, EngineError> {
        let arrangement_id = arrangement_id.trim();
        if arrangement_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "arrangement_id cannot be empty".to_string(),
            });
        }
        let arrangement = self
            .state
            .container_state
            .arrangements
            .get(arrangement_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Arrangement '{arrangement_id}' not found"),
            })?;
        let payload = self.arrangement_rack_payload(arrangement_id)?;
        let profile_kind = match profile {
            Some(profile) => {
                if payload.len() > profile.capacity() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Rack profile '{}' cannot fit arrangement '{}' with {} occupied positions",
                            profile.as_str(),
                            arrangement_id,
                            payload.len()
                        ),
                    });
                }
                profile
            }
            None => Self::choose_smallest_rack_profile_for_slots(payload.len())?,
        };
        let rack_id = rack_id
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
            .unwrap_or_else(|| self.next_rack_id());
        if self.state.container_state.racks.contains_key(&rack_id) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Rack '{rack_id}' already exists"),
            });
        }
        let profile_snapshot = RackProfileSnapshot::from_kind(profile_kind);
        let placements = payload
            .into_iter()
            .enumerate()
            .map(|(order_index, (occupant, role_label))| {
                Ok(RackPlacementEntry {
                    coordinate: Self::rack_coordinate_from_index(&profile_snapshot, order_index)?,
                    occupant: Some(occupant),
                    arrangement_id: arrangement_id.to_string(),
                    order_index,
                    role_label,
                })
            })
            .collect::<Result<Vec<_>, EngineError>>()?;
        let rack_name = name
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
            .or_else(|| {
                arrangement.name.as_ref().map(|name| {
                    format!("Rack draft: {}", name.trim())
                })
            })
            .unwrap_or_else(|| format!("Rack draft: {arrangement_id}"));
        self.state.container_state.racks.insert(
            rack_id.clone(),
            Rack {
                rack_id: rack_id.clone(),
                name: rack_name,
                profile: profile_snapshot,
                placements,
                created_by_op: created_by_op.map(ToString::to_string),
                created_at_unix_ms: Self::now_unix_ms(),
            },
        );
        if set_as_default {
            if let Some(arrangement) = self.state.container_state.arrangements.get_mut(arrangement_id)
            {
                arrangement.default_rack_id = Some(rack_id.clone());
            }
        }
        Ok(rack_id)
    }

    pub(super) fn place_arrangement_on_rack(
        &mut self,
        arrangement_id: &str,
        rack_id: &str,
    ) -> Result<usize, EngineError> {
        let arrangement_id = arrangement_id.trim();
        let rack_id = rack_id.trim();
        if arrangement_id.is_empty() || rack_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "arrangement_id and rack_id must be non-empty".to_string(),
            });
        }
        let payload = self.arrangement_rack_payload(arrangement_id)?;
        let rack = self
            .state
            .container_state
            .racks
            .get_mut(rack_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Rack '{rack_id}' not found"),
            })?;
        if rack
            .placements
            .iter()
            .any(|entry| entry.arrangement_id == arrangement_id)
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack '{}' already contains arrangement '{}'",
                    rack_id, arrangement_id
                ),
            });
        }
        let start_index = rack.placements.len();
        if start_index + payload.len() > rack.profile.capacity() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack '{}' does not have enough free positions for arrangement '{}' (need {}, remaining {})",
                    rack_id,
                    arrangement_id,
                    payload.len(),
                    rack.profile.capacity().saturating_sub(start_index)
                ),
            });
        }
        for (offset, (occupant, role_label)) in payload.into_iter().enumerate() {
            rack.placements.push(RackPlacementEntry {
                coordinate: Self::rack_coordinate_from_index(&rack.profile, start_index + offset)?,
                occupant: Some(occupant),
                arrangement_id: arrangement_id.to_string(),
                order_index: offset,
                role_label,
            });
        }
        if let Some(arrangement) = self.state.container_state.arrangements.get_mut(arrangement_id) {
            arrangement
                .default_rack_id
                .get_or_insert_with(|| rack_id.to_string());
        }
        Ok(start_index)
    }

    pub(super) fn move_rack_placement(
        &mut self,
        rack_id: &str,
        from_coordinate: &str,
        to_coordinate: &str,
        move_block: bool,
    ) -> Result<(), EngineError> {
        let rack_id = rack_id.trim();
        let from_coordinate = from_coordinate.trim();
        let to_coordinate = to_coordinate.trim();
        let rack = self
            .state
            .container_state
            .racks
            .get(rack_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Rack '{rack_id}' not found"),
            })?;
        let from_index = Self::rack_index_from_coordinate(&rack.profile, from_coordinate)?;
        let to_index = Self::rack_index_from_coordinate(&rack.profile, to_coordinate)?;
        let mut ordered = self.sorted_rack_placements(&rack)?;
        let from_pos = ordered
            .iter()
            .position(|(idx, _)| *idx == from_index)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Rack '{}' has no occupied position at '{}'",
                    rack_id, from_coordinate
                ),
            })?;
        if move_block {
            let arrangement_id = ordered[from_pos].1.arrangement_id.clone();
            let block_positions = ordered
                .iter()
                .enumerate()
                .filter_map(|(pos, (_, entry))| {
                    if entry.arrangement_id == arrangement_id {
                        Some(pos)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            let block_start = *block_positions.first().unwrap_or(&from_pos);
            let block_len = block_positions.len();
            let block_entries = ordered
                .iter()
                .filter(|(_, entry)| entry.arrangement_id == arrangement_id)
                .map(|(_, entry)| entry.clone())
                .collect::<Vec<_>>();
            ordered.retain(|(_, entry)| entry.arrangement_id != arrangement_id);
            let insertion_index = to_index.min(ordered.len());
            if ordered.len() + block_len > rack.profile.capacity() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Rack '{}' profile '{}' cannot fit moved arrangement block",
                        rack_id,
                        rack.profile.kind.as_str()
                    ),
                });
            }
            for (offset, entry) in block_entries.into_iter().enumerate() {
                ordered.insert(insertion_index + offset, (block_start + offset, entry));
            }
        } else {
            let arrangement_id = ordered[from_pos].1.arrangement_id.clone();
            let block_positions = ordered
                .iter()
                .enumerate()
                .filter_map(|(pos, (_, entry))| {
                    if entry.arrangement_id == arrangement_id {
                        Some(pos)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            let block_start = *block_positions.first().unwrap_or(&from_pos);
            let block_end = *block_positions.last().unwrap_or(&from_pos);
            if to_index < block_start || to_index > block_end {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Sample moves must stay within arrangement block '{}'",
                        arrangement_id
                    ),
                });
            }
            let mut block_entries = ordered[block_start..=block_end]
                .iter()
                .map(|(_, entry)| entry.clone())
                .collect::<Vec<_>>();
            let local_from = from_pos - block_start;
            let local_to = to_index - block_start;
            let moved = block_entries.remove(local_from);
            block_entries.insert(local_to, moved);
            for (offset, entry) in block_entries.into_iter().enumerate() {
                ordered[block_start + offset].1 = entry;
            }
        }
        let mut reflowed = ordered
            .into_iter()
            .map(|(_, entry)| entry)
            .collect::<Vec<_>>();
        let mut order_by_arrangement: HashMap<String, usize> = HashMap::new();
        for entry in &mut reflowed {
            let next = order_by_arrangement
                .entry(entry.arrangement_id.clone())
                .and_modify(|idx| *idx += 1)
                .or_insert(0usize);
            entry.order_index = *next;
        }
        Self::reflow_rack_placements(&rack.profile, &mut reflowed)?;
        if let Some(target_rack) = self.state.container_state.racks.get_mut(rack_id) {
            target_rack.placements = reflowed;
        }
        Ok(())
    }

    pub(super) fn set_rack_profile(
        &mut self,
        rack_id: &str,
        profile: RackProfileKind,
    ) -> Result<(), EngineError> {
        let existing_rack = self
            .state
            .container_state
            .racks
            .get(rack_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Rack '{rack_id}' not found"),
            })?;
        let mut placements = self
            .sorted_rack_placements(&existing_rack)?
            .into_iter()
            .map(|(_, entry)| entry)
            .collect::<Vec<_>>();
        let rack = self
            .state
            .container_state
            .racks
            .get_mut(rack_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Rack '{rack_id}' not found"),
            })?;
        let new_profile = RackProfileSnapshot::from_kind(profile);
        if rack.placements.len() > new_profile.capacity() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack '{}' cannot switch to profile '{}' because {} positions are occupied",
                    rack_id,
                    profile.as_str(),
                    rack.placements.len()
                ),
            });
        }
        Self::reflow_rack_placements(&new_profile, &mut placements)?;
        rack.profile = new_profile;
        rack.placements = placements;
        Ok(())
    }

    fn xml_escape(raw: &str) -> String {
        raw.replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;")
            .replace('"', "&quot;")
            .replace('\'', "&apos;")
    }

    pub(super) fn export_rack_labels_svg(
        &self,
        rack_id: &str,
        arrangement_id: Option<&str>,
        path: &str,
    ) -> Result<usize, EngineError> {
        let rack = self
            .state
            .container_state
            .racks
            .get(rack_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Rack '{rack_id}' not found"),
            })?;
        let arrangement_filter = arrangement_id.map(str::trim).filter(|value| !value.is_empty());
        let rows = self
            .sorted_rack_placements(rack)?
            .into_iter()
            .map(|(_, entry)| entry)
            .filter(|entry| {
                arrangement_filter
                    .map(|arrangement_id| entry.arrangement_id == arrangement_id)
                    .unwrap_or(true)
            })
            .collect::<Vec<_>>();
        if rows.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack '{}' has no placements matching the requested arrangement scope",
                    rack_id
                ),
            });
        }
        let label_width = 280.0f32;
        let label_height = 120.0f32;
        let columns = 2usize;
        let margin = 18.0f32;
        let svg_width = margin * 2.0 + columns as f32 * label_width;
        let svg_height =
            margin * 2.0 + ((rows.len() + columns - 1) / columns) as f32 * label_height;
        let mut svg = format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{svg_width}\" height=\"{svg_height}\" viewBox=\"0 0 {svg_width} {svg_height}\">"
        );
        svg.push_str(
            "<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>",
        );
        for (idx, entry) in rows.iter().enumerate() {
            let col = idx % columns;
            let row = idx / columns;
            let x = margin + col as f32 * label_width;
            let y = margin + row as f32 * label_height;
            let mut lines = vec![
                format!("{}  {}", rack.rack_id, entry.coordinate),
                format!("role: {}", entry.role_label),
            ];
            match entry.occupant.as_ref() {
                Some(RackOccupant::Container { container_id }) => {
                    lines.push(format!("container: {container_id}"));
                    if let Some(container) = self.state.container_state.containers.get(container_id) {
                        if let Some(name) = container.name.as_ref().filter(|name| !name.trim().is_empty()) {
                            lines.push(name.clone());
                        }
                        if let Some(seq_id) = container.members.first() {
                            lines.push(format!("seq: {seq_id}"));
                            if let Some(dna) = self.state.sequences.get(seq_id) {
                                let topology = if dna.is_circular() { "circular" } else { "linear" };
                                lines.push(format!("{} bp | {}", dna.len(), topology));
                            }
                        }
                        if let Some(op_id) = container.created_by_op.as_ref() {
                            lines.push(format!("op: {op_id}"));
                        }
                    }
                }
                Some(RackOccupant::LadderReference { ladder_name }) => {
                    lines.push(format!("ladder: {ladder_name}"));
                }
                None => lines.push("empty".to_string()),
            }
            if let Some(arrangement) = self.state.container_state.arrangements.get(&entry.arrangement_id) {
                lines.push(format!(
                    "arrangement: {}",
                    arrangement
                        .name
                        .as_ref()
                        .filter(|name| !name.trim().is_empty())
                        .cloned()
                        .unwrap_or_else(|| arrangement.arrangement_id.clone())
                ));
                if let Some(op_id) = arrangement.created_by_op.as_ref() {
                    lines.push(format!("origin: {op_id}"));
                }
            }
            svg.push_str(&format!(
                "<rect x=\"{x}\" y=\"{y}\" width=\"{}\" height=\"{}\" rx=\"8\" ry=\"8\" fill=\"#f7fbfd\" stroke=\"#5d8aa8\" stroke-width=\"1.2\"/>",
                label_width - 12.0,
                label_height - 12.0
            ));
            for (line_idx, line) in lines.iter().enumerate() {
                let font_size = if line_idx == 0 { 15 } else { 13 };
                let font_weight = if line_idx == 0 { "700" } else { "400" };
                let line_y = y + 24.0 + line_idx as f32 * 15.0;
                svg.push_str(&format!(
                    "<text x=\"{}\" y=\"{}\" font-family=\"monospace\" font-size=\"{}\" font-weight=\"{}\" fill=\"#173042\">{}</text>",
                    x + 12.0,
                    line_y,
                    font_size,
                    font_weight,
                    Self::xml_escape(line)
                ));
            }
        }
        svg.push_str("</svg>");
        fs::write(path, svg).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write rack labels SVG '{}': {e}", path),
        })?;
        Ok(rows.len())
    }

    pub(super) fn set_arrangement_ladders(
        &mut self,
        arrangement_id: &str,
        ladders: Option<Vec<String>>,
    ) -> Result<Vec<String>, EngineError> {
        let arrangement_id = arrangement_id.trim();
        if arrangement_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "arrangement_id cannot be empty".to_string(),
            });
        }
        let arrangement = self
            .state
            .container_state
            .arrangements
            .get_mut(arrangement_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Arrangement '{arrangement_id}' not found"),
            })?;
        if arrangement.mode != ArrangementMode::Serial {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Arrangement '{}' is mode '{:?}', only serial arrangements can update gel ladders",
                    arrangement_id, arrangement.mode
                ),
            });
        }
        let normalized = Self::normalize_serial_gel_ladders_owned(ladders);
        arrangement.ladders = normalized.clone();
        Ok(normalized)
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
