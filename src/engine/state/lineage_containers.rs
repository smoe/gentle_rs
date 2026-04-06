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

    fn choose_smallest_rack_profile_for_slots(
        slot_count: usize,
    ) -> Result<RackProfileKind, EngineError> {
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

    fn validate_custom_rack_profile_dimensions(
        rows: usize,
        columns: usize,
    ) -> Result<(), EngineError> {
        if rows == 0 || columns == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Custom rack profiles require rows >= 1 and columns >= 1".to_string(),
            });
        }
        Ok(())
    }

    pub(crate) fn rack_coordinate_from_index(
        profile: &RackProfileSnapshot,
        index: usize,
    ) -> Result<String, EngineError> {
        let available = Self::rack_available_coordinates_in_fill_order(profile)?;
        if index >= available.len() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack slot index {index} is outside available rack capacity {}",
                    available.len()
                ),
            });
        }
        Ok(available[index].clone())
    }

    pub(crate) fn rack_row_label_from_index(mut row: usize) -> String {
        let mut out = String::new();
        row += 1;
        while row > 0 {
            let rem = ((row - 1) % 26) as u8;
            out.push((b'A' + rem) as char);
            row = (row - 1) / 26;
        }
        out.chars().rev().collect()
    }

    fn rack_row_index_from_label(label: &str) -> Result<usize, EngineError> {
        let trimmed = label.trim().to_ascii_uppercase();
        if trimmed.is_empty() || !trimmed.chars().all(|ch| ch.is_ascii_uppercase()) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Invalid rack row label '{label}'"),
            });
        }
        let mut value = 0usize;
        for ch in trimmed.chars() {
            value = value
                .checked_mul(26)
                .and_then(|acc| acc.checked_add((ch as u8 - b'A' + 1) as usize))
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Rack row label '{}' is too large", label),
                })?;
        }
        Ok(value.saturating_sub(1))
    }

    fn rack_parse_row_column_label(coordinate: &str) -> Result<(usize, usize), EngineError> {
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
        if letters.is_empty() || digits.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Invalid rack coordinate '{coordinate}'"),
            });
        }
        let row = Self::rack_row_index_from_label(&letters)?;
        let column = digits.parse::<usize>().map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Invalid rack coordinate '{coordinate}': {e}"),
        })?;
        if column == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Invalid rack coordinate '{coordinate}'"),
            });
        }
        Ok((row, column - 1))
    }

    pub(crate) fn rack_coordinate_from_row_column(
        profile: &RackProfileSnapshot,
        row: usize,
        column: usize,
    ) -> Result<String, EngineError> {
        if row >= profile.rows || column >= profile.columns {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack position ({}, {}) is outside profile {}x{}",
                    row, column, profile.rows, profile.columns
                ),
            });
        }
        Ok(format!(
            "{}{}",
            Self::rack_row_label_from_index(row),
            column + 1
        ))
    }

    fn rack_fill_ordinal(profile: &RackProfileSnapshot, row: usize, column: usize) -> usize {
        match profile.fill_direction {
            RackFillDirection::RowMajor => row * profile.columns + column,
            RackFillDirection::ColumnMajor => column * profile.rows + row,
        }
    }

    fn normalized_blocked_coordinates_for_profile(
        profile: &RackProfileSnapshot,
        blocked_coordinates: &[String],
        drop_out_of_bounds: bool,
    ) -> Result<Vec<String>, EngineError> {
        let mut rows = Vec::new();
        let mut seen = HashSet::new();
        for raw in blocked_coordinates {
            let trimmed = raw.trim();
            if trimmed.is_empty() {
                continue;
            }
            let (row, column) = Self::rack_parse_row_column_label(trimmed)?;
            if row >= profile.rows || column >= profile.columns {
                if drop_out_of_bounds {
                    continue;
                }
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Blocked rack coordinate '{}' is outside profile {}x{}",
                        trimmed, profile.rows, profile.columns
                    ),
                });
            }
            let canonical = Self::rack_coordinate_from_row_column(profile, row, column)?;
            if seen.insert(canonical.clone()) {
                rows.push((Self::rack_fill_ordinal(profile, row, column), canonical));
            }
        }
        rows.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        Ok(rows.into_iter().map(|(_, coordinate)| coordinate).collect())
    }

    pub(crate) fn rack_available_coordinates_in_fill_order(
        profile: &RackProfileSnapshot,
    ) -> Result<Vec<String>, EngineError> {
        let blocked = Self::normalized_blocked_coordinates_for_profile(
            profile,
            &profile.blocked_coordinates,
            false,
        )?
        .into_iter()
        .collect::<HashSet<_>>();
        let mut out = Vec::with_capacity(profile.capacity().saturating_sub(blocked.len()));
        match profile.fill_direction {
            RackFillDirection::RowMajor => {
                for row in 0..profile.rows {
                    for column in 0..profile.columns {
                        let coordinate =
                            Self::rack_coordinate_from_row_column(profile, row, column)?;
                        if !blocked.contains(&coordinate) {
                            out.push(coordinate);
                        }
                    }
                }
            }
            RackFillDirection::ColumnMajor => {
                for column in 0..profile.columns {
                    for row in 0..profile.rows {
                        let coordinate =
                            Self::rack_coordinate_from_row_column(profile, row, column)?;
                        if !blocked.contains(&coordinate) {
                            out.push(coordinate);
                        }
                    }
                }
            }
        }
        Ok(out)
    }

    pub(crate) fn rack_index_from_coordinate(
        profile: &RackProfileSnapshot,
        coordinate: &str,
    ) -> Result<usize, EngineError> {
        let (row, column) = Self::rack_parse_row_column_label(coordinate)?;
        if row >= profile.rows || column >= profile.columns {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack coordinate '{}' is outside profile {}x{}",
                    coordinate, profile.rows, profile.columns
                ),
            });
        }
        let coordinate = Self::rack_coordinate_from_row_column(profile, row, column)?;
        let available = Self::rack_available_coordinates_in_fill_order(profile)?;
        available
            .iter()
            .position(|candidate| *candidate == coordinate)
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack coordinate '{}' is blocked or unavailable in profile {}x{}",
                    coordinate, profile.rows, profile.columns
                ),
            })
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
        let available = Self::rack_available_coordinates_in_fill_order(profile)?;
        if entries.len() > available.len() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Profile {}x{} with current blocked coordinates can fit only {} occupied positions",
                    profile.rows,
                    profile.columns,
                    available.len()
                ),
            });
        }
        for (idx, entry) in entries.iter_mut().enumerate() {
            entry.coordinate = available[idx].clone();
        }
        Ok(())
    }

    fn reproject_rack_with_profile(
        &mut self,
        rack_id: &str,
        new_profile: RackProfileSnapshot,
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
        let available_capacity =
            Self::rack_available_coordinates_in_fill_order(&new_profile)?.len();
        if placements.len() > available_capacity {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack '{}' cannot switch to the requested profile because {} positions are occupied but only {} remain available",
                    rack_id,
                    placements.len(),
                    available_capacity
                ),
            });
        }
        Self::reflow_rack_placements(&new_profile, &mut placements)?;
        let rack = self
            .state
            .container_state
            .racks
            .get_mut(rack_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Rack '{rack_id}' not found"),
            })?;
        rack.profile = new_profile;
        rack.placements = placements;
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
                &container_ids
                    .iter()
                    .map(|id| id.trim().to_string())
                    .filter(|id| !id.is_empty())
                    .collect::<Vec<_>>(),
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
        if let Err(err) = self.create_rack_from_arrangement(
            &arrangement_id,
            None,
            None,
            None,
            created_by_op,
            true,
        ) {
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
                arrangement
                    .name
                    .as_ref()
                    .map(|name| format!("Rack draft: {}", name.trim()))
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
            if let Some(arrangement) = self
                .state
                .container_state
                .arrangements
                .get_mut(arrangement_id)
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
        let available_capacity =
            Self::rack_available_coordinates_in_fill_order(&rack.profile)?.len();
        if start_index + payload.len() > available_capacity {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack '{}' does not have enough free positions for arrangement '{}' (need {}, remaining {})",
                    rack_id,
                    arrangement_id,
                    payload.len(),
                    available_capacity.saturating_sub(start_index)
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
        if let Some(arrangement) = self
            .state
            .container_state
            .arrangements
            .get_mut(arrangement_id)
        {
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
            let selected_entries = ordered
                .iter()
                .filter(|(_, entry)| entry.arrangement_id == arrangement_id)
                .map(|(_, entry)| entry.clone())
                .collect::<Vec<_>>();
            let remaining = ordered
                .iter()
                .filter(|(_, entry)| entry.arrangement_id != arrangement_id)
                .map(|(idx, entry)| (*idx, entry.clone()))
                .collect::<Vec<_>>();
            let insertion_index =
                Self::rack_block_insertion_index_after_target(&remaining, to_index);
            let mut reflowed = Vec::with_capacity(selected_entries.len() + remaining.len());
            reflowed.extend(
                remaining[..insertion_index]
                    .iter()
                    .map(|(_, entry)| entry.clone()),
            );
            reflowed.extend(selected_entries);
            reflowed.extend(
                remaining[insertion_index..]
                    .iter()
                    .map(|(_, entry)| entry.clone()),
            );
            let ordered = reflowed
                .into_iter()
                .enumerate()
                .map(|(index, entry)| (index, entry))
                .collect::<Vec<_>>();
            return self.reflow_ordered_rack_entries(rack_id, &rack.profile, ordered);
        } else {
            let arrangement_id = ordered[from_pos].1.arrangement_id.clone();
            let target_pos = ordered
                .iter()
                .position(|(idx, _)| *idx == to_index)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Sample moves must target an occupied coordinate within arrangement block '{}'",
                        arrangement_id
                    ),
                })?;
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
            if target_pos < block_start || target_pos > block_end {
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
            let local_to = target_pos - block_start;
            let moved = block_entries.remove(local_from);
            block_entries.insert(local_to, moved);
            for (offset, entry) in block_entries.into_iter().enumerate() {
                ordered[block_start + offset].1 = entry;
            }
        }
        self.reflow_ordered_rack_entries(rack_id, &rack.profile, ordered)
    }

    fn rack_selected_arrangement_ids_in_order(
        ordered: &[(usize, RackPlacementEntry)],
        requested_arrangement_ids: &[String],
    ) -> Vec<String> {
        let requested = requested_arrangement_ids
            .iter()
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .collect::<HashSet<_>>();
        let mut seen = HashSet::new();
        let mut out = Vec::new();
        for (_, entry) in ordered {
            if requested.contains(entry.arrangement_id.as_str())
                && seen.insert(entry.arrangement_id.clone())
            {
                out.push(entry.arrangement_id.clone());
            }
        }
        out
    }

    fn rack_block_insertion_index(
        remaining: &[(usize, RackPlacementEntry)],
        target_index: usize,
    ) -> usize {
        let mut idx = 0usize;
        while idx < remaining.len() {
            let block_start_index = remaining[idx].0;
            if target_index <= block_start_index {
                return idx;
            }
            let arrangement_id = remaining[idx].1.arrangement_id.clone();
            while idx < remaining.len() && remaining[idx].1.arrangement_id == arrangement_id {
                idx += 1;
            }
        }
        remaining.len()
    }

    fn rack_block_insertion_index_after_target(
        remaining: &[(usize, RackPlacementEntry)],
        target_index: usize,
    ) -> usize {
        let mut idx = 0usize;
        while idx < remaining.len() {
            let block_start_index = remaining[idx].0;
            if target_index < block_start_index {
                return idx;
            }
            let arrangement_id = remaining[idx].1.arrangement_id.clone();
            let mut block_end_index = block_start_index;
            while idx < remaining.len() && remaining[idx].1.arrangement_id == arrangement_id {
                block_end_index = remaining[idx].0;
                idx += 1;
            }
            if target_index <= block_end_index {
                return idx;
            }
        }
        remaining.len()
    }

    fn rack_group_insertion_index_within_block(local_to: usize, remaining_len: usize) -> usize {
        local_to.min(remaining_len)
    }

    fn reflow_ordered_rack_entries(
        &mut self,
        rack_id: &str,
        profile: &RackProfileSnapshot,
        ordered: Vec<(usize, RackPlacementEntry)>,
    ) -> Result<(), EngineError> {
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
        Self::reflow_rack_placements(profile, &mut reflowed)?;
        if let Some(target_rack) = self.state.container_state.racks.get_mut(rack_id) {
            target_rack.placements = reflowed;
        }
        Ok(())
    }

    pub(super) fn move_rack_samples(
        &mut self,
        rack_id: &str,
        from_coordinates: &[String],
        to_coordinate: &str,
    ) -> Result<(), EngineError> {
        let rack_id = rack_id.trim();
        let to_coordinate = to_coordinate.trim();
        if from_coordinates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Rack sample move requires at least one source coordinate".to_string(),
            });
        }
        let requested_coordinates = from_coordinates
            .iter()
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
            .collect::<BTreeSet<_>>();
        if requested_coordinates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Rack sample move requires non-empty source coordinates".to_string(),
            });
        }
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
        let ordered = self.sorted_rack_placements(&rack)?;
        let selected = ordered
            .iter()
            .filter_map(|(index, entry)| {
                if requested_coordinates.contains(entry.coordinate.as_str()) {
                    Some((*index, entry.clone()))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        if selected.len() != requested_coordinates.len() {
            let found = selected
                .iter()
                .map(|(_, entry)| entry.coordinate.clone())
                .collect::<BTreeSet<_>>();
            let missing = requested_coordinates
                .difference(&found)
                .cloned()
                .collect::<Vec<_>>();
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Rack '{}' does not contain requested sample coordinate(s): {}",
                    rack_id,
                    missing.join(", ")
                ),
            });
        }
        let arrangement_id = selected
            .first()
            .map(|(_, entry)| entry.arrangement_id.clone())
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Rack '{}' does not contain requested sample coordinates",
                    rack_id
                ),
            })?;
        if selected
            .iter()
            .any(|(_, entry)| entry.arrangement_id != arrangement_id)
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Multi-sample rack moves must stay within one arrangement block"
                    .to_string(),
            });
        }
        let to_index = Self::rack_index_from_coordinate(&rack.profile, to_coordinate)?;
        let target_pos = ordered
            .iter()
            .position(|(index, _)| *index == to_index)
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Multi-sample moves must target an occupied coordinate within arrangement block '{}'",
                    arrangement_id
                ),
            })?;
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
        let from_pos = ordered
            .iter()
            .position(|(index, _)| *index == selected[0].0)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Rack '{}' does not contain requested sample coordinates",
                    rack_id
                ),
            })?;
        let block_start = *block_positions.first().unwrap_or(&from_pos);
        let block_end = *block_positions.last().unwrap_or(&from_pos);
        if target_pos < block_start || target_pos > block_end {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Multi-sample moves must stay within arrangement block '{}'",
                    arrangement_id
                ),
            });
        }
        let mut block_entries = ordered[block_start..=block_end]
            .iter()
            .map(|(_, entry)| entry.clone())
            .collect::<Vec<_>>();
        let selected_local_positions = selected
            .iter()
            .filter_map(|(index, _)| {
                ordered[block_start..=block_end]
                    .iter()
                    .position(|(candidate_index, _)| candidate_index == index)
            })
            .collect::<Vec<_>>();
        let selected_local_set = selected_local_positions
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        let selected_entries = block_entries
            .iter()
            .enumerate()
            .filter_map(|(idx, entry)| {
                if selected_local_set.contains(&idx) {
                    Some(entry.clone())
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        let remaining_entries = block_entries
            .drain(..)
            .enumerate()
            .filter_map(|(idx, entry)| {
                if selected_local_set.contains(&idx) {
                    None
                } else {
                    Some(entry)
                }
            })
            .collect::<Vec<_>>();
        let local_to = target_pos - block_start;
        let insertion_index =
            Self::rack_group_insertion_index_within_block(local_to, remaining_entries.len());
        let mut reordered_block = Vec::with_capacity(block_end - block_start + 1);
        reordered_block.extend(remaining_entries[..insertion_index].iter().cloned());
        reordered_block.extend(selected_entries);
        reordered_block.extend(remaining_entries[insertion_index..].iter().cloned());
        let mut reordered = ordered;
        for (offset, entry) in reordered_block.into_iter().enumerate() {
            reordered[block_start + offset].1 = entry;
        }
        self.reflow_ordered_rack_entries(rack_id, &rack.profile, reordered)
    }

    pub(super) fn move_rack_arrangement_blocks(
        &mut self,
        rack_id: &str,
        arrangement_ids: &[String],
        to_coordinate: &str,
    ) -> Result<(), EngineError> {
        let rack_id = rack_id.trim();
        let to_coordinate = to_coordinate.trim();
        if arrangement_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Rack block move requires at least one arrangement id".to_string(),
            });
        }
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
        let to_index = Self::rack_index_from_coordinate(&rack.profile, to_coordinate)?;
        let ordered = self.sorted_rack_placements(&rack)?;
        let selected_arrangement_ids =
            Self::rack_selected_arrangement_ids_in_order(&ordered, arrangement_ids);
        if selected_arrangement_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Rack '{}' does not contain any of the requested arrangement blocks",
                    rack_id
                ),
            });
        }
        let selected_id_set = selected_arrangement_ids
            .iter()
            .cloned()
            .collect::<HashSet<_>>();
        let selected_entries = ordered
            .iter()
            .filter(|(_, entry)| selected_id_set.contains(&entry.arrangement_id))
            .map(|(_, entry)| entry.clone())
            .collect::<Vec<_>>();
        let remaining = ordered
            .iter()
            .filter(|(_, entry)| !selected_id_set.contains(&entry.arrangement_id))
            .map(|(idx, entry)| (*idx, entry.clone()))
            .collect::<Vec<_>>();
        let insertion_index = Self::rack_block_insertion_index(&remaining, to_index);
        let mut reflowed = Vec::with_capacity(selected_entries.len() + remaining.len());
        reflowed.extend(
            remaining[..insertion_index]
                .iter()
                .map(|(_, entry)| entry.clone()),
        );
        reflowed.extend(selected_entries);
        reflowed.extend(
            remaining[insertion_index..]
                .iter()
                .map(|(_, entry)| entry.clone()),
        );
        let ordered = reflowed
            .into_iter()
            .enumerate()
            .map(|(index, entry)| (index, entry))
            .collect::<Vec<_>>();
        self.reflow_ordered_rack_entries(rack_id, &rack.profile, ordered)
    }

    pub(super) fn set_rack_profile(
        &mut self,
        rack_id: &str,
        profile: RackProfileKind,
    ) -> Result<(), EngineError> {
        if profile == RackProfileKind::Custom {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Use SetRackProfileCustom for custom row/column rack dimensions"
                    .to_string(),
            });
        }
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
        let mut new_profile = RackProfileSnapshot::from_kind(profile);
        new_profile.fill_direction = existing_rack.profile.fill_direction;
        new_profile.blocked_coordinates = Self::normalized_blocked_coordinates_for_profile(
            &new_profile,
            &existing_rack.profile.blocked_coordinates,
            true,
        )?;
        self.reproject_rack_with_profile(rack_id, new_profile)
    }

    pub(super) fn set_rack_profile_custom(
        &mut self,
        rack_id: &str,
        rows: usize,
        columns: usize,
    ) -> Result<(), EngineError> {
        Self::validate_custom_rack_profile_dimensions(rows, columns)?;
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
        let mut new_profile = RackProfileSnapshot::custom(rows, columns);
        new_profile.fill_direction = existing_rack.profile.fill_direction;
        new_profile.blocked_coordinates = Self::normalized_blocked_coordinates_for_profile(
            &new_profile,
            &existing_rack.profile.blocked_coordinates,
            true,
        )?;
        self.reproject_rack_with_profile(rack_id, new_profile)
    }

    fn rack_template_blocked_coordinates(
        profile: &RackProfileSnapshot,
        template: RackAuthoringTemplate,
    ) -> Result<Vec<String>, EngineError> {
        match template {
            RackAuthoringTemplate::BenchRows | RackAuthoringTemplate::PlateColumns => Ok(vec![]),
            RackAuthoringTemplate::PlateEdgeAvoidance => {
                if profile.rows < 3 || profile.columns < 3 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Rack template '{}' requires at least a 3x3 profile, found {}x{}",
                            template.as_str(),
                            profile.rows,
                            profile.columns
                        ),
                    });
                }
                let mut blocked = Vec::new();
                match profile.fill_direction {
                    RackFillDirection::RowMajor => {
                        for row in 0..profile.rows {
                            for column in 0..profile.columns {
                                if row == 0
                                    || row + 1 == profile.rows
                                    || column == 0
                                    || column + 1 == profile.columns
                                {
                                    blocked.push(Self::rack_coordinate_from_row_column(
                                        profile, row, column,
                                    )?);
                                }
                            }
                        }
                    }
                    RackFillDirection::ColumnMajor => {
                        for column in 0..profile.columns {
                            for row in 0..profile.rows {
                                if row == 0
                                    || row + 1 == profile.rows
                                    || column == 0
                                    || column + 1 == profile.columns
                                {
                                    blocked.push(Self::rack_coordinate_from_row_column(
                                        profile, row, column,
                                    )?);
                                }
                            }
                        }
                    }
                }
                Self::normalized_blocked_coordinates_for_profile(profile, &blocked, false)
            }
        }
    }

    pub(super) fn apply_rack_template(
        &mut self,
        rack_id: &str,
        template: RackAuthoringTemplate,
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
        let mut new_profile = existing_rack.profile.clone();
        new_profile.fill_direction = match template {
            RackAuthoringTemplate::BenchRows => RackFillDirection::RowMajor,
            RackAuthoringTemplate::PlateColumns | RackAuthoringTemplate::PlateEdgeAvoidance => {
                RackFillDirection::ColumnMajor
            }
        };
        new_profile.blocked_coordinates =
            Self::rack_template_blocked_coordinates(&new_profile, template)?;
        self.reproject_rack_with_profile(rack_id, new_profile)
    }

    pub(super) fn set_rack_fill_direction(
        &mut self,
        rack_id: &str,
        fill_direction: RackFillDirection,
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
        let mut new_profile = existing_rack.profile.clone();
        new_profile.fill_direction = fill_direction;
        new_profile.blocked_coordinates = Self::normalized_blocked_coordinates_for_profile(
            &new_profile,
            &new_profile.blocked_coordinates,
            false,
        )?;
        self.reproject_rack_with_profile(rack_id, new_profile)
    }

    pub(super) fn set_rack_blocked_coordinates(
        &mut self,
        rack_id: &str,
        blocked_coordinates: Vec<String>,
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
        let mut new_profile = existing_rack.profile.clone();
        new_profile.blocked_coordinates = Self::normalized_blocked_coordinates_for_profile(
            &new_profile,
            &blocked_coordinates,
            false,
        )?;
        self.reproject_rack_with_profile(rack_id, new_profile)
    }

    fn xml_escape(raw: &str) -> String {
        raw.replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;")
            .replace('"', "&quot;")
            .replace('\'', "&apos;")
    }

    fn rack_arrangement_display_name(&self, arrangement_id: &str) -> String {
        self.state
            .container_state
            .arrangements
            .get(arrangement_id)
            .and_then(|arrangement| arrangement.name.clone())
            .filter(|name| !name.trim().is_empty())
            .unwrap_or_else(|| arrangement_id.to_string())
    }

    fn scoped_rack_entries(
        &self,
        rack: &Rack,
        arrangement_id: Option<&str>,
    ) -> Result<Vec<RackPlacementEntry>, EngineError> {
        let arrangement_filter = arrangement_id
            .map(str::trim)
            .filter(|value| !value.is_empty());
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
        Ok(rows)
    }

    fn rack_role_summary(role_labels: &[String]) -> String {
        const MAX_ROLES: usize = 4;
        if role_labels.is_empty() {
            return "-".to_string();
        }
        let mut labels = role_labels
            .iter()
            .map(|label| label.trim())
            .filter(|label| !label.is_empty())
            .map(|label| label.to_string())
            .collect::<Vec<_>>();
        labels.dedup();
        if labels.len() <= MAX_ROLES {
            labels.join(", ")
        } else {
            format!(
                "{}, +{} more",
                labels[..MAX_ROLES].join(", "),
                labels.len() - MAX_ROLES
            )
        }
    }

    pub(super) fn export_rack_labels_svg(
        &self,
        rack_id: &str,
        arrangement_id: Option<&str>,
        preset: RackLabelSheetPreset,
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
        let arrangement_filter = arrangement_id
            .map(str::trim)
            .filter(|value| !value.is_empty());
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
        struct RackLabelSheetLayout {
            label_width: f32,
            label_height: f32,
            columns: usize,
            margin: f32,
            title_font_size: usize,
            body_font_size: usize,
            line_step: f32,
            corner_radius: f32,
        }
        let layout = match preset {
            RackLabelSheetPreset::CompactCards => RackLabelSheetLayout {
                label_width: 280.0,
                label_height: 120.0,
                columns: 2,
                margin: 18.0,
                title_font_size: 15,
                body_font_size: 13,
                line_step: 15.0,
                corner_radius: 8.0,
            },
            RackLabelSheetPreset::PrintA4 => RackLabelSheetLayout {
                label_width: 178.0,
                label_height: 104.0,
                columns: 3,
                margin: 14.0,
                title_font_size: 12,
                body_font_size: 10,
                line_step: 12.0,
                corner_radius: 6.0,
            },
            RackLabelSheetPreset::WideCards => RackLabelSheetLayout {
                label_width: 520.0,
                label_height: 88.0,
                columns: 1,
                margin: 18.0,
                title_font_size: 14,
                body_font_size: 11,
                line_step: 12.0,
                corner_radius: 8.0,
            },
        };
        let label_width = layout.label_width;
        let label_height = layout.label_height;
        let columns = layout.columns;
        let margin = layout.margin;
        let svg_width = margin * 2.0 + columns as f32 * label_width;
        let svg_height =
            margin * 2.0 + ((rows.len() + columns - 1) / columns) as f32 * label_height;
        let mut svg = format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{svg_width}\" height=\"{svg_height}\" viewBox=\"0 0 {svg_width} {svg_height}\" data-label-preset=\"{}\">",
            preset.as_str()
        );
        svg.push_str("<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>");
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
                    if let Some(container) = self.state.container_state.containers.get(container_id)
                    {
                        if let Some(name) = container
                            .name
                            .as_ref()
                            .filter(|name| !name.trim().is_empty())
                        {
                            lines.push(name.clone());
                        }
                        if let Some(seq_id) = container.members.first() {
                            lines.push(format!("seq: {seq_id}"));
                            if let Some(dna) = self.state.sequences.get(seq_id) {
                                let topology = if dna.is_circular() {
                                    "circular"
                                } else {
                                    "linear"
                                };
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
            if let Some(arrangement) = self
                .state
                .container_state
                .arrangements
                .get(&entry.arrangement_id)
            {
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
                "<rect x=\"{x}\" y=\"{y}\" width=\"{}\" height=\"{}\" rx=\"{}\" ry=\"{}\" fill=\"#f7fbfd\" stroke=\"#5d8aa8\" stroke-width=\"1.2\"/>",
                label_width - 12.0,
                label_height - 12.0,
                layout.corner_radius,
                layout.corner_radius
            ));
            for (line_idx, line) in lines.iter().enumerate() {
                let font_size = if line_idx == 0 {
                    layout.title_font_size
                } else {
                    layout.body_font_size
                };
                let font_weight = if line_idx == 0 { "700" } else { "400" };
                let line_y = y + 24.0 + line_idx as f32 * layout.line_step;
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

    pub(super) fn export_rack_carrier_labels_svg(
        &self,
        rack_id: &str,
        arrangement_id: Option<&str>,
        template: RackPhysicalTemplateKind,
        preset: RackCarrierLabelPreset,
        path: &str,
    ) -> Result<(usize, RackPhysicalTemplateSpec), EngineError> {
        let rack = self
            .state
            .container_state
            .racks
            .get(rack_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Rack '{rack_id}' not found"),
            })?;
        let spec = Self::rack_physical_template_spec(template, &rack.profile);
        let rows = self.scoped_rack_entries(rack, arrangement_id)?;
        if rows.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Rack '{}' has no placements matching the requested arrangement scope",
                    rack_id
                ),
            });
        }

        let mut arrangement_ids = Vec::new();
        let mut seen_arrangements = HashSet::new();
        for entry in &rows {
            if seen_arrangements.insert(entry.arrangement_id.clone()) {
                arrangement_ids.push(entry.arrangement_id.clone());
            }
        }
        let mut arrangement_summaries = Vec::new();
        for arrangement_id in &arrangement_ids {
            let arrangement_rows = rows
                .iter()
                .filter(|entry| entry.arrangement_id == *arrangement_id)
                .collect::<Vec<_>>();
            let start_coordinate = arrangement_rows
                .first()
                .map(|entry| entry.coordinate.clone())
                .unwrap_or_else(|| "-".to_string());
            let end_coordinate = arrangement_rows
                .last()
                .map(|entry| entry.coordinate.clone())
                .unwrap_or_else(|| "-".to_string());
            let role_labels = arrangement_rows
                .iter()
                .map(|entry| entry.role_label.clone())
                .collect::<Vec<_>>();
            let arrangement = self.state.container_state.arrangements.get(arrangement_id);
            arrangement_summaries.push((
                arrangement_id.clone(),
                self.rack_arrangement_display_name(arrangement_id),
                start_coordinate,
                end_coordinate,
                Self::rack_role_summary(&role_labels),
                arrangement.and_then(|arrangement| arrangement.created_by_op.clone()),
                arrangement_rows.len(),
            ));
        }

        let margin_mm = 4.0;
        let front_strip_width_mm = spec.overall_width_mm.max(48.0);
        let front_strip_height_mm = spec.front_label_strip_depth_mm.max(6.0);
        let module_card_width_mm = front_strip_width_mm.max(72.0);
        let module_card_height_mm = 22.0;
        let module_gap_mm = 4.0;
        let page_width_mm = module_card_width_mm + margin_mm * 2.0;
        let include_front_strip = !matches!(preset, RackCarrierLabelPreset::ModuleCardsOnly);
        let include_module_cards = !matches!(preset, RackCarrierLabelPreset::FrontStripOnly);
        let strip_section_height_mm = if include_front_strip {
            front_strip_height_mm + 6.0
        } else {
            0.0
        };
        let cards_section_height_mm = if include_module_cards {
            arrangement_summaries.len() as f32 * (module_card_height_mm + module_gap_mm)
        } else {
            0.0
        };
        let page_height_mm = margin_mm * 2.0 + strip_section_height_mm + cards_section_height_mm;
        let scope_label = if arrangement_summaries.len() == 1 {
            arrangement_summaries[0].1.clone()
        } else {
            format!("{} modules", arrangement_summaries.len())
        };
        let rack_title = if rack.name.trim().is_empty() {
            rack.rack_id.clone()
        } else {
            format!("{} ({})", rack.name.trim(), rack.rack_id.trim())
        };

        let mut svg = format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.1}mm\" height=\"{:.1}mm\" viewBox=\"0 0 {:.1} {:.1}\" data-rack-carrier-template=\"{}\" data-rack-carrier-preset=\"{}\" data-rack-id=\"{}\">",
            page_width_mm,
            page_height_mm,
            page_width_mm,
            page_height_mm,
            spec.kind.as_str(),
            preset.as_str(),
            Self::xml_escape(&rack.rack_id)
        );
        svg.push_str("<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>");
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"3.2\" font-weight=\"700\" fill=\"#173042\">GENtle carrier labels</text>",
            margin_mm,
            margin_mm - 0.6 + 3.2
        ));

        let strip_x = margin_mm;
        let strip_y = margin_mm + 4.8;
        if include_front_strip {
            svg.push_str(&format!(
                "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"1.5\" ry=\"1.5\" fill=\"#f8fafc\" stroke=\"#5d8aa8\" stroke-width=\"0.35\"/>",
                strip_x,
                strip_y,
                front_strip_width_mm,
                front_strip_height_mm
            ));
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"2.5\" font-weight=\"700\" fill=\"#173042\">{} | {}</text>",
                strip_x + 1.8,
                strip_y + front_strip_height_mm * 0.62,
                Self::xml_escape(&rack_title),
                Self::xml_escape(&scope_label)
            ));
            let total_scope_slots = rows.len().max(1) as f32;
            let mut chip_x = strip_x + 0.7;
            let chip_y = strip_y + front_strip_height_mm - 1.6;
            let chip_height = 1.0;
            for (arrangement_id, _, _, _, _, _, count) in &arrangement_summaries {
                let chip_width =
                    ((front_strip_width_mm - 1.4) * (*count as f32 / total_scope_slots)).max(2.0);
                svg.push_str(&format!(
                    "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"0.5\" ry=\"0.5\" fill=\"{}\" fill-opacity=\"0.88\"/>",
                    chip_x,
                    chip_y,
                    chip_width,
                    chip_height,
                    Self::rack_arrangement_color(arrangement_id)
                ));
                chip_x += chip_width;
            }
        }

        let mut card_y = strip_y + strip_section_height_mm;
        if include_module_cards {
            for (arrangement_id, label, start, end, role_summary, created_by_op, count) in
                &arrangement_summaries
            {
                svg.push_str(&format!(
                    "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"2.0\" ry=\"2.0\" fill=\"#f7fbfd\" stroke=\"{}\" stroke-width=\"0.6\"/>",
                    margin_mm,
                    card_y,
                    module_card_width_mm,
                    module_card_height_mm,
                    Self::rack_arrangement_color(arrangement_id)
                ));
                svg.push_str(&format!(
                    "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"3.4\" height=\"{:.1}\" rx=\"1.1\" ry=\"1.1\" fill=\"{}\" fill-opacity=\"0.92\"/>",
                    margin_mm + 1.2,
                    card_y + 1.2,
                    module_card_height_mm - 2.4,
                    Self::rack_arrangement_color(arrangement_id)
                ));
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"3.0\" font-weight=\"700\" fill=\"#173042\">{}</text>",
                    margin_mm + 6.0,
                    card_y + 5.2,
                    Self::xml_escape(label)
                ));
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"2.35\" fill=\"#334155\">rack={} | span={}..{} | slots={}</text>",
                    margin_mm + 6.0,
                    card_y + 9.6,
                    Self::xml_escape(&rack.rack_id),
                    Self::xml_escape(start),
                    Self::xml_escape(end),
                    count
                ));
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"2.25\" fill=\"#334155\">roles: {}</text>",
                    margin_mm + 6.0,
                    card_y + 13.8,
                    Self::xml_escape(role_summary)
                ));
                let origin_text = created_by_op
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("-");
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"2.15\" fill=\"#475569\">template={} | preset={} | origin={}</text>",
                    margin_mm + 6.0,
                    card_y + 18.0,
                    spec.kind.as_str(),
                    preset.as_str(),
                    Self::xml_escape(origin_text)
                ));
                card_y += module_card_height_mm + module_gap_mm;
            }
        }
        svg.push_str("</svg>");
        fs::write(path, svg).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write rack carrier labels SVG '{}': {e}", path),
        })?;
        Ok((arrangement_summaries.len() + 1, spec))
    }

    fn rack_physical_template_spec(
        template: RackPhysicalTemplateKind,
        profile: &RackProfileSnapshot,
    ) -> RackPhysicalTemplateSpec {
        let (
            family,
            pitch_x_mm,
            pitch_y_mm,
            opening_diameter_mm,
            inner_wall_mm,
            outer_wall_mm,
            floor_thickness_mm,
            rack_height_mm,
            edge_margin_mm,
            corner_radius_mm,
            front_top_clearance_mm,
            front_label_strip_depth_mm,
            front_label_strip_recess_mm,
        ) = match template {
            RackPhysicalTemplateKind::StoragePcrTubeRack => (
                RackPhysicalTemplateFamily::Storage,
                10.0,
                10.0,
                9.2,
                0.8,
                1.2,
                1.2,
                12.0,
                1.8,
                1.4,
                3.0,
                6.0,
                0.8,
            ),
            RackPhysicalTemplateKind::PipettingPcrTubeRack => (
                RackPhysicalTemplateFamily::Pipetting,
                14.0,
                14.0,
                9.2,
                4.8,
                2.4,
                2.0,
                18.0,
                4.0,
                2.4,
                5.0,
                8.0,
                1.2,
            ),
        };
        let overall_width_mm = if profile.columns == 0 {
            0.0
        } else {
            edge_margin_mm * 2.0
                + opening_diameter_mm
                + profile.columns.saturating_sub(1) as f32 * pitch_x_mm
        };
        let overall_depth_mm = if profile.rows == 0 {
            front_top_clearance_mm
        } else {
            front_top_clearance_mm
                + edge_margin_mm
                + opening_diameter_mm
                + profile.rows.saturating_sub(1) as f32 * pitch_y_mm
        };
        RackPhysicalTemplateSpec {
            kind: template,
            family,
            container_format: "pcr_tube_0_2ml".to_string(),
            rows: profile.rows,
            columns: profile.columns,
            pitch_x_mm,
            pitch_y_mm,
            opening_diameter_mm,
            inner_wall_mm,
            outer_wall_mm,
            floor_thickness_mm,
            rack_height_mm,
            edge_margin_mm,
            corner_radius_mm,
            front_top_clearance_mm,
            front_label_strip_depth_mm,
            front_label_strip_recess_mm,
            overall_width_mm,
            overall_depth_mm,
        }
    }

    fn rack_arrangement_color(arrangement_id: &str) -> &'static str {
        const PALETTE: [&str; 8] = [
            "#0f766e", "#2563eb", "#9333ea", "#ca8a04", "#c2410c", "#be123c", "#0369a1", "#4d7c0f",
        ];
        let mut hasher = DefaultHasher::new();
        arrangement_id.hash(&mut hasher);
        let idx = (hasher.finish() as usize) % PALETTE.len();
        PALETTE[idx]
    }

    fn rack_physical_hole_center_mm(
        spec: &RackPhysicalTemplateSpec,
        row: usize,
        column: usize,
    ) -> (f32, f32) {
        (
            spec.edge_margin_mm + spec.opening_diameter_mm * 0.5 + column as f32 * spec.pitch_x_mm,
            spec.front_top_clearance_mm
                + spec.opening_diameter_mm * 0.5
                + row as f32 * spec.pitch_y_mm,
        )
    }

    fn rack_compact_display_label(label: &str, max_chars: usize) -> String {
        let trimmed = label.trim();
        if trimmed.is_empty() {
            return "-".to_string();
        }
        let char_count = trimmed.chars().count();
        if char_count <= max_chars || max_chars <= 1 {
            return trimmed.to_string();
        }
        let keep = max_chars.saturating_sub(1);
        format!("{}…", trimmed.chars().take(keep).collect::<String>())
    }

    fn export_rack_physical_svg_string(
        &self,
        rack: &Rack,
        spec: &RackPhysicalTemplateSpec,
    ) -> Result<String, EngineError> {
        let placements = self.sorted_rack_placements(rack)?;
        let mut placements_by_coordinate = HashMap::new();
        for (_, entry) in &placements {
            placements_by_coordinate.insert(entry.coordinate.clone(), entry.clone());
        }
        let blocked = rack
            .profile
            .blocked_coordinates
            .iter()
            .cloned()
            .collect::<HashSet<_>>();
        let mut arrangement_legend = Vec::new();
        let mut seen_arrangements = HashSet::new();
        for (_, entry) in &placements {
            if seen_arrangements.insert(entry.arrangement_id.clone()) {
                let label = self
                    .state
                    .container_state
                    .arrangements
                    .get(&entry.arrangement_id)
                    .and_then(|arrangement| arrangement.name.clone())
                    .filter(|name| !name.trim().is_empty())
                    .unwrap_or_else(|| entry.arrangement_id.clone());
                arrangement_legend.push((
                    entry.arrangement_id.clone(),
                    label,
                    Self::rack_arrangement_color(&entry.arrangement_id).to_string(),
                ));
            }
        }
        let legend_height = if arrangement_legend.is_empty() {
            0.0
        } else {
            5.0 + arrangement_legend.len() as f32 * 4.5
        };
        let header_height = 16.0 + legend_height;
        let footer_height = 8.0;
        let rack_x = 2.0;
        let rack_y = header_height;
        let svg_width = spec.overall_width_mm + 4.0;
        let svg_height = header_height + spec.overall_depth_mm + footer_height;
        let mut svg = format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.1}mm\" height=\"{:.1}mm\" viewBox=\"0 0 {:.1} {:.1}\" data-rack-physical-template=\"{}\" data-rack-id=\"{}\">",
            svg_width,
            svg_height,
            svg_width,
            svg_height,
            spec.kind.as_str(),
            Self::xml_escape(&rack.rack_id)
        );
        svg.push_str("<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>");
        svg.push_str(&format!(
            "<text x=\"2\" y=\"6\" font-family=\"monospace\" font-size=\"4.2\" font-weight=\"700\" fill=\"#0f172a\">GENtle rack fabrication sketch</text><text x=\"2\" y=\"10.5\" font-family=\"monospace\" font-size=\"3.2\" fill=\"#334155\">rack={} | template={} | family={} | format={} | profile={}x{}</text>",
            Self::xml_escape(&rack.rack_id),
            spec.kind.as_str(),
            spec.family.as_str(),
            Self::xml_escape(&spec.container_format),
            spec.rows,
            spec.columns
        ));
        let mut legend_y = 15.0;
        if !arrangement_legend.is_empty() {
            for (_, label, color) in &arrangement_legend {
                svg.push_str(&format!(
                    "<rect x=\"2\" y=\"{:.1}\" width=\"3.2\" height=\"3.2\" rx=\"0.8\" ry=\"0.8\" fill=\"{}\"/><text x=\"6.4\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"2.8\" fill=\"#334155\">{}</text>",
                    legend_y - 2.6,
                    color,
                    legend_y,
                    Self::xml_escape(label)
                ));
                legend_y += 4.5;
            }
            if !blocked.is_empty() {
                svg.push_str(&format!(
                    "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"3.2\" height=\"3.2\" rx=\"0.8\" ry=\"0.8\" fill=\"#cbd5e1\" stroke=\"#64748b\" stroke-width=\"0.35\"/><text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"2.8\" fill=\"#334155\">blocked / reserved</text>",
                    svg_width - 56.0,
                    12.4,
                    svg_width - 51.0,
                    15.0,
                ));
            }
        }
        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"{:.1}\" ry=\"{:.1}\" fill=\"#f8fafc\" stroke=\"#0f172a\" stroke-width=\"0.6\"/>",
            rack_x,
            rack_y,
            spec.overall_width_mm,
            spec.overall_depth_mm,
            spec.corner_radius_mm,
            spec.corner_radius_mm
        ));
        if spec.front_top_clearance_mm > 0.8 {
            svg.push_str(&format!(
                "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#e2e8f0\" stroke=\"#94a3b8\" stroke-width=\"0.3\"/>",
                rack_x + 0.6,
                rack_y + 0.6,
                spec.overall_width_mm - 1.2,
                (spec.front_top_clearance_mm - 1.2).max(0.6),
            ));
        }
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"2.6\" fill=\"#334155\">front lip {:.1} mm | front label face {:.1} mm</text>",
            rack_x + 1.0,
            rack_y + 3.4,
            spec.front_top_clearance_mm,
            spec.front_label_strip_depth_mm
        ));
        for column in 0..spec.columns {
            let (cx, _) = Self::rack_physical_hole_center_mm(spec, 0, column);
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"2.6\" fill=\"#475569\">{}</text>",
                rack_x + cx,
                rack_y + spec.front_top_clearance_mm - 0.8,
                column + 1
            ));
        }
        for row in 0..spec.rows {
            let (_, cy) = Self::rack_physical_hole_center_mm(spec, row, 0);
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" dominant-baseline=\"middle\" font-family=\"monospace\" font-size=\"2.6\" fill=\"#475569\">{}</text>",
                rack_x + spec.edge_margin_mm - 0.6,
                rack_y + cy,
                Self::rack_row_label_from_index(row)
            ));
            for column in 0..spec.columns {
                let coordinate = Self::rack_coordinate_from_row_column(&rack.profile, row, column)?;
                let (cx, cy) = Self::rack_physical_hole_center_mm(spec, row, column);
                let x = rack_x + cx;
                let y = rack_y + cy;
                if blocked.contains(&coordinate) {
                    svg.push_str(&format!(
                        "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"#e2e8f0\" stroke=\"#64748b\" stroke-width=\"0.35\"/><line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#64748b\" stroke-width=\"0.35\"/><line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#64748b\" stroke-width=\"0.35\"/>",
                        x,
                        y,
                        spec.opening_diameter_mm * 0.5,
                        x - 2.2,
                        y - 2.2,
                        x + 2.2,
                        y + 2.2,
                        x - 2.2,
                        y + 2.2,
                        x + 2.2,
                        y - 2.2
                    ));
                    continue;
                }
                svg.push_str(&format!(
                    "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"#ffffff\" stroke=\"#334155\" stroke-width=\"0.35\"/>",
                    x,
                    y,
                    spec.opening_diameter_mm * 0.5
                ));
                if let Some(entry) = placements_by_coordinate.get(&coordinate) {
                    match entry.occupant.as_ref() {
                        Some(RackOccupant::LadderReference { .. }) => {
                            svg.push_str(&format!(
                                "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"#d97706\" fill-opacity=\"0.85\"/>",
                                x,
                                y,
                                spec.opening_diameter_mm * 0.18
                            ));
                        }
                        Some(RackOccupant::Container { .. }) => {
                            let fill = Self::rack_arrangement_color(&entry.arrangement_id);
                            svg.push_str(&format!(
                                "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"{}\" fill-opacity=\"0.28\" stroke=\"{}\" stroke-width=\"0.25\"/>",
                                x,
                                y,
                                spec.opening_diameter_mm * 0.33,
                                fill,
                                fill
                            ));
                        }
                        None => {}
                    }
                }
            }
        }
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"2.6\" fill=\"#64748b\">GENtle fabrication SVG</text>",
            svg_width - 1.4,
            svg_height - 1.6
        ));
        svg.push_str("</svg>");
        Ok(svg)
    }

    fn export_rack_isometric_svg_string(
        &self,
        rack: &Rack,
        spec: &RackPhysicalTemplateSpec,
    ) -> Result<String, EngineError> {
        let placements = self.sorted_rack_placements(rack)?;
        let mut placements_by_coordinate = HashMap::new();
        for (_, entry) in &placements {
            placements_by_coordinate.insert(entry.coordinate.clone(), entry.clone());
        }
        let blocked = rack
            .profile
            .blocked_coordinates
            .iter()
            .cloned()
            .collect::<HashSet<_>>();
        let mut arrangement_legend = Vec::new();
        let mut seen_arrangements = HashSet::new();
        for (_, entry) in &placements {
            if seen_arrangements.insert(entry.arrangement_id.clone()) {
                let label = self
                    .state
                    .container_state
                    .arrangements
                    .get(&entry.arrangement_id)
                    .and_then(|arrangement| arrangement.name.clone())
                    .filter(|name| !name.trim().is_empty())
                    .unwrap_or_else(|| entry.arrangement_id.clone());
                arrangement_legend.push((
                    entry.arrangement_id.clone(),
                    label,
                    Self::rack_arrangement_color(&entry.arrangement_id).to_string(),
                ));
            }
        }

        let (top_fill, front_fill, side_fill, strip_fill, strip_text_fill, edge_stroke) =
            match spec.family {
                RackPhysicalTemplateFamily::Storage => (
                    "#f8fafc", "#e2e8f0", "#cbd5e1", "#dbeafe", "#173042", "#475569",
                ),
                RackPhysicalTemplateFamily::Pipetting => (
                    "#e7f6f2", "#b8dfd5", "#8ecaba", "#0f766e", "#f8fafc", "#115e59",
                ),
            };
        let skew_x = 0.48;
        let skew_y = 0.32;
        let height_scale = 0.78;
        let margin = 8.0;
        let top_left_x = margin;
        let top_left_y = margin + spec.overall_depth_mm * skew_y;
        let top_front_left = (top_left_x, top_left_y);
        let top_front_right = (top_left_x + spec.overall_width_mm, top_left_y);
        let top_back_left = (
            top_left_x + spec.overall_depth_mm * skew_x,
            top_left_y - spec.overall_depth_mm * skew_y,
        );
        let top_back_right = (top_back_left.0 + spec.overall_width_mm, top_back_left.1);
        let front_bottom_left = (
            top_front_left.0,
            top_front_left.1 + spec.rack_height_mm * height_scale,
        );
        let front_bottom_right = (
            top_front_right.0,
            top_front_right.1 + spec.rack_height_mm * height_scale,
        );
        let back_bottom_right = (
            top_back_right.0,
            top_back_right.1 + spec.rack_height_mm * height_scale,
        );
        let legend_width = if arrangement_legend.is_empty() {
            0.0
        } else {
            70.0
        };
        let svg_width = top_back_right.0 + margin + legend_width;
        let svg_height = front_bottom_left.1 + margin + 16.0;
        let occupied_count = placements.len();
        let total_slots =
            Self::rack_available_coordinates_in_fill_order(&rack.profile)?.len() + blocked.len();
        let rack_title = rack.rack_id.trim().to_string();
        let scope_label = if arrangement_legend.len() == 1 {
            "1 arrangement".to_string()
        } else {
            format!("{} arrangements", arrangement_legend.len())
        };

        let mut svg = format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.1}mm\" height=\"{:.1}mm\" viewBox=\"0 0 {:.1} {:.1}\" data-rack-isometric-template=\"{}\" data-rack-id=\"{}\" data-rack-front-top-clearance-mm=\"{:.1}\">",
            svg_width,
            svg_height,
            svg_width,
            svg_height,
            spec.kind.as_str(),
            Self::xml_escape(&rack.rack_id),
            spec.front_top_clearance_mm
        );
        svg.push_str("<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>");
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"4.0\" font-weight=\"700\" fill=\"#173042\">GENtle rack isometric sketch</text>",
            margin,
            4.8
        ));

        svg.push_str(&format!(
            "<polygon points=\"{:.2},{:.2} {:.2},{:.2} {:.2},{:.2} {:.2},{:.2}\" fill=\"{}\" stroke=\"{}\" stroke-width=\"0.6\"/>",
            top_front_left.0,
            top_front_left.1,
            top_front_right.0,
            top_front_right.1,
            top_back_right.0,
            top_back_right.1,
            top_back_left.0,
            top_back_left.1,
            top_fill,
            edge_stroke
        ));
        svg.push_str(&format!(
            "<polygon points=\"{:.2},{:.2} {:.2},{:.2} {:.2},{:.2} {:.2},{:.2}\" fill=\"{}\" stroke=\"{}\" stroke-width=\"0.6\"/>",
            top_front_left.0,
            top_front_left.1,
            top_front_right.0,
            top_front_right.1,
            front_bottom_right.0,
            front_bottom_right.1,
            front_bottom_left.0,
            front_bottom_left.1,
            front_fill,
            edge_stroke
        ));
        svg.push_str(&format!(
            "<polygon points=\"{:.2},{:.2} {:.2},{:.2} {:.2},{:.2} {:.2},{:.2}\" fill=\"{}\" stroke=\"{}\" stroke-width=\"0.6\"/>",
            top_front_right.0,
            top_front_right.1,
            top_back_right.0,
            top_back_right.1,
            back_bottom_right.0,
            back_bottom_right.1,
            front_bottom_right.0,
            front_bottom_right.1,
            side_fill,
            edge_stroke
        ));

        let strip_margin = 1.4;
        let strip_x = top_front_left.0 + strip_margin;
        let strip_y = top_front_left.1 + 1.2;
        let strip_w = (spec.overall_width_mm - strip_margin * 2.0).max(16.0);
        let strip_h = (spec.front_label_strip_depth_mm * 0.82).max(7.2);
        svg.push_str(&format!(
            "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" rx=\"1.2\" ry=\"1.2\" fill=\"{}\" fill-opacity=\"0.94\" stroke=\"{}\" stroke-width=\"0.35\"/>",
            strip_x,
            strip_y,
            strip_w,
            strip_h,
            strip_fill,
            edge_stroke
        ));
        svg.push_str(&format!(
            "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"monospace\" font-size=\"2.75\" font-weight=\"700\" fill=\"{}\">{}</text>",
            strip_x + 1.4,
            strip_y + 3.0,
            strip_text_fill,
            Self::xml_escape(&rack_title)
        ));
        svg.push_str(&format!(
            "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"monospace\" font-size=\"2.15\" fill=\"{}\">{} | {} / {} occupied</text>",
            strip_x + 1.4,
            strip_y + 5.9,
            strip_text_fill,
            Self::xml_escape(&scope_label),
            occupied_count,
            total_slots
        ));

        for column in 0..spec.columns {
            let (cx, _) = Self::rack_physical_hole_center_mm(spec, 0, column);
            let x = top_left_x + cx;
            let y = top_left_y + 1.4;
            svg.push_str(&format!(
                "<text x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"2.4\" fill=\"#475569\">{}</text>",
                x,
                y,
                column + 1
            ));
        }
        for row in 0..spec.rows {
            let (_, cy) = Self::rack_physical_hole_center_mm(spec, row, 0);
            let x = top_left_x + cy * skew_x - 2.0;
            let y = top_left_y - cy * skew_y + 0.8;
            svg.push_str(&format!(
                "<text x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"2.4\" fill=\"#475569\">{}</text>",
                x,
                y,
                Self::rack_row_label_from_index(row)
            ));
        }

        let hole_rx = spec.opening_diameter_mm * 0.42;
        let hole_ry = spec.opening_diameter_mm * 0.24;
        let cap_height = match spec.family {
            RackPhysicalTemplateFamily::Storage => 3.0,
            RackPhysicalTemplateFamily::Pipetting => 4.0,
        };
        for row in 0..spec.rows {
            for column in 0..spec.columns {
                let coordinate = Self::rack_coordinate_from_row_column(&rack.profile, row, column)?;
                let (cx_mm, cy_mm) = Self::rack_physical_hole_center_mm(spec, row, column);
                let x = top_left_x + cx_mm + cy_mm * skew_x;
                let y = top_left_y - cy_mm * skew_y;
                if blocked.contains(&coordinate) {
                    svg.push_str(&format!(
                        "<ellipse cx=\"{:.2}\" cy=\"{:.2}\" rx=\"{:.2}\" ry=\"{:.2}\" fill=\"#e2e8f0\" stroke=\"#64748b\" stroke-width=\"0.35\"/><line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#64748b\" stroke-width=\"0.35\"/><line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#64748b\" stroke-width=\"0.35\"/>",
                        x,
                        y,
                        hole_rx,
                        hole_ry,
                        x - hole_rx * 0.75,
                        y - hole_ry * 0.75,
                        x + hole_rx * 0.75,
                        y + hole_ry * 0.75,
                        x - hole_rx * 0.75,
                        y + hole_ry * 0.75,
                        x + hole_rx * 0.75,
                        y - hole_ry * 0.75
                    ));
                    continue;
                }
                svg.push_str(&format!(
                    "<ellipse cx=\"{:.2}\" cy=\"{:.2}\" rx=\"{:.2}\" ry=\"{:.2}\" fill=\"#ffffff\" stroke=\"#475569\" stroke-width=\"0.35\"/>",
                    x,
                    y,
                    hole_rx,
                    hole_ry
                ));
                if let Some(entry) = placements_by_coordinate.get(&coordinate) {
                    let (fill, top_fill, stroke) = match entry.occupant.as_ref() {
                        Some(RackOccupant::LadderReference { .. }) => {
                            ("#b45309", "#f59e0b", "#92400e")
                        }
                        Some(RackOccupant::Container { .. }) => {
                            let color = Self::rack_arrangement_color(&entry.arrangement_id);
                            (color, color, "#173042")
                        }
                        None => continue,
                    };
                    let cap_rx = hole_rx * 0.84;
                    let cap_ry = hole_ry * 0.92;
                    let cap_top_y = y - cap_height;
                    let shell_left_x = x - cap_rx;
                    let shell_right_x = x + cap_rx;
                    let shell_interface_y = y - cap_ry * 0.06;
                    let shell_curve_y = y + cap_ry * 0.88;
                    svg.push_str(&format!(
                        "<ellipse cx=\"{:.2}\" cy=\"{:.2}\" rx=\"{:.2}\" ry=\"{:.2}\" fill=\"#0f172a\" fill-opacity=\"0.12\"/>",
                        x,
                        y + 0.4,
                        cap_rx,
                        cap_ry
                    ));
                    svg.push_str(&format!(
                        "<path data-rack-tube-shell=\"1\" d=\"M {:.2} {:.2} L {:.2} {:.2} Q {:.2} {:.2} {:.2} {:.2} L {:.2} {:.2} Z\" fill=\"{}\" fill-opacity=\"0.92\"/>",
                        shell_left_x,
                        cap_top_y,
                        shell_left_x,
                        shell_interface_y,
                        x,
                        shell_curve_y,
                        shell_right_x,
                        shell_interface_y,
                        shell_right_x,
                        cap_top_y,
                        fill
                    ));
                    svg.push_str(&format!(
                        "<path data-rack-tube-interface=\"1\" d=\"M {:.2} {:.2} Q {:.2} {:.2} {:.2} {:.2}\" fill=\"none\" stroke=\"{}\" stroke-opacity=\"0.35\" stroke-width=\"0.28\" stroke-linecap=\"round\"/>",
                        shell_left_x,
                        shell_interface_y,
                        x,
                        shell_curve_y,
                        shell_right_x,
                        shell_interface_y,
                        stroke
                    ));
                    svg.push_str(&format!(
                        "<ellipse cx=\"{:.2}\" cy=\"{:.2}\" rx=\"{:.2}\" ry=\"{:.2}\" fill=\"{}\" stroke=\"{}\" stroke-width=\"0.25\"/>",
                        x,
                        cap_top_y,
                        cap_rx,
                        cap_ry,
                        top_fill,
                        stroke
                    ));
                    svg.push_str(&format!(
                        "<path d=\"M {:.2} {:.2} C {:.2} {:.2}, {:.2} {:.2}, {:.2} {:.2}\" fill=\"none\" stroke=\"#ffffff\" stroke-opacity=\"0.35\" stroke-width=\"0.25\"/>",
                        x - cap_rx * 0.55,
                        cap_top_y + 0.05,
                        x - cap_rx * 0.2,
                        cap_top_y - cap_ry * 0.7,
                        x + cap_rx * 0.15,
                        cap_top_y - cap_ry * 0.7,
                        x + cap_rx * 0.4,
                        cap_top_y + 0.05
                    ));
                }
            }
        }

        if !arrangement_legend.is_empty() {
            let legend_x = top_back_right.0 + 8.0;
            let mut legend_y = top_back_left.1 + 8.0;
            svg.push_str(&format!(
                "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"monospace\" font-size=\"3.0\" font-weight=\"700\" fill=\"#173042\">Occupancy</text>",
                legend_x,
                legend_y
            ));
            legend_y += 4.0;
            for (_, label, color) in &arrangement_legend {
                svg.push_str(&format!(
                    "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"4.0\" height=\"4.0\" rx=\"1.0\" ry=\"1.0\" fill=\"{}\"/><text x=\"{:.2}\" y=\"{:.2}\" font-family=\"monospace\" font-size=\"2.5\" fill=\"#334155\">{}</text>",
                    legend_x,
                    legend_y - 3.0,
                    color,
                    legend_x + 6.0,
                    legend_y,
                    Self::xml_escape(&Self::rack_compact_display_label(label, 26))
                ));
                legend_y += 5.0;
            }
            svg.push_str(&format!(
                "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"4.0\" height=\"4.0\" rx=\"1.0\" ry=\"1.0\" fill=\"#f59e0b\"/><text x=\"{:.2}\" y=\"{:.2}\" font-family=\"monospace\" font-size=\"2.5\" fill=\"#334155\">ladder reference</text>",
                legend_x,
                legend_y - 3.0,
                legend_x + 6.0,
                legend_y
            ));
            legend_y += 5.0;
            if !blocked.is_empty() {
                svg.push_str(&format!(
                    "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"4.0\" height=\"4.0\" rx=\"1.0\" ry=\"1.0\" fill=\"#e2e8f0\" stroke=\"#64748b\" stroke-width=\"0.25\"/><text x=\"{:.2}\" y=\"{:.2}\" font-family=\"monospace\" font-size=\"2.5\" fill=\"#334155\">blocked / reserved</text>",
                    legend_x,
                    legend_y - 3.0,
                    legend_x + 6.0,
                    legend_y
                ));
            }
        }

        svg.push_str(&format!(
            "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"monospace\" font-size=\"2.4\" fill=\"#64748b\">template={} | family={} | positions={}</text>",
            margin,
            svg_height - 4.0,
            spec.kind.as_str(),
            spec.family.as_str(),
            occupied_count
        ));
        svg.push_str("</svg>");
        Ok(svg)
    }

    pub(super) fn export_rack_fabrication_svg(
        &self,
        rack_id: &str,
        template: RackPhysicalTemplateKind,
        path: &str,
    ) -> Result<RackPhysicalTemplateSpec, EngineError> {
        let rack = self
            .state
            .container_state
            .racks
            .get(rack_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Rack '{rack_id}' not found"),
            })?;
        let spec = Self::rack_physical_template_spec(template, &rack.profile);
        let svg = self.export_rack_physical_svg_string(rack, &spec)?;
        fs::write(path, svg).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write rack fabrication SVG '{}': {e}", path),
        })?;
        Ok(spec)
    }

    pub(super) fn export_rack_isometric_svg(
        &self,
        rack_id: &str,
        template: RackPhysicalTemplateKind,
        path: &str,
    ) -> Result<RackPhysicalTemplateSpec, EngineError> {
        let rack = self
            .state
            .container_state
            .racks
            .get(rack_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Rack '{rack_id}' not found"),
            })?;
        let spec = Self::rack_physical_template_spec(template, &rack.profile);
        let svg = self.export_rack_isometric_svg_string(rack, &spec)?;
        fs::write(path, svg).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write rack isometric SVG '{}': {e}", path),
        })?;
        Ok(spec)
    }

    pub(super) fn export_rack_openscad(
        &self,
        rack_id: &str,
        template: RackPhysicalTemplateKind,
        path: &str,
    ) -> Result<RackPhysicalTemplateSpec, EngineError> {
        let rack = self
            .state
            .container_state
            .racks
            .get(rack_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Rack '{rack_id}' not found"),
            })?;
        let spec = Self::rack_physical_template_spec(template, &rack.profile);
        let blocked = rack
            .profile
            .blocked_coordinates
            .iter()
            .cloned()
            .collect::<HashSet<_>>();
        let mut scad = format!(
            "// GENtle rack OpenSCAD export\n// rack_id={}\n// template={}\n// family={}\n// format={}\n$fn = 56;\n\nouter_width = {:.3};\nouter_depth = {:.3};\nrack_height = {:.3};\nopening_diameter = {:.3};\nfloor_thickness = {:.3};\nfront_top_clearance = {:.3};\nfront_label_strip_depth = {:.3};\nfront_label_strip_recess = {:.3};\ncorner_radius = {:.3};\n\nmodule gentle_rack() {{\n    difference() {{\n        cube([outer_width, outer_depth, rack_height], false);\n",
            rack.rack_id,
            spec.kind.as_str(),
            spec.family.as_str(),
            spec.container_format,
            spec.overall_width_mm,
            spec.overall_depth_mm,
            spec.rack_height_mm,
            spec.opening_diameter_mm,
            spec.floor_thickness_mm,
            spec.front_top_clearance_mm,
            spec.front_label_strip_depth_mm,
            spec.front_label_strip_recess_mm,
            spec.corner_radius_mm
        );
        for row in 0..spec.rows {
            for column in 0..spec.columns {
                let coordinate = Self::rack_coordinate_from_row_column(&rack.profile, row, column)?;
                if blocked.contains(&coordinate) {
                    continue;
                }
                let (cx, cy) = Self::rack_physical_hole_center_mm(&spec, row, column);
                scad.push_str(&format!(
                    "        translate([{:.3}, {:.3}, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);\n",
                    cx, cy
                ));
            }
        }
        scad.push_str(&format!(
            "        translate([{:.3}, {:.3}, rack_height - front_label_strip_recess]) cube([{:.3}, {:.3}, front_label_strip_recess + 0.05], false);\n",
            spec.outer_wall_mm * 0.5,
            0.8,
            (spec.overall_width_mm - spec.outer_wall_mm).max(0.0),
            (spec.front_top_clearance_mm - 1.6).max(0.0)
        ));
        scad.push_str("    }\n}\n\ngentle_rack();\n");
        fs::write(path, scad).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write rack OpenSCAD '{}': {e}", path),
        })?;
        Ok(spec)
    }

    pub(super) fn export_rack_simulation_json(
        &self,
        rack_id: &str,
        template: RackPhysicalTemplateKind,
        path: &str,
    ) -> Result<RackPhysicalTemplateSpec, EngineError> {
        let rack = self
            .state
            .container_state
            .racks
            .get(rack_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Rack '{rack_id}' not found"),
            })?;
        let spec = Self::rack_physical_template_spec(template, &rack.profile);
        let placements = self.sorted_rack_placements(rack)?;
        let blocked = rack
            .profile
            .blocked_coordinates
            .iter()
            .cloned()
            .collect::<HashSet<_>>();
        let mut placements_by_coordinate = HashMap::new();
        for (_, entry) in &placements {
            placements_by_coordinate.insert(entry.coordinate.clone(), entry.clone());
        }

        let mut arrangement_ids = Vec::new();
        let mut seen_arrangements = HashSet::new();
        for (_, entry) in &placements {
            if seen_arrangements.insert(entry.arrangement_id.clone()) {
                arrangement_ids.push(entry.arrangement_id.clone());
            }
        }
        let arrangement_blocks = arrangement_ids
            .into_iter()
            .filter_map(|arrangement_id| {
                let entries = placements
                    .iter()
                    .filter(|(_, entry)| entry.arrangement_id == arrangement_id)
                    .map(|(_, entry)| entry)
                    .collect::<Vec<_>>();
                let first = entries.first()?;
                let last = entries.last()?;
                let arrangement = self.state.container_state.arrangements.get(&arrangement_id);
                let role_labels = entries
                    .iter()
                    .map(|entry| entry.role_label.clone())
                    .collect::<Vec<_>>();
                Some(serde_json::json!({
                    "arrangement_id": arrangement_id,
                    "arrangement_name": arrangement
                        .and_then(|arrangement| arrangement.name.clone())
                        .filter(|name| !name.trim().is_empty())
                        .unwrap_or_else(|| first.arrangement_id.clone()),
                    "coordinate_start": first.coordinate,
                    "coordinate_end": last.coordinate,
                    "placement_count": entries.len(),
                    "role_labels": role_labels,
                    "created_by_op": arrangement.and_then(|arrangement| arrangement.created_by_op.clone()),
                }))
            })
            .collect::<Vec<_>>();

        let mut slots = Vec::with_capacity(spec.rows.saturating_mul(spec.columns));
        for row in 0..spec.rows {
            for column in 0..spec.columns {
                let coordinate = Self::rack_coordinate_from_row_column(&rack.profile, row, column)?;
                let (center_x_mm, center_y_mm) =
                    Self::rack_physical_hole_center_mm(&spec, row, column);
                let entry = placements_by_coordinate.get(&coordinate);
                let occupant = match entry.and_then(|entry| entry.occupant.as_ref()) {
                    Some(RackOccupant::Container { container_id }) => {
                        let container = self.state.container_state.containers.get(container_id);
                        let sequence_id = container
                            .and_then(|container| container.members.first())
                            .cloned();
                        let sequence = sequence_id
                            .as_ref()
                            .and_then(|seq_id| self.state.sequences.get(seq_id));
                        serde_json::json!({
                            "kind": "container",
                            "container_id": container_id,
                            "container_name": container.and_then(|container| container.name.clone()),
                            "sequence_id": sequence_id,
                            "length_bp": sequence.map(|dna| dna.len()),
                            "topology": sequence.map(|dna| if dna.is_circular() { "circular" } else { "linear" }),
                            "created_by_op": container.and_then(|container| container.created_by_op.clone()),
                        })
                    }
                    Some(RackOccupant::LadderReference { ladder_name }) => {
                        serde_json::json!({
                            "kind": "ladder_reference",
                            "ladder_name": ladder_name,
                        })
                    }
                    None => serde_json::Value::Null,
                };
                slots.push(serde_json::json!({
                    "coordinate": coordinate,
                    "row": row,
                    "column": column,
                    "fill_ordinal": Self::rack_fill_ordinal(&rack.profile, row, column),
                    "blocked": blocked.contains(&coordinate),
                    "center_mm": {
                        "x": center_x_mm,
                        "y": center_y_mm,
                    },
                    "arrangement_id": entry.map(|entry| entry.arrangement_id.clone()),
                    "arrangement_name": entry.map(|entry| self.rack_arrangement_display_name(&entry.arrangement_id)),
                    "order_index": entry.map(|entry| entry.order_index),
                    "role_label": entry.map(|entry| entry.role_label.clone()),
                    "occupant": occupant,
                }));
            }
        }

        let payload = serde_json::json!({
            "schema": "gentle.rack_simulation_export.v1",
            "rack": {
                "rack_id": &rack.rack_id,
                "name": &rack.name,
                "created_by_op": &rack.created_by_op,
                "created_at_unix_ms": rack.created_at_unix_ms,
                "profile": &rack.profile,
            },
            "physical_template": &spec,
            "arrangement_blocks": arrangement_blocks,
            "slots": slots,
        });
        let payload_text = serde_json::to_string_pretty(&payload).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not serialize rack simulation JSON '{}': {e}", path),
        })?;
        fs::write(path, payload_text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write rack simulation JSON '{}': {e}", path),
        })?;
        Ok(spec)
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
