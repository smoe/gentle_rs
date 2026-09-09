//! Thin editing and file-binding helpers for the shared homology request.

use eframe::egui;
use gentle_protocol as gp;

pub(super) fn parse_bound_conservation_request(
    bytes: &[u8],
    current: &gp::GenomicRegionHomologyScreenRequest,
) -> Result<gp::GenomicRegionHomologyScreenRequest, String> {
    let mut request: gp::GenomicRegionHomologyScreenRequest =
        serde_json::from_slice(bytes).map_err(|error| error.to_string())?;
    if request.set_id != current.set_id
        || request.region_id != current.region_id
        || request
            .expected_region_content_sha256
            .as_ref()
            .is_some_and(|digest| Some(digest) != current.expected_region_content_sha256.as_ref())
    {
        return Err("Request belongs to a different or changed saved region".to_string());
    }
    request
        .expected_region_content_sha256
        .clone_from(&current.expected_region_content_sha256);
    Ok(request)
}

fn text_field(ui: &mut egui::Ui, label: &str, value: &mut String) {
    ui.label(label);
    ui.add(egui::TextEdit::singleline(value).desired_width(230.0));
    ui.end_row();
}

fn optional_field(ui: &mut egui::Ui, label: &str, value: &mut Option<String>) {
    let mut text = value.clone().unwrap_or_default();
    ui.label(label);
    if ui
        .add(egui::TextEdit::singleline(&mut text).desired_width(230.0))
        .changed()
    {
        *value = (!text.trim().is_empty()).then_some(text);
    }
    ui.end_row();
}

pub(super) fn render_conservation_request(
    ui: &mut egui::Ui,
    request: &mut gp::GenomicRegionHomologyScreenRequest,
) {
    egui::Grid::new("conservation_sources")
        .num_columns(2)
        .show(ui, |ui| {
            optional_field(
                ui,
                "Query genome ID (optional)",
                &mut request.query_genome_id,
            );
            optional_field(ui, "Genome catalog (optional)", &mut request.catalog_path);
            optional_field(ui, "Genome cache (optional)", &mut request.cache_dir);
        });
    ui.horizontal(|ui| {
        ui.strong("Target genomes");
        if ui
            .small_button("+")
            .on_hover_text("Add an explicit target genome")
            .clicked()
        {
            request.targets.push(Default::default());
        }
    });
    if request.targets.is_empty() {
        ui.label("All validated local genomic indexes");
    }
    let mut remove = None;
    for (index, target) in request.targets.iter_mut().enumerate() {
        ui.push_id(index, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label("Genome ID");
                ui.add(egui::TextEdit::singleline(&mut target.genome_id).desired_width(180.0));
                ui.checkbox(&mut target.required, "Required");
                egui::ComboBox::from_id_salt("role")
                    .selected_text(target.role.as_str())
                    .show_ui(ui, |ui| {
                        for role in [
                            gp::GenomicRegionHomologyTargetRole::SameGenome,
                            gp::GenomicRegionHomologyTargetRole::ExpectedOrtholog,
                            gp::GenomicRegionHomologyTargetRole::CrossSpeciesUnassigned,
                        ] {
                            ui.selectable_value(&mut target.role, role, role.as_str());
                        }
                    });
                if ui
                    .small_button("x")
                    .on_hover_text("Remove this target")
                    .clicked()
                {
                    remove = Some(index);
                }
            });
            if target.role == gp::GenomicRegionHomologyTargetRole::ExpectedOrtholog
                || !target.expected_loci.is_empty()
            {
                egui::CollapsingHeader::new("Expected loci and orthology evidence").show(
                    ui,
                    |ui| {
                        let mut remove_locus = None;
                        for (locus_index, locus) in target.expected_loci.iter_mut().enumerate() {
                            egui::Grid::new(("expected_locus", locus_index))
                                .num_columns(2)
                                .show(ui, |ui| {
                                    text_field(ui, "Locus ID", &mut locus.expected_locus_id);
                                    text_field(ui, "Assembly", &mut locus.reference.assembly_name);
                                    text_field(ui, "Contig", &mut locus.reference.contig_name);
                                    ui.label("Start (0-based)");
                                    ui.add(egui::DragValue::new(&mut locus.start_0based));
                                    ui.end_row();
                                    ui.label("End (exclusive)");
                                    ui.add(egui::DragValue::new(&mut locus.end_0based_exclusive));
                                    ui.end_row();
                                    ui.label("Strand");
                                    egui::ComboBox::from_id_salt(("strand", locus_index))
                                        .selected_text(format!("{:?}", locus.strand))
                                        .show_ui(ui, |ui| {
                                            for strand in [
                                                gp::GenomicRegionStrand::Plus,
                                                gp::GenomicRegionStrand::Minus,
                                                gp::GenomicRegionStrand::Unstranded,
                                            ] {
                                                ui.selectable_value(
                                                    &mut locus.strand,
                                                    strand,
                                                    format!("{strand:?}"),
                                                );
                                            }
                                        });
                                    ui.end_row();
                                    text_field(ui, "Orthology evidence ID", &mut locus.evidence_id);
                                    text_field(ui, "Evidence source ID", &mut locus.source_id);
                                    optional_field(
                                        ui,
                                        "Source SHA-256 (optional)",
                                        &mut locus.source_sha256,
                                    );
                                });
                            if ui
                                .small_button("x")
                                .on_hover_text("Remove this expected locus")
                                .clicked()
                            {
                                remove_locus = Some(locus_index);
                            }
                            ui.separator();
                        }
                        if let Some(index) = remove_locus {
                            target.expected_loci.remove(index);
                        }
                        if ui
                            .small_button("+")
                            .on_hover_text("Add an expected ortholog locus")
                            .clicked()
                        {
                            target.expected_loci.push(Default::default());
                        }
                    },
                );
            }
            ui.separator();
        });
    }
    if let Some(index) = remove {
        request.targets.remove(index);
    }
    ui.strong("Search policy");
    let policy = &mut request.policy;
    egui::Grid::new("conservation_policy")
        .num_columns(2)
        .show(ui, |ui| {
            ui.label("Minimum identity (%)");
            ui.add(egui::DragValue::new(&mut policy.min_identity_percent).range(0.0..=100.0));
            ui.end_row();
            ui.label("Maximum E-value");
            ui.add(
                egui::DragValue::new(&mut policy.max_evalue)
                    .speed(0.0001)
                    .range(0.0..=f64::MAX),
            );
            ui.end_row();
            for (label, value) in [
                ("Minimum aligned bases", &mut policy.min_alignment_length_bp),
                ("Maximum chain gap (bp)", &mut policy.max_chain_gap_bp),
                ("Retained loci per target", &mut policy.max_loci_per_target),
                ("HSP processing budget", &mut policy.max_hsps_per_target),
                (
                    "Minimum exact block (bp)",
                    &mut policy.min_conserved_block_bp,
                ),
            ] {
                ui.label(label);
                ui.add(egui::DragValue::new(value));
                ui.end_row();
            }
        });
    let mut promoter_matrix_enabled = policy.promoter_similarity_matrix.is_some();
    if ui
        .checkbox(
            &mut promoter_matrix_enabled,
            "Annotate same-genome hits as transcript-promoter matrix",
        )
        .changed()
    {
        policy.promoter_similarity_matrix =
            promoter_matrix_enabled.then(gp::PromoterSimilarityMatrixPolicy::default);
    }
    if let Some(matrix) = policy.promoter_similarity_matrix.as_mut() {
        egui::Grid::new("promoter_similarity_policy")
            .num_columns(2)
            .show(ui, |ui| {
                for (label, value) in [
                    ("Promoter upstream (bp)", &mut matrix.upstream_bp),
                    ("Promoter downstream (bp)", &mut matrix.downstream_bp),
                    ("Displayed promoter rows", &mut matrix.max_rows),
                ] {
                    ui.label(label);
                    ui.add(egui::DragValue::new(value));
                    ui.end_row();
                }
            });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conservation_request_import_preserves_targets_policy_and_binding() {
        let current = gp::GenomicRegionHomologyScreenRequest {
            set_id: "synthetic".into(),
            region_id: "roi".into(),
            expected_region_content_sha256: Some("sha256:current".into()),
            ..Default::default()
        };
        let mut input = current.clone();
        input.catalog_path = Some("custom.json".into());
        input.policy.max_hsps_per_target = 123;
        input.targets.push(gp::GenomicRegionHomologyTargetRequest {
            genome_id: "synthetic_ortholog".into(),
            required: true,
            role: gp::GenomicRegionHomologyTargetRole::ExpectedOrtholog,
            expected_loci: vec![gp::GenomicRegionHomologyExpectedLocus {
                evidence_id: "orthology:1".into(),
                ..Default::default()
            }],
        });
        let bytes = serde_json::to_vec(&input).expect("serialize");
        assert_eq!(
            parse_bound_conservation_request(&bytes, &current).expect("import"),
            input
        );
        input.expected_region_content_sha256 = Some("sha256:old".into());
        assert!(
            parse_bound_conservation_request(
                &serde_json::to_vec(&input).expect("serialize"),
                &current
            )
            .is_err()
        );
        input.expected_region_content_sha256 = None;
        assert_eq!(
            parse_bound_conservation_request(
                &serde_json::to_vec(&input).expect("serialize"),
                &current
            )
            .expect("bind")
            .expected_region_content_sha256,
            current.expected_region_content_sha256
        );
        input.region_id = "another".into();
        assert!(
            parse_bound_conservation_request(
                &serde_json::to_vec(&input).expect("serialize"),
                &current
            )
            .is_err()
        );
    }
}
