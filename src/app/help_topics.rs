//! Stable, presentation-only window help destinations. Manual navigation never
//! captures a window, sends an agent request, or changes the engine project.

use eframe::egui;

macro_rules! help_topics {
    ($($name:ident => ($topic:literal, $fragment:literal)),+ $(,)?) => {
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub(crate) enum HelpTopic { $($name),+ }

        impl HelpTopic {
            pub(super) fn destination(self) -> (&'static str, &'static str) {
                match self { $(Self::$name => ($topic, $fragment)),+ }
            }

            #[cfg(test)]
            pub(super) const ALL: &[Self] = &[$(Self::$name),+];
        }
    };
}

help_topics! {
    MainWindow => ("main-window-layout", ""),
    DnaViewer => ("dna-window-panels", ""),
    Configuration => ("configuration-window", ""),
    Help => ("help-manuals", ""),
    About => ("about-gentle", ""),
    CommandPalette => ("command-palette-and-historyjobs-panels", ""),
    History => ("status-history-and-background-jobs", ""),
    Jobs => ("status-history-and-background-jobs", ""),
    GenomePrepare => ("reference-genomes-main-menu", "prepare-reference-or-helper-genome"),
    GenomeRetrieve => ("reference-genomes-main-menu", "retrieve-genomic-or-helper-sequence"),
    GenomeBlast => ("reference-genomes-main-menu", "blast-genome-or-helper-sequence"),
    GenomeInspector => ("reference-genomes-main-menu", "prepared-references"),
    GenomeTracks => ("reference-genomes-main-menu", "import-genome-tracks"),
    GenomeCaches => ("reference-genomes-main-menu", "clear-caches"),
    Gibson => ("gibson-workspace", ""),
    Rack => ("rack-view", ""),
    RackLabels => ("rack-view", ""),
    SerialGel => ("serial-gel-ladder-view", ""),
    GelImage => ("measured-gel-images-11-development", ""),
    PcrDesign => ("primer-qpcr-and-sequence-specific-rt-design", ""),
    SequencingConfirmation => ("sequencing-confirmation", ""),
    Planning => ("experiment-planning", ""),
    RoutineAssistant => ("routine-assistant", ""),
    AgentAssistant => ("agent-assistant", ""),
    ExternalServices => ("external-services-window", ""),
    GeneSets => ("gene-set-inspector-and-collection-operations", ""),
    EvidencePreparation => ("array-tracks-in-dna-viewer", ""),
    MicroRna => ("microrna-target-scan", ""),
    ClawBio => ("gentle-shell-gui", "clawbio-window"),
    NewSequence => ("loading-and-creating-sequences", ""),
    Genbank => ("loading-and-creating-sequences", "genbank-and-dbsnp-fetch"),
    ProteinEvidence => ("loading-and-creating-sequences", "protein-evidence-window"),
    Jaspar => ("loading-and-creating-sequences", "jaspar-expert-window"),
    ProteinExpression => ("protein-expression-handoff", ""),
    RnaPanel => ("rna-read-panel", ""),
    Dotplot => ("map-modes-and-dotplot-workspace", ""),
    Splicing => ("splicing-expert-workspace", ""),
    RnaReads => ("rna-read-mapping-workspace", ""),
    PromoterDesign => ("loading-and-creating-sequences", "promoter-design-window"),
    Isoforms => ("isoform-architecture-panels-engine-ops", ""),
}

fn request_id() -> egui::Id {
    egui::Id::new("gentle.manual_help_request")
}

pub(crate) fn request(ctx: &egui::Context, topic: HelpTopic) {
    ctx.data_mut(|data| data.insert_temp(request_id(), topic));
    ctx.request_repaint_of(egui::ViewportId::ROOT);
}

pub(super) fn take_request(ctx: &egui::Context) -> Option<HelpTopic> {
    ctx.data_mut(|data| {
        let topic = data.get_temp::<HelpTopic>(request_id());
        data.remove::<HelpTopic>(request_id());
        topic
    })
}

pub(crate) fn button(ui: &mut egui::Ui, topic: HelpTopic) -> egui::Response {
    let response = ui.button(crate::i18n::tr("button.help"));
    if response.clicked() {
        request(ui.ctx(), topic);
    }
    response
}
