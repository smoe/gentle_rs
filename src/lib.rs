//! GENtle shared library crate.
//!
//! This crate exposes the shared biology engine, rendering paths, and adapter
//! bridges used by GUI, CLI, and optional JavaScript/Lua frontends.

use amino_acids::AminoAcids;
use dna_ladder::{LadderCatalog, default_dna_ladders, default_rna_ladders};
use enzymes::Enzymes;
use std::sync::LazyLock;

/// About/help metadata and version presentation helpers.
pub mod about;
/// Agent-assistant bridge models, transports, and execution guardrails.
pub mod agent_bridge;
/// Amino-acid lookup tables and codon translation helpers.
pub mod amino_acids;
/// Top-level GUI application wiring and event loop state.
pub mod app;
/// Shared DNA display configuration and visibility policies.
pub mod dna_display;
/// DNA/RNA ladder catalogs and migration helper utilities.
pub mod dna_ladder;
/// Core DNA sequence model and sequence-level biological operations.
pub mod dna_sequence;
/// Internal egui compatibility helpers used during GUI API migrations.
pub(crate) mod egui_compat;
/// Shared deterministic operation engine and state model.
pub mod engine;
/// Shared shell command parser/executor reused by GUI shell and CLI shell mode.
pub mod engine_shell;
/// Ensembl protein-entry parsing contracts and REST normalization helpers.
pub mod ensembl_protein;
/// Restriction-enzyme catalog loading and convenience selection helpers.
pub mod enzymes;
/// Expert-view data contracts for feature-centric deep-inspection UIs.
pub mod feature_expert;
/// Feature location/strand utilities used across render and engine code.
pub mod feature_location;
/// GC-content computations and display helpers.
pub mod gc_contents;
/// Genome catalog, preparation, indexing, extraction, and BLAST integration.
pub mod genomes;
/// Gibson assembly planning contracts and deterministic preview derivation.
pub mod gibson_planning;
/// Embedded icon/resource helpers for GUI rendering.
pub mod icons;
/// IUPAC nucleotide-code conversion and validation helpers.
pub mod iupac_code;
/// JavaScript adapter wrappers over shared engine contracts.
#[cfg(feature = "js-interface")]
pub mod js_interface;
/// Lineage graph export and serialization utilities.
pub mod lineage_export;
/// Adaptive routing decisions for linear DNA base-letter rendering.
pub mod linear_base_routing;
/// Lua adapter wrappers over shared engine contracts.
#[cfg(feature = "lua-interface")]
pub mod lua_interface;
/// Main sequence-window GUI controller and interaction orchestration.
pub mod main_area_dna;
/// MCP stdio server adapter exposing tool execution and capability discovery.
pub mod mcp_server;
/// Methylation-site detection and state helpers.
pub mod methylation_sites;
/// NCBI GenBank XML (`GBSet/GBSeq`) parsing and dialect detection helpers.
pub mod ncbi_genbank_xml;
/// Open-reading-frame detection logic.
pub mod open_reading_frame;
/// Virtual pool gel model and rendering primitives.
pub mod pool_gel;
/// Protease digest definitions and helpers.
pub mod protease;
/// Protocol-cartoon catalog and deterministic SVG rendering helpers.
pub mod protocol_cartoon;
/// Position-specific scoring matrix (motif) primitives.
pub mod pssm;
/// Shared DNA rendering entry points.
pub mod render_dna;
/// Circular-map DNA renderer.
pub mod render_dna_circular;
/// Linear-map DNA renderer.
pub mod render_dna_linear;
/// Shared export surfaces (SVG and snapshot pathways).
pub mod render_export;
/// Feature expert-view SVG renderer.
pub mod render_feature_expert;
/// Sequence export helpers and render-side formatting.
pub mod render_sequence;
/// Runtime status helpers for built-in and overrideable external resources.
pub mod resource_status;
/// Resource synchronization (REBASE/JASPAR) parsing and snapshot writing.
pub mod resource_sync;
/// Restriction-enzyme site model and cut geometry utilities.
pub mod restriction_enzyme;
/// RNA structure wrappers and tool integration glue.
pub mod rna_structure;
/// Shared wheel/key/cursor interaction policy for scroll, pan, and zoom.
pub mod scroll_input_policy;
/// Sequence-row abstraction shared by specialized row renderers.
pub mod sequence_rows;
/// Blank-row renderer implementation.
pub mod sequence_rows_blank;
/// DNA base-row renderer implementation.
pub mod sequence_rows_dna;
/// Restriction-enzyme sequence-row renderer implementation.
pub mod sequence_rows_restriction_enzymes;
/// Combined readiness summary for prepared references, helpers, and resources.
pub mod service_readiness;
/// Generated/derived shell-help documentation helpers.
pub mod shell_docs;
/// Hidden shared support helpers for internal tests across library and binaries.
#[doc(hidden)]
pub mod test_support;
/// TF-motif registry and matching support.
pub mod tf_motifs;
/// Process-local tool-path override registry.
pub mod tool_overrides;
/// UniProt/SWISS-PROT parsing contracts and projection payload models.
pub mod uniprot;
/// Generic GUI window abstraction.
pub mod window;
/// Window backdrop configuration and rendering helpers.
pub mod window_backdrop;
/// DNA sequence-window wrapper and per-window controls.
pub mod window_dna;
/// Curated workflow example payloads and templates.
pub mod workflow_examples;

// Restriction enzymes and proteases
pub static ENZYMES: LazyLock<Enzymes> = LazyLock::new(Enzymes::default);

// Amino acids
pub static AMINO_ACIDS: LazyLock<AminoAcids> = LazyLock::new(AminoAcids::default);

// DNA ladders
pub static DNA_LADDERS: LazyLock<LadderCatalog> = LazyLock::new(default_dna_ladders);

// RNA ladders
pub static RNA_LADDERS: LazyLock<LadderCatalog> = LazyLock::new(default_rna_ladders);
