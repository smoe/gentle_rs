"""Fast source/wiring checks; Rust and native GUI acceptance require Glen/CI.

The Rust tests create synthetic DNA/report records entirely in memory. These
guards do not execute that code or establish PCR/scientific correctness.
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


class SimplePcrHandoffSourceTests(unittest.TestCase):
    def test_shared_preparation_is_read_only_and_preserves_exact_selection(self):
        source = (ROOT / "src/engine/analysis/simple_pcr_handoff.rs").read_text()
        production, tests = source.split("#[cfg(test)]", 1)
        for required in (
            "get_primer_design_report(&displayed_report.report_id)",
            "report_bytes != encoded(displayed_report)?",
            "pair.rank == pair_rank",
            "pair_rank_ambiguous",
            "non_annealing_5prime_tail_bp != 0",
            "primer.anneal_length_bp != primer.sequence.len()",
            "template_binding_changed",
            "reverse_primer: pair.reverse.sequence.clone()",
        ):
            self.assertIn(required, production)
        for forbidden in (".apply(", "state_mut(", "std::fs", "write_primer_design_store"):
            self.assertNotIn(forbidden, production)
        self.assertIn("simple_pcr_handoff_engine_refuses_changed_report_and_binding", tests)

    def test_ui_preparation_never_executes_or_saves_and_confirmation_is_bound(self):
        source = (ROOT / "src/main_area_dna/simple_pcr_handoff.rs").read_text()
        production, tests = source.split("#[cfg(test)]", 1)
        self.assertRegex(production, r"engine\s*\.\s*prepare_simple_pcr_pair")
        for required in (
            "engine.instance_id() != pending.engine_instance",
            "self.pcr_forward != pending.previous_forward",
            "current != pending.prepared",
            "self.pcr_reverse = pending.prepared.reverse_primer.clone()",
            'ui.button("Replace PCR primers")',
            'ui.button("Cancel primer replacement")',
        ):
            self.assertIn(required, production)
        for forbidden in ("apply_operation_with_feedback", "save_engine_ops_state", "state_mut(",
                          "self.pcr_adv_forward =", "self.pcr_mut_position ="):
            self.assertNotIn(forbidden, production)
        self.assertIn("simple_pcr_handoff_gui_report_preview_paints_without_copying_or_running_pcr", tests)

    def test_ui_traits_and_paint_follow_current_egui_contract(self):
        source = (ROOT / "src/main_area_dna/simple_pcr_handoff.rs").read_text()
        production, tests = source.split("#[cfg(test)]", 1)
        self.assertRegex(
            production,
            r"#\[derive\(Clone, Debug, Default\)\]\s+pub\(super\) struct SimplePcrHandoffUiState",
        )
        self.assertRegex(
            production,
            r"#\[derive\(Clone, Debug\)\]\s+struct PendingSimplePcrHandoff",
        )
        self.assertIn("ctx.run_ui(egui::RawInput::default(), |ui|", tests)
        self.assertNotIn("ctx.run(", tests)
        self.assertIn("output.drop_without_applying_deltas()", tests)
        self.assertIn("require_traits::<SimplePcrHandoffUiState>()", tests)
        self.assertIn("require_traits::<PendingSimplePcrHandoff>()", tests)
        engine_tests = (ROOT / "src/engine/analysis/simple_pcr_handoff.rs").read_text().split(
            "#[cfg(test)]", 1
        )[1]
        self.assertNotIn("use super::*;", engine_tests)

    def test_existing_report_and_form_are_wired_without_new_operations(self):
        preview = (ROOT / "src/main_area_dna/primer_design_ui.rs").read_text()
        start = preview.index("fn render_primer_design_report_preview")
        end = preview.index("fn render_qpcr_design_report_preview", start)
        preview = preview[start:end]
        self.assertIn('egui::Button::new("Use selected pair in simple PCR")', preview)
        self.assertIn("request_simple_pcr_handoff(Arc::clone(&report), rank)", preview)
        self.assertIn("render_simple_pcr_handoff_confirmation(ui)", preview)
        window = (ROOT / "src/main_area_dna.rs").read_text()
        self.assertIn(".open(self.simple_pcr_handoff.focus_requested.then_some(true))", window)
        self.assertIn("ui.scroll_to_rect(simple_pcr_form.response.rect", window)
        state = window[window.index("struct EngineOpsUiState"):window.index("pub struct MainAreaDna")]
        self.assertNotIn("simple_pcr_handoff", state)


if __name__ == "__main__":
    unittest.main()
