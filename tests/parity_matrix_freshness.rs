use std::{fs, path::Path};

#[test]
fn gui_cli_mcp_parity_matrix_is_fresh() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/gui_cli_mcp_parity.md");
    let actual =
        fs::read_to_string(&path).unwrap_or_else(|err| panic!("read {}: {err}", path.display()));
    let expected = gentle_protocol::render_gui_cli_mcp_parity_matrix_markdown();
    assert_eq!(
        actual, expected,
        "docs/gui_cli_mcp_parity.md is stale; run scripts/regenerate_parity_matrix.sh"
    );
}
