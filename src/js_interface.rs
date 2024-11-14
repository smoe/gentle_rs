use crate::{
    app::GENtleApp, dna_sequence::DNAsequence, methylation_sites::MethylationMode, ENZYMES,
};
use deno_core::*;

#[op2]
#[serde]
fn load_dna(#[string] path: &str) -> Result<DNAsequence, deno_core::anyhow::Error> {
    let mut dna = GENtleApp::load_from_file(path)?;

    // Add default enzymes and stuff
    ENZYMES
        .restriction_enzymes()
        .clone_into(dna.restriction_enzymes_mut());
    dna.set_max_restriction_enzyme_sites(Some(2));
    dna.set_methylation_mode(MethylationMode::both());
    dna.update_computed_features();
    Ok(dna)
}

#[op2]
fn write_gb(
    #[serde] seq: DNAsequence,
    #[string] path: &str,
) -> Result<(), deno_core::anyhow::Error> {
    seq.write_genbank_file(path)?;
    Ok(())
}

#[op2]
#[serde]
fn digest(
    #[serde] seq: DNAsequence,
    #[string] enzymes: &str,
) -> Result<Vec<DNAsequence>, deno_core::anyhow::Error> {
    let enzymes = enzymes.split(',').map(|s| s.trim()).collect::<Vec<_>>();
    let res = ENZYMES.restriction_enzymes_by_name(&enzymes);
    let seqs = seq.restriction_enzymes_full_digest(res);
    Ok(seqs)
}

pub struct JavaScriptInterface {
    runtime: JsRuntime,
}

impl JavaScriptInterface {
    pub fn new() -> Self {
        // Build a deno_core::Extension providing custom ops
        const LOAD_DNA: OpDecl = load_dna();
        const WRITE_GB: OpDecl = write_gb();
        const DIGEST: OpDecl = digest();
        let ext = Extension {
            name: "my_ext",
            ops: std::borrow::Cow::Borrowed(&[LOAD_DNA, WRITE_GB, DIGEST]),
            ..Default::default()
        };

        let mut ret = Self {
            runtime: JsRuntime::new(RuntimeOptions {
                extensions: vec![ext],
                ..Default::default()
            }),
        };
        let init_code = r#"
        	function load_dna(path) {return Deno.core.ops.load_dna(path)}
         	function write_gb(seq,path) {return Deno.core.ops.write_gb(seq,path)}
          	function digest(seq,enzymes) {return Deno.core.ops.digest(seq,enzymes)}
        "#
        .to_string();
        ret.run(init_code);
        println!("Interactive JavaScript Shell (type 'exit' to quit)");
        ret
    }

    pub fn run(&mut self, code: String) {
        match self.runtime.execute_script("<usage>", code) {
            Ok(_) => {}
            Err(e) => eprintln!("{}", e),
        }
    }
}

impl Default for JavaScriptInterface {
    fn default() -> Self {
        Self::new()
    }
}
