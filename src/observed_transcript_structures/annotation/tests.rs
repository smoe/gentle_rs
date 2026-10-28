//! Inline invented GTF, reference metadata and loci; no biological/private input.
//!
//! Recreate via the `observed_transcript_structures::annotation` library filter.
//! Inputs exercise the pure line consumer only, not file verification/publication.

use super::*;
use crate::digest_utils::sha256_prefixed_bytes;

fn sha(c: char) -> String {
    format!("sha256:{}", c.to_string().repeat(64))
}

fn context() -> Context {
    Context {
        reference_id: "invented_reference".into(),
        reference: ReferenceBinding {
            assembly: "invented_assembly".into(),
            fasta_sha256: sha('1'),
        },
        selector: ReferenceGene {
            source_id: "source_a".into(),
            gene_id: "gene_A.7".into(),
        },
        provider: Some("declared_source_label".into()),
        declared_provider_release: Some("caller_release".into()),
        scope: InputScope::ProducerExcerpt,
        policy: Policy {
            transcript_feature_rows: TranscriptRows::Required,
            minimum_reference_gap_bases: 5,
        },
    }
}

fn dictionary(c: &Context) -> Dictionary {
    Dictionary {
        reference_id: c.reference_id.clone(),
        reference: c.reference.clone(),
        contigs: vec![
            Contig {
                name: "ctg_a".into(),
                length: 1000,
            },
            Contig {
                name: "ctg_b".into(),
                length: 2000,
            },
        ],
    }
}

fn limits() -> Limits {
    Limits {
        max_input_bytes: 100_000,
        max_line_bytes: 10_000,
        max_lines: 1000,
        max_dictionary_entries: 10,
        max_indexed_genes: 100,
        max_indexed_transcripts: 100,
        max_exon_rows: 1000,
        max_retained_bytes: 2_000_000,
        max_work: 1_000_000,
        max_diagnostics: 100,
    }
}

fn line(
    feature: &str,
    contig: &str,
    strand: char,
    start: usize,
    end: usize,
    attrs: &str,
) -> String {
    format!("{contig}\tinvented\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")
}

fn exon(gene: &str, tx: &str, contig: &str, strand: char, start: usize, end: usize) -> String {
    line(
        "exon",
        contig,
        strand,
        start,
        end,
        &format!("gene_id \"{gene}\"; transcript_id \"{tx}\";"),
    )
}

fn source() -> String {
    let a = "gene_id \"gene_A.7\"; transcript_id \"tx_A.2\"; gene_version \"7\"; transcript_version \"2\";";
    let b =
        "gene_id \"gene_B\"; transcript_id \"tx_B\"; gene_version \"3\"; transcript_version \"9\";";
    [
        "#!genome-version not_an_authority\n".into(),
        line("exon", "ctg_b", '-', 221, 230, b),
        line(
            "gene",
            "ctg_a",
            '+',
            101,
            130,
            "gene_id \"gene_A.7\"; gene_version \"7\";",
        ),
        line("transcript", "ctg_a", '+', 101, 130, a),
        line("exon", "ctg_a", '+', 121, 130, a),
        line("transcript", "ctg_b", '-', 201, 230, b),
        line("exon", "ctg_a", '+', 101, 110, a),
        line(
            "gene",
            "ctg_b",
            '-',
            201,
            230,
            "gene_id \"gene_B\"; gene_version \"3\";",
        ),
        line("exon", "ctg_b", '-', 201, 210, b),
    ]
    .concat()
}

fn parse(
    c: &Context,
    d: &Dictionary,
    b: &Limits,
    bytes: &[u8],
) -> Result<Preparation, EngineError> {
    let mut p = Consumer::new(c, d, b, &mut || Ok(()))?;
    for raw in bytes.split_inclusive(|&v| v == b'\n') {
        p.push_line(raw, &mut || Ok(()))?;
    }
    p.finish(&mut || Ok(()))
}

fn run(c: &Context, raw: &str) -> Preparation {
    parse(c, &dictionary(c), &limits(), raw.as_bytes()).unwrap()
}

fn unrequired() -> Context {
    let mut c = context();
    c.policy.transcript_feature_rows = TranscriptRows::NotRequired;
    c
}

fn geometry(a: &PrimerPayloadAnnotation) -> serde_json::Value {
    serde_json::to_value(&a.transcripts).unwrap()
}

#[test]
fn gtf_preparation_two_configurable_genes_preserve_geometry_and_source_claims() {
    let c = context();
    let raw = source();
    let a = run(&c, &raw);
    assert_eq!(a.outcome, Outcome::Completed);
    assert_eq!(a.grammar, "gtf_strict.v1");
    assert_eq!(a.coordinate_authority, "declared");
    assert_eq!(a.raw_sha256, sha256_prefixed_bytes(raw.as_bytes()));
    assert_eq!(a.raw_bytes, raw.len());
    assert_eq!(a.counts.lines, 9);
    assert_eq!(
        (a.counts.comments, a.counts.data_rows, a.counts.exon_rows),
        (1, 8, 4)
    );
    assert_eq!((a.indexed_genes, a.indexed_transcripts), (2, 2));
    assert_eq!(a.counts.selected_rows, 4);
    assert_eq!(a.counts.selected_exon_rows, 2);
    assert_eq!(a.provider, c.provider);
    assert_eq!(a.declared_provider_release, c.declared_provider_release);
    assert_eq!(a.scope, InputScope::ProducerExcerpt);
    assert_eq!(a.policy, c.policy);
    assert_eq!(a.selector.gene_id, "gene_A.7");
    assert_eq!(a.gene_version.as_deref(), Some("7"));
    assert_eq!(a.transcript_versions["tx_A.2"].as_deref(), Some("2"));
    let annotation = a.annotation.as_ref().unwrap();
    assert_eq!(annotation.release, format!("snapshot-{}", a.raw_sha256));
    assert_eq!(annotation.schema, PRIMER_PAYLOAD_ANNOTATION_SCHEMA);
    assert_eq!(annotation.transcripts.len(), 1);
    assert_eq!(
        annotation.transcripts[0].genomic_exons,
        [(100, 110), (120, 130)]
    );
    assert_eq!(annotation.transcripts[0].target_length, 1000);
    assert_eq!(annotation.transcripts[0].gene_ids, ["gene_A.7"]);
    assert_eq!(a.included_contigs, ["ctg_a"]);
    let mut other = c.clone();
    other.selector.gene_id = "gene_B".into();
    other.scope = InputScope::WholeSourceFile;
    other.declared_provider_release = None;
    let b = run(&other, &raw);
    assert_eq!(b.outcome, Outcome::Completed);
    assert_eq!(b.raw_sha256, a.raw_sha256);
    assert_eq!(b.annotation.as_ref().unwrap().release, annotation.release);
    let model = &b.annotation.as_ref().unwrap().transcripts[0];
    assert_eq!(model.genomic_exons, [(200, 210), (220, 230)]);
    assert_eq!(model.genomic_strand, Some('-'));
    assert_eq!(model.target_length, 2000);
    assert_eq!(b.declared_provider_release, None);
    assert!(
        b.nonclaims
            .iter()
            .any(|s| s.contains("provider completeness"))
    );
    assert_ne!(
        geometry(b.annotation.as_ref().unwrap()),
        geometry(annotation)
    );
    let shuffled = raw
        .lines()
        .rev()
        .map(|l| format!("{l}\n"))
        .collect::<String>();
    assert_eq!(
        geometry(run(&c, &shuffled).annotation.as_ref().unwrap()),
        geometry(annotation)
    );
}

#[test]
fn gtf_preparation_lf_crlf_and_final_line_keep_geometry_but_not_raw_identity() {
    let c = context();
    let raw = source();
    let a = run(&c, &raw);
    let b = run(&c, &raw.replace('\n', "\r\n"));
    let z = run(&c, raw.trim_end_matches('\n'));
    for p in [&b, &z] {
        assert_eq!(
            geometry(p.annotation.as_ref().unwrap()),
            geometry(a.annotation.as_ref().unwrap())
        );
        assert_ne!(p.raw_sha256, a.raw_sha256);
        assert_ne!(
            p.annotation.as_ref().unwrap().release,
            a.annotation.as_ref().unwrap().release
        );
    }
    let raw = exon("gene_A.7", "one_base", "ctg_a", '+', 101, 101);
    assert_eq!(
        run(&unrequired(), &raw).annotation.unwrap().transcripts[0].genomic_exons,
        [(100, 101)]
    );
}

#[test]
fn gtf_preparation_exact_ids_missing_and_empty_sources_are_not_empty_annotations() {
    for id in ["gene_A", "gene_a.7", "gene_A.8", "unseen"] {
        let mut c = context();
        c.selector.gene_id = id.into();
        let p = run(&c, &source());
        assert_eq!(p.outcome, Outcome::MissingGene, "{id}");
        assert!(p.annotation.is_none());
        assert!(p.included_contigs.is_empty());
        assert!(p.nonclaims.iter().any(|s| s.contains("exact gene_id")));
    }
    let p = run(&context(), "");
    assert_eq!(p.outcome, Outcome::MissingGene);
    assert_eq!(p.raw_sha256, sha256_prefixed_bytes(b""));
    assert_eq!(p.counts.lines, 0);
    let raw = line("gene", "ctg_a", '+', 101, 130, "gene_id \"gene_A.7\";");
    let p = run(&context(), &raw);
    assert_eq!(p.outcome, Outcome::TranscriptModelConflict);
    assert!(p.annotation.is_none());
    assert!(p.selected_reasons.contains_key(&Reason::NoExons));
}

#[test]
fn gtf_preparation_strict_attributes_preserve_unicode_quotes_backslashes_and_unknown_values() {
    let mut c = unrequired();
    c.selector.gene_id = "Géne".into();
    let raw = line(
        "exon",
        "ctg_a",
        '+',
        101,
        110,
        r#"gene_id "Géne"; transcript_id "tx\"quoted\\slash"; tag "a;b"; tag repeated; transcript_version " 9 ";"#,
    );
    let p = run(&c, &raw);
    assert_eq!(p.outcome, Outcome::Completed);
    assert_eq!(
        p.annotation.unwrap().transcripts[0].transcript_id,
        "tx\"quoted\\slash"
    );
    assert_eq!(
        p.transcript_versions["tx\"quoted\\slash"].as_deref(),
        Some(" 9 ")
    );
}

#[test]
fn gtf_preparation_grammar_errors_poison_consumer_even_outside_selected_gene() {
    let base = exon("other_gene", "other_tx", "ctg_b", '-', 101, 110);
    let mut cases = vec![
        "\n".into(),
        "   \n".into(),
        base.replace("\texon\t", "\t\t"),
        base.replace("\tinvented\t", "\t\t"),
        base.replace("\tinvented\t", "\t \t"),
        base.replace("\texon\t", "\t exon\t"),
        base.replacen('\t', " ", 1),
        format!("{}\textra\n", base.trim_end()),
        base.replace("\t101\t", "\t0\t"),
        base.replace("\t101\t", "\t111\t"),
        base.replace("\t101\t", "\t+101\t"),
        base.replace("\t101\t", "\t1e2\t"),
        base.replace("\t110\t", "\t2001\t"),
        base.replace("ctg_b", "missing_contig"),
        base.replace("\t-\t", "\t?\t"),
        base.replace("\t.\t-\t", "\tNaN\t-\t"),
        base.replace("\t-\t.\t", "\t-\t3\t"),
        base.replace("gene_id", "Gene_id"),
        base.replace("gene_id \"other_gene\";", "gene_id other_gene;"),
        base.replace(
            "gene_id \"other_gene\";",
            "gene_id \"other_gene\"; gene_id \"other_gene\";",
        ),
        base.replace("transcript_id", "transcript-id"),
        base.replace("other_gene", "bad id"),
        base.replace("other_gene", ""),
        base.replace("other_gene", &"a".repeat(513)),
        base.replace("transcript_id \"other_tx\";", ""),
        base.replace("\texon\t", "\tgene\t"),
        base.replace("other_gene", r"bad\q"),
        base.replace(";\n", "\n"),
        base.replace(";\n", "; trailing"),
        base.replace("\"other_tx\";", "\"unclosed;"),
        base.replace(";\n", "; tag \"\";\n"),
        base.replace(";\n", "; tag \"a\"b;\n"),
    ];
    cases.push(base.replace("101", &(usize::MAX as u128 + 1).to_string()));
    for (i, raw) in cases.iter().enumerate() {
        let c = context();
        let d = dictionary(&c);
        let b = limits();
        let mut p = Consumer::new(&c, &d, &b, &mut || Ok(())).unwrap();
        assert!(
            p.push_line(raw.as_bytes(), &mut || Ok(())).is_err(),
            "grammar {i}: {raw}"
        );
        assert_eq!(
            p.finish(&mut || Ok(())).unwrap_err().message,
            "observed_gtf_poisoned"
        );
    }
    for raw in [
        b"# invalid \xff\n".as_slice(),
        b"# embedded\0control\n",
        b"# a\nb\n",
        b"# interior\rc\n",
    ] {
        assert!(parse(&context(), &dictionary(&context()), &limits(), raw).is_err());
    }
}

#[test]
fn gtf_preparation_required_transcript_rows_and_hulls_are_explicit_not_completeness_proof() {
    let raw = exon("gene_A.7", "tx", "ctg_a", '+', 101, 110);
    let p = run(&context(), &raw);
    assert_eq!(p.outcome, Outcome::TranscriptModelConflict);
    assert!(
        p.selected_reasons
            .contains_key(&Reason::MissingTranscriptFeature)
    );
    let p = run(&unrequired(), &raw);
    assert_eq!(p.outcome, Outcome::Completed);
    assert!(p.nonclaims.iter().any(|s| s.contains("same-locus")));
    let matching = line(
        "transcript",
        "ctg_a",
        '+',
        101,
        110,
        "gene_id \"gene_A.7\"; transcript_id \"tx\";",
    );
    assert_eq!(
        run(&context(), &(matching.clone() + &raw)).outcome,
        Outcome::Completed
    );
    let divergent = matching.replace("\t110\t", "\t115\t");
    let p = run(&context(), &(divergent + &raw));
    assert_eq!(p.outcome, Outcome::TranscriptModelConflict);
    assert!(
        p.selected_reasons
            .contains_key(&Reason::TranscriptSpanMismatch)
    );
    let p = run(&context(), &(matching.clone() + &matching + &raw));
    assert_eq!(p.outcome, Outcome::DuplicateTranscriptAccession);
    assert!(
        p.selected_reasons
            .contains_key(&Reason::MultipleTranscriptRows)
    );
}

#[test]
fn gtf_preparation_global_identity_index_cannot_hide_selected_collisions() {
    for (gene, ctg, strand, expected) in [
        ("other", "ctg_a", '+', Reason::CrossGeneLink),
        ("gene_A.7", "ctg_b", '+', Reason::MultipleLoci),
        ("gene_A.7", "ctg_a", '-', Reason::MultipleLoci),
    ] {
        let a = exon("gene_A.7", "same_tx", "ctg_a", '+', 101, 110);
        let b = exon(gene, "same_tx", ctg, strand, 221, 230);
        for raw in [a.clone() + &b, b.clone() + &a] {
            let p = run(&unrequired(), &raw);
            assert_eq!(p.outcome, Outcome::DuplicateTranscriptAccession);
            assert!(p.selected_reasons.contains_key(&expected));
            assert!(p.annotation.is_none());
            assert_eq!(p.indexed_transcripts, 1);
        }
    }
    let raw = exon("gene_A.7", "tx_1", "ctg_a", '+', 101, 110)
        + &exon("gene_A.7", "tx_2", "ctg_b", '-', 201, 210);
    assert_eq!(
        run(&unrequired(), &raw).outcome,
        Outcome::AmbiguousGeneLocus
    );
    let raw = exon("gene_A.7", "tx", "ctg_a", '.', 101, 110);
    assert_eq!(
        run(&unrequired(), &raw).outcome,
        Outcome::TranscriptModelConflict
    );
}

#[test]
fn gtf_preparation_duplicate_and_unrelated_model_conflicts_are_counted_not_pruned() {
    let raw = exon("gene_A.7", "tx", "ctg_a", '+', 101, 110);
    let p = run(&unrequired(), &(raw.clone() + &raw));
    assert_eq!(p.outcome, Outcome::IdenticalDuplicateRow);
    assert!(p.annotation.is_none());
    assert_eq!(
        (p.counts.duplicate_rows, p.counts.selected_duplicate_rows),
        (1, 1)
    );
    let unrelated = exon("other", "bad_model", "ctg_b", '.', 201, 210);
    let p = run(&unrequired(), &(raw + &unrelated.clone() + &unrelated));
    assert_eq!(p.outcome, Outcome::Completed);
    assert_eq!(p.counts.duplicate_rows, 1);
    assert_eq!(p.counts.selected_duplicate_rows, 0);
    assert!(p.unrelated_reasons.contains_key(&Reason::UnsupportedStrand));
    assert!(
        p.unrelated_reasons
            .contains_key(&Reason::IdenticalDuplicateRow)
    );
    assert_eq!(p.indexed_transcripts, 2);
    assert_eq!(p.annotation.unwrap().transcripts.len(), 1);
}

#[test]
fn gtf_preparation_version_presence_conflicts_are_order_independent_and_scoped() {
    let rows = [
        exon("gene_A.7", "tx", "ctg_a", '+', 101, 110),
        exon("gene_A.7", "tx", "ctg_a", '+', 121, 130)
            .replace(";\n", "; gene_version \"1\"; transcript_version \"2\";\n"),
        exon("gene_A.7", "tx", "ctg_a", '+', 141, 150)
            .replace(";\n", "; gene_version \"3\"; transcript_version \"4\";\n"),
    ];
    for order in [[0, 1, 2], [2, 0, 1], [1, 2, 0]] {
        let raw = order.iter().map(|&i| rows[i].clone()).collect::<String>();
        let p = run(&unrequired(), &raw);
        assert_eq!(p.outcome, Outcome::TranscriptModelConflict);
        assert!(p.selected_reasons.contains_key(&Reason::VersionConflict));
        assert!(
            p.selected_reasons
                .contains_key(&Reason::VersionPresenceMismatch)
        );
        assert_eq!(p.gene_version, None);
        assert_eq!(p.transcript_versions["tx"], None);
    }
    let mut c = context();
    c.declared_provider_release = Some("different declaration".into());
    let p = run(&c, &source());
    assert_eq!(p.outcome, Outcome::Completed);
    assert_eq!(
        geometry(p.annotation.as_ref().unwrap()),
        geometry(run(&context(), &source()).annotation.as_ref().unwrap())
    );
}

#[test]
fn gtf_preparation_gap_and_overlap_refuse_entire_gene_without_touching_legacy_rules() {
    for (start, outcome, reason) in [
        (
            110,
            Outcome::TranscriptModelConflict,
            Reason::ExonsOverlapOrTouch,
        ),
        (
            111,
            Outcome::TranscriptModelConflict,
            Reason::ExonsOverlapOrTouch,
        ),
        (
            115,
            Outcome::ReferenceGapPolicyRefused,
            Reason::GapBelowPolicy,
        ),
    ] {
        // First interval [100,110); 111 is touching, 115 leaves a 4-base gap.
        let raw = exon("gene_A.7", "tx", "ctg_a", '+', 101, 110)
            + &exon("gene_A.7", "tx", "ctg_a", '+', start, 130)
            + &exon("gene_A.7", "otherwise_valid", "ctg_a", '+', 201, 210);
        let p = run(&unrequired(), &raw);
        assert_eq!(p.outcome, outcome);
        assert!(p.selected_reasons.contains_key(&reason));
        assert!(p.annotation.is_none());
    }
    let raw = exon("gene_A.7", "tx", "ctg_a", '+', 101, 110)
        + &exon("gene_A.7", "tx", "ctg_a", '+', 116, 130);
    assert_eq!(run(&unrequired(), &raw).outcome, Outcome::Completed);
    let mut c = unrequired();
    c.policy.minimum_reference_gap_bases = 6;
    assert_eq!(run(&c, &raw).outcome, Outcome::ReferenceGapPolicyRefused);
    c.policy.minimum_reference_gap_bases = 1;
    let a = run(&c, &raw);
    c.policy.minimum_reference_gap_bases = 5;
    let b = run(&c, &raw);
    assert_eq!(
        serde_json::to_vec(&a.annotation).unwrap(),
        serde_json::to_vec(&b.annotation).unwrap()
    );
    assert_ne!(a.policy, b.policy); // S2.1b must bind this difference in its receipt.
}

#[test]
fn gtf_preparation_dictionary_and_context_admission_are_checked_without_aliases() {
    for case in 0..15 {
        let mut c = context();
        let mut d = dictionary(&c);
        match case {
            0 => d.contigs.clear(),
            1 => d.contigs[0].length = 0,
            2 => d.contigs.push(d.contigs[0].clone()),
            3 => d.contigs[0].name = "*".into(),
            4 => d.contigs[0].name = "contains space".into(),
            5 => d.contigs[0].name = "nul\0name".into(),
            6 => d.contigs[0].name = "a".repeat(4097),
            7 => d.reference.assembly = "other".into(),
            8 => d.reference_id = "other".into(),
            9 => d.reference.fasta_sha256 = sha('2'),
            10 => c.reference.fasta_sha256 = "unknown".into(),
            11 => c.selector.gene_id.clear(),
            12 => c.policy.minimum_reference_gap_bases = 0,
            13 => {
                c.reference.assembly = " ".into();
                d.reference.assembly = " ".into();
            }
            _ => {
                c.reference_id = " ".into();
                d.reference_id = " ".into();
            }
        }
        assert!(
            Consumer::new(&c, &d, &limits(), &mut || Ok(())).is_err(),
            "admission {case}"
        );
    }
    let c = context();
    let mut d = dictionary(&c);
    d.contigs.reverse();
    let p = parse(&c, &d, &limits(), source().as_bytes()).unwrap();
    assert_eq!(p.outcome, Outcome::Completed);
    assert_eq!(d.contigs[0].name, "ctg_b");
    assert_eq!(p.annotation.unwrap().transcripts[0].target_name, "ctg_a");
}

#[test]
fn gtf_preparation_unknown_strand_is_not_a_second_locus_and_id_limits_are_bytes() {
    let raw = exon("gene_A.7", "tx", "ctg_a", '+', 101, 110)
        + &line(
            "CDS",
            "ctg_a",
            '.',
            101,
            105,
            "gene_id \"gene_A.7\"; transcript_id \"tx\";",
        );
    let p = run(&unrequired(), &raw);
    assert_eq!(p.outcome, Outcome::TranscriptModelConflict);
    assert!(p.selected_reasons.contains_key(&Reason::UnsupportedStrand));
    assert!(!p.selected_reasons.contains_key(&Reason::MultipleLoci));
    assert!(!p.selected_reasons.contains_key(&Reason::AmbiguousGeneLocus));
    assert_eq!(p.counts.features["CDS"], 1);
    let mut c = unrequired();
    c.selector.gene_id = "é".repeat(256);
    let raw = exon(&c.selector.gene_id, "tx", "ctg_a", '+', 101, 110);
    assert_eq!(run(&c, &raw).outcome, Outcome::Completed);
    let raw = exon(&"é".repeat(257), "tx", "ctg_a", '+', 101, 110);
    assert!(parse(&c, &dictionary(&c), &limits(), raw.as_bytes()).is_err());
}

#[test]
fn gtf_preparation_every_resource_bound_has_exact_and_one_over_coverage() {
    let c = context();
    let d = dictionary(&c);
    let raw = source();
    let baseline = run(&c, &raw);
    for case in 0..9 {
        let mut b = limits();
        let exact = match case {
            0 => raw.len(),
            1 => raw.split_inclusive('\n').map(str::len).max().unwrap(),
            2 => baseline.counts.lines,
            3 => d.contigs.len(),
            4 => baseline.indexed_genes,
            5 => baseline.indexed_transcripts,
            6 => baseline.counts.exon_rows,
            7 => baseline.retained_peak_bytes,
            _ => baseline.work,
        };
        let set = |b: &mut Limits, n| match case {
            0 => b.max_input_bytes = n,
            1 => b.max_line_bytes = n,
            2 => b.max_lines = n,
            3 => b.max_dictionary_entries = n,
            4 => b.max_indexed_genes = n,
            5 => b.max_indexed_transcripts = n,
            6 => b.max_exon_rows = n,
            7 => b.max_retained_bytes = n,
            _ => b.max_work = n,
        };
        set(&mut b, exact);
        assert_eq!(
            parse(&c, &d, &b, raw.as_bytes()).unwrap().outcome,
            Outcome::Completed,
            "exact {case}"
        );
        set(&mut b, exact - 1);
        assert!(
            parse(&c, &d, &b, raw.as_bytes()).is_err(),
            "one over {case}"
        );
    }
    for bad in [0, usize::MAX] {
        let mut b = limits();
        b.max_work = bad;
        assert!(Consumer::new(&c, &d, &b, &mut || Ok(())).is_err());
    }
}

#[test]
fn gtf_preparation_whole_file_limits_include_comments_and_unselected_models() {
    let c = unrequired();
    let raw = exon("unselected", "x", "ctg_b", '-', 201, 210)
        + &exon("unselected2", "y", "ctg_b", '-', 221, 230);
    let mut b = limits();
    b.max_indexed_genes = 1;
    assert!(parse(&c, &dictionary(&c), &b, raw.as_bytes()).is_err());
    b = limits();
    b.max_indexed_transcripts = 1;
    assert!(parse(&c, &dictionary(&c), &b, raw.as_bytes()).is_err());
    b = limits();
    b.max_exon_rows = 1;
    assert!(parse(&c, &dictionary(&c), &b, raw.as_bytes()).is_err());
    let mut b = limits();
    b.max_lines = 1;
    assert!(parse(&c, &dictionary(&c), &b, b"# a\n# b\n").is_err());
}

#[test]
fn gtf_preparation_diagnostic_cap_does_not_drop_counts_or_strengthen_outcomes() {
    let c = unrequired();
    let raw = exon("gene_A.7", "tx", "ctg_a", '.', 101, 110)
        + &exon("gene_A.7", "tx", "ctg_a", '.', 105, 115);
    let baseline = run(&c, &raw);
    assert!(baseline.diagnostic_count > 1);
    for max in [0, 1, baseline.diagnostic_count] {
        let mut b = limits();
        b.max_diagnostics = max;
        let p = parse(&c, &dictionary(&c), &b, raw.as_bytes()).unwrap();
        assert_eq!(p.outcome, baseline.outcome);
        assert_eq!(p.selected_reasons, baseline.selected_reasons);
        assert_eq!(p.diagnostic_count, baseline.diagnostic_count);
        assert_eq!(p.diagnostics.len(), max);
        assert_eq!(p.omitted_diagnostics, p.diagnostic_count - max);
    }
}

#[test]
fn gtf_preparation_control_error_propagates_at_every_bounded_poll_without_partial_result() {
    let c = context();
    let d = dictionary(&c);
    let b = limits();
    let raw = source();
    let execute = |check: &mut dyn FnMut() -> Result<(), EngineError>| {
        let mut p = Consumer::new(&c, &d, &b, check)?;
        for line in raw.split_inclusive('\n') {
            p.push_line(line.as_bytes(), check)?;
        }
        p.finish(check)
    };
    let mut polls = 0;
    execute(&mut || {
        polls += 1;
        Ok(())
    })
    .unwrap();
    for stop in 1..=polls {
        let mut seen = 0;
        let result = execute(&mut || {
            seen += 1;
            if seen == stop {
                Err(bad("invented_cancel_or_control_error"))
            } else {
                Ok(())
            }
        });
        assert_eq!(
            result.unwrap_err().message,
            "invented_cancel_or_control_error",
            "poll {stop}"
        );
    }
    let mut p = Consumer::new(&c, &d, &b, &mut || Ok(())).unwrap();
    assert!(
        p.push_line(b"# no terminal newline", &mut || Ok(()))
            .is_ok()
    );
    assert_eq!(
        p.push_line(b"# another\n", &mut || Ok(()))
            .unwrap_err()
            .message,
        "observed_gtf_line_after_unterminated"
    );
    assert!(p.finish(&mut || Ok(())).is_err());
}
