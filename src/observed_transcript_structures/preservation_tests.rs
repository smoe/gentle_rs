//! S2.0a frozen S0/v1 boundaries, using only the parent module's hand-crafted toy loci.
//!
//! Recreate with the `s2a_` library-test filter. Fingerprints are SHA-256 over
//! actual serialized bytes at the pre-S2.0b production baseline, not new biology.
//! No external fixture, private sequence or expectation-regeneration mode exists.

use super::*;
use crate::digest_utils::sha256_prefixed_bytes;

fn bytes<T: Serialize>(value: &T) -> Vec<u8> {
    serde_json::to_vec(value).unwrap()
}

fn cases() -> Vec<Input> {
    let a = input();
    let mut reverse = a.clone();
    reverse.queries[0].masks = Some(vec![(8, 12)]);
    reverse.queries[0].segments[0] = segment(
        100,
        vec![('M', 10), ('N', 10), ('M', 10), ('N', 10), ('M', 10)],
        '-',
        &[(8, 12)],
    );
    let mut unavailable = a.clone();
    unavailable.queries[0].search_complete = false;
    unavailable.queries[0].masks = None;
    unavailable
        .unassessed_records
        .push(unavailable.queries[0].record.clone());
    let mut alternatives = a.clone();
    let mut other = segment(
        200,
        vec![('M', 10), ('N', 15), ('M', 10), ('N', 15), ('M', 10)],
        '-',
        &[],
    );
    other.target_name = "contig_b".into();
    alternatives.queries[0].segments.push(other);
    let mut repeated = alternatives.clone();
    let mut q = repeated.queries[0].clone();
    q.query_id = "query_a".into();
    q.record.record_index = 7;
    repeated.queries.insert(0, q);
    vec![a, reverse, unavailable, alternatives, repeated]
}

#[test]
fn observed_structures_s2a_frozen_bytes_ids_work_and_order() {
    let fingerprints = cases()
        .iter()
        .enumerate()
        .map(|(i, input)| {
            let s = summarize(input).unwrap();
            format!(
                "{i} {} {} {} {} {} {} {} {} {}",
                sha256_prefixed_bytes(&bytes(input)),
                sha256_prefixed_bytes(&bytes(&s)),
                sha256_prefixed_bytes(&serde_json::to_vec_pretty(&s).unwrap()),
                s.context_sha256,
                s.comparison_work,
                sha256_prefixed_bytes(&bytes(&s.observations)),
                sha256_prefixed_bytes(&bytes(&s.records)),
                sha256_prefixed_bytes(&bytes(&s.families)),
                serde_json::to_string(&(
                    s.observations
                        .iter()
                        .map(|o| &o.observation_id)
                        .collect::<Vec<_>>(),
                    s.families.iter().map(|f| &f.chain_id).collect::<Vec<_>>(),
                ))
                .unwrap(),
            )
        })
        .collect::<Vec<_>>()
        .join("\n");
    // Frozen at e84b3f73 with unchanged runtime code. Columns: input/compact/pretty
    // SHA, context SHA, work, observation/record/family SHAs, then ordered IDs.
    assert_eq!(
        fingerprints,
        r#"0 sha256:430b4529ec73f86700ffb6e3cd3d37860496df7cf2cd0247d6f5fdbd58f9f8fb sha256:cff3791f12b0a3b0276a5cbdad0425855899eecc0a1132fd72e9707c0fedbbcc sha256:051f9c47e131c701896f379f990d2bed0597417e071576a7fbf379812f5e7cda sha256:dad4f1832db387aa7b62c7418bedb03c38fcb355021ccb61e2d54aaa95a89377 54 sha256:d93b9e3dbe02544cd0993455273f2f76aac353d2879a682125c6c6720dbb3603 sha256:bdaa06f16e0bae823afe3dbeb9799cbb108dd59cef2dca1b194271d187121607 sha256:fd9bcf182bd78d61f7e595553b0d1a10c21150f4064e1d3a39d1e6b7e3d7b607 [["sha256:08d86f8f51ca860fcc23e0d8d43c7075c5eff5904785f51616ce6036b0cd7040"],["sha256:e691d25883d50f3dc1f1568fb977db69b2aebfe2fda3f50a5a47abc401f61c23"]]
1 sha256:56e4afc68db8f954124f1e23f832c357e0942fbaafe8f15c7ebc0d66e1273188 sha256:6b59ef36ad80cfe3207ca2a2b8010c72c9f67e16c1706abff25c8c064f717992 sha256:b25cae6be66c7b975b1e9dc8331a4c1e9a5c0909799a80cd55b900fdcd05eed9 sha256:d3ee6bc5476df9c10707a134bc9940d397e0441a4e4ec5de0033a2b681f0914b 51 sha256:1b5f35bef7c9323176962f92c404e7b8293c3a6ce5ffff992c9e2128740137ae sha256:1add76a4beb63848d36de72c4ba50c149c89466cbb015851a70b13ebc729b1a3 sha256:fd9bcf182bd78d61f7e595553b0d1a10c21150f4064e1d3a39d1e6b7e3d7b607 [["sha256:6b9a07c8ee5c92fd73654cbb2385e6d1abfc71bcb9da4ff133a8a9baf4e94029"],["sha256:e691d25883d50f3dc1f1568fb977db69b2aebfe2fda3f50a5a47abc401f61c23"]]
2 sha256:e168538fb11fe2395388604fe12c907ffecf4dfd9e68fd186ada0b20cec095ae sha256:b7259f2b9b56c4731a30dfd4301c40ad7800ff8735bac3db86f7611d4b4ac3e7 sha256:5df38dd40582170e47ab2df84a466366543d8b6c217e7b9b52c51fd0feb2043f sha256:b0d1453f2e24849b84db32b6af0b30c4813f7cfc221d5975d66f18fe529551ef 24 sha256:5b09a3dd08c9aa7d927da70e9809456a566a3a130b98b10cdb6cf0144a373605 sha256:cac842a28b8c55767ace13c0a63bc7f43ff4549832c377137fca3c1418b8983d sha256:fcdad0c341a248fe92c93b8c0ee7f7d3430739157d425bd1736e695d65232e25 [["sha256:49806399b75c5dbd72ca477c38d14997f427b8e91f735ae2d7b5b5495a6ed684"],["sha256:e691d25883d50f3dc1f1568fb977db69b2aebfe2fda3f50a5a47abc401f61c23"]]
3 sha256:9c5de758ac6da898916a41a8838693b6069297b2a5aea9e273ade4c6f03364b1 sha256:7cbbb76ee53ad9a96da777d6624f6ca5ba0be90e424c523f1a0cd1677693068d sha256:337785d2e23b92d5a39e6dd1c3297cf0eeb8d6083013fea51de3f7e197aa0d48 sha256:3a3e794d7cc0fc7ac100fb00e3bc20600372b6760a49665ce50feed50d2d8d4b 71 sha256:1b384fe2c5dc8dcc5dbac65d5f32b2c4d0b922039ca5d9d68b8823b63d220abb sha256:c4299d438fd60a00857d268542f30509d9854947a30fe957e2c8e9dc9c0ce04e sha256:62a3695257497a24fb0e18734c75eb16f91ca02d7c6cfcd804c7f8511eea1465 [["sha256:08d86f8f51ca860fcc23e0d8d43c7075c5eff5904785f51616ce6036b0cd7040","sha256:56f070240da12406c83e3dad5608f4ad29f1e7457fbdb572e065cdab60a9758b"],["sha256:20379db973743fada55c73b574664137aa174404b0053c33af8373fe96aa7a54","sha256:e691d25883d50f3dc1f1568fb977db69b2aebfe2fda3f50a5a47abc401f61c23"]]
4 sha256:438b3ab1ce3e251d1dd329cbb0a151c1aecd2d19ac0a350f30a977965e148973 sha256:20edf3721cf6d90655c58184fcbdf3b829f78f20d322517c08ed8db7be79ef5e sha256:d27e5f75a20f26c5b0b2242112142c980ef08f31f8b0614c4406350359f369b0 sha256:28e63ec694b9c7367a1044644d607ab839d29cd0fc94d5edfceab65ca16af40d 142 sha256:0777068118e2fb5a5d2424b2b259b8cbaf09de6479b3ef8caa8e131876b6c0c2 sha256:a4dd5c831f7a3ebf884f2225f2c8699d83fda992587e7c39e70ef87d04b4798e sha256:89c5e57117aeea9ed67f84ac61360b2a4df4d120529e56e9e14363d358003f30 [["sha256:08d86f8f51ca860fcc23e0d8d43c7075c5eff5904785f51616ce6036b0cd7040","sha256:56f070240da12406c83e3dad5608f4ad29f1e7457fbdb572e065cdab60a9758b","sha256:c5001703196b68d58c46ddf56544ff31b9543bb0f91ff45b17e73e817b8ae662","sha256:dd3941f686cd2ad2ff1e4863d8d4f8142213f4e697c6d11e87a9bb5aa6af541c"],["sha256:20379db973743fada55c73b574664137aa174404b0053c33af8373fe96aa7a54","sha256:e691d25883d50f3dc1f1568fb977db69b2aebfe2fda3f50a5a47abc401f61c23"]]"#
    );
}

#[test]
fn observed_structures_s2a_query_and_node_order_is_not_input_order() {
    let mut a = cases().pop().unwrap();
    let expected = summarize(&a).unwrap();
    a.queries.reverse();
    for q in &mut a.queries {
        q.segments.reverse();
    }
    let mut reordered = summarize(&a).unwrap();
    assert_ne!(expected.context_sha256, reordered.context_sha256);
    // Only input-order provenance changes; compare all ordered result fields.
    reordered.context_sha256 = expected.context_sha256.clone();
    assert_eq!(bytes(&expected), bytes(&reordered));
}

#[test]
fn observed_structures_s2a_all_segment_fields_remain_in_the_context() {
    let a = input();
    let base = summarize(&a).unwrap();
    let fields = serde_json::to_value(&a.queries[0].segments[0]).unwrap();
    assert_eq!(
        fields
            .as_object()
            .unwrap()
            .keys()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        [
            "aligned_query_intervals",
            "alignment_detail_verified",
            "alignment_operations",
            "edits",
            "gene_ids",
            "mapq",
            "query_interval",
            "reported_identity_per_mille",
            "reported_role",
            "splice_edges",
            "strand",
            "target_interval",
            "target_name",
            "transcript_ids",
        ]
    );
    // These fields do not determine observed geometry, but must still bind context.
    for (name, value) in [
        ("mapq", serde_json::json!(0)),
        ("reported_role", serde_json::json!("secondary")),
        ("reported_identity_per_mille", serde_json::json!(123)),
        ("edits", serde_json::json!(9)),
        ("gene_ids", serde_json::json!(["another_untrusted_label"])),
        (
            "transcript_ids",
            serde_json::json!(["another_untrusted_tx"]),
        ),
        (
            "splice_edges",
            serde_json::json!([{
                "contig":"contig_a", "donor":110, "acceptor":120,
                "strand":"+", "transcript_ids":["untrusted_edge"]
            }]),
        ),
    ] {
        let mut changed = a.clone();
        let mut seg = fields.clone();
        seg[name] = value;
        changed.queries[0].segments[0] = serde_json::from_value(seg).unwrap();
        let s = summarize(&changed).unwrap();
        assert_eq!(s.context_sha256, digest(&changed).unwrap(), "{name}");
        assert_ne!(s.context_sha256, base.context_sha256, "{name}");
        assert_eq!(bytes(&s.families), bytes(&base.families), "{name}");
        assert_eq!(s.comparison_work, base.comparison_work, "{name}");
    }
    let mut absent = serde_json::to_value(&a).unwrap();
    absent.as_object_mut().unwrap().remove("unassessed_records");
    let decoded: Input = serde_json::from_value(absent).unwrap();
    assert_eq!(bytes(&decoded), bytes(&a));
    assert_eq!(bytes(&summarize(&decoded).unwrap()), bytes(&base));
}

#[test]
fn observed_structures_s2a_error_strings_and_stored_support_are_preserved() {
    let errors = (0..10)
        .map(|case| {
            let mut a = input();
            match case {
                0 => a.queries[0].segments[0].aligned_query_intervals = vec![(0, 1)],
                1 => a.queries[0].segments[0].target_interval.1 = 149,
                2 => a.annotation_sources[0].annotation.assembly = "other".into(),
                3 => a.budgets.max_total_work = 1,
                4 => a.queries[0].masks = Some(vec![(8, 12), (11, 15)]),
                5 => a.queries[0].segments[0].alignment_operations[0].0 = 'S',
                6 => a.queries[0].segments[0].alignment_operations[0].0 = '\u{14d}',
                7 => a.queries[0].segments[0].alignment_operations[0].1 = 0,
                8 => a.queries[0].target_lengths.clear(),
                _ => a.annotation_sources[0].annotation.transcripts[0].target_length = 999,
            }
            let error = summarize(&a).unwrap_err();
            format!("{case}: {}", error.message)
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert_eq!(
        errors,
        r#"0: observed_structure_stored_support_mismatch
1: primer_payload_cigar_span
2: observed_structure_annotation_binding
3: primer_payload_chain_comparison_budget
4: observed_structure_record_or_query
5: primer_payload_cigar_operation
6: observed_structure_segment_geometry
7: primer_payload_cigar_operation
8: observed_structure_segment_geometry
9: observed_structure_reference_length_mismatch"#
    );
    for strand in ['+', '-'] {
        let mut a = input();
        a.queries[0].segments[0].strand = strand;
        assert!(summarize(&a).is_ok());
        a.queries[0].segments[0].aligned_query_intervals = vec![(0, 29)];
        assert_eq!(
            summarize(&a).unwrap_err().message,
            "observed_structure_stored_support_mismatch"
        );
    }
}
