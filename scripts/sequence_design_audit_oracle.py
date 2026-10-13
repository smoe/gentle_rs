"""Literal, bounded independent oracles for sequence-design acceptance.

All inputs are hand-crafted synthetic inserts, not genes. Recreate them with
audit_cases(). This test-only module does not import GENtle, approximate its
search, or serve as a production translator/designer. Exhaustive spaces have
at most 36 words; GC checks use cross multiplication, not the engine's helpers.
"""

import copy
import hashlib
import itertools
import json


SYNONYMS = {
    "M": ("ATG",), "F": ("TTT", "TTC"), "K": ("AAA", "AAG"),
    "E": ("GAA", "GAG"),
    "L": ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
}
IUPAC = dict(zip("ACGTRYSWKMBDHVN", (
    "A", "C", "G", "T", "AG", "CT", "CG", "AT", "GT", "AC",
    "CGT", "AGT", "ACT", "ACG", "ACGT")))


def bounds(minimum, maximum):
    return {"min_basis_points": minimum, "max_basis_points": maximum}


def request(sequence="ATGTTTAAATAA", protein="MFK", start=0):
    return {
        "schema": "gentle.dna_sequence_design_request.v1",
        "target": {"kind": "inline_sequence", "sequence": sequence},
        "purpose": "synthetic_coding_insert",
        "cds": {"start_0based": start, "end_0based_exclusive": start + 3 * (len(protein) + 1)},
        "protein_sequence": protein, "genetic_code": 1,
        "protected_intervals": [], "avoid_motifs": [],
        "max_evaluations": 4096, "output_seq_id": "synthetic_audit_output",
    }


def audit_cases():
    """Fixed cases exercise global/local constraints separately and together."""
    cases = []
    for name, motifs, gc, window in (
        ("motif_only", ["TTTAAA"], None, None),
        ("coupled_motifs", ["TTTAAA", "TTCAAA", "TTTAAG"], None, None),
        ("motif_infeasible", ["TTTAAA", "TTCAAA", "TTTAAG", "TTCAAG"], None, None),
        ("gc_only", [], bounds(1666, 1667), None),
        ("coupled_motif_gc", ["GAATTC"], bounds(2500, 2500), None),
        ("window_only", [], None, dict(bounds(1666, 1667), window_bp=6)),
        ("joint", ["TTTAAA"], bounds(1666, 1667), dict(bounds(1666, 1667), window_bp=6)),
        ("empty_count_window", [], bounds(1667, 1667), None),
        ("frozen_gc_infeasible", [], bounds(0, 0), None),
    ):
        item = (request("ATGGAATTCTAA", "MEF") if name == "coupled_motif_gc" else request())
        item["avoid_motifs"] = [{"pattern": p, "strand": "both"} for p in motifs]
        if gc is not None:
            item["gc_content"] = gc
        if window is not None:
            item["gc_window"] = window
        cases.append((name, item))
    for size in (2, 5, 7, 12):
        item = request("CATGCTACGTTAAG", "MLR", 1)
        item["avoid_motifs"] = [{"pattern": "CTN", "strand": "reverse"}]
        item["gc_content"] = bounds(2500, 3334)
        item["gc_window"] = dict(bounds(0, 5000), window_bp=size)
        cases.append((f"partial_window_{size}", item))
        protected = copy.deepcopy(item)
        protected["protected_intervals"] = [{"start_0based": 6, "end_0based_exclusive": 7}]
        cases.append((f"protected_window_{size}", protected))
    item = request("GATGTTTAAATAAC", "MFK", 1)
    item["avoid_motifs"] = [{"pattern": "GATGTTT", "strand": "forward"}]
    cases.append(("flank_boundary", item))
    return cases


def words(item):
    source = item["target"]["sequence"]
    start, end = item["cds"]["start_0based"], item["cds"]["end_0based_exclusive"]
    protein = item["protein_sequence"]
    choices = [SYNONYMS[residue] for residue in protein]
    count = 1
    for choice in choices:
        count *= len(choice)
    if (count > 36 or len(source) > 20 or protein[0] != "M"
            or end - start != 3 * (len(protein) + 1) or source[end - 3:end] != "TAA"):
        raise ValueError("Oracle admits only its tiny synthetic spaces")
    return [source[:start] + "".join(codons) + "TAA" + source[end:]
            for codons in itertools.product(*choices)]


def motif_present(sequence, pattern, strand):
    masks = [IUPAC[letter] for letter in pattern]
    reverse = ["".join(base.translate(str.maketrans("ACGT", "TGCA")) for base in mask)
               for mask in reversed(masks)]
    directions = ([masks] if strand == "forward" else [reverse] if strand == "reverse"
                  else [masks, reverse])
    return any(all(base in mask for base, mask in zip(sequence[start:start + len(pattern)], row))
               for row in directions for start in range(len(sequence) - len(pattern) + 1))


def within_gc(sequence, limits):
    count = sum(base in "GC" for base in sequence)
    return (limits["min_basis_points"] * len(sequence) <= count * 10000
            <= limits["max_basis_points"] * len(sequence))


def satisfies(item, sequence):
    source = item["target"]["sequence"]
    start, end = item["cds"]["start_0based"], item["cds"]["end_0based_exclusive"]
    if sequence not in words(item):
        return False
    for span in item["protected_intervals"]:
        a, b = span["start_0based"], span["end_0based_exclusive"]
        if sequence[a:b] != source[a:b]:
            return False
    if any(motif_present(sequence, motif["pattern"], motif["strand"])
           for motif in item["avoid_motifs"]):
        return False
    cds = sequence[start:end]
    if "gc_content" in item and not within_gc(cds, item["gc_content"]):
        return False
    if "gc_window" in item:
        policy = item["gc_window"]
        size = policy["window_bp"]
        if not all(within_gc(cds[pos:pos + size], policy) for pos in range(len(cds) - size + 1)):
            return False
    return True


def expected_output(item):
    source = item["target"]["sequence"]
    candidates = [word for word in words(item) if satisfies(item, word)]
    return min(candidates, key=lambda word: (sum(a != b for a, b in zip(source, word)), word),
               default=None)


def assert_report(item, report):
    """Verify complete tiny searches; counts/history need not agree by strategy."""
    def require(condition, label):
        if not condition:
            raise AssertionError(label)

    expected = expected_output(item)
    normalized = copy.deepcopy(item)
    if normalized.get("search_strategy") == "full_enumeration":
        del normalized["search_strategy"]
    require(report["request"] == normalized, "exact request policy")
    require(report["source_sequence"] == item["target"]["sequence"], "exact source")
    require(report["status"] == ("feasible" if expected is not None else "proven_infeasible"), "status")
    require(report["output_sequence"] == expected, "independent minimum-edit DNA")
    if expected is None:
        require(report["approval_digest"] is None, "infeasible approval")
        require(report["output_sha256"] is None, "missing output hash")
        require(not report["minimum_edits_proven"], "infeasible minimum claim")
    else:
        require(report["optimization_complete"] and report["minimum_edits_proven"], "complete tiny search")
        source = item["target"]["sequence"]
        require(report["output_sha256"] == "sha256:" + hashlib.sha256(expected.encode()).hexdigest(), "output hash")
        require(report["edits"] == [{"position_0based": i, "before": a, "after": b}
                                    for i, (a, b) in enumerate(zip(source, expected)) if a != b], "exact edits")
    start, end = item["cds"]["start_0based"], item["cds"]["end_0based_exclusive"]
    for policy, field in (("gc_content", "gc_content"), ("gc_window", "gc_window")):
        if policy not in item:
            require(field not in report, "unrequested GC is absent")
            continue
        facts = report[field]
        denominator = (end - start if policy == "gc_content" else item[policy]["window_bp"])
        permitted = [count for count in range(denominator + 1)
                     if item[policy]["min_basis_points"] * denominator <= count * 10000]
        maximums = [count for count in range(denominator + 1)
                    if count * 10000 <= item[policy]["max_basis_points"] * denominator]
        require(facts["minimum_gc_bases"] == min(permitted, default=denominator + 1), "minimum count")
        require(facts["maximum_gc_bases"] == max(maximums, default=-1), "maximum count")
        if policy == "gc_window":
            require(facts["window_bp"] == denominator, "window length")
        for side, dna in (("input", item["target"]["sequence"]), ("output", expected)):
            if dna is None:
                require(facts[side] is None, "missing output facts")
            elif policy == "gc_content":
                require(facts["denominator_bases"] == end - start, "CDS denominator")
                require(facts[side] == {"gc_bases": sum(b in "GC" for b in dna[start:end]),
                                        "satisfies_bounds": within_gc(dna[start:end], item[policy])}, "global facts")
            else:
                size = item[policy]["window_bp"]
                rows = [{"interval": {"start_0based": pos, "end_0based_exclusive": pos + size},
                         "gc_bases": sum(b in "GC" for b in dna[pos:pos + size]),
                         "satisfies_bounds": within_gc(dna[pos:pos + size], item[policy])}
                        for pos in range(start, end - size + 1)]
                require(facts["window_count"] == len(rows), "complete window inventory")
                require(facts[side] == rows, "window coordinates/counts")


# Wire ordering is explicit because serde serializes structs in declaration order,
# not alphabetically. This is an adversarial test encoder, never a production API.
REPORT_ORDER = "schema request source_sequence source_sha256 source_features_sha256 omitted_source_feature_count genetic_code_mapping genetic_code_mapping_sha256 algorithm status reason output_sequence output_sha256 edits initial_matches evaluated_candidates effective_evaluation_budget search_space optimization_complete minimum_edits_proven approval_digest nonclaims gc_content gc_window".split()
REQUEST_ORDER = "schema target purpose cds protein_sequence genetic_code protected_intervals avoid_motifs gc_content gc_window max_evaluations search_strategy output_seq_id".split()


def ordered(value, fields):
    if set(value) - set(fields):
        raise ValueError("Unknown audit wire fields")
    return {field: value[field] for field in fields if field in value}


def approval_digest(report):
    """Rehash adversarial reports honestly, including full declared constraints."""
    item = copy.deepcopy(report)
    item["approval_digest"] = None
    req = ordered(item["request"], REQUEST_ORDER)
    req["target"] = ordered(req["target"], ("kind", "sequence", "seq_id"))
    req["cds"] = ordered(req["cds"], ("start_0based", "end_0based_exclusive"))
    req["protected_intervals"] = [ordered(row, ("start_0based", "end_0based_exclusive"))
                                  for row in req["protected_intervals"]]
    req["avoid_motifs"] = [ordered(row, ("pattern", "strand")) for row in req["avoid_motifs"]]
    for name, fields in (("gc_content", ("min_basis_points", "max_basis_points")),
                         ("gc_window", ("window_bp", "min_basis_points", "max_basis_points"))):
        if name in req:
            req[name] = ordered(req[name], fields)
    item["request"] = req
    item["genetic_code_mapping"] = [ordered(row, ("codon", "residue")) for row in item["genetic_code_mapping"]]
    item["edits"] = [ordered(row, ("position_0based", "before", "after")) for row in item["edits"]]
    item["initial_matches"] = [ordered(row, ("motif_index", "interval", "strand")) for row in item["initial_matches"]]
    for row in item["initial_matches"]:
        row["interval"] = ordered(row["interval"], ("start_0based", "end_0based_exclusive"))
    if "gc_content" in item:
        facts = ordered(item["gc_content"], "denominator_bases minimum_gc_bases maximum_gc_bases reachable_minimum_gc_bases reachable_maximum_gc_bases input output".split())
        for side in ("input", "output"):
            if facts[side] is not None:
                facts[side] = ordered(facts[side], ("gc_bases", "satisfies_bounds"))
        item["gc_content"] = facts
    if "gc_window" in item:
        facts = ordered(item["gc_window"], "window_bp window_count minimum_gc_bases maximum_gc_bases input output".split())
        for side in ("input", "output"):
            if facts[side] is not None:
                facts[side] = [ordered(row, ("interval", "gc_bases", "satisfies_bounds")) for row in facts[side]]
                for row in facts[side]:
                    row["interval"] = ordered(row["interval"], ("start_0based", "end_0based_exclusive"))
        item["gc_window"] = facts
    raw = json.dumps(ordered(item, REPORT_ORDER), ensure_ascii=False, separators=(",", ":")).encode()
    return "sha256:" + hashlib.sha256(raw).hexdigest()
