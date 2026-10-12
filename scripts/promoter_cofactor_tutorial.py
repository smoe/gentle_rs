"""Offline synthetic promoter-cofactor tutorial preparation and shared-CLI replay.

Only repository-authored fixture.sql is executed to create temporary Parquet.
Queries and capture validation belong to GENtle, not this teaching harness.
No downloads, production package or external publication. The tiny synthetic
reference is prepared locally to verify explicit DNA-feature attachment.
"""

import argparse
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "test_files/fixtures/promoter_cofactors"
TABLES = {"anchors": "anchors", "feature": "feature_MA9000.1",
          "anchor_promoter": "anchor_promoter", "promoter": "promoter",
          "promoter_gene": "promoter_gene",
          "cofactor_distance_isoform_comparison": "cofactor_distance_isoform_comparison"}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def inventory(path):
    raw = path.read_bytes()
    return {"path": path.name, "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def prepare(package, duckdb):
    package.mkdir(parents=True, exist_ok=False)
    sql = (FIXTURE / "fixture.sql").read_text(encoding="utf-8")
    for table, name in TABLES.items():
        path = str(package / f"{name}.parquet").replace("'", "''")
        sql += f"\nCOPY {table} TO '{path}' (FORMAT PARQUET);"
    subprocess.run([str(duckdb), "-no-init", ":memory:", "-c", sql],
                   check=True, capture_output=True, text=True, encoding="utf-8", timeout=120)
    manifest = json.loads((FIXTURE / "manifest.template.json").read_text(encoding="utf-8"))
    manifest["files"] = [inventory(package / f"{name}.parquet") for name in TABLES.values()]
    write_json(package / "manifest.json", manifest)
    write_json(package / "complete.json", {"identity": manifest["identity"],
               "files": manifest["files"] + [inventory(package / "manifest.json")]})


def replay(output, cli, duckdb):
    output.mkdir(parents=True, exist_ok=False)
    package = output / "package"
    prepare(package, duckdb)
    commands = []
    state = str(output / "tutorial.state.json")

    def run(*args):
        command = [str(cli), "--state", state, *args]
        commands.append(command)
        proc = subprocess.run(command, cwd=ROOT, check=True, capture_output=True,
                              text=True, encoding="utf-8", timeout=180)
        return json.loads(proc.stdout)

    def query(name, **fields):
        request = dict(package_path=str(package), assembly="GRCh38", timeout_seconds=120,
                       duckdb_executable=str(duckdb), **fields)
        path = output / f"{name}.request.json"
        write_json(path, request)
        result = run("features", "promoter-cofactors", "@" + str(path))
        report = result["report"]
        if report["availability"] != "available":
            raise RuntimeError(f"{name}: {report['availability']}: {report.get('diagnostic')}")
        write_json(output / f"{name}.report.json", report)
        return report

    query("inspect", query="inspect")
    rankings = query("rankings", query="rankings", distance_band="gap_6_20")
    assert rankings["rankings"][0]["positive_anchor_fraction"] == 0.5
    anchors = query("anchors", query="anchors", region=dict(chromosome="1", start_0based=100, end_0based_exclusive=116))
    anchor = next(a for a in anchors["anchors"] if (a["chrom"], a["anchor_start"], a["anchor_end"]) == ("1", 100, 116))
    detail = query("detail", query="anchor_detail", anchor_id=anchor["anchor_id"], motif="MA9000.1")
    hit = next(d for d in detail["details"] if d["distance_band"] == "gap_6_20")
    assert (hit["hit_start"], hit["hit_end"], hit["best_score"], hit["interval_distance_bp"]) == (133, 148, 4.25, 17)
    higher = query("higher_threshold", query="anchor_detail", anchor_id=anchor["anchor_id"], motif="MA9000.1", presence_threshold=5)
    high_hit = next(d for d in higher["details"] if d["distance_band"] == "gap_6_20")
    assert not high_hit["present_at_requested_threshold"]
    assert high_hit["n_source_loci"] == hit["n_source_loci"]
    assert high_hit["n_score_zero_loci"] == hit["n_score_zero_loci"]
    capture = dict(set_id="promoter_cofactors", label="Synthetic retained cofactor hit",
                   source=dict(source_kind="promoter_cofactor", report=detail,
                               target=dict(kind="hit", anchor_id=anchor["anchor_id"],
                                           motif_id="MA9000.1", distance_band="gap_6_20")))
    write_json(output / "capture.request.json", capture)
    run("shell", "regions capture " + shlex.quote("@" + str(output / "capture.request.json")))
    export = dict(set_id="promoter_cofactors", json_path=str(output / "regions.json"))
    write_json(output / "export.request.json", export)
    run("shell", "regions export " + shlex.quote("@" + str(output / "export.request.json")))
    regions = json.loads((output / "regions.json").read_text(encoding="utf-8"))
    saved = regions["regions"][0]
    assert (saved["interval"]["start_0based"], saved["interval"]["end_0based_exclusive"]) == (133, 148)
    assert saved["evidence"][0]["source_record"] == detail
    write_json(output / "inspect_region.request.json", dict(set_id="promoter_cofactors", region_id=saved["region_id"]))
    run("shell", "regions inspect " + shlex.quote("@" + str(output / "inspect_region.request.json")))
    # A non-palindromic synthetic reference tests orientation, not motif biology.
    reference = "AACCGT" * 40
    (output / "toy.fa").write_text(">1\n" + reference + "\n", encoding="utf-8")
    (output / "toy.gtf").write_text('1\tsynthetic\tgene\t1\t240\t.\t+\t.\tgene_id "GENE-A"; gene_name "GENE-A";\n', encoding="utf-8")
    catalog = output / "catalog.json"
    write_json(catalog, {"Cofactor synthetic reference": dict(
        ncbi_taxonomy_id=9606, ncbi_assembly_name="GRCh38", ncbi_assembly_accession="GCA_000000000.1", sequence_local="toy.fa",
        annotations_local="toy.gtf", cache_dir=str(output / "cache"))})
    run("genomes", "prepare", "Cofactor synthetic reference", "--catalog", str(catalog))
    for seq_id, expected_start, expected_end in [("demo_plus", 33, 48)]:
        run("genomes", "extract-region", "Cofactor synthetic reference", "1", "101", "200",
            "--output-id", seq_id, "--annotation-scope", "none", "--catalog", str(catalog))
        feature_request = dict(set_id="promoter_cofactors", region_id=saved["region_id"],
                               seq_id=seq_id, expected_region_content_sha256=saved["content_sha256"])
        preview_path = output / f"{seq_id}.preview.request.json"
        write_json(preview_path, feature_request)
        preview = run("shell", "regions preview-feature " + shlex.quote("@" + str(preview_path)))
        write_json(output / f"{seq_id}.preview.json", preview)
        feature = preview["feature_materialization"]
        projection = feature["projection"]
        assert (projection["local_start_0based"], projection["local_end_0based_exclusive"]) == (expected_start, expected_end)
        assert projection["local_strand"] == "plus"
        assert not feature["curation"]["applied"]
        shutil.copyfile(state, output / f"{seq_id}.preview.state.json")
        feature_request["expected_approval_sha256"] = feature["approval_sha256"]
        apply_path = output / f"{seq_id}.apply.request.json"
        write_json(apply_path, feature_request)
        applied = run("shell", "regions materialize-feature " + shlex.quote("@" + str(apply_path)))
        assert applied["feature_materialization"]["curation"]["applied"]
        write_json(output / f"{seq_id}.applied.json", applied)
    write_json(output / "replay.json", {"status": "passed", "synthetic": True,
               "commands": commands, "feature_sequences": ["demo_plus"],
               "non_claim": "Offline CLI replay, not live GUI or scientific acceptance"})
    return saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New directory only")
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--duckdb", type=Path, required=True)
    args = parser.parse_args()
    replay(args.output.resolve(), args.cli.resolve(strict=True), args.duckdb.resolve(strict=True))
    print(f"Synthetic CLI replay passed: {args.output.resolve()}")


if __name__ == "__main__":
    main()
