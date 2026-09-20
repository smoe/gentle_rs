#!/usr/bin/env python3
"""Recreate a tiny synthetic producer-layout regulatory TFBS package offline.

Uses only Python's standard library and the DuckDB CLI. All coordinates,
scores, owners and annotations are hand-crafted, not biological measurements.
The layout follows jaspar-mapping manage_regulatory_tfbs.py at 77acb384.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def quote(value):
    return "'" + str(value).replace("'", "''") + "'"


def build(root, duckdb):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=False)
    annotation = root / "annotation"
    annotation.mkdir()

    def save(path, value):
        path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")

    def sql(query, database=":memory:"):
        result = subprocess.run([duckdb, "-no-init", "-batch", "-bail", "-json", str(database), "-c", query], capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(result.stderr)
        return json.loads(result.stdout or "[]")

    def record(path, base):
        return dict(path=path.relative_to(base).as_posix(), bytes=path.stat().st_size,
                    sha256=hashlib.sha256(path.read_bytes()).hexdigest())

    definition = "tss_upstream_700_downstream_300_v1"
    promoter_sql = f"""SELECT 'fixture-genome' AS genome_id, 'fixture-gtf' AS annotation_release,
        {quote(definition)} AS promoter_definition_id,700::BIGINT AS upstream_bp,300::BIGINT AS downstream_bp,
        *,tss_start+1 AS tss_end FROM (VALUES
        ('p1','tss1','1','+',300::BIGINT,1301::BIGINT,1000::BIGINT),
        ('p2','tss2','1','-',1100,2101,1400),
        ('p3','tss3','2','+',300,1301,1000),
        ('p4','tss4','MT','+',300,1301,1000))
        p(promoter_id,tss_id,chrom,strand,promoter_start,promoter_end,tss_start)"""
    feature_sql = """SELECT 'fixture-assembly' AS assembly,*,NULL::BIGINT AS extended_start,NULL::BIGINT AS extended_end,
        NULL::BIGINT AS core_start,NULL::BIGINT AS core_end FROM (VALUES
        ('reg1','1','enhancer',305::BIGINT,350::BIGINT),
        ('reg2','1','open_chromatin_region',1200,1220),
        ('reg3','1','enhancer',1600,1620),
        ('reg4','2','enhancer',305,350)) f(regulatory_feature_id,chrom,feature_type,start,"end")"""
    owners_sql = """SELECT 'fixture-genome' AS genome_id,'fixture-gtf' AS annotation_release,* FROM
        (VALUES ('tss1','tx1','geneA','TOY'),('tss1','tx2','geneA','TOY'),
        ('tss1','tx3','geneB','NEIGHBOR'),('tss2','tx4','geneA','TOY'),
        ('tss3','tx5','geneC','EMPTY'),('tss4','tx6','geneM','MITO'))
        o(tss_id,transcript_id,gene_id,gene_name)"""

    def copy(query, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        sql(f"COPY ({query}) TO {quote(path)} (FORMAT PARQUET, COMPRESSION ZSTD);")

    copy(promoter_sql, annotation / "promoter.parquet")
    copy(owners_sql, annotation / "transcript_tss.parquet")
    copy(f'SELECT genome_id,annotation_release,tss_id,chrom,strand,tss_start AS start,tss_end AS "end" FROM ({promoter_sql})', annotation / "transcription_start_site.parquet")
    copy(feature_sql, annotation / "regulatory_feature.parquet")
    copy("SELECT 'reg1' AS regulatory_feature_id,'geneB' AS gene_id", annotation / "regulatory_feature_gene.parquet")
    coverage = []
    for chrom in ("1", "2", "MT"):
        f = annotation / "chromosomes" / chrom / "features.parquet"
        p = annotation / "chromosomes" / chrom / "promoters.parquet"
        copy(f"SELECT * FROM ({feature_sql}) WHERE chrom={quote(chrom)}", f)
        copy(f"SELECT * FROM ({promoter_sql}) WHERE chrom={quote(chrom)}", p)
        coverage.append(dict(chrom=chrom, length=3000, feature_count={"1": 3, "2": 1, "MT": 0}[chrom],
                             promoter_count=2 if chrom == "1" else 1,
                             coverage="known_empty_intersection" if chrom == "MT" else "annotated_intersection",
                             features=record(f, annotation), promoters=record(p, annotation)))
    am = dict(schema_version=1, kind="regulatory_tfbs_annotation", state="complete",
              genome_id="fixture-genome", assembly="fixture-assembly", annotation_release="fixture-gtf",
              regulatory_release="fixture-regulatory", promoter_definition_id=definition, upstream_bp=700, downstream_bp=300,
              coordinate_rule="BED_half_open_offsets_include_TSS_base_clamped_to_sequence", regulatory_gff_sha256="a" * 64,
              coordinate_audit=dict(status="passed", assembly="fixture-assembly", gff_sha256="a" * 64, features_compared=4,
                                    note="Synthetic expected geometry only; not a real Ensembl audit"),
              chromosomes=coverage, inputs=[], files=[record(p, annotation) for p in sorted(annotation.rglob("*.parquet"))])
    save(annotation / "manifest.json", am)

    inventory = []
    source_inventory = []
    for chrom in ("1", "2", "MT"):
        for motif, length, name, floor in [("MA0861.2", 16, "TP73", -5), ("MA1961.2", 11, "PATZ1", -1)]:
            for strand in ("+", "-"):
                source_inventory.append(dict(run_id="synthetic-source-run", genome_id="fixture-genome", motif_set_id="synthetic-motifs",
                    chrom=chrom,motif_id=motif,strand=strand,state="complete",coordinate_mode="bed",minimum_score=floor,
                    minimum_pwm_relative_score=None,maximum_pwm_relative_score=None,emitted_hits=999999,
                    output_relative_path=f"/unavailable/original-atlas/{chrom}/{motif}/{strand}.parquet",sha256="b"*64,
                    score_mode="log2_relative_risk",pseudocount=1.0,pseudocount_scheme="additive_per_base",
                    background_model_id="uniform_acgt_v1",n_policy="skip"))
            if chrom == "MT":
                inventory.append(dict(chrom=chrom,motif_id=motif,rows=0,bytes=0,path=None,sha256=None,state="known_empty_intersection"))
                continue
            data = [(295, "+", -4.0, 132), (295, "-", -3.5, 132), (1200, "+", 0.0, 136), (1600, "-", 6.0, 132)] if motif == "MA0861.2" else [(320, "-", -0.5, 132), (320, "+", 1.25, 132), (1200, "+", 0.0, 136)]
            rows = ",".join(f"({s}::BIGINT,{quote(strand)},{score}::FLOAT,{tags}::USMALLINT)" for s,strand,score,tags in data)
            query = f"""SELECT 'fixture-genome' AS genome_id,'synthetic-motifs' AS motif_set_id,
                {quote(chrom)} AS chrom,start,start+{length} AS "end",{quote(motif)} AS motif_id,{quote(name)} AS motif_name,
                strand,score,0.5::FLOAT AS pwm_relative_score,'log2_relative_risk' AS score_mode,
                1.0::DOUBLE AS pseudocount,'uniform_acgt_v1' AS background_model_id,'additive_per_base' AS pseudocount_scheme,
                {floor}::DOUBLE AS source_minimum_score,'skip' AS n_policy,regulation_tags,
                true AS overlaps_regulatory_tss_intersection,
                false AS overlaps_promoter_core,false AS overlaps_promoter_extended,
                (regulation_tags & 4)<>0 AS overlaps_enhancer,(regulation_tags & 8)<>0 AS overlaps_open_chromatin,
                false AS overlaps_ctcf,false AS overlaps_emar,false AS overlaps_other_regulatory,
                true AS overlaps_tss_window,true AS overlaps_regulatory,true AS regulation_candidate
                FROM (VALUES {rows}) h(start,strand,score,regulation_tags) {"WHERE false" if chrom == "2" else ""}"""
            path = root / "hits" / f"chrom={chrom}" / f"{motif}.parquet"
            copy(query, path)
            inventory.append(dict(chrom=chrom,motif_id=motif,rows=len(data) if chrom == "1" else 0,
                                  state="complete",**record(path,root)))
    save(root / "file_inventory.json", inventory)
    save(root / "source-inventory.tmp.json", source_inventory)
    db = root / "regulatory_tfbs.duckdb"
    sql(f"""CREATE TABLE genome AS SELECT 'fixture-genome' AS genome_id,'fixture-assembly' AS assembly_name,
        'fixture-accession.1' AS assembly_accession,'source-fasta-release' AS ensembl_release;
        CREATE TABLE motif_metadata AS SELECT 'synthetic-motifs' AS motif_set_id,* FROM
          (VALUES ('MA0861.2','TP73',16,'fixture-human'),('MA1961.2','PATZ1',11,'fixture-other-taxon'))
          m(motif_id,motif_name,motif_length,species);
        CREATE TABLE sequence_region AS SELECT 'fixture-genome' AS genome_id,chrom,3000::BIGINT AS length,true AS included_in_scan
          FROM (VALUES ('1'),('2'),('MT')) c(chrom);
        CREATE TABLE scan_file_inventory AS SELECT * FROM read_json_auto({quote(root / 'source-inventory.tmp.json')});
        CREATE TABLE file_inventory AS SELECT * FROM read_json_auto({quote(root / 'file_inventory.json')});
        CREATE MACRO motif_hits(files) AS TABLE SELECT * FROM read_parquet(files,hive_partitioning=false);
        """, db)
    (root / "source-inventory.tmp.json").unlink()
    (root / "schema.sql").write_text("-- Query aid only; GENtle must not execute this file.\n", encoding="utf-8", newline="\n")
    manifest = dict(schema_version=1,kind="genome_regulatory_tfbs_subset",state="complete",genome_id=am["genome_id"],assembly=am["assembly"],
        annotation_release=am["annotation_release"],regulatory_release=am["regulatory_release"],promoter_definition_id=definition,
        production_plan_sha256="c"*64,source_commit="synthetic-fixture-not-a-production-commit",database=db.name,
        coordinate_mode="bed_0based_half_open",scope="regulatory_and_tss",score_selection="source_retention",upstream_bp=700,downstream_bp=300,
        chromosomes=["1","2","MT"],motif_count=2,requires_tp73=False,rows=sum(r["rows"] for r in inventory),
        parquet_bytes=sum(r["bytes"] for r in inventory),complete_for_declared_annotation_and_source_floors=True,complete_genome_scan=False,
        files=[record(root / p, root) for p in ("file_inventory.json",db.name,"schema.sql","annotation/manifest.json")])
    save(root / "manifest.json", manifest)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--duckdb", default="duckdb")
    args = parser.parse_args()
    build(args.output, args.duckdb)
