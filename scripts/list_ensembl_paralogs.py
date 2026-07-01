#!/usr/bin/env python3
"""List Ensembl same-species paralog gene IDs for RNA-read controls.

The script intentionally uses Ensembl's homology endpoint before any BLAST-style
search: curated same-species paralogues are better first-pass negative controls
than ad hoc sequence-similarity hits. Output defaults to one target Ensembl gene
ID per line so callers can pass the IDs directly to GENtle's Ensembl fixture
retrieval helper.
"""

from __future__ import annotations

import argparse
import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

DEFAULT_SERVER = "https://rest.ensembl.org"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="List Ensembl paralogues for a gene symbol or stable id.")
    parser.add_argument("query", help="Gene symbol or Ensembl gene stable ID, e.g. PATZ1 or ENSG...")
    parser.add_argument("--species", default="homo_sapiens", help="Source Ensembl species; default homo_sapiens.")
    parser.add_argument("--target-species", default="homo_sapiens", help="Target species; default homo_sapiens.")
    parser.add_argument("--server", default=DEFAULT_SERVER, help="Ensembl REST base URL.")
    parser.add_argument("--compara", default="vertebrates", help="Ensembl compara database; default vertebrates.")
    parser.add_argument("--from-json", type=Path, default=None, help="Read an existing homology JSON payload instead of contacting Ensembl.")
    parser.add_argument("--json-output", type=Path, default=None, help="Write raw Ensembl homology JSON payload for provenance.")
    parser.add_argument("--format", choices=["ids", "tsv", "json"], default="ids", help="Output format; default ids.")
    return parser.parse_args()


def is_stable_gene_id(query: str) -> bool:
    return query.upper().startswith("ENSG")


def fetch_payload(args: argparse.Namespace) -> dict[str, Any]:
    if args.from_json is not None:
        return json.loads(args.from_json.read_text(encoding="utf-8"))

    quoted_species = urllib.parse.quote(args.species.strip())
    quoted_query = urllib.parse.quote(args.query.strip())
    endpoint = "id" if is_stable_gene_id(args.query.strip()) else "symbol"
    path = f"/homology/{endpoint}/{quoted_species}/{quoted_query}"
    query = urllib.parse.urlencode(
        {
            "type": "paralogues",
            "target_species": args.target_species,
            "format": "full",
            "sequence": "none",
            "compara": args.compara,
            "content-type": "application/json",
        }
    )
    url = args.server.rstrip("/") + path + "?" + query
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(f"Ensembl homology request failed HTTP {exc.code}: {body[:1000]}") from exc
    except urllib.error.URLError as exc:
        raise SystemExit(f"Ensembl homology request failed: {exc}") from exc


def first_text(mapping: dict[str, Any], *keys: str) -> str:
    for key in keys:
        value = mapping.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def collect_paralogs(payload: dict[str, Any], source_query: str, target_species: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    source_query_lower = source_query.strip().lower()
    for data_row in payload.get("data", []) if isinstance(payload, dict) else []:
        if not isinstance(data_row, dict):
            continue
        source_id = first_text(data_row, "id", "gene_id", "stable_id")
        for homology in data_row.get("homologies", []) or []:
            if not isinstance(homology, dict):
                continue
            homology_type = first_text(homology, "type", "description")
            if homology_type and "paralog" not in homology_type.lower():
                continue
            target = homology.get("target")
            if not isinstance(target, dict):
                continue
            target_id = first_text(target, "id", "gene_id", "stable_id")
            if not target_id:
                continue
            target_species_value = first_text(target, "species", "taxon")
            if target_species and target_species_value and target_species_value != target_species:
                continue
            target_symbol = first_text(
                target,
                "display_id",
                "display_label",
                "external_gene_name",
                "gene_symbol",
                "name",
            )
            if target_id == source_id or target_id.lower() == source_query_lower or target_symbol.lower() == source_query_lower:
                continue
            if target_id in seen:
                continue
            seen.add(target_id)
            rows.append(
                {
                    "target_id": target_id,
                    "target_symbol": target_symbol,
                    "target_species": target_species_value,
                    "homology_type": homology_type,
                    "perc_id": first_text(target, "perc_id", "percentage_identity"),
                    "perc_pos": first_text(target, "perc_pos", "percentage_positive"),
                    "dn_ds": first_text(homology, "dn_ds"),
                }
            )
    rows.sort(key=lambda row: (row["target_symbol"].lower(), row["target_id"]))
    return rows


def main() -> int:
    args = parse_args()
    payload = fetch_payload(args)
    if args.json_output is not None:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    rows = collect_paralogs(payload, args.query, args.target_species)
    if args.format == "json":
        print(json.dumps({"query": args.query, "species": args.species, "target_species": args.target_species, "paralogs": rows}, indent=2))
    elif args.format == "tsv":
        print("target_id\ttarget_symbol\ttarget_species\thomology_type\tperc_id\tperc_pos\tdn_ds")
        for row in rows:
            print("\t".join(row[key] for key in ["target_id", "target_symbol", "target_species", "homology_type", "perc_id", "perc_pos", "dn_ds"]))
    else:
        for row in rows:
            print(row["target_id"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
