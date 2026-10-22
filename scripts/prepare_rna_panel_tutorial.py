"""Optional terminal JSON convenience for tutorial 10.04; no read bytes are opened.

The inner agent instead constructs the same typed request for panel-screen.
Synthetic data preparation and all scientific decisions belong to GENtle.
"""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["request"])
    parser.add_argument("--panel", required=True)
    parser.add_argument("--md5-manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--input", nargs=4, metavar=("SAMPLE", "BARCODE", "PRIMING", "FILE"), action="append", required=True)
    parser.add_argument("--declared-complete", action="store_true", required=True,
                        help="Attest copying each explicitly named file has finished; GENtle still verifies checksums at EOF")
    parser.add_argument("--write-locus-fastas", action="store_true")
    parser.add_argument("--cancel-path")
    args = parser.parse_args()
    files = [{"sample_id": sample, "barcode": barcode, "priming": priming,
              "path": str(Path(path).resolve()), "declared_complete": args.declared_complete}
             for sample, barcode, priming, path in args.input]
    request = {"schema": "gentle.rna_read_panel_screen_request.v1",
               "panel_path": str(Path(args.panel).resolve()), "files": files,
               "md5_manifest_path": str(Path(args.md5_manifest).resolve()),
               "output_dir": str(Path(args.output_dir).resolve()),
               "policy": {"write_locus_fastas": args.write_locus_fastas},
               "cancel_path": str(Path(args.cancel_path).resolve()) if args.cancel_path else None}
    with Path(args.output).open("x", encoding="utf-8", newline="\n") as out:
        json.dump(request, out, indent=2)
        out.write("\n")


if __name__ == "__main__":
    main()
