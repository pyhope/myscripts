#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from pathlib import Path
import re

def parse_last_float(tokens):
    """Return the last float-parsable token from a list, or None if not found."""
    for tok in reversed(tokens):
        try:
            float(tok)
            return tok
        except ValueError:
            continue
    return None

def main():
    ap = argparse.ArgumentParser(description="Extract last numbers after marker blocks.")
    ap.add_argument("-i", "--input", default="OUTCAR", help="Input OUTCAR file (default: OUTCAR)")
    ap.add_argument("-o", "--output", default="mag", help="Output file prefix (default: mag)")
    ap.add_argument(
        "-n", "--index", type=int, nargs="+", default=[1,2],
        help="Leading integer index/indices to match after each marker, e.g. -n 2 5 7"
    )
    ap.add_argument(
        "-m", "--marker", default="magnetization (x)",
        help="Marker to search for (default: 'magnetization (x)')"
    )
    ap.add_argument(
        "-c", "--charge", action="store_true",
        help="If set, use marker 'total charge\\n' (overrides -m)."
    )
    ap.add_argument("-w", "--write", action="store_true", help="Write output to file (default: False)")
    args = ap.parse_args()

    infile = Path(args.input)

    # Marker selection
    if args.charge:
        marker = "total charge\n"
    else:
        marker = args.marker

    # Compile regexes for all requested indices
    indices = args.index[:]  # keep user order
    starts_with = {
        idx: re.compile(rf"^\s*{re.escape(str(idx))}(\s|$)")
        for idx in indices
    }

    results = []  # lines to write/print
    marker_idx = 0
    looking = False
    pending_vals = None  # dict idx->value (string token)

    def flush_pending():
        """Emit a result line if we have any marker pending."""
        nonlocal pending_vals, marker_idx
        if pending_vals is None:
            return
        # Build line in the user-specified order; missing indices -> leave空或用空白
        row = [str(marker_idx)]
        for idx in indices:
            row.append(pending_vals.get(idx, ""))
        results.append(" ".join(row))
        pending_vals = None

    def line_has_marker(line: str) -> bool:
        """Robust marker match. If marker contains a trailing newline, compare stripped forms too."""
        if marker in line:
            return True
        # handle the requested 'total charge\n' vs file's 'total charge'
        if marker.endswith("\n") and line.strip() == marker.strip():
            return True
        return False

    with infile.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            # If we hit a new marker and were still collecting for previous one, flush it first
            if line_has_marker(raw):
                if looking:
                    flush_pending()
                marker_idx += 1
                looking = True
                pending_vals = {}
                continue

            if looking:
                # Try match each still-missing index
                for idx in indices:
                    if idx in pending_vals:
                        continue
                    if starts_with[idx].match(raw):
                        tokens = raw.strip().split()
                        if not tokens:
                            continue
                        last_float_tok = parse_last_float(tokens)
                        if last_float_tok is None:
                            print(
                                f"Warning: could not parse float for marker {marker_idx}, index {idx}: {raw.strip()}",
                                file=sys.stderr
                            )
                            # leave it missing; continue scanning
                        else:
                            pending_vals[idx] = last_float_tok
                # If all gathered, we can flush now
                if pending_vals is not None and len(pending_vals) == len(indices):
                    flush_pending()
                    looking = False

    # End-of-file: if a marker block was open, flush what we have
    if looking:
        flush_pending()

    # Print to stdout
    for r in results:
        print(r)

    # Optional write to file
    if args.write:
        idx_tag = "-".join(str(i) for i in indices)
        outpath = Path(f"{args.output}-{idx_tag}.txt")
        with outpath.open("w", encoding="utf-8") as fw:
            for r in results:
                fw.write(r + "\n")
        print(f"# Wrote {len(results)} lines to {outpath}", file=sys.stderr)

if __name__ == "__main__":
    main()
