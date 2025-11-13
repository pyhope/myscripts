#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Task:
1) Validate formats of eqp1.dat and vhub.dat:
   - 1281 lines total each
   - First line is a header/meaningless; data start from line 2
   - 1280 data rows
   - Each data row has 4 columns separated by whitespace
     * col1: int in {1, 2}
     * col2: int in [1, 640]
     * col3: float
     * col4: float
2) For each row in eqp1.dat, find the row in vhub.dat with the same (col1, col2).
   - Keep eqp1.col3 unchanged
   - Replace eqp1.col4 with: eqp1.col4 - vhub.col3 (note: subtract vhub's 3rd column)
3) Write rows in the original order to eqp2.dat, copying eqp1.dat's first line as-is.

Usage:
    python process_eqp_vhub.py \
        --eqp1 eqp1.dat --vhub vhub.dat --out eqp2.dat
"""

import argparse
from typing import List, Tuple, Dict

def parse_args():
    ap = argparse.ArgumentParser(description="Process eqp1.dat and vhub.dat into eqp2.dat with validation.")
    ap.add_argument("--eqp1", '-e', default="eqp1.dat", help="Path to eqp1.dat (input).")
    ap.add_argument("--vhub", '-v', default="vhub.dat", help="Path to vhub.dat (input).")
    ap.add_argument("--out", '-o', default="eqp2.dat", help="Path to output eqp2.dat (output).")
    return ap.parse_args()

def read_and_validate(path: str) -> Tuple[str, List[Tuple[int, int, float, float]]]:
    """
    Returns:
        header_line (str), rows (list of tuples)
    Raises:
        ValueError upon any validation failure.
    """
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    if len(lines) != 1281:
        raise ValueError(f"{path}: expected 1281 lines, got {len(lines)}")

    header = lines[0].rstrip("\n")

    rows: List[Tuple[int, int, float, float]] = []
    for i, raw in enumerate(lines[1:], start=2):  # 1-based line number; data starts at line 2
        s = raw.strip()
        if not s:
            raise ValueError(f"{path}: empty line at #{i}")
        parts = s.split()
        if len(parts) != 4:
            raise ValueError(f"{path}: line #{i} should have 4 columns, got {len(parts)}: {parts}")

        try:
            c1 = int(parts[0])
            c2 = int(parts[1])
        except Exception:
            raise ValueError(f"{path}: line #{i} first two columns must be integers; got {parts[:2]}")

        if c1 not in (1, 2):
            raise ValueError(f"{path}: line #{i} col1 must be 1 or 2; got {c1}")
        if not (1 <= c2 <= 640):
            raise ValueError(f"{path}: line #{i} col2 must be in [1,640]; got {c2}")

        try:
            c3 = float(parts[2])
            c4 = float(parts[3])
        except Exception:
            raise ValueError(f"{path}: line #{i} col3/col4 must be floats; got {parts[2:4]}")

        rows.append((c1, c2, c3, c4))

    if len(rows) != 1280:
        raise ValueError(f"{path}: expected 1280 data rows, got {len(rows)}")

    return header, rows

def build_index(rows: List[Tuple[int, int, float, float]], source_name: str) -> Dict[Tuple[int,int], Tuple[float,float]]:
    """
    Build dict: (col1, col2) -> (col3, col4)
    """
    idx: Dict[Tuple[int,int], Tuple[float,float]] = {}
    for j, (c1, c2, c3, c4) in enumerate(rows, start=1):
        key = (c1, c2)
        if key in idx:
            raise ValueError(f"{source_name}: duplicate key (col1,col2)={key} at row #{j}")
        idx[key] = (c3, c4)
    # Optional sanity check: should contain exactly 2*640 keys
    if len(idx) != 1280:
        raise ValueError(f"{source_name}: index size {len(idx)} != 1280 (expected 1280 unique (col1,col2) pairs)")
    return idx

def main():
    args = parse_args()

    # Read and validate both files
    header_eqp1, rows_eqp1 = read_and_validate(args.eqp1)
    header_vhub, rows_vhub = read_and_validate(args.vhub)

    # Build (col1,col2) index for vhub
    vhub_index = build_index(rows_vhub, args.vhub)

    # Compute new rows for eqp2: col3 unchanged, col4 := eqp1.col4 - vhub.col3 matched by key
    out_lines: List[str] = []
    out_lines.append(header_eqp1 + "\n")  # copy header as-is

    SEP = "  "
    for k, (c1, c2, c3_eqp, c4_eqp) in enumerate(rows_eqp1, start=2):  # line numbers if needed
        key = (c1, c2)
        if key not in vhub_index:
            raise ValueError(f"No matching (col1,col2)={key} found in {args.vhub} for eqp1 line #{k}")
        vhub_c3, _vhub_c4 = vhub_index[key]

        new_c3 = c3_eqp
        new_c4 = c4_eqp - vhub_c3  # the required operation
        out_lines.append(f"{c1:8d}{c2:8d}{SEP}{new_c3:13.9f}{SEP}{new_c4:13.9f}\n")

    with open(args.out, "w", encoding="utf-8") as f:
        f.writelines(out_lines)

    print(f"Done. Wrote {len(out_lines)} lines to {args.out} (including header).")

if __name__ == "__main__":
    main()
