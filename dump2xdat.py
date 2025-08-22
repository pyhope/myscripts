#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from collections import defaultdict

def parse_box_bounds(header_line, lines3):
    """
    Parse BOX BOUNDS lines. Support both triclinic ('xy xz yz ...') and orthogonal.
    Return (lx, ly, lz, xy, xz, yz).
    """
    tokens = header_line.strip().split()
    triclinic = ("xy" in tokens and "xz" in tokens and "yz" in tokens)

    def parse_three(vals):
        parts = vals.strip().split()
        if len(parts) < 2:
            raise ValueError("BOX BOUNDS line has <2 values.")
        if len(parts) == 2:
            lo, hi = float(parts[0]), float(parts[1])
            tilt = 0.0
        else:
            lo, hi, tilt = float(parts[0]), float(parts[1]), float(parts[2])
        return lo, hi, tilt

    (xlo_b, xhi_b, xy) = parse_three(lines3[0])
    (ylo_b, yhi_b, xz) = parse_three(lines3[1])
    (zlo_b, zhi_b, yz) = parse_three(lines3[2])

    # Recover lx, ly, lz from bounding-box + tilts (LAMMPS Howto_triclinic)
    # xlo = ?, xhi = ?  but lx = (xhi_b - xlo_b) - (max - min) over {0,xy,xz,xy+xz}
    if triclinic:
        max_x = max(0.0, xy, xz, xy + xz)
        min_x = min(0.0, xy, xz, xy + xz)
        lx = (xhi_b - xlo_b) - (max_x - min_x)
        ly = (yhi_b - ylo_b) - abs(yz)
        lz = (zhi_b - zlo_b)
    else:
        # orthogonal box
        lx = (xhi_b - xlo_b)
        ly = (yhi_b - ylo_b)
        lz = (zhi_b - zlo_b)
        xy = xz = yz = 0.0

    # sanity checks
    if lx <= 0 or ly <= 0 or lz <= 0:
        raise ValueError("Non-positive box length recovered. Check dump format.")

    return lx, ly, lz, xy, xz, yz


def cart_to_frac_triclinic(x, y, z, lx, ly, lz, xy, xz, yz):
    """
    Inverse mapping of the LAMMPS upper-triangular cell:
      x = s*lx + t*xy + u*xz
      y = t*ly + u*yz
      z = u*lz
    Solve:
      u = z/lz
      t = (y - u*yz)/ly
      s = (x - t*xy - u*xz)/lx
    """
    u = z / lz
    t = (y - u * yz) / ly
    s = (x - t * xy - u * xz) / lx
    return s, t, u


def read_dump(filename):
    """
    Read a LAMMPS dump (id type x y z) with BOX BOUNDS header.
    Return:
      meta: dict with lx,ly,lz,xy,xz,yz, natoms, types_set (sorted), first_frame_order (type,id order)
      frames: list of per-frame lists of tuples (id, type, x, y, z)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    n = len(lines)
    frames = []
    box_params = None
    natoms_expected = None
    columns = None

    while i < n:
        if not lines[i].startswith("ITEM: TIMESTEP"):
            i += 1
            continue
        # TIMESTEP
        if i + 1 >= n:
            break
        timestep = int(lines[i+1].strip())
        i += 2

        # NUMBER OF ATOMS
        if i >= n or not lines[i].startswith("ITEM: NUMBER OF ATOMS"):
            raise ValueError("Expected 'ITEM: NUMBER OF ATOMS'")
        if i + 1 >= n:
            raise ValueError("Missing natoms value")
        natoms = int(lines[i+1].strip())
        if natoms_expected is None:
            natoms_expected = natoms
        elif natoms != natoms_expected:
            raise ValueError("natoms changed across frames, not supported.")
        i += 2

        # BOX BOUNDS
        if i >= n or not lines[i].startswith("ITEM: BOX BOUNDS"):
            raise ValueError("Expected 'ITEM: BOX BOUNDS ...'")
        box_header = lines[i]
        if i + 3 >= n:
            raise ValueError("BOX BOUNDS lines incomplete")
        box_lines = [lines[i+1], lines[i+2], lines[i+3]]
        lx, ly, lz, xy, xz, yz = parse_box_bounds(box_header, box_lines)
        if box_params is None:
            box_params = (lx, ly, lz, xy, xz, yz)
        else:
            # assert fixed box (NVT)
            EPS = 1e-8
            if any(abs(a-b) > 1e-6 for a, b in zip(box_params, (lx, ly, lz, xy, xz, yz))):
                raise ValueError("Box changed across frames; only fixed-shape (NVT) supported.")
        i += 4

        # ATOMS header
        if i >= n or not lines[i].startswith("ITEM: ATOMS"):
            raise ValueError("Expected 'ITEM: ATOMS ...'")
        atom_head = lines[i].strip().split()[2:]  # after 'ITEM:' 'ATOMS'
        # Build column map; require id,type,x,y,z exist
        col_idx = {}
        for idx, name in enumerate(atom_head):
            col_idx[name] = idx
        required = ["id", "type", "x", "y", "z"]
        for r in required:
            if r not in col_idx:
                raise ValueError(f"'ITEM: ATOMS' must contain column '{r}'.")
        columns = col_idx
        i += 1

        # Read natoms lines
        frame = []
        for _ in range(natoms):
            if i >= n:
                raise ValueError("Unexpected EOF while reading atoms.")
            parts = lines[i].split()
            # tolerate extra columns
            try:
                aid = int(parts[columns["id"]])
                typ = int(parts[columns["type"]])
                x = float(parts[columns["x"]])
                y = float(parts[columns["y"]])
                z = float(parts[columns["z"]])
            except Exception:
                raise ValueError(f"Malformed ATOMS line: {lines[i].rstrip()}")
            frame.append((aid, typ, x, y, z))
            i += 1

        frames.append(frame)

    if not frames:
        raise ValueError("No frames parsed from dump.")

    # types present
    types_set = sorted({typ for (aid, typ, *_rest) in frames[0]})
    return {
        "natoms": natoms_expected,
        "box": box_params,           # (lx, ly, lz, xy, xz, yz)
        "types_set": types_set,
    }, frames


def write_xdatcar(output_filename, title, species, box, frames, wrap_frac=True):
    """
    Write frames (list of per-frame atom tuples (id, typ, x,y,z)) to XDATCAR.
    - species: list like ["Fe","Mg","Si","O"], corresponds to LAMMPS types 1..N.
    - Atoms are ordered by species (type) ascending and then by id ascending in every frame.
    - Coordinates are converted to fractional (Direct).
    """
    lx, ly, lz, xy, xz, yz = box
    natoms = len(frames[0])

    # count per type from first frame
    counts = defaultdict(int)
    for (aid, typ, *_xyz) in frames[0]:
        counts[typ] += 1
    # Build counts in species order
    counts_list = []
    for t_id in range(1, len(species)+1):
        cnt = counts.get(t_id, 0)
        if cnt == 0:
            raise ValueError(f"No atoms of type {t_id} found, but '--types' includes {len(species)} entries.")
        counts_list.append(cnt)

    # prepare constant header
    with open(output_filename, 'w') as f:
        # line 1: title
        f.write(f"{title}\n")
        # line 2: scale
        f.write("           1\n")  # use 1.0 scale
        # lines 3-5: lattice vectors (upper-triangular form)
        f.write(f" {lx:12.6f} {0.0:12.6f} {0.0:12.6f}\n")
        f.write(f" {xy:12.6f} {ly:12.6f} {0.0:12.6f}\n")
        f.write(f" {xz:12.6f} {yz:12.6f} {lz:12.6f}\n")
        # line 6: species
        f.write("  " + "   ".join(species) + " \n")
        # line 7: counts
        f.write("  " + "   ".join(f"{c:4d}" for c in counts_list) + "\n")

        # frames: re-order by species then id; convert to fractional
        for i_frame, frame in enumerate(frames, start=1):
            # build per-type list
            by_type = defaultdict(list)
            for (aid, typ, x, y, z) in frame:
                by_type[typ].append((aid, x, y, z))
            # sort within each type by id
            ordered_atoms = []
            for t_id in range(1, len(species)+1):
                if t_id not in by_type:
                    raise ValueError(f"Type {t_id} missing in a later frame.")
                ordered_atoms.extend(sorted(by_type[t_id], key=lambda r: r[0]))  # (id,x,y,z)

            # write frame label
            f.write(f"Direct configuration=     {i_frame}\n")

            # coords
            for (_aid, x, y, z) in ordered_atoms:
                s, t, u = cart_to_frac_triclinic(x, y, z, lx, ly, lz, xy, xz, yz)
                if wrap_frac:
                    # wrap to [0,1)
                    s = s - float(int(s))
                    t = t - float(int(t))
                    u = u - float(int(u))
                f.write(f"  {s: .9f}  {t: .9f}  {u: .9f}\n")


def main():
    ap = argparse.ArgumentParser(description="Convert LAMMPS dump (id type x y z) to VASP XDATCAR (NVT, fixed box).")
    ap.add_argument("-i", "--input_dump", type=str, default="nvt.dump", help="Input LAMMPS dump file")
    ap.add_argument("-o", "--output_xdatcar", type=str, default="XDATCAR", help="Output XDATCAR filename")
    ap.add_argument("-t", "--types", type=str, required=True,
                    help='Space-separated element symbols mapping to LAMMPS types 1..N, e.g. --types "Fe Mg Si O"')
    ap.add_argument("--title", type=str, default="unknown system", help="Title line for XDATCAR")
    ap.add_argument("-s", "--start_step", type=int, default=1, help="Start frame index (1-based, inclusive) in file order")
    ap.add_argument("-e", "--end_step", type=int, default=10**12, help="End frame index (1-based, inclusive)")
    ap.add_argument("--no_wrap", action="store_true", help="Do not wrap fractional coords into [0,1)")
    args = ap.parse_args()

    try:
        meta, frames = read_dump(args.input_dump)
    except Exception as exc:
        print(f"[ERROR] Failed to read dump: {exc}", file=sys.stderr)
        sys.exit(1)

    species = args.types.strip().split()
    n_types_seen = len(meta["types_set"])
    if len(species) != n_types_seen:
        print(f"[ERROR] '--types' count ({len(species)}) != number of LAMMPS types seen ({n_types_seen}).",
              file=sys.stderr)
        print(f"       Types seen in dump: {meta['types_set']}", file=sys.stderr)
        sys.exit(2)

    # clip frames in file order using 1-based indices
    total = len(frames)
    start_1b = max(1, args.start_step)
    end_1b = min(args.end_step, total)
    if start_1b > end_1b:
        print("[ERROR] Invalid clip range.", file=sys.stderr)
        sys.exit(3)
    frames_clip = frames[start_1b-1:end_1b]

    try:
        write_xdatcar(
            args.output_xdatcar,
            title=args.title,
            species=species,
            box=meta["box"],
            frames=frames_clip,
            wrap_frac=(not args.no_wrap),
        )
    except Exception as exc:
        print(f"[ERROR] Writing XDATCAR failed: {exc}", file=sys.stderr)
        sys.exit(4)

    print(f"[OK] Wrote XDATCAR '{args.output_xdatcar}' with frames {start_1b}..{end_1b} "
          f"(total {end_1b - start_1b + 1}), natoms={meta['natoms']}.")


if __name__ == "__main__":
    main()
