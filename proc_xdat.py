#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys

def read_xdatcar(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    if len(lines) < 7:
        raise ValueError("XDATCAR file seems too short (<7 lines).")

    header = [ln if ln.endswith('\n') else ln + '\n' for ln in lines[:7]]
    title = header[0].rstrip('\n')
    scale_raw = header[1].strip()
    try:
        scale = float(scale_raw)
    except Exception:
        raise ValueError(f"Invalid scale factor on line 2: {scale_raw!r}")

    # lattice vectors
    a = list(map(float, header[2].split()))
    b = list(map(float, header[3].split()))
    c = list(map(float, header[4].split()))

    elements = header[5].split()
    counts = list(map(int, header[6].split()))
    natoms = sum(counts)

    if scale > 0:
        a = [scale * x for x in a]
        b = [scale * x for x in b]
        c = [scale * x for x in c]
    # scale <= 0: assume Angstrom already

    coords = []
    i = 7
    while i < len(lines):
        if lines[i].lstrip().startswith('Direct configuration'):
            i += 1
            frame = []
            for _ in range(natoms):
                if i >= len(lines):
                    raise ValueError("Unexpected end of file while reading coordinates.")
                parts = lines[i].split()
                if len(parts) < 3:
                    raise ValueError(f"Malformed coordinate line at line {i+1}: '{lines[i].rstrip()}'")
                s, t, u = map(float, parts[:3])
                frame.append([s, t, u])
                i += 1
            coords.append(frame)
        else:
            i += 1

    return {
        "header": header,
        "title": title,
        "a": a, "b": b, "c": c,
        "elements": elements, "counts": counts,
        "natoms": natoms,
        "frames_frac": coords
    }

def build_types(elements, counts):
    types = []
    for t_id, (_, cnt) in enumerate(zip(elements, counts), start=1):
        types.extend([t_id] * cnt)
    return types

def write_xdatcar_clip(output_filename, header, frames_frac, start_1b, end_1b):
    if start_1b < 1 or end_1b > len(frames_frac) or start_1b > end_1b:
        raise ValueError("Invalid clipping range.")

    with open(output_filename, 'w') as f:
        for line in header:
            f.write(line)
        new_counter = 0
        for idx in range(start_1b - 1, end_1b):
            new_counter += 1
            f.write(f"Direct configuration=     {new_counter}\n")
            for s, t, u in frames_frac[idx]:
                f.write(f"  {s: .9f}  {t: .9f}  {u: .9f}\n")

def triclinic_params(a, b, c):
    """Return (lx, ly, lz, xy, xz, yz) from arbitrary cell vectors a,b,c."""
    import numpy as np
    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)
    c = np.array(c, dtype=float)

    lx = float(np.linalg.norm(a))
    if lx <= 0:
        raise ValueError("Lattice vector a has zero length.")
    ex = a / lx

    b_x = float(np.dot(b, ex))
    b_perp = b - b_x * ex
    ly = float(np.linalg.norm(b_perp))
    if ly <= 0:
        raise ValueError("Lattice vectors a and b are colinear.")
    ey = b_perp / ly

    ez = np.cross(ex, ey)
    c_x = float(np.dot(c, ex))
    c_y = float(np.dot(c, ey))
    c_z = float(np.dot(c, ez))
    if c_z < 0:
        ez = -ez
        c_z = -c_z
        c_y = -c_y

    xy, xz, yz = b_x, c_x, c_y
    return lx, ly, c_z, xy, xz, yz

def dump_bounds_lines_lammps(lx, ly, lz, xy, xz, yz):
    """
    LAMMPS triclinic bounding box (Howto_triclinic):
      xlo_bound = xlo + min(0, xy, xz, xy+xz)
      xhi_bound = xhi + max(0, xy, xz, xy+xz)
      ylo_bound = ylo + min(0, yz)
      yhi_bound = yhi + max(0, yz)
      zlo_bound = zlo
      zhi_bound = zhi
    We use xlo=ylo=zlo=0; xhi=lx, yhi=ly, zhi=lz.
    """
    xlo = 0.0
    xhi = lx
    ylo = 0.0
    yhi = ly
    zlo = 0.0
    zhi = lz
    min_x = min(0.0, xy, xz, xy + xz)
    max_x = max(0.0, xy, xz, xy + xz)
    min_y = min(0.0, yz)
    max_y = max(0.0, yz)

    xlo_b = xlo + min_x
    xhi_b = xhi + max_x
    ylo_b = ylo + min_y
    yhi_b = yhi + max_y
    zlo_b = zlo
    zhi_b = zhi

    line1 = f"{xlo_b:.16e} {xhi_b:.16e} {xy:.16e}\n"
    line2 = f"{ylo_b:.16e} {yhi_b:.16e} {xz:.16e}\n"
    line3 = f"{zlo_b:.16e} {zhi_b:.16e} {yz:.16e}\n"
    return line1, line2, line3

def frac_to_triclinic_cart(frac, lx, ly, lz, xy, xz, yz):
    s, t, u = frac
    x = s * lx + t * xy + u * xz
    y = t * ly + u * yz
    z = u * lz
    return x, y, z

def write_dump(output_filename, a, b, c, frames_frac, elements, counts):
    lx, ly, lz, xy, xz, yz = triclinic_params(a, b, c)
    natoms = sum(counts)
    types = build_types(elements, counts)

    with open(output_filename, 'w') as f:
        for step_idx, frame in enumerate(frames_frac, start=1):  # 1-based timestep
            f.write("ITEM: TIMESTEP\n")
            f.write(f"{step_idx}\n")
            f.write("ITEM: NUMBER OF ATOMS\n")
            f.write(f"{natoms}\n")
            f.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
            l1, l2, l3 = dump_bounds_lines_lammps(lx, ly, lz, xy, xz, yz)
            f.write(l1); f.write(l2); f.write(l3)
            f.write("ITEM: ATOMS id type x y z\n")
            atom_id = 1
            for (t_id), (s, t, u) in zip(types, frame):
                x, y, z = frac_to_triclinic_cart((s, t, u), lx, ly, lz, xy, xz, yz)
                f.write(f"{atom_id} {t_id} {x:.5f} {y:.5f} {z:.5f}\n")
                atom_id += 1

def write_poscar_for_indices(header, frames_frac, indices_0b, poscar_prefix="POSCAR"):
    for idx in indices_0b:
        if idx < 0 or idx >= len(frames_frac):
            raise ValueError(f"Index {idx} out of range (0..{len(frames_frac)-1}).")
        filename = f"{poscar_prefix}.{idx}"
        with open(filename, 'w') as f:
            for ln in header:
                f.write(ln)
            f.write("Direct\n")
            for s, t, u in frames_frac[idx]:
                f.write(f"  {s: .9f}  {t: .9f}  {u: .9f}\n")

def parse_indices_file(path):
    indices = []
    with open(path, 'r') as f:
        for lineno, raw in enumerate(f, start=1):
            s = raw.strip()
            if not s:
                continue
            try:
                val = int(s)
            except ValueError:
                raise ValueError(f"Line {lineno} in index file is not an integer: {s}")
            indices.append(val)
    if not indices:
        raise ValueError("Index file is empty.")
    return indices

def main():
    parser = argparse.ArgumentParser(description="XDATCAR utility: clip or export (POSCARs or LAMMPS dump).")
    parser.add_argument("--input_filename", "-i", type=str, default="XDATCAR", help="XDATCAR input filename")
    parser.add_argument("--output_filename", "-o", type=str, default="XDATCAR_new", help="Output filename for clipped XDATCAR or dump")
    parser.add_argument("--start_step", "-s", type=int, default=2500, help="Starting frame index (1-based, inclusive) for clipping")
    parser.add_argument("--end_step", "-e", type=int, default=500000, help="Ending frame index (1-based, inclusive) for clipping")
    parser.add_argument("--index_file", "-f", type=str, help="Path to index file (0-based indices per line) for export")
    parser.add_argument("--poscar_prefix", "-p", type=str, default="POSCAR", help="Prefix for exported POSCAR files when not dumping")
    parser.add_argument("--dump", action="store_true", help="Output in LAMMPS dump format")
    args = parser.parse_args()

    try:
        data = read_xdatcar(args.input_filename)
    except Exception as exc:
        print(f"[ERROR] Failed to read XDATCAR '{args.input_filename}': {exc}", file=sys.stderr)
        sys.exit(1)

    header     = data["header"]
    a, b, c    = data["a"], data["b"], data["c"]
    elements   = data["elements"]
    counts     = data["counts"]
    frames_frac= data["frames_frac"]

    if args.index_file:
        try:
            indices = parse_indices_file(args.index_file)
        except Exception as exc:
            print(f"[ERROR] Reading indices: {exc}", file=sys.stderr)
            sys.exit(2)

        selected_frames = []
        for idx in indices:
            if idx < 0 or idx >= len(frames_frac):
                print(f"[ERROR] Index {idx} out of range (0..{len(frames_frac)-1}).", file=sys.stderr)
                sys.exit(3)
            selected_frames.append(frames_frac[idx])

        if args.dump:
            out = args.output_filename if args.output_filename else "dump.lammpstrj"
            try:
                write_dump(out, a, b, c, selected_frames, elements, counts)
            except Exception as exc:
                print(f"[ERROR] Writing dump: {exc}", file=sys.stderr)
                sys.exit(4)
            print(f"[OK] Wrote dump with {len(selected_frames)} selected frames to '{out}'.")
        else:
            try:
                write_poscar_for_indices(header, frames_frac, indices, poscar_prefix=args.poscar_prefix)
            except Exception as exc:
                print(f"[ERROR] Writing POSCARs failed: {exc}", file=sys.stderr)
                sys.exit(5)
            print(f"[OK] Exported {len(indices)} POSCAR files with prefix '{args.poscar_prefix}'.")
    else:
        start_1b = args.start_step
        end_1b   = min(args.end_step, len(frames_frac))
        if start_1b < 1 or start_1b > end_1b:
            print("[ERROR] Invalid clip range.", file=sys.stderr)
            sys.exit(6)

        frames_clip = frames_frac[start_1b-1:end_1b]

        if args.dump:
            out = args.output_filename if args.output_filename else "dump.lammpstrj"
            try:
                write_dump(out, a, b, c, frames_clip, elements, counts)
            except Exception as exc:
                print(f"[ERROR] Writing dump: {exc}", file=sys.stderr)
                sys.exit(7)
            print(f"[OK] Wrote dump for clipped frames {start_1b}..{end_1b} to '{out}'.")
        else:
            try:
                write_xdatcar_clip(args.output_filename, header, frames_frac, start_1b, end_1b)
            except Exception as exc:
                print(f"[ERROR] Clipping XDATCAR failed: {exc}", file=sys.stderr)
                sys.exit(8)
            print(f"[OK] Clipped XDATCAR written to '{args.output_filename}' "
                  f"(frames {start_1b}..{end_1b}, total {end_1b - start_1b + 1}).")

if __name__ == "__main__":
    main()
