#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Generate Monkhorst-Pack q-point grid.")

parser.add_argument("nq1", type=int, help="Number of q-points along direction 1")
parser.add_argument("nq2", type=int, nargs="?", help="Number of q-points along direction 2 (default = nq1)")
parser.add_argument("nq3", type=int, nargs="?", help="Number of q-points along direction 3 (default = nq1)")

group = parser.add_mutually_exclusive_group()
group.add_argument("-ct", "--cartesian", action="store_true", help="Output Cartesian coordinates (Å⁻¹)")
group.add_argument("-qe", action="store_true", help="Output QE-style coordinates (Cartesian, unit = 2π/a)")

parser.add_argument("-a", type=float, default=4.376208478, help="Lattice constant a (Å)")
parser.add_argument("-b", type=float, default=4.636170450, help="Lattice constant b (Å)")
parser.add_argument("-c", type=float, default=6.382051746, help="Lattice constant c (Å)")

parser.add_argument("--ba", type=float, help="b/a ratio (used only in -qe mode)", default=1.059403)
parser.add_argument("--ca", type=float, help="c/a ratio (used only in -qe mode)", default=1.458352)

parser.add_argument("-ir","--irreducible", action="store_true", help="Output irreducible q-points and weights using spglib")

args = parser.parse_args()

def wrap_reduced_coords(q):
    return ((q + 0.5) % 1.0) - 0.5

nq1 = args.nq1
nq2 = args.nq2 if args.nq2 is not None else nq1
nq3 = args.nq3 if args.nq3 is not None else nq1

if args.cartesian:
    mode = "cartesian"
elif args.qe:
    mode = "qe"
else:
    mode = "reduced"

a = args.a

if mode == "qe":
    b = a * args.ba
    c = a * args.ca
else:
    b = args.b
    c = args.c

if args.irreducible:
    import spglib

    # Build dummy crystal structure for spglib
    lattice = np.array([[a, 0, 0],
                        [0, b, 0],
                        [0, 0, c]])
    # Use P1 symmetry, minimal atoms
    positions = [[0, 0, 0]]
    numbers = [1]
    cell = (lattice, positions, numbers)

    mesh = [nq1, nq2, nq3]
    mapping, grid = spglib.get_ir_reciprocal_mesh(mesh=mesh, cell=cell, is_shift=[0, 0, 0])
    ir_points, weights = np.unique(grid[mapping], axis=0, return_counts=True)

    total_weight = np.sum(weights)
    ir_points = ir_points / np.array(mesh)
    ir_points = wrap_reduced_coords(ir_points)

    print(len(ir_points))
    for q, w in zip(ir_points, weights):
        if mode == "reduced":
            print(f"  {q[0]:.6f}  {q[1]:.6f}  {q[2]:.6f}  {w/total_weight:.6f}")
        elif mode == "cartesian":
            q_cart = np.array([q[0] * 2 * np.pi / a,
                               q[1] * 2 * np.pi / b,
                               q[2] * 2 * np.pi / c])
            print(f"  {q_cart[0]:.6f}  {q_cart[1]:.6f}  {q_cart[2]:.6f}  {w/total_weight:.6f}")
        elif mode == "qe":
            qx_qe = q[0]
            qy_qe = (a / b) * q[1]
            qz_qe = (a / c) * q[2]
            print(f"  {qx_qe:.6f}  {qy_qe:.6f}  {qz_qe:.6f}  {w/total_weight:.6f}")

else:
    q_points = []

    for i in range(nq1):
        for j in range(nq2):
            for k in range(nq3):
                q1 = wrap_reduced_coords(i / nq1)
                q2 = wrap_reduced_coords(j / nq2)
                q3 = wrap_reduced_coords(k / nq3)

                if mode == "reduced":
                    q = (q1, q2, q3)

                elif mode == "cartesian":
                    bx = 2 * np.pi / a
                    by = 2 * np.pi / b
                    bz = 2 * np.pi / c
                    q = (q1 * bx, q2 * by, q3 * bz)

                elif mode == "qe":
                    qx_qe = q1
                    qy_qe = (a / b) * q2
                    qz_qe = (a / c) * q3
                    q = (qx_qe, qy_qe, qz_qe)

                q_points.append(q)

    print(len(q_points))
    for q in q_points:
        print(f"  {q[0]:.6f}  {q[1]:.6f}  {q[2]:.6f}")
