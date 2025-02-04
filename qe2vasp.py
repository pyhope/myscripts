#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Convert QE coordinates to VASP input format
"""

import click
from pathlib import Path
import sys
from typing import DefaultDict

@click.command(help=__doc__)
@click.argument("fname", type=click.Path(exists=True))
def main(fname: str):
    with open(fname) as fp:

        sys.stdout.write("{}\n".format(Path(fname).name))
        sys.stdout.write("1.0\n")

        for line in fp:
            if not line:
                raise RuntimeError("Unexpected EOF")
            if line.startswith("CELL_PARAMETERS"):
                cell_unit = line.split()[1][1:-1]
                if cell_unit != "angstrom":
                    raise NotImplementedError(f"{cell_unit} units not supported")
                break

        for _ in range(3):
            sys.stdout.write(fp.readline())

        for line in fp:
            if not line:
                raise RuntimeError("Unexpected EOF")
            if line.startswith("ATOMIC_POSITIONS"):
                coords_unit = line.split()[1][1:-1]
                if coords_unit == "crystal":
                    out_coords_unit = "Direct"
                elif coords_unit == "angstrom":
                    out_coords_unit = "Cartesian"
                else:
                    raise NotImplementedError(f"{coords_unit} units not supported")
                break

        coords = DefaultDict(list)
        while True:
            line = fp.readline()
            if not line:
                break
            atype, acoord = line.split(maxsplit=1)
            coords[atype].append(acoord)

        for atype, _ in coords.items():
            sys.stdout.write(f"{atype:>4}")

        sys.stdout.write("\n")

        for _, acoords in coords.items():
            sys.stdout.write(f"{len(acoords):>4}")

        sys.stdout.write("\n")

        sys.stdout.write(f"{out_coords_unit}\n")

        for _, acoords in coords.items():
            for acoord in acoords:
                sys.stdout.write(acoord)

        sys.stdout.write("\n")


if __name__ == "__main__":
    main()
