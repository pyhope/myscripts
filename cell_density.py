#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute density (g/cm^3) for an orthorhombic unit cell given a, b, c (Å) and
the counts of Mg, Si, O, Fe atoms in the cell.

Example:
    python cell_density.py -a 4.75 -b 5.20 -c 6.10 -Mg 2 -Si 1 -O 3 -Fe 0
"""
import argparse

AVOGADRO = 6.02214076e23  # mol^-1
ANGSTROM3_TO_CM3 = 1e-24  # 1 Å^3 = 1e-24 cm^3

# Standard atomic weights (IUPAC 2019 conventional values)
ATOMIC_MASS = {
    "Mg": 24.305,    # g/mol
    "Si": 28.0855,   # g/mol
    "O":  15.999,    # g/mol
    "Fe": 55.845,    # g/mol
}

def positive_float(value: str) -> float:
    x = float(value)
    if x <= 0:
        raise argparse.ArgumentTypeError("Value must be > 0.")
    return x

def nonnegative_int(value: str) -> int:
    x = int(value)
    if x < 0:
        raise argparse.ArgumentTypeError("Count must be >= 0.")
    return x

def compute_density(a, b, c, n_Mg, n_Si, n_O, n_Fe):
    # Volume in cm^3
    volume_A3 = a * b * c
    volume_cm3 = volume_A3 * ANGSTROM3_TO_CM3

    # Mass of the unit cell in grams
    molar_mass = (
        n_Mg * ATOMIC_MASS["Mg"]
        + n_Si * ATOMIC_MASS["Si"]
        + n_O  * ATOMIC_MASS["O"]
        + n_Fe * ATOMIC_MASS["Fe"]
    )  # g/mol for the cell contents
    mass_g = molar_mass / AVOGADRO  # grams per cell

    density = mass_g / volume_cm3  # g/cm^3
    return density, volume_cm3, mass_g, molar_mass

def main():
    parser = argparse.ArgumentParser(
        description="Compute density (g/cm^3) of an orthorhombic unit cell from a,b,c (Å) and atom counts."
    )
    parser.add_argument("-a", type=positive_float, required=True, help="a (Å)")
    parser.add_argument("-b", type=positive_float, required=True, help="b (Å)")
    parser.add_argument("-c", type=positive_float, required=True, help="c (Å)")
    parser.add_argument("-Mg", type=nonnegative_int, required=True, help="Number of Mg atoms in the cell")
    parser.add_argument("-Si", type=nonnegative_int, required=True, help="Number of Si atoms in the cell")
    parser.add_argument("-O",  type=nonnegative_int, required=True, help="Number of O atoms in the cell")
    parser.add_argument("-Fe", type=nonnegative_int, required=True, help="Number of Fe atoms in the cell")
    args = parser.parse_args()

    density, volume_cm3, mass_g, molar_mass = compute_density(
        args.a, args.b, args.c, args.Mg, args.Si, args.O, args.Fe
    )

    print(f"Cell parameters (Å): a={args.a:.6f}, b={args.b:.6f}, c={args.c:.6f}")
    print(f"Unit cell volume: {args.a*args.b*args.c:.6f} Å^3  =  {volume_cm3:.6e} cm^3")
    print(f"Cell composition: Mg={args.Mg}, Si={args.Si}, O={args.O}, Fe={args.Fe}")
    print(f"Molar mass of cell contents: {molar_mass:.6f} g/mol")
    print(f"Mass per cell: {mass_g:.6e} g")
    print(f"Density: {density:.6f} g/cm^3")

if __name__ == "__main__":
    main()
