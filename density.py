#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute crystal density (g/cm^3) from unit-cell volume and composition.

Only elements considered: Mg, Si, O, Fe.

Density = (mass of atoms in unit cell) / (unit-cell volume)
- mass in grams: sum_i (n_i * atomic_weight_i) / N_A
- volume in cm^3: V(Å^3) * 1e-24  OR  V(cm^3) as provided
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

AVOGADRO = 6.02214076e23  # 1/mol, exact by SI definition

# Standard atomic weights (g/mol), representative IUPAC values
ATOMIC_WEIGHTS = {
    "Mg": 24.305,
    "Si": 28.085,
    "O": 15.999,
    "Fe": 55.845,
}


@dataclass(frozen=True)
class Composition:
    Mg: int = 0
    Si: int = 0
    O: int = 0
    Fe: int = 0

    def validate(self) -> None:
        for el in ("Mg", "Si", "O", "Fe"):
            val = getattr(self, el)
            if val < 0:
                raise ValueError(f"Atom count for {el} must be non-negative (got {val}).")
        if (self.Mg + self.Si + self.O + self.Fe) == 0:
            raise ValueError("Composition has zero atoms; provide at least one atom count.")


def density_g_cm3(volume: float, volume_unit: str, comp: Composition) -> float:
    """
    Parameters
    ----------
    volume : float
        Unit-cell volume.
    volume_unit : str
        Either 'A3' (angstrom^3) or 'cm3' (cm^3).
    comp : Composition
        Counts of Mg, Si, O, Fe atoms in the unit cell.

    Returns
    -------
    float
        Density in g/cm^3.
    """
    comp.validate()

    if volume <= 0.0:
        raise ValueError(f"Volume must be > 0 (got {volume}).")

    vu = volume_unit.strip().lower()
    if vu in ("a3", "ang3", "angstrom3", "angstrom^3", "å3"):
        volume_cm3 = volume * 1.0e-24  # 1 Å^3 = 1e-24 cm^3
    elif vu in ("cm3", "cm^3"):
        volume_cm3 = volume
    else:
        raise ValueError("volume_unit must be 'A3' or 'cm3'.")

    molar_mass_unitcell = (
        comp.Mg * ATOMIC_WEIGHTS["Mg"]
        + comp.Si * ATOMIC_WEIGHTS["Si"]
        + comp.O * ATOMIC_WEIGHTS["O"]
        + comp.Fe * ATOMIC_WEIGHTS["Fe"]
    )  # g/mol for the whole unit cell "formula" (counts as given)

    mass_g_unitcell = molar_mass_unitcell / AVOGADRO  # g
    rho = mass_g_unitcell / volume_cm3  # g/cm^3
    return rho


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Calculate density (g/cm^3) from unit-cell volume and Mg/Si/O/Fe counts."
    )
    p.add_argument(
        "--volume",
        '-v',
        type=float,
        required=True,
        help="Unit-cell volume (use --unit to specify A3 or cm3).",
    )
    p.add_argument(
        "--unit",
        type=str,
        default="A3",
        help="Volume unit: A3 (angstrom^3) or cm3. Default: A3",
    )
    p.add_argument("-Mg", type=int, default=0, help="Number of Mg atoms in unit cell.")
    p.add_argument("-Si", type=int, default=0, help="Number of Si atoms in unit cell.")
    p.add_argument("-O", type=int, default=0, help="Number of O atoms in unit cell.")
    p.add_argument("-Fe", type=int, default=0, help="Number of Fe atoms in unit cell.")
    return p


def main() -> None:
    args = build_argparser().parse_args()
    comp = Composition(Mg=args.Mg, Si=args.Si, O=args.O, Fe=args.Fe)
    rho = density_g_cm3(args.volume, args.unit, comp)
    print(f"Density = {rho:.6f} g/cm^3")


if __name__ == "__main__":
    main()
