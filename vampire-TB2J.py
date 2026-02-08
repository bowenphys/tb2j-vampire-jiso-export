#!/usr/bin/env python3
"""
   Export TB2J results to a Vampire isotropic UCF file (Jiso only)
   - Eastern Institute For Advanced Study 
   - Bowen
"""

from __future__ import annotations

import argparse
import os
import warnings

import numpy as np
from ase.units import J

from TB2J.io_exchange.edit import load


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Export a Vampire UCF file with isotropic interactions only from "
            "TB2J results (TB2J.pickle)."
        )
    )
    parser.add_argument(
        "--inpath",
        default="TB2J_results",
        help=(
            "Path to TB2J results directory or TB2J.pickle file. "
            "Default: TB2J_results"
        ),
    )
    parser.add_argument(
        "--out",
        default="TB2J_results/Vampire/vampire_isotropic.UCF",
        help=(
            "Output Vampire UCF file path. "
            "Default: TB2J_results/Vampire/vampire_isotropic.UCF"
        ),
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=None,
        help=(
            "Optional absolute cutoff for filtering interactions. "
            "Only interactions with |Jij| >= cutoff are written."
        ),
    )
    parser.add_argument(
        "--cutoff-unit",
        choices=["joule", "eV"],
        default="joule",
        help=(
            "Unit for --cutoff. 'joule' applies cutoff on exported Jij (default). "
            "'eV' applies cutoff on 2*Jiso in eV."
        ),
    )
    return parser.parse_args()


def _warn_if_non_orthogonal(cell: np.ndarray, atol: float = 1e-8) -> None:
    v1, v2, v3 = cell
    if (
        abs(float(np.dot(v1, v2))) > atol
        or abs(float(np.dot(v1, v3))) > atol
        or abs(float(np.dot(v2, v3))) > atol
    ):
        warnings.warn(
            "Non-orthogonal lattice vectors detected. Vampire may require "
            "orthogonal vectors depending on your setup.",
            RuntimeWarning,
            stacklevel=2,
        )


def _prepare_isotropic_interactions(exchange_dict, cutoff, cutoff_unit):
    interactions = []
    for key, jiso in exchange_dict.items():
        jij_joule = float(np.real(jiso * 2.0 / J))
        if abs(jij_joule) < 1e-30:
            jij_joule = 0.0

        if cutoff is not None:
            if cutoff_unit == "joule":
                metric = abs(jij_joule)
            else:  # cutoff_unit == "eV"
                metric = abs(float(np.real(jiso * 2.0)))
            if metric < cutoff:
                continue

        interactions.append((key, jij_joule))
    return interactions


def export_vampire_isotropic_ucf(
    spinio, out_path: str, cutoff: float | None = None, cutoff_unit: str = "joule"
) -> tuple[int, int]:
    exchange_dict = getattr(spinio, "exchange_Jdict", None)
    if not exchange_dict:
        raise ValueError("exchange_Jdict is empty or missing in TB2J results.")

    nspins = sum(1 for idx in spinio.index_spin if idx >= 0)
    if nspins == 0:
        raise ValueError("No magnetic atoms found (index_spin has no values >= 0).")

    cell = spinio.atoms.get_cell().array
    _warn_if_non_orthogonal(cell)
    lattice_parameters = np.linalg.norm(cell, axis=-1)
    if np.any(np.isclose(lattice_parameters, 0.0)):
        raise ValueError("Invalid cell: zero lattice parameter detected.")

    scaled_positions = spinio.atoms.get_scaled_positions()
    interactions = _prepare_isotropic_interactions(exchange_dict, cutoff, cutoff_unit)

    with open(out_path, "w", encoding="utf-8") as myfile:
        myfile.write("# Unit cell size (Angstrom):\n")
        np.savetxt(myfile, [lattice_parameters], fmt="%f")

        myfile.write("# Unit cell lattice vectors:\n")
        # Keep the same normalization style as the existing TB2J writer.
        np.savetxt(myfile, cell / lattice_parameters, fmt="%f")

        myfile.write("# Atoms\n")
        myfile.write(f"{nspins} {nspins}\n")
        for i, id_spin in enumerate(spinio.index_spin):
            if id_spin >= 0:
                pos = scaled_positions[i]
                myfile.write(
                    f"{id_spin} {pos[0]} {pos[1]} {pos[2]} {id_spin}\n"
                )

        myfile.write("# Interactions\n")
        myfile.write(f"{len(interactions)} isotropic\n")

        for iid, (key, jval) in enumerate(interactions):
            R, ispin, jspin = key
            myfile.write(
                f"{iid:5d} {ispin:3d} {jspin:3d} "
                f"{R[0]:3d} {R[1]:3d} {R[2]:3d} {jval:< 12.5e}\n"
            )

    return len(interactions), len(exchange_dict)


def main() -> int:
    args = parse_args()

    try:
        spinio = load(args.inpath)
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc
    except Exception as exc:  # pragma: no cover - guard for corrupted inputs
        raise SystemExit(
            f"Failed to load TB2J results from '{args.inpath}': {exc}"
        ) from exc

    out_path = os.path.abspath(args.out)
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    try:
        n_written, n_total = export_vampire_isotropic_ucf(
            spinio,
            out_path,
            cutoff=args.cutoff,
            cutoff_unit=args.cutoff_unit,
        )
    except ValueError as exc:
        raise SystemExit(f"Error: {exc}") from exc

    print(f"Wrote isotropic Vampire UCF to: {out_path}")
    print(f"Interactions written: {n_written}/{n_total}")
    if args.cutoff is not None:
        print(f"Cutoff: {args.cutoff} ({args.cutoff_unit})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

