"""
Compute eigenvalue-based shape descriptors for protein binding pockets from PDB files.

Usage:
    python pocket_eigen_ratios.py <pockets_dir> <output_csv>
"""

from Bio.PDB import PDBParser
import numpy as np
import pandas as pd
import os
import sys


def extract_coords(pdb_path, parser):
    """Extract all atom (x,y,z) coordinates from a pocket PDB file."""
    structure = parser.get_structure("pocket", pdb_path)

    coords = []

    for atom in structure.get_atoms():
        coords.append(atom.coord)
    return np.array(coords)

def compute_eigen_ratios(coords):
    """Compute sorted eigenvalues (λ1≥λ2≥λ3) of the coordinate covariance matrix and their ratios."""
    cov_matrix = np.cov(coords, rowvar=False)
    eigvals = np.linalg.eigvalsh(cov_matrix)
    eigvals = np.sort(eigvals)[::-1]

    ev1, ev2, ev3 = eigvals[0], eigvals[1], eigvals[2]
    rEV32 = ev3 / ev2
    rEV31 = ev3 / ev1
    rEV21 = ev2 / ev1

    return ev1, ev2, ev3, rEV32, rEV31, rEV21


def main(pockets_dir, output_csv):
    parser = PDBParser(QUIET=True)

    pdb_files = sorted(
        f for f in os.listdir(pockets_dir)
        if f.lower().endswith(".pdb")
    )

    rows = []
    for fname in pdb_files:
        pocket_id = os.path.splitext(fname)[0]
        pdb_path = os.path.join(pockets_dir, fname)

        coords = extract_coords(pdb_path, parser)
        ev1, ev2, ev3, rEV32, rEV31, rEV21 = compute_eigen_ratios(coords)

        rows.append({
            "PocketID": pocket_id,
            "Eigval1": ev1,
            "Eigval2": ev2,
            "Eigval3": ev3,
            "rEV32": rEV32,
            "rEV31": rEV31,
            "rEV21": rEV21,
            "NumAtoms": len(coords),
        })

    df = pd.DataFrame(rows)
    df[["Eigval1", "Eigval2", "Eigval3", "rEV32", "rEV31", "rEV21"]] = df[
        ["Eigval1", "Eigval2", "Eigval3", "rEV32", "rEV31", "rEV21"]
    ].round(4)
    df.to_csv(output_csv, index=False)
    print(f'CSV successfully generated: "{output_csv}"')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python pocket_eigen_ratios.py <pockets_dir> <output_csv>")
    main(sys.argv[1], sys.argv[2])