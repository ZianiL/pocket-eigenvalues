# pocket-eigenvalues
Python script for extracting atomic coordinates from PDB pocket structures and computing eigenvalue-based geometric descriptors (covariance eigenvalues and ratios) for RNA &amp; protein binding pockets.

## REQUIREMENTS

- Python ≥ 3.8
- Biopython
- NumPy
- Pandas

Dependencies can be installed using:

```bash
pip install -r requirements.txt

## INPUT

- A directory containing pocket PDB file(s) (one pocket per file).
- Each pocket is identified by its filename (without the .pdb extension).

## USAGE

python pocket_eigen_ratios.py data/ output.csv

## OUTPUT

The script produces a CSV file containing the following columns:

- PocketID: pocket identifier (filename without .pdb)
- Eigval1, Eigval2, Eigval3: eigenvalues of the coordinate covariance matrix
- rEV32, rEV31, rEV21: eigenvalue ratios
- NumAtoms: number of atoms in the pocket

## EXAMPLE DATA

The data/ directory contains two example pocket PDB files provided for testing the pipeline.
