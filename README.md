## pocket-eigenvalues
Python script for extracting atomic coordinates from PDB pocket structures and computing eigenvalue-based geometric descriptors (covariance eigenvalues and ratios) for RNA &amp; protein binding pockets.

#### REQUIREMENTS

- Python ≥ 3.8
- Biopython
- NumPy
- Pandas

Dependencies can be installed using:
```bash
pip install -r requirements.txt
```
#### INPUT

- A directory containing pocket PDB file(s) (one pocket per file).
- Each pocket is identified by its filename.

#### USAGE
```bash
python pocket_xyz_eigenvalues.py data/ results.csv
```
Where:
- <input_folder> is the path to a directory containing pocket PDB file(s) (e.g. `data/`)
- <output_csv> is the name of the CSV file to be generated (e.g. `results.csv`)

#### OUTPUT

The script produces a CSV file containing the following columns:

- `PocketID`: pocket identifier (filename without .pdb)
- `Eigval1, Eigval2, Eigval3`: eigenvalues of the coordinate covariance matrix
- `rEV32, rEV31, rEV21`: eigenvalue ratios
- `NumAtoms`: number of atoms in the pocket

#### EXAMPLE DATA

The `data/` directory contains two example pocket PDB files provided for testing the pipeline.
