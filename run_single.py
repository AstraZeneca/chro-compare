"""
run_single.py
Runs the ChroCompare algorithm directly on a single pair of input SMILES
strings, without requiring a pre-built dataset or fingerprint file.

For each molecule, the normalised physicochemical descriptor fingerprint is
calculated on the fly (using the descriptor set and normalisation constants
in physico_descs.csv), and then passed to ChroCompare.compare() to predict,
for every chromatographic column listed in columns.csv, the retention factor
of each molecule and the probability that the pair will be adequately
separated.

Requires the trained models and calibration data produced by
generate_calibration.py to already exist in ./Models/ and
./CalibrationData/.

Usage:
    python run_single.py "<SMILES_1>" "<SMILES_2>"

Functions:
  get_fingerprint(smiles, names, means, stds, descriptor_fns): generates the
     normalised physicochemical fingerprint for a single SMILES string.
  run(smiles1, smiles2): generates fingerprints for both molecules and
     returns the ChroCompare prediction for the pair.
"""

import sys

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

from ChroCompare import compare

COLUMN_FILE = 'columns.csv'
DESCRIPTOR_FILE = 'physico_descs.csv'

REQUIRED_RANGE = [0.05, 0.95]
REQUIRED_SEPARATION = 0.1


def get_fingerprint(smiles, names, means, stds, descriptor_fns):
    """Generates the normalised physicochemical fingerprint for one SMILES string.

    This mirrors get_fing() in generate_fingerprints.py, but returns a numpy
    array directly (rather than a string for saving to csv) and validates
    the SMILES string, since here there is no upstream data-cleaning step.

    Args:
        smiles: SMILES string of the molecule.
        names: list of RDKit descriptor names to calculate, in order.
        means: list of per-descriptor means used for normalisation.
        stds: list of per-descriptor standard deviations used for
            normalisation.
        descriptor_fns: dict mapping descriptor name to its RDKit function,
            e.g. dict(Descriptors._descList).

    Returns:
        A 1 x n_descriptors numpy array of normalised descriptor values,
        matching the shape expected by the trained models in compare().

    Raises:
        ValueError: if the SMILES string cannot be parsed by RDKit.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('Could not parse SMILES string: ' + smiles)

    fingerprint = [
        (descriptor_fns[names[i]](mol) - means[i]) / stds[i]
        for i in range(len(names))
    ]
    return np.array(fingerprint).reshape(1, -1)


def run(smiles1, smiles2):
    """Runs the ChroCompare algorithm for a single pair of SMILES strings.

    Args:
        smiles1: SMILES string of the first molecule.
        smiles2: SMILES string of the second molecule.

    Returns:
        The tuple (preds1, preds2, probs, best_column) as returned by
        ChroCompare.compare():
          preds1, preds2: predicted (normalised) retention factor of
            molecule 1 and 2 respectively, for each column in columns.csv.
          probs: predicted probability of adequate separation, for each
            column in columns.csv.
          best_column: name of the column with the highest predicted
            separation probability.
    """
    columns = pd.read_csv(COLUMN_FILE)
    cols = list(columns['name'])

    descriptors = pd.read_csv(DESCRIPTOR_FILE)
    names = list(descriptors['name'])
    means = list(descriptors['mean'])
    stds = list(descriptors['std'])
    descriptor_fns = dict(Descriptors._descList)

    fing1 = get_fingerprint(smiles1, names, means, stds, descriptor_fns)
    fing2 = get_fingerprint(smiles2, names, means, stds, descriptor_fns)

    params = REQUIRED_RANGE + [REQUIRED_SEPARATION]
    return compare(fing1, fing2, cols, params)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python run_single.py "<SMILES_1>" "<SMILES_2>"')
        sys.exit(1)

    run(sys.argv[1], sys.argv[2])
