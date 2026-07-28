![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)

# ChroCompare

ChroCompare predicts whether two compounds can be adequately separated on a
given chromatographic method, and helps identify which of several available
methods is most likely to succeed. For each candidate column it predicts the
retention factor of both compounds and, from a calibrated distribution of
model errors, a probability that the pair will be separated within a
required region and by a required margin. This supports choosing the most
promising chromatographic method for a separation, as well as assigning
likely elution order and peak identity.

This repository contains the code used to train the retention-factor models,
calibrate the separation-probability predictions, and reproduce the pairwise
test-set evaluation described in the accompanying manuscript
(`DRAFT18_paper.pdf`).

## Software requirements

- Python 3.9+ (developed and tested with CPython; no OS-specific dependencies)
- [RDKit](https://www.rdkit.org/) (descriptor calculation)
- pandas, numpy, scipy
- scikit-learn (`RandomForestRegressor`)
- cloudpickle
- tqdm

No `requirements.txt` / dependency pinning is included yet — see the note in
`CONTRIBUTING.md` about the project's current maturity level.

## Repository contents

| File | Purpose |
|---|---|
| `generate_fingerprints.py` | Calculates normalised physicochemical descriptor fingerprints for a set of SMILES strings, using the descriptors and normalisation constants in `physico_descs.csv`. |
| `generate_calibration.py` | Trains a random forest retention-factor model per chromatographic column (as listed in `columns.csv`) and builds the calibration PDF/CDF used to turn a pair of predictions into a separation probability. Saves models to `./Models/` and calibration data to `./CalibrationData/`. |
| `ChroCompare.py` | Core prediction function (`compare`): given two fingerprints, returns the predicted retention factor of each on every column and the probability of separation, plus the best column overall. |
| `testing.py` | Runs `compare()` for every pair of compounds in a held-out test set and writes a detailed per-pair results file (`results.csv`), including comparisons against two "naive" (non-probabilistic) baseline strategies. |
| `run_single.py` | Runs the full pipeline for a single pair of SMILES strings, without needing a pre-built dataset file. See "Usage" below. |
| `analysis.py`, `percentage_plot.py`, `AUC_curve.py` | Analyse and plot the results produced by `testing.py` (calibration curves, ROC/AUC, method-recommendation accuracy). |
| `columns.csv` | List of chromatographic columns/methods used, with their run length (used to normalise retention time to retention factor). |
| `physico_descs.csv` | List of RDKit physicochemical descriptors used as model input, with the mean/standard deviation used to normalise each one. |

## Usage

1. **Train models and calibration data** (one-off, requires a labelled
   dataset with a `Fingerprint` column, as produced by
   `generate_fingerprints.py`, and one retention-time column per
   chromatographic method):

   ```bash
   python generate_calibration.py
   ```

   This populates `./Models/` and `./CalibrationData/`, which are required
   by everything below.

2. **Predict for a single pair of compounds:**

   ```bash
   python run_single.py "<SMILES_1>" "<SMILES_2>"
   ```

3. **Evaluate on a full test set** (a CSV with `SMILES` and `Fingerprint`
   columns, plus one retention-time column per method, as `TEST_SET.csv`):

   ```bash
   python testing.py
   ```

   followed by `analysis.py`, `percentage_plot.py`, and `AUC_curve.py` to
   analyse and plot the resulting `results.csv`.

## Citing

If you use this code, please cite the accompanying manuscript (will attach).

## License

Apache License 2.0 — see `LICENSE.md`.
