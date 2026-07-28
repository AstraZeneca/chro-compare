"""
generate_fingerprints.py
Generates required normalised physicochemical descriptors from the physico_descs.txt
file for each SMILES string in the csv file specified in FILENAME (here SMILES.csv)

Functions:
  get_fing(smiles): generates the normalised physicochemical descriptor for the set
     specified in physico_descs.csv for the given smiles
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def get_fing(smiles):
    """Generates physicochemical fingerprint"""
    mol = Chem.MolFromSmiles(smiles)
    fp = []
    for i in range(len(names)):
        fp.append((fns[names[i]](mol) - means[i])/stds[i])
    return str(list(fp))

df = pd.read_csv('SMILES.csv')

descriptors = pd.read_csv('physico_descs.csv')
names = list(descriptors['name'])
means = list(descriptors['mean'])
stds = list(descriptors['std'])

fns = dict(Descriptors._descList)

print('Calculating fingerprints:')
df['Fingerprint'] = df['SMILES'].apply(get_fing)

print('Writing file:')
df.to_csv('data_fingerprints.csv')
