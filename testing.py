"""
testing.py
Performs the ChroCompare algorithm for all pairs of compounds the specified
test set. The results file is written line by line and is initialised in
each instance of running this script.

Functions:
  get_list(list_in): returns the list of floats from the string present
     in the saved csv.
  in_range(pred1, pred2): checks if both predictions are in the required
     range.
  good_sep(pred1, pred2): checks if two retentions meet the requirements
     for good separation in this analysis that both are in the required
     region and that they are at least the required separation apart.
  naive_model1(preds1, preds2): returns the index of the method with the
     largest predicted separation and whether the molecule listed first
     elutes second.
  naive_model2(preds1, preds2): returns the index of the method with the
     largest predicted separation where both predictions are within the
     required range and whether the molecule listed first elutes second.
"""

import pandas as pd
import numpy as np
from tqdm import tqdm
from ChroCompare import compare

REQUIRED_RANGE = [0.05, 0.95]
REQUIRED_SEPARATION = 0.1

def get_list(list_in):
    """Converts from string in csv to list of floats"""
    result = list_in.strip('[]').split(', ')
    result = [float(item) for item in result]
    return result

def in_range(val1, val2):
    """Checks if both values are in required range"""
    return (REQUIRED_RANGE[0] <= val1 <= REQUIRED_RANGE[1] and
            REQUIRED_RANGE[0] <= val2 <= REQUIRED_RANGE[1])

def good_sep(val1, val2):
    """Checks if values fulfill separation criteria"""
    if in_range(val1, val2) and abs(val1 - val2) >= REQUIRED_SEPARATION:
        return True, val1 > val2
    return False, val1 > val2

def naive_model1(preds1, preds2):
    """Returns the index of the method with the largest predicted separation
       and whether the molecule listed first elutes second"""
    comparison = np.array(preds1) - np.array(preds2)
    index = np.argmax(abs(comparison))
    return index, preds1[index] > preds2[index]

def naive_model2(preds1, preds2):
    """Returns the index of the method with the largest perdicted separation
       where both are in the required range and whether the molecules listed
       first elutes second"""
    comparison = []
    for i in range(len(preds1)):
        if in_range(preds1[i], preds2[i]):
            comparison.append(abs(preds1[i]-preds2[i]))
        else:
            comparison.append(0)
    index = np.argmax(comparison)
    return index, preds1[index] > preds2[index]


df = pd.read_csv('TEST_SET.csv')

columns = pd.read_csv('columns.csv')
cols = list(columns['name'])
lengths = list(columns['length'])

PARAMS = REQUIRED_RANGE + [REQUIRED_SEPARATION]

line = 'SMILES1,SMILES2,'
for col in cols:
    add = col+'_true1,'+col+'_true2,'+col+'_pred1,'+col+'_pred2,'+col+'_sep,'
    line += add

line += 'col_choice,choice_first,choice_sep,pred_naive1,pred_1_first,true_naive1,true_1_first,'
line += 'pred_naive2,pred_2_first,true_naive2,true_2_first\n'

FILE_OUT = 'results.csv'
with open(FILE_OUT, 'w', encoding='UTF-8') as file:
    file.write(line)

N = len(list(df['SMILES']))
print('Running test - saving results to: results.csv')
print('(Note the progress bar proceeds triangularly)')
for i in tqdm(range(N)):
    for j in range(i+1, N):
        print([i,j])
        two = df[df['Unnamed: 0'] == i]
        two = pd.concat([two,df[df['Unnamed: 0'] == j]])

        fings = list(two['Fingerprint'].apply(get_list))

        row = str(list(two['SMILES'])).strip('[]')

        fing1 = np.array(list(fings[0])).reshape(1,-1)
        fing2 = np.array(list(fings[1])).reshape(1,-1)
        row += ','

        pred1, pred2, probs, prediction = compare(fing1, fing2, cols, PARAMS)

        true1 = []
        true2 = []
        for k in range(len(cols)):
            row += str(list(np.array(list(two[cols[k]]))/lengths[k])).strip('[]')
            row += ','+str(pred1[k])+','+str(pred2[k])+','+str(probs[k])+','
            true1.append(list(two[cols[k]])[0]/lengths[k])
            true2.append(list(two[cols[k]])[1]/lengths[k])

        row += prediction + ','

        ind = cols.index(prediction)
        row += str(pred1[ind] > pred2[ind])+','
        separated, first = good_sep(true1[ind],true2[ind])
        if (pred1[ind] > pred2[ind]) == first:
            row += str(separated)+','
        else:
            row += 'False,'

        naive1_pred, naive1_pred_first = naive_model1(pred1,pred2)
        naive1_true, naive1_true_first = naive_model1(true1,true2)
        naive2_pred, naive2_pred_first = naive_model2(pred1,pred2)
        naive2_true, naive2_true_first = naive_model2(true1,true2)

        row += cols[naive1_pred]+','+str(naive1_pred_first)+','
        row += cols[naive1_true]+','+str(naive1_true_first)+','
        row += cols[naive2_pred]+','+str(naive2_pred_first)+','
        row += cols[naive2_true]+','+str(naive2_true_first)+'\n'

        with open(FILE_OUT, 'a',encoding='UTF-8') as file:
            file.write(row)
