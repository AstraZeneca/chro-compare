"""
ChroCompare.py
Defines the function used to perform the ChroCompare analysis.

Functions:
  integ(rf1, rf2, params, pdf, cdf): approximates the integral of the
    precalulated PDF at the predicted rf1 and rf2 given the criteria
    defined by the list params [lower bound, upper bound, required separation]
  compare(fing1, fing2, cols, params): Returns the predicted retention
    factors, and the probability of separation for each of the specified
    columns for two molecular fingerprints given specified separation
    criteria. It also returns the name of the column with the highest
    predicted probability of separation.
"""
import numpy as np
import cloudpickle

def integ(rf1, rf2, params, pdf, cdf):
    """Approximates required double integral"""
    running_sum = 0
    for x in range(int(round(params[0]*1000)),int(round((params[1]-params[2])*1000))):
        running_sum += pdf[int(round(rf1*1000)),x]*(cdf[int(round(rf2*1000)),
                           int(round(x+params[2]*1000))]-cdf[int(round(rf2*1000)),
                           int(round(params[1]*1000))])/1000
    return running_sum

def compare(fing1, fing2, cols, params):
    """Predicts and compares retention factors"""
    preds1 = []
    preds2 = []
    probs = []

    for col in cols:
        with open('./Models/'+col+'_model.cp.pkl', 'rb') as f:
            rf = cloudpickle.load(f)
        with open('./CalibrationData/'+col+'_PDF.cp.pkl', 'rb') as f:
            pdf = cloudpickle.load(f)
        with open('./CalibrationData/'+col+'_CDF.cp.pkl', 'rb') as f:
            cdf = cloudpickle.load(f)

        preds1.append(rf.predict(fing1)[0])
        preds2.append(rf.predict(fing2)[0])

        if preds1[-1] < preds2[-1]:
            prob = integ(preds1[-1], preds2[-1], params, pdf, cdf)
        else:
            prob = integ(preds2[-1], preds1[-1], params, pdf, cdf)

        probs.append(prob)
        print('Tested '+col + ' (' + str(round(preds1[-1],4)) + ',' +
              str(round(preds2[-1],4)) + ',' + str(round(prob,4)) + ')')

    print('Best column: '+cols[np.argmax(probs)]+' with probability '+str(max(probs)))
    return preds1, preds2, probs, cols[np.argmax(probs)]
