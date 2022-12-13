import pandas as pd
import pickle
import ast
import os

os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab") #insert your own pathname here.

### As the name of this script implies - it prepared prediction data for evaluation, reducing overhead.
### Takes output of kinase_prediction.py and prepares it for use with evaluate_predictions.py


def parsePredictions(prediction_filename):
    with open(prediction_filename) as f:
        lines = f.readlines()
    f.close()
    
    for id, line in enumerate(lines):
        lines[id] = line.split(";")
        lines[id][1] = ast.literal_eval(lines[id][1])
    return lines


def parseKinaseSubstrateDS(filename):
    with open(filename) as f:
        df = pd.read_table(f, header=2)
    f.close()
    return df


def dataPrep(predictions, dataset):
    preds = parsePredictions(predictions)
    data = parseKinaseSubstrateDS(dataset)
    
    data['SITE_+/-7_AA']=data['SITE_+/-7_AA'].str.upper()
    data['SITE_+/-7_AA']=data['SITE_+/-7_AA'].str.strip("_")

    with open(fr'Phosphoproteomics/steps/predictions/{predictions}.pickle', 'wb') as f:
        pickle.dump(preds, f)
    f.close()
    with open(fr'Phosphoproteomics/steps/predictions/{dataset}.pickle', 'wb') as f:
        pickle.dump(data, f)
    f.close()    







dataPrep(r"Phosphoproteomics/steps/predictions/predictions_sugiyama_filtered.txt", r"Phosphoproteomics/data/kin_sub_ds.txt")
