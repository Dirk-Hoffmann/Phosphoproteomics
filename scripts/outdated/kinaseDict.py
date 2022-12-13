# import pickle

# with open(r"/Users/dirk/Documents/UniBas/Zavolab/Phosphoproteomics/selbach_pwms.pickle", 'rb') as f:
#     pwm = pickle.load(f)

# print(pwm)

import requests as r
import pandas as pd
import os

os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab") #insert your own pathname here.


def parseKinaseSubstrateDS(filename):
    with open(filename) as f:
        df = pd.read_table(f, header=2)
    f.close()
    return df

df = parseKinaseSubstrateDS(r"Phosphoproteomics/data/kin_sub_ds.txt")


def nameFromAccID(Acc_ID):
    response = r.post("https://www.uniprot.org/uniprot/{}.fasta".format(Acc_ID))
    data=''.join([str(lines, 'utf-8') for lines in response.iter_lines()][:1])

    for text in data.split(" "):
        if "GN=" in text:
            return text.strip("GN=")

AccIDToName = {}

for AccID in df["KIN_ACC_ID"]:
    if AccID not in AccIDToName:
        None



#print(df)
