import numpy as np
import math
import pickle
import os
import requests as r
import json
import pandas as pd

os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab") #insert your own pathname here.

# protToAcc = pd.read_csv(r"Phosphoproteomics/data/prot_to_acc.csv", sep=";", na_filter=False)
# protList = protToAcc["Prot"].tolist()
# acclist = protToAcc["Acc"].tolist()
# protAccDict = {}
# for count, n in enumerate(protList):
#     protAccDict[n] = acclist[count]

def parseKinaseSubstrateDS(filename):
    with open(filename) as f:
        df = pd.read_table(f, header=2)
    f.close()
    return df

kin_sub_DS = parseKinaseSubstrateDS(r"Phosphoproteomics/data/kin_sub_ds.txt")

def unpickle(filename):
    with open(fr'Phosphoproteomics/steps/pwms/{filename}.pickle', 'rb') as f:
        out = pickle.load(f)
    f.close()
    return out

def average_kinase(pwm_dict):
    counter = 0
    for subkinase in pwm_dict:
        if counter == 0:
            pwm_total = pwm_dict[subkinase]
            counter =1
        else:
            pwm_total += pwm_dict[subkinase]
    return pwm_total/len(pwm_dict)

def cross_entropy_calculation(pwm_1, pwm_2):
    output = []
    
    pwm_1 = average_kinase(pwm_1)+1e-3
    pwm_2 = average_kinase(pwm_2)+1e-3

    #calculate cross entropy loss
    loss = np.sum(pwm_1 * np.log(pwm_1/pwm_2))

    #normalize by length of pwm
    loss = loss/len(pwm_1)

    #apply sigmoid function to make sure the output is between 0 and 1
    return 1/1*np.exp(-loss)

def get_subcellular_specificity(UniProt_AccID):
    url = "https://www.uniprot.org/uniprotkb/"+UniProt_AccID+".json"
    response = r.get(url)
    subcellular_specificity = []
    try:
        for comment in response.json()["comments"]:
            if comment["commentType"] == "SUBCELLULAR LOCATION":
                for subcellularLocation in comment["subcellularLocations"]:
                    subcellular_specificity.append(subcellularLocation["location"]["value"])
                break
    finally:
        return(subcellular_specificity)


def pwm_dict_cross_entropy(pwm_dict_1, pwm_dict_2):
    #Make sure that pwm_dict_1 is <= pwm_dict_2
    pwm_dict_1 = unpickle(pwm_dict_1)
    pwm_dict_2 = unpickle(pwm_dict_2)
    cross_entropy_output = []
    len_1 = len(pwm_dict_1)
    len_2 = len(pwm_dict_2)
    

    shared_pwms = set(pwm_dict_1.keys()).intersection(pwm_dict_2.keys())
    print(len(shared_pwms), shared_pwms)
    #first we calculate cross entropy scores for the pwms that are found in shared pwms
    for kinase in shared_pwms:
        cross_entropy_output.append([cross_entropy_calculation(pwm_dict_1.pop(kinase), pwm_dict_2.pop(kinase)), kinase])

    shared_kinases = dict()
    ### then we try to pick up stragglers that might've dodged out first intersection search.
    for kinase_1 in pwm_dict_1.keys():
        for kinase_2 in pwm_dict_2.keys():
            if kinase_1.upper() in kinase_2.upper() or kinase_2.upper() in kinase_1.upper():
                if kinase_1 not in shared_kinases:
                    shared_kinases[kinase_1]=kinase_2
       
    for kinase in shared_kinases.keys():
        if kinase != "LCK": ###Excluding this because it's not in the PPS+ database but gets added due to being a substring of MLCK
            cross_entropy_output.append([cross_entropy_calculation(pwm_dict_1.pop(kinase), pwm_dict_2.pop(shared_kinases[kinase])), kinase, shared_kinases[kinase]])

    cross_entropy_output.sort()

    #### THIS BIT BELOW IS MADE TO FIND SUBCELLULAR SPECIFITIES FOR EACH KINASE
    for cross_entropy in cross_entropy_output:
        ACC_IDs = set()
        subcellular_specifities = set()
        for ACC_ID in kin_sub_DS[kin_sub_DS['KINASE'].str.match(cross_entropy[1], na=False)]['KIN_ACC_ID']:
            ACC_IDs.add(ACC_ID)
        for ACC_ID in ACC_IDs:
            for subcellular_location in get_subcellular_specificity(str(ACC_ID)):
                subcellular_specifities.add(subcellular_location)
        cross_entropy.append(subcellular_specifities)
        print(cross_entropy)

    print("outputlength: ",len(cross_entropy_output), " input_1 length: ",len_1," input_2 length: ", len_2)
    return cross_entropy_output
        



print( get_subcellular_specificity("Q05655"))

with open(fr"Phosphoproteomics/steps/sugiyama_ppsplus_cross_entropy.pickle", "wb") as f:
    pickle.dump(pwm_dict_cross_entropy("sugiyama_pwms", "ppsplus_pwms"), f)
f.close()




