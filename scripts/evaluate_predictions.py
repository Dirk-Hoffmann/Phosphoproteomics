import ast
import enum
import pandas as pd
import pickle
import pylcs
import re
import os

os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab") #insert your own pathname here.


prediction_filename = "predictions_selbach"
path_to_predictions = fr"Phosphoproteomics/steps/predictions/{prediction_filename}.txt"


with open(path_to_predictions) as f:
    lines = f.readlines()
f.close()

protToAcc = pd.read_csv(r"Phosphoproteomics/data/prot_to_acc.csv", sep=";", na_filter=False)
protList = protToAcc["Prot"].tolist()
acclist = protToAcc["Acc"].tolist()
protAccDict = {}
for count, n in enumerate(protList):
    protAccDict[n] = acclist[count]

def parsePredictions(prediction_filename):
    with open(prediction_filename) as f:
        lines = f.readlines()
    f.close()
    
    for id, line in enumerate(lines):
        lines[id] = line.split(";")
        lines[id][1] = ast.literal_eval(lines[id][1])
    return lines

def isomer(known_kinase,predicted):
    #Used to determine whether or not a miss is an isomer or not. 
    if pylcs.lcs2(known_kinase,predicted) >= len(re.sub(r'[0-9]+', '', known_kinase)):
        #print(re.sub(r'[0-9]+', '', known_kinase))
        return True
    return False

def parseKinaseSubstrateDS(filename):
    with open(filename) as f:
        df = pd.read_table(f, header=2)
    f.close()
    return df


def makeKinaseList(dataset, phosphosite):
    
    ### dataset is a pandas dataframe containing kinase-substrate pairings from phosphosite+
    kinaseNameList = [] #### LIST OF KNOWN KINASES TARGETING THIS SUBSTRATE
    for kinase in dataset[dataset['SITE_+/-7_AA'].str.match(phosphosite,na=False)]['KINASE']:
            kinase = kinase.upper()
            if kinase not in kinaseNameList:
                kinaseNameList.append(kinase)
    return kinaseNameList    

def makePredictionSet(predictions, length):
    topPredictions = sorted(predictions)[::-1]
    predictionSet = set()
    for prediction in topPredictions:
        if len(predictionSet)==length:
            break
        predictionSet.add(prediction[1]) 
    return(predictionSet)
    
def evaluation(predictions, dataset):

    with open('{}.pickle'.format(predictions), 'rb') as f:
        preds = pickle.load(f)
    f.close()
    data = parseKinaseSubstrateDS(dataset)
    data['SITE_+/-7_AA']=data['SITE_+/-7_AA'].str.upper()
    data['SITE_+/-7_AA']=data['SITE_+/-7_AA'].str.strip("_")

    ratios = [0,0]

    #acc_id_to_kinase = dict(zip(data.KIN_ACC_ID, data.KINASE))
    
    f = open(fr'Phosphoproteomics/results/{prediction_filename}_known_vs_predicted.txt', "a")

    hitsPerKinase = {}
    sitesPerKinase = {}
    isomersPerKinase = {}

    for line in preds:
        phosphosite = line[0]
        predictions = line[1]
        tracker = {}
        counter = 0

        kinaseNameList = makeKinaseList(data, phosphosite)
#####################Metadata##################################
        for kinase in kinaseNameList:
            if kinase not in sitesPerKinase:
                isomersPerKinase[kinase]= {"sum":0}
                hitsPerKinase[kinase] = 0
                sitesPerKinase[kinase] =1
            else:
                sitesPerKinase[kinase]+=1
################################################################
        predictionSet = makePredictionSet(predictions, len(kinaseNameList))
        for knownKinase in kinaseNameList:
            ratios[1]+=1
            if str(predictionSet).find(knownKinase)>-1:
                    print(str(predictionSet).find(knownKinase))
                    ratios[0]+=1
                    hitsPerKinase[knownKinase]+=1
                    counter += 1
######################metadata###########################################
            else:
                for predictedKinase in predictionSet:
                    if isomer(knownKinase, predictedKinase):
                        if predictedKinase not in isomersPerKinase[knownKinase]:
                            isomersPerKinase[knownKinase][predictedKinase]=1
                            isomersPerKinase[knownKinase]["sum"] +=1
                        else:
                            isomersPerKinase[knownKinase][predictedKinase]+=1
                            isomersPerKinase[knownKinase]["sum"] +=1    
############################################################################        
        keys = list(tracker.keys())
        l = "Correctly identified {} kinases | ".format(counter)
        if counter == 0:
            l = ''
        if len(tracker)>0:
            values = tracker[keys[0]]
            f.write(line[0]+ " | " +l +"Known Kinases = "+str(keys) +" | "+ "Predicted Kinases = " +  str(values)+"\n")
    
    f.close()

    with open(fr"Phosphoproteomics/results/hits_sites_per_kinase_{prediction_filename}.txt", "w") as L:

        for key in sitesPerKinase:
            #L.write(str(key) + " hits/sites: "+ str(hitsPerKinase[key])+"/" + str(sitesPerKinase[key]) + " Isomers = "+ str(isomersPerKinase[key])+"\n")
            L.write(F"{key}   hits/sites:  {hitsPerKinase[key]}/  {sitesPerKinase[key]}   Isomers =  {isomersPerKinase[key]}\n")
    L.close()
    print(ratios)


evaluation(path_to_predictions, r"Phosphoproteomics/data/kin_sub_ds.txt")
