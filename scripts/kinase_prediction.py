import pandas as pd
import pickle
import numpy as np
import math
import requests as r
from background import bg_rates
import os

os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab") #insert your own pathname here.

### This script calculates the probability of a phosphorylated peptide being explained by a kinase based on the 2022 paper by Ghosh, Souvik et. al. https://pubmed.ncbi.nlm.nih.gov/35234914/
### Specifically from the methods section "Inference of kinase activity from phosphoproteome data"


pwms_filename = "pps+_pwms"


def import_FastaDB(path:str):
    ### Imports a dictionary of ACCID:FASTA as created by ConstructFastaDB.py
    with open(path, "rb") as f:
        fastaDB = pickle.load(f)
    f.close()
    return fastaDB

### df is the phosphosites from the 2022 paper by Ghosh, Souvik et. al.
df = pd.read_csv(r"Phosphoproteomics/data/CFIm_KD_Phospho_Peptides.csv", sep=";")

def idToUrl(id):
    currentUrl = "https://www.uniprot.org/uniprot/"+id+".fasta"
    return currentUrl
    

phosphosites = df.values.tolist()

def indexConverter(sequence1, sequence2, index1):
    #converts index of phosphosite to index in Fasta sequence
    for i in range(len(sequence2)):
        x = 0
        for n in range(len(sequence1)):
            x+=1
            if sequence2[i+n]!=sequence1[n]:
                break
        if x == len(sequence1):
            return index1+i

with open(fr'Phosphoproteomics/steps/pwms/{pwms_filename}.pickle', 'rb') as f:
    pwms_loaded = pickle.load(f) 
f.close()


def centre(phosphosite):
    ### Finds the centre amino acid of a phosphosite.
    idx = []
    for i in range(len(phosphosite[1])):
        if phosphosite[1][i] == "[":
            it = i
            id = ""
            while phosphosite[1][it+1] != "]":
                it+=1
                id+= phosphosite[1][it]
            #-1 because the phosphosite data indexes from 1 for some reason >:(
            if phosphosite[1][it+3] != "O":
                idx.append(int(id)-1)
    return idx

#relativeIndexes = [i for i in range(-7,8) if i!= 0]
aas="ARNDCQEGHILKMFPSTWYV"

#{key: value for (key, value) in iterable}
aaToIndex = {key:value for (key, value) in zip(aas, range(len(aas)))}

def likelihoodOfKinase(background, pwm, peptideSequence):
    #calculates the likelihood of a kinase phosphorylating a phosphosite

    l_k = pwm.shape[1] # length of PWM
    l_i = len(peptideSequence) # length of peptide

    likelihood = 0

    ### This catches cases where the phosphosite is at the very end of a peptide. 
    ###In this case we just predict as much as possible
    if l_i-l_k <= 0: ### checks that the PWM is longer than the peptideSequence
        fullyPredictedProbability = 0
        for index, AA in enumerate(peptideSequence):
            if pwm[aaToIndex[AA], index] == 0:
                fullyPredictedProbability+= math.log(0.1*10**-10)
            else:
                fullyPredictedProbability += math.log(pwm[aaToIndex[AA], index])
        likelihood += fullyPredictedProbability
        return likelihood

    ### From here on were continuing as planned.

    for j in range(l_i-l_k): ### the sum        
        subSum = 0
        ### part 1
        for n in range(0, j-1):
            subSum += math.log(background[aaToIndex[peptideSequence[n]]])
        
        ### part 2
        for n in range(j, j+l_k-1):
            if pwm[aaToIndex[peptideSequence[n]], n-j] == 0:
                subSum+= math.log(0.1*10**-10)
            else:
                subSum += math.log(pwm[aaToIndex[peptideSequence[n]], n-j])
        ### part 3
        for n in range(j+l_k,l_i-1):
            subSum += math.log(background[aaToIndex[peptideSequence[n]]])
        likelihood+= math.e**subSum 
    return likelihood

def slicer(phosphosite:str, proteinSeq:str, relativeIndexes:list):
    ### slice returns the peptide + the number of aa's on either side specified by relativeIndexes to enable "full prediction" of the site.
    if indexConverter(phosphosite,proteinSeq, 0)+relativeIndexes[0] < 0:
        left = 0
    else:
        left = indexConverter(phosphosite,proteinSeq,0)+relativeIndexes[0]
    if indexConverter(phosphosite,proteinSeq,len(phosphosite))+relativeIndexes[-1] > len(proteinSeq):
        right = len(proteinSeq)-1
    else:
        right = indexConverter(phosphosite,proteinSeq,len(phosphosite))+relativeIndexes[-1]

    return proteinSeq[left:right]

def likelihoodDictionary(pwmDictionary:dict, peptideSequence:str, bgRates):
    #creates a dictionary of likelihoods for each PWM to the phosphosite.
    likelihoodDict = {}
    for kinase in pwmDictionary:
        likelihoodDict[kinase] = 0
        for subKinase in pwmDictionary[kinase]:
            likelihoodDict[kinase]+=likelihoodOfKinase(bgRates, pwmDictionary[kinase][subKinase], peptideSequence)
    return likelihoodDict

fastaDB = import_FastaDB(r"Phosphoproteomics/data/CFImKDFastaDB.pickle")

def probabilities(pwmDictionary:dict, phosphosite:list, relativeIndexes:list):
    AccID = phosphosite[2]
    data = fastaDB[AccID]
    phosphositeSeq = slicer(phosphosite[0], data, relativeIndexes)     ### slice returns the peptide + the number of aa's on either side specified by relativeIndexes to enable "full prediction" of the site.
    print(phosphositeSeq, "here", phosphosite[0])
    #pwmDictionary is a dict of dicts, because the kinases have seperate pwms for the different phosphorylation targets.
    ### first off to calculate the bottom part of the probability //// basically the sum of the likelihoods multiplied by the prior distribution which we assume to be even
    likelihoodDict = likelihoodDictionary(pwmDictionary, phosphositeSeq, bg_rates)
    out = []
    p_k = 1/(len(likelihoodDict)+1)
    bottom = 0
    
    for kinase in likelihoodDict:
        bottom += likelihoodDict[kinase]*p_k
    for kinase in likelihoodDict:
        out.append([(likelihoodDict[kinase]*p_k)/bottom, kinase])

    return out


f = open(fr"Phosphoproteomics/steps/predictions/{pwms_filename}_predictions.txt", "a")

relativeIndexes = [i for i in range(-7,8)]

for phosphosite in phosphosites:
    print(phosphosite)
    AccID  = phosphosite[2]
    data=fastaDB[AccID]
    centre_idx = centre(phosphosite)

    updated_idx=[]
    for id in centre_idx:
        updated_idx.append(indexConverter(phosphosite[0],data,int(id)))
    
    for id in updated_idx:
        name_of_site = ""

        for offset in relativeIndexes:
            try:
                name_of_site += data[id+offset]
            
            except:
                print("An index was out of range", phosphosite[1])
                break  
        print(name_of_site)  
        f.write("{name};{predictions}\n".format(name = name_of_site, predictions = probabilities(pwms_loaded, phosphosite, relativeIndexes)))
    
f.close()