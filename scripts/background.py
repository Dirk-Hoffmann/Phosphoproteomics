import requests as r
import os

os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab") #insert your own pathname here.

### This script calculates the background rates of each amino acid across a set of phosphosites.
### The phosphosites have to be extended by 7 amino acids upstream and downstream, as this is the 
### peptide fragment which is used for the sliding window likelihood calculation which we calculate the background for.

def id_to_url(id):
    currentUrl = "https://www.uniprot.org/uniprot/"+id+".fasta"
    return currentUrl
    
with open(r"Phosphoproteomics/data/CFIm_KD_Phospho_Peptides.csv") as f:
    lines = [line for line in f]
f.close()

def calculate_background(phosphopeptides):
    aas="ARNDCQEGHILKMFPSTWYV"
    bound = 7 ### Use to define the length in either direction that we need to look in uniprot
    bg_list = [0 for i in aas]

    for phosphosite in phosphopeptides:
        phosphosite = phosphosite.split(';')
        phosphosite_sequence = phosphosite[0]
        response = r.post(id_to_url(phosphosite[2].replace("\n","")))
        data=''.join(response.text).replace("\n","")
        for count,n in enumerate(data):
            if data[count:count+len(phosphosite_sequence)] == phosphosite_sequence:
                phosphosite_with_bounds = data[count-bound:count+len(phosphosite_sequence)+bound]
                break
        if phosphosite_with_bounds:
            for aa in phosphosite_with_bounds:
                for count, i in enumerate(aas):
                    if aa == i:
                        bg_list[count]+= 1
    return bg_list

#calculate_background(lines) ## output below:

background_list = [17510, 21996, 6196, 10363, 1533, 10499, 17847, 16042, 4666, 5965, 14455, 21018, 2174, 4161, 25658, 35493, 14478, 887, 3201, 10704]

bg_rates = []

for n in background_list:
    bg_rates.append(n/sum(background_list))
