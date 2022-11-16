import pandas as pd
import pickle
import numpy as np
import math
import requests as r
from background import bg_rates


df = pd.read_csv(r"C:\Users\dirkj\Documents\Aarhus University\UniBas\Zavolab\Phosphoproteomics\CFIm_KD_Phospho_Peptides.csv", sep=";")

def id_to_url(id):
    currentUrl = "https://www.uniprot.org/uniprot/"+id+".fasta"
    return currentUrl
    

phosphosites = df.values.tolist()

def id_converter(seq_1, seq_2, id_1):
    #converts index of phosphosite to index in sequence
    for i in range(len(seq_2)):
        x = 0
        for n in range(len(seq_1)):
            x+=1
            if seq_2[i+n]!=seq_1[n]:
                break
        if x == len(seq_1):
            return id_1+i

with open('selbach_pwms.pickle', 'rb') as f:
    pwms_loaded = pickle.load(f) 
f.close()

def centre(phosphosite):
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

relative_idx = [i for i in range(-7,8) if i!= 0]
aas="ARNDCQEGHILKMFPSTWYV"

#{key: value for (key, value) in iterable}
aa_to_idx = {key:value for (key, value) in zip(aas, range(len(aas)))}

def predict_kinase(pwm_dict, phosphosite):
    seqUrl = id_to_url(phosphosite[2])
    response = r.post(seqUrl)
    data=''.join(response.text).replace("\n","")
    centre_idx = centre(phosphosite)
    updated_idx=[]
    for id in centre_idx:
        updated_idx.append(id_converter(phosphosite[0],data,int(id)))

    probs = []

    for kinase in pwm_dict:
        for target in pwm_dict[kinase]:
            for n in updated_idx:
                check = 0
                cumprob = 0
                for id in range(len(relative_idx)):
                    indent = relative_idx[id]+7
                    AA = aa_to_idx[str(data[n+relative_idx[id]])]

                        ### we want to catch cases, where the probability would be set  to 0 and just remove them.
                    if pwm_dict[kinase][target][AA,indent] == 0:
                        cumprob += math.log(0.000000000001)
                        break
                    cumprob+=math.log(pwm_dict[kinase][target][AA,indent])
            if check != 1:
                probs.append([math.e**cumprob,target])
    out=sorted(probs)[::-1][0:10]
    print(out)
    return out

#predict_kinase(pwms_loaded, phosphosites) 

def likelihood_of_kinase(background, pwm, peptide):
    peptide_seq = peptide
    l_k = pwm.shape[1]
    l_i = len(peptide_seq)

    out = 0

    ### This catches cases where the phosphosite is at the very end of a peptide.
    if l_i-l_k <= 0:
        subsum = 0
        for n in range(len(peptide_seq)):
            if pwm[aa_to_idx[peptide_seq[n]], n] == 0:
                subsum+= math.log(0.1*10**-100)
            else:
                subsum += math.log(pwm[aa_to_idx[peptide_seq[n]], n])
        out += subsum

    for j in range(l_i-l_k):
        
        subsum = 0
        ### part 1
        for n in range(0, j-1):
            subsum += math.log(background[aa_to_idx[peptide_seq[n]]])
        
        ### part 2
        for n in range(j, j+l_k-1):
            if pwm[aa_to_idx[peptide_seq[n]], n-j] == 0:
                subsum+= math.log(0.1*10**-100)
            else:
                subsum += math.log(pwm[aa_to_idx[peptide_seq[n]], n-j])
        ### part 3
        for n in range(j+l_k,l_i-1):
            subsum += math.log(background[aa_to_idx[peptide_seq[n]]])
        out+= math.e**subsum 
    return out




#likelihood_of_kinase(bg_rates, pwms_loaded["ZAK-i"]["ZAK-S"], phosphosites)

def probabilities(pwm_dict, phosphosite):
    seqUrl = id_to_url(phosphosite[2])
    response = r.post(seqUrl)
    data=''.join([str(lines, 'utf-8') for lines in response.iter_lines()][1:])
    ### slice below returns the peptide + the number of aa's on either side specified by relative_idx to enable "full prediction" of the site.
    if id_converter(phosphosite[0],data,0)+relative_idx[0] < 0:
        left = 0
    else: 
        left = id_converter(phosphosite[0],data,0)+relative_idx[0]
    if id_converter(phosphosite[0],data,len(phosphosite[0]))+relative_idx[-1] > len(data):
        right = len(data)-1
        #print(data)
    else: 
        right = id_converter(phosphosite[0],data,len(phosphosite[0]))+relative_idx[-1]

    phosphosite_seq = data[left:right]
    #pwm_dict is a dict of dicts, because the kinases have seperate pwms for the different phosphorylation targets.
    #honestly idk about this ### first off to calculate the bottom part of the probability //// basically the sum of the likelihoods multiplied by the prior distribution which we assume to be even
    ### LETS GO construct a dict of likelihoods for the kinases.
    likelihood_dict = {}
    out = []
    p_k = 1/(len(pwm_dict)+1)

    for kinase in pwm_dict:
        likelihood_dict[kinase]={}
        for subkinase in pwm_dict[kinase]:
            likelihood_dict[kinase][subkinase]=likelihood_of_kinase(bg_rates, pwm_dict[kinase][subkinase], phosphosite_seq)
    
    #uncomment the stuff above when you're ready, for now its nice to have the dict calculated

    bottom = 0

    for kinase in likelihood_dict:
        for subkinase in likelihood_dict[kinase]:
            #print(math.e**likelihood_dict[kinase][subkinase], math.e)
            bottom += likelihood_dict[kinase][subkinase]*p_k
    for kinase in likelihood_dict:
        for subkinase in likelihood_dict[kinase]:
            out.append([(likelihood_dict[kinase][subkinase]*p_k)/bottom, kinase])
    

    #print(bottom)
    #print(likelihood_dict)
    #out.sort()
    #print(out[::-1][0:10])
    return out


#probabilities(pwms_loaded, phosphosites[29])



f = open("predictions_2.txt", "a")

relative_idx = [i for i in range(-7,8)]


for site in phosphosites:
    print(site)
    seqUrl = id_to_url(site[2])
    response = r.post(seqUrl)
    data=''.join([str(lines, 'utf-8') for lines in response.iter_lines()][1:])
    centre_idx = centre(site)

    updated_idx=[]
    for id in centre_idx:
        updated_idx.append(id_converter(site[0],data,int(id)))
    
    for id in updated_idx:
        name_of_site = ""

        for offset in relative_idx:
            try:
                name_of_site += data[id+offset]
            
            except:
                print("An index was out of range", site[1])
                break  
        print(name_of_site)  
        f.write("{name};{predictions}\n".format(name = name_of_site, predictions = probabilities(pwms_loaded, site)))
    
f.close()