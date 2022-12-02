import numpy as np
import requests as r
import pickle
import pandas as pd
import os
#Goal is to create a PWM constructor that relies on relative indexing from the phosphosite.

os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab") #insert your own pathname here.

with open(r"Phosphoproteomics/data/Sugiyama_phosphosites_updated.txt") as f:
    lines = [line.strip("\n").split("\t") for line in f]
f.close()

# def parseKinaseSubstrateDS(filename):
#     with open(filename) as f:
#         df = pd.read_table(f, header=2)
#     f.close()
#     return df

# kin_sub_DS = parseKinaseSubstrateDS(r"Phosphoproteomics/data/kin_sub_ds.txt")

# kin_sub_DS_list = [[kinase_name, [pps.upper() for pps in kin_sub_DS[kin_sub_DS['KINASE'].str.match(kinase_name, na=False)]['SITE_+/-7_AA']]] for kinase_name in kin_sub_DS['KINASE']]

# print(kin_sub_DS_list[:10])

relativeIndex = [i for i in range(-7,8)]


def accessionIdToUrl(Acc):
    currentUrl = "https://www.uniprot.org/uniprot/"+Acc+".fasta"
    return currentUrl
    
def pwmFromList(phosphositeList):
    ### Input should be a list of length 2
    ### index 0 being the name of the kinase  
    ### index 1 being a list of phosphosites.
    aas="ARNDCQEGHILKMFPSTWYV"
    
    seqList = phosphositeList[1]

    out = {}
    for phosphosite in seqList:
        name = phosphositeList[0] + "-" + phosphosite[7]
        if name not in out:
            out[name]=np.zeros((len(aas),len(phosphosite)))
        
        for id_1, aa in enumerate(aas):
            for id_2, aa_pps in enumerate(phosphosite):
                if aa==aa_pps:
                    out[name][id_1,id_2]+=1

    #Convert freqs to rates
    for table in out:
        sums = np.sum(out[table],axis=0)
        #print("sums =",sums)
        for i in range(len(sums)):
            #print("here",len(sums))
            for aa in range(len(aas)):
                if out[table][aa,i] != 0:
                    out[table][aa,i]=out[table][aa,i]/sums[i]    
    return out


pwms = {}

for i in lines:
    i[1] = i[1].strip('"').split(";")
    pwms[i[0]] = pwmFromList(i)

# for i in kin_sub_DS_list:

#     pwms[i[0].upper()] = pwmFromList(i)

output_filename = "sugiyama_pwms"

with open(fr'Phosphoproteomics/steps/pwms/{output_filename}.pickle', 'wb') as f:
    pickle.dump(pwms, f)
f.close()
with open(fr'Phosphoproteomics/steps/pwms/{output_filename}.pickle', 'rb') as f:
    pwms_loaded = pickle.load(f)
f.close()
print(len(pwms_loaded.keys()), pwms_loaded.keys())



#Function below is outdated, made for use with .gmt files
def pwm_from_line(gmt_line):
    aas="ARNDCQEGHILKMFPSTWYV"
    out={}
    gmt_list = gmt_line.split()
    #print(gmt_line)
    #print(gmt_list)

    seqList = []

    acc_id_list = gmt_list[1].split("|")
    #print(acc_id_list)
    for i in acc_id_list[1:]:
        seq = ""
        li = i.split('.')
        seqUrl = accessionIdToUrl(li[0])
        response = r.post(seqUrl)
        data=''.join([str(lines, 'utf-8') for lines in response.iter_lines()][1:])
        #print(data[int(li[2])-1], li)
        for n in relativeIndex:
            try:
                seq+= data[int(li[2])-1+n]
            except:
                break
        seqList.append(seq)

    print(seqList)
    counter = 2

    for i in seqList:
        #Add PWM to output if it isn't in output yet.
        nam = gmt_list[0][12:]+"-"+i[7]
        if nam not in out:
            out[nam]=np.zeros((len(aas),len(i)))
        
        #Count AA freqs in each position
        for aa in range(len(aas)):
            switch =0
            for n in range(len(i)):

                if aas[aa]==i[n]:
                    out[nam][aa,n]+=1
        counter+=1
   
    #Convert freqs to rates
    for table in out:
        sums = np.sum(out[table],axis=0)
        #print("sums =",sums)
        for i in range(len(sums)):
            #print("here",len(sums))
            for aa in range(len(aas)):
                if out[table][aa,i] != 0:
                    out[table][aa,i]=out[table][aa,i]/sums[i]

    return out

