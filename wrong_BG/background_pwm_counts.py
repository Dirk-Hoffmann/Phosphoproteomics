import pandas as pd
import numpy as np
import requests as r
import pickle
from fetch_from_uniprot import id_to_url
from kinase_prediction import centre


### we want to construct at background PWM based on the average rate of all AA's in our phosphosites.
### I'm trying to write this function to be adaptable to a number of intervals, with the intention of experimenting on what the ideal amount of surrounding peptides to look at is (terrible phrasing but you get it)

df = pd.read_csv(r"C:\Users\dirkj\Documents\Aarhus University\UniBas\Zavolab\Phosphoproteomics\CFIm_KD_Phospho_Peptides.csv", sep=";")

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


#Indexes to count from uniprot sequence
relative_idx = [i for i in range(-6,7)]
aas="ARNDCQEGHILKMFPSTWYV"
bg_pwms = {}
zero_matrix = np.zeros((len(aas),len(relative_idx)))

counter = 0
listen = df.values.tolist()

for i in listen:
    print(counter,"/",len(listen))
    centre_idx = centre(i)
    seqUrl = id_to_url(i[2])
    response = r.post(seqUrl)
    data=''.join(response.text).replace("\n","")
    updated_idx=[]
    for id in centre_idx:
        updated_idx.append(id_converter(i[0],data,int(id)))

    for id in updated_idx:
        if data[id] not in bg_pwms:
            bg_pwms[data[id]]= zero_matrix
        for offset in range(len(relative_idx)):
            for aa in range(len(aas)):
                try:
                    if data[id+relative_idx[offset]]==aas[aa]:
                        bg_pwms[data[id]][aa,offset]+=1
                        break
                except:
                    print("An index was out of range", i[1])
                    break
    counter +=1

print(bg_pwms)


with open('background_pwms.pickle', 'wb') as f:
    pickle.dump(bg_pwms, f)