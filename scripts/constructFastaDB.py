import requests as r
import pickle
import os

os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab") #insert your own pathname here.


### This program downloads a list of FASTA files and constructs a dictionary for pulling sequences based on Acc-ID ###
### Used to reduce the overhead on a lot of the calculations in this project, as posting requests to uniprot takes quite a bit of time ###

def import_FastaDB(path:str):
    with open(path, "rb") as f:
        fastaDB = pickle.load(f)
    f.close()
    return fastaDB

def idToUrl(id):
    currentUrl = "https://www.uniprot.org/uniprot/"+id+".fasta"
    return currentUrl

with open(r"Phosphoproteomics/data/CFIm_KD_Phospho_Peptides.csv") as f:
    lines = f.readlines()
f.close

def accToSeq(phosphosite):
    seqUrl = idToUrl(phosphosite.split(";")[2].strip("\n"))
    response = r.post(seqUrl)
    data=''.join([str(lines, 'utf-8') for lines in response.iter_lines()][1:])
    return data

fastaDB = dict()

for n, line in enumerate(lines):
    AccID = line.split(";")[2].strip("\n")
    fastaDB[AccID] = accToSeq(line)
    print(AccID," ",fastaDB[AccID])

with open(r"Phosphoproteomics/data/CFImKDFastaDB.pickle", "wb") as f:
    pickle.dump(fastaDB, f)
f.close()

