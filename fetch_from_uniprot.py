from ast import Break
from io import StringIO
from urllib import response
import requests as r
from Bio import SeqIO

###This program is built to load sequence data from uniprot and try and fill in blanks in sequences, but it turned out that usually the blanks are just there because of sequence ends. womp womp
### bad luck

def id_to_url(id):
    currentUrl = "https://www.uniprot.org/uniprot/"+id+".fasta"
    return currentUrl
    
def exact_seq_finder(seq, seq_fragment):
    None

def fetcher(id, seq_fragment):
    aas="ARNDCQEGHILKMFPSTWYV_"

    response = r.post(id_to_url(id))
    data=''.join(response.text).replace("\n","")
    print(data, response)
    Seq = ''
    out_idxs = []
    for n in range(len(seq_fragment)):
        if seq_fragment[n] == "_":
            out_idxs.append(n)
    for i in range(len(data)):
        out = ''
        n=0
        for n in range(len(seq_fragment)):
            out+=data[i+n]
            if n not in out_idxs:
                if data[i+n] != seq_fragment[n]:
                    break
        if len(out) == len(seq_fragment):
            break
    for i in out_idxs:
        if out[i] not in aas:
            out[i]=="_"
        #Seq += data[i]
    print(out, "yeha", out)


fetcher("P05387", "_MRYVASYLLAALGG")