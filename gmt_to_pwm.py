from string import whitespace
import numpy as np
import pandas as pd
import pickle
from background import bg_rates

#Load PWMs from .gmt

with open(r"C:\Users\dirkj\Documents\Aarhus University\UniBas\Zavolab\Phosphoproteomics\Selbach_PWMs\iKiP-DB.gmt") as f:
    lines = [line for line in f]
f.close()

def underscores(gmt_list, index_1, name):
    #This function is supposed to catch _ in sequences and write sequence + uniprot ID to a file.
    #index 1 == index in gmt_list that corresponds to seq with _
    #index 2 == index in the seq that has a _ 
    uniprot_ids = gmt_list[1].split("|")
    out="Name:|" + name + "| uniprot ID:|"+uniprot_ids[index_1-1] + "| and sequence:|" + gmt_list[index_1].strip("-p;u")+"\n"
    f=open("underscores.txt","a")
    f.write(out)
    f.close()

def pwm_from_line(gmt_line):
    aas="ARNDCQEGHILKMFPSTWYV"
    out={}
    gmt_list = gmt_line.split()
    counter = 2

    for i in gmt_list[2:]:
        i=i.strip("-p;u")
        #Add PWM to output if it isn't in output yet.
        nam = gmt_list[0][12:]+"-"+i[7]
        if nam not in out:
            out[nam]=np.zeros((len(aas),len(i)))
        
        #Count AA freqs in each position
        for aa in range(len(aas)):
            switch =0
            for n in range(len(i)):

                if aas[aa]==i[n]:
                    #This if statement is just here to catch undefined AA's and record them so we can check uniprot for the seq.
                    if aa == 20 and switch == 0:
                        switch =1
                        underscores(gmt_list,counter, nam)
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

pwms = {}

#print(pwm_from_line(lines[0]))

for i in lines:
    pwms[i.split()[0][12:]+"-"+i[7]] = pwm_from_line(i)

with open('selbach_pwms.pickle', 'wb') as f:
    pickle.dump(pwms, f)

with open('selbach_pwms.pickle', 'rb') as f:
    pwms_loaded = pickle.load(f)

print(len(pwms_loaded.keys()))