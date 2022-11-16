import ast
import enum
import numpy as np
import pandas as pd
import cProfile

with open('predictions.txt') as f:
    lines = f.readlines()
f.close()

prot_to_ACC = pd.read_csv(r"C:\Users\dirkj\Documents\Aarhus University\UniBas\Zavolab\Phosphoproteomics\predictions\bbabab (version 1).xlsb.csv", sep=";")
protlist = prot_to_ACC["Prot"].tolist()
acclist = prot_to_ACC["Acc"].tolist()
Prot_ACC_dict = {}
for count, n in enumerate(protlist):
    Prot_ACC_dict[n] = acclist[count]


#lines = lines.split(";") 
#print(L)
def parse_predictions(prediction_filename):
    with open(prediction_filename) as f:
        lines = f.readlines()
    f.close()
    for id, line in enumerate(lines):
        lines[id] = line.split(";")
        lines[id][1] = ast.literal_eval(lines[id][1])
    return lines




def parse_kinase_substrate(filename):
    with open(filename) as f:
        df = pd.read_table(f, header=2)
    f.close()
    return df
    print(df[df['KIN_ACC_ID'].str.match('Q05655', na=False)])

def evaluation(predictions, dataset):
    preds = parse_predictions(predictions)
    data = parse_kinase_substrate(dataset)

    #mismatches = {}

    ratios = [0,0]
    #jabadabadu = {}

    data['SITE_+/-7_AA']=data['SITE_+/-7_AA'].str.upper()
    data['SITE_+/-7_AA']=data['SITE_+/-7_AA'].str.strip("_")
    for line in preds:
        #print(line[0])
        top_preds = sorted(line[1])[::-1]
        kinase_list = [] #### LIST OF KNOWN KINASES TARGETING THIS SUBSTRATE
        temp_ =[] ### DOES NOT MATTER.
        for n in data[data['SITE_+/-7_AA'].str.match(line[0],na=False)]['KIN_ACC_ID']:
            if n not in kinase_list:
                kinase_list.append(n)
     #   for i in data[data['SUB_ACC_ID'].str.match(line[0],na=False)]['KIN_ACC_ID']:
     #       if i not in temp_:
     #           temp_.append(i)
     #   for n, o in enumerate(kinase_list):
     #       jabadabadu[o]=temp_[n]
        for count, kinase in enumerate(kinase_list):
            top_preds[count][1]=top_preds[count][1].replace("-i","")
            ratios[1]+=1 
            for n in kinase_list:
     #           if n not in jabadabadu:
     #               jabadabadu[n]={}
                if n in Prot_ACC_dict[top_preds[count][1]]:
                    ratios[0]+=1
                else:
                    print(n, Prot_ACC_dict[top_preds[count][1]])

                # else:
                #     if top_preds[count][1] not in jabadabadu[n]:
                #         jabadabadu[n][top_preds[count][1]] = 1
                #     else:
                #         jabadabadu[n][top_preds[count][1]]+=1
               
            # else:
            #     print(top_preds[count], kinase_list)
            #     mismatches[top_preds[count][1]]= kinase_list
    #print(mismatches)
    print(ratios)

    # for n in jabadabadu:
    #     print(n, jabadabadu[n])

#[191, 14472]

#parse_kinase_substrate('Kinase_Substrate_Dataset')
cProfile.run('evaluation("predictions.txt", "Kinase_Substrate_Dataset")')
#evaluation("predictions.txt", "Kinase_Substrate_Dataset")

#L[1] = ast.literal_eval(L[1])

