import pandas as pd


with open(r"/Users/dirk/Documents/UniBas/Zavolab/Phosphoproteomics/Selbach_PWMs/Sugiyama_phosphosites.txt") as f:
    lines =[line.split("\t")[0] for line in f.readlines()]
f.close()

prot_to_ACC = pd.read_csv(r"/Users/dirk/Documents/UniBas/Zavolab/Phosphoproteomics/predictions/bbabab (version 1).xlsb.csv", sep=";")
protlist = prot_to_ACC["Prot"].tolist()
acclist = prot_to_ACC["Acc"].tolist()
Prot_ACC_dict = {}
for count, n in enumerate(protlist):
    Prot_ACC_dict[n] = acclist[count]

lmao = {}

for n in lines:
    print(n)
    if n in Prot_ACC_dict:
        lmao[n]=Prot_ACC_dict[n]
        print(n, Prot_ACC_dict[n])
    else:
        lmao[n]=None

ff = pd.DataFrame.from_dict([lmao])

ff = pd.melt(ff)

print(ff)


ff.to_excel('out.xlsx', index=False)