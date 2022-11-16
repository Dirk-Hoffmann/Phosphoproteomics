import requests as r

def id_to_url(id):
    currentUrl = "https://www.uniprot.org/uniprot/"+id+".fasta"
    return currentUrl
    

with open(r"C:\Users\dirkj\Documents\Aarhus University\UniBas\Zavolab\Phosphoproteomics\CFIm_KD_Phospho_Peptides.csv") as f:
    lines = [line for line in f]
f.close()

#print(lines[0])

def calculate_background(phosphopeptides):
    aas="ARNDCQEGHILKMFPSTWYV"
    bound = 6 ### Use to define the length in either direction that we need to look in uniprot
    bg_list = [0 for i in aas]

    for line in phosphopeptides:
        line = line.split(';')
        response = r.post(id_to_url(line[2].replace("\n","")))
        data=''.join(response.text).replace("\n","")
        id = 0
        for count,n in enumerate(data):
            seq = False
            f = len(line[0])
            # print(count, f)
            # print(data[count:count+f])
            if data[count:count+len(line[0])] == line[0]:
                seq = data[count-bound:count+len(line[0])+bound]
                print("got'em", seq, line[0])
                break
        if seq:
            for n in seq:
                for count, i in enumerate(aas):
                    if n == i:
                        bg_list[count]+= 1
        print(bg_list)
    return bg_list

        #print(line)
    #print(bg_list)
    #print(phosphopeptides[1])



#calculate_background(lines) ## output below:

background_list = [17510, 21996, 6196, 10363, 1533, 10499, 17847, 16042, 4666, 5965, 14455, 21018, 2174, 4161, 25658, 35493, 14478, 887, 3201, 10704]

bg_rates = []

for n in background_list:
    bg_rates.append(n/sum(background_list))

#print(bg_rates, sum(bg_rates))