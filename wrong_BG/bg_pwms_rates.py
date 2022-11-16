import pickle
from re import X
import numpy as np

with open('background_pwms.pickle', 'rb') as f:
    counts = pickle.load(f)
#f.close()

print(counts)

aas="ARNDCQEGHILKMFPSTWYV"

#Separated site-specific bg_calculation

# for table in counts:
#     for i in range(13):
#         for aa in range(len(aas)):
#             counts[table][aa,i]=counts[table][aa,i]/sum(counts[table][i:len(aas),0])


### Actual BG calculation eg. just the chance of any amino acid to appear in that position.


for n in counts:
    counts[n]=np.delete(arr=counts[n], obj=6, axis=1)

bg = []

for i in range(len(aas)):
    count= 0
    for array in counts:
        print(array)
        for n in np.nditer(counts[array][i]):
            #print(n)
            count += n
    bg.append(count)

print(bg)

bg_rates = []

for i in bg:
    bg_rates.append(i/sum(bg))

print(bg_rates)