import matplotlib.pyplot as plt
import pickle
import pandas as pd
import os
import seaborn as sns

sns.set()

os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab") #insert your own pathname here.

def unpickle(relative_path):
    with open(fr'{relative_path}', 'rb') as f:
        out = pickle.load(f)
    f.close()
    return out
    
cross_entropy = unpickle("Phosphoproteomics/steps/sugiyama_ppsplus_cross_entropy.pickle")
CE_updated = []

for list in cross_entropy:
    if len(list) == 4:
        list = [float(list[0]),[list[1],list[2]],list[3]]
    for locale in list[2]:
        for index in locale.split(","):
            list = [list[0], list[1], index.strip().upper()]
            CE_updated.append(list)


df = pd.DataFrame(CE_updated, columns=['Cross Entropy Loss', 'Kinase Name', 'Subcellular Localization'])

# with open(fr"Phosphoproteomics/plots/temp_files/sugiyama_ppsplus_CE_df.pickle", "wb") as f:
#     pickle.dump(df, f)
# f.close()

# print(df)

df_wide = unpickle("Phosphoproteomics/plots/temp_files/sugiyama_ppsplus_CE_df.pickle")

print(df)



#### first off we produce a histogram of the cross entropy loss:


df_no_subcellular = df_wide.loc[df["Subcellular Localization"]==set()]

df_only_subcellular = df_wide.loc[df["Subcellular Localization"]!=set()]

# print(df_no_subcellular)

#df_wide = df.pivot(columns='Subcellular Localization')

histogram = df_wide.plot.hist(y='Cross Entropy Loss', bins=22 , title = "Full Histogram", legend = False)
histogram_no_subcellular = df_no_subcellular.plot.hist(y='Cross Entropy Loss', bins=22 , title = "Histogram No subcellular")
histogram_only_subcellular = df_only_subcellular.plot.hist(y='Cross Entropy Loss', bins=22 , title = "Histogram Only subcellular")

### Calculate median CE_loss for each subcellular localization

# subcellular_locales = dict()

# for set in df['Subcellular Localization'].values.tolist():
#     for locale in set:
#         if locale not in subcellular_locales:
#             subcellular_locales[locale]=1
#         else:
#             subcellular_locales[locale]+=1

# print(df.melt(id_vars = ('Kinase Name', 'Cross Entropy Loss'), value_vars = 'Subcellular Localization'))

medians = df.groupby('Subcellular Localization').median()
print(df.groupby('Subcellular Localization').count())
medians['Count'] = df.groupby('Subcellular Localization').count()['Kinase Name']


#medians.to_excel(excel_writer='Phosphoproteomics/plots/medians_per_subcel_location.xlsx')

print(medians.sort_values(by='Cross Entropy Loss'))

plt.show()


