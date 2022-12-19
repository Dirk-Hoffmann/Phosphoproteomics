import matplotlib.pyplot as plt
import pickle
import pandas as pd
import os
import seaborn as sns
from scipy import stats

sns.set()

os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab") #insert your own pathname here.

def unpickle(relative_path):
    with open(fr'{relative_path}', 'rb') as f:
        out = pickle.load(f)
    f.close()
    return out
    
cross_entropy = unpickle("Phosphoproteomics/steps/sugiyama_ppsplus_cross_entropy.pickle")
CE_updated = []
CE_updated_wide = []


for list in cross_entropy:
    if len(list) == 4:
        list = [float(list[0]),list[2],list[3]]
        CE_updated_wide.append(list)
    for locale in list[2]:
        for index in locale.split(","):
            list = [list[0], list[1], index.strip().upper()]
            CE_updated.append(list)


df = pd.DataFrame(CE_updated, columns=['Cross Entropy Loss', 'Kinase', 'Subcellular Localization'])

# with open(fr"Phosphoproteomics/plots/temp_files/sugiyama_ppsplus_CE_df.pickle", "wb") as f:
#     pickle.dump(df, f)
# f.close()

# print(df)

df_wide = pd.DataFrame(CE_updated_wide, columns=['Cross Entropy Loss', 'Kinase', 'Subcellular Localization'])

# print(df)
# print(df_wide)


#### first off we produce a histogram of the cross entropy loss:


df_no_subcellular = df_wide.loc[df["Subcellular Localization"]==set()]

df_only_subcellular = df_wide.loc[df["Subcellular Localization"]!=set()]

# print(df_no_subcellular)

#df_wide = df.pivot(columns='Subcellular Localization')

# histogram = df_wide.plot.hist(y='Cross Entropy Loss', bins=22 , title = "Full Histogram", legend = False)
# histogram_no_subcellular = df_no_subcellular.plot.hist(y='Cross Entropy Loss', bins=22 , title = "Histogram No subcellular")
# histogram_only_subcellular = df_only_subcellular.plot.hist(y='Cross Entropy Loss', bins=22 , title = "Histogram Only subcellular")

### Calculate median CE_loss for each subcellular localization

# subcellular_locales = dict()

# for set in df['Subcellular Localization'].values.tolist():
#     for locale in set:
#         if locale not in subcellular_locales:
#             subcellular_locales[locale]=1
#         else:
#             subcellular_locales[locale]+=1

# print(df.melt(id_vars = ('Kinase', 'Cross Entropy Loss'), value_vars = 'Subcellular Localization'))



medians = df.groupby('Subcellular Localization').median()
print(df.groupby('Subcellular Localization').count())
medians['Count'] = df.groupby('Subcellular Localization').count()['Kinase']


medians.to_excel(excel_writer='Phosphoproteomics/plots/medians_per_subcel_location.xlsx')


print(medians.sort_values(by='Cross Entropy Loss'))

filtered_df = df.groupby('Subcellular Localization').filter(lambda x : len(x)> 4)


index_sort = filtered_df.groupby('Subcellular Localization')['Cross Entropy Loss'].median().sort_values().index

print(index_sort.values)

#df_sorted = filtered_df[index_sort]

print(df_wide["Cross Entropy Loss"])

print(df_wide)

# l = stats.ranksums(x = df_wide["Cross Entropy Loss"], y=filtered_df.loc[filtered_df["Subcellular Localization"]=="Cytoplasm"]["Cross Entropy Loss"])

print(df[df['Kinase']=='CDK2'])

df.groupby('Kinase').count().to_csv('grouped.csv')



print("HERE",filtered_df[filtered_df['Subcellular Localization']!= 'CYTOPLASM'].drop_duplicates(subset='Kinase'))

ranksums = dict()

######## The rank sums or Mann-Whitney U-test calculates whether or not the distribution of cross entropy scores for a specific subcellular localization is significantly different
######## from the rest of the measurements. As this test only assumes that the samples are unrelated we should be in the clear.
######## output of for loop is dict with key:value (subcellular localization: pvalue)

sns.set(rc={"figure.figsize":(15, 7)})



for subcellular_loc in index_sort.values:
    os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab/Phosphoproteomics/plots/distribution_kdeplots")

    temp_df = filtered_df[filtered_df['Subcellular Localization']== subcellular_loc].drop_duplicates(subset='Kinase')

    other_df = filtered_df[~filtered_df['Kinase'].isin(temp_df['Kinase'])].drop_duplicates(subset='Kinase')

    temp_df['Localization'] = subcellular_loc
    other_df['Localization'] = "Not "+subcellular_loc

    # print("temp_df",temp_df)
    # print("other_df",other_df)
    # print(len(other_df)+len(temp_df))
    # for i in temp_df['Kinase']:
    #     print(i)
    #     other_df['Kinase'].drop(i, axis=0)
    #     print(other_df['Kinase'])

    temp_df['hue']=(temp_df['Subcellular Localization']==subcellular_loc).map({False:'Not {}'.format(subcellular_loc), True:subcellular_loc})
    ranksums[subcellular_loc] = stats.ranksums(x = other_df['Cross Entropy Loss'], y=temp_df["Cross Entropy Loss"]).pvalue
    print(subcellular_loc,ranksums[subcellular_loc])
    joined = pd.concat([temp_df,other_df])
    if ranksums[subcellular_loc] <= 0.1:
            os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab/Phosphoproteomics/plots/distribution_kdeplots/significant")
    # s = sns.FacetGrid(data=joined, row= 'Localization', aspect=3)
    # s.map(sns.kdeplot, 'Cross Entropy Loss')
    # s.set(title='{subloc} - p-value: {p_value}'.format(subloc = subcellular_loc, p_value= ranksums[subcellular_loc]))
    sns.kdeplot(data=joined, x = 'Cross Entropy Loss', hue='Localization', clip=0).set(title='Subcellular Localization: {subloc} \n p-value: {p_value} \n n: {n}'.format(subloc = subcellular_loc, p_value= ranksums[subcellular_loc], n=len(temp_df['Kinase'])))
    plt.savefig('{}_kdeplot.pdf'.format(subcellular_loc))
    plt.clf()

##############

# sns.set(rc={"figure.figsize":(28, 15)})

# subcel_boxplot = sns.boxplot(data= filtered_df, order=index_sort.values, y = 'Subcellular Localization', x = "Cross Entropy Loss")
# subcel_boxplot = sns.swarmplot(data= filtered_df, order=index_sort.values, y = 'Subcellular Localization', x = "Cross Entropy Loss", alpha = 0.5, linewidth=1)

# plt.savefig('subcel_boxplot.pdf')
# plt.clf()

# sns.kdeplot(data=filtered_df,x= 'Cross Entropy Loss', hue="Subcellular Localization", cumulative=True)

# plt.show()

#################