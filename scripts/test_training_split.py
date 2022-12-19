import numpy as np
import pandas as pd
import pickle
from sklearn.model_selection import StratifiedShuffleSplit
import os

os.chdir(r"/Users/dirk/Documents/UniBas/Zavolab") #insert your own pathname here.

### This script is meant to split data into test/training sets whilst keeping the class proportions consistent between the splits.
### This is called stratified splitting, and I'm using the sklearn model_selection package to do so
### I want to take one pandas dataframe as as input and output two data frames called training and testing

def import_dataframe_pickle(path):
    with open(path, "rb") as f:
        out = pickle.load(f)
    f.close
    return out

kin_sub_ds = import_dataframe_pickle("Phosphoproteomics/data/kin_sub_ds.txt.pickle")

kin_sub_ds_filtered = kin_sub_ds[kin_sub_ds['KINASE'].map(kin_sub_ds['KINASE'].value_counts())> 5]

print(kin_sub_ds_filtered)

def stratified_shuffle_split(dataset, ratio):
    sss = StratifiedShuffleSplit(n_splits=1, test_size=ratio, random_state= 0)
    
    for i, (train_index, test_index) in enumerate(sss.split(dataset['SITE_+/-7_AA'],dataset['KINASE'])):
        print(f"Fold {i}:")
        print(f"  Train: index={train_index}")
        print(f"  Test:  index={test_index}")


stratified_shuffle_split(kin_sub_ds_filtered, .20)