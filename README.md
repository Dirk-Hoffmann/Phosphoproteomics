# Phosphoproteomics /// Sarcopenia
 Git repository of my project with phosphoproteomics data for Zavolan Lab anno 2022

First off - the structure of this repository is supposed to be more or less intuitive but I'll explain it here regardless.

//data

    This folder stores raw data. 
    
    //CFIm_KD_Phospho_peptides.csv 
    
        Phosphosites from the Zavolan Lab Sarcopenia project.
        
    //CFImKDFastaDB.pickle
    
        Pickle of a dictionary with {Accession ID : FASTA} for each phosphopeptide in the phosphosites.
        
    //kin_sub_ds.txt.pickle
    
        The Kinase - Substrate Dataset from Phosphosite+ (from November 2022)
        
    //Sugiyama_phosphosites.csv 
    
        Unfiltered in vitro kinase-substrate data from Sugiyama et. al. 
        
    //Sugiyama_phosphosites_updated.tx
    
        Sugiyama phosphosites, that have had their nomenclature updated so they fit the phosphosite+ Kinase-substrate dataset

//plots

    A folder made for storing and generating plots from data. This folder has some messy code and some plots and is in general more there to segment those away from the rest of the code. I might eventually clean it up, but in the end that's not really the purpose of this folder.

//results

     contains final results/plots

//scripts

    contains code in general

//steps

    contains intermediate results (like pwms or predictions)





Phosphoproteomics
For this project, my goal was to calculate position-weight-matrices for 385 human kinases (354 wt, 21 mutant, 10 lipid) as described by Sugiyama et al. 2019. These updated PWMs are used to predict kinases that act on riboseq phosphosite data from aging mice.

folder structure should be intuitive but here goes. 

data// 
 contains raw/semi-raw data

plots//
 contains plots + code for producing plots (more messy in a lot of cases)

results//
 contains final results/plots

scripts//
 contains code in general

steps//
 contains intermediate results (like pwms or predictions)
 
//Packages: versions

    Python: 3.9.12 (main, Apr 5 2022, 01:53:17)

    requests: 2.27.1

    numpy: 1.21.5

    pandas: 1.4.2

    pylcs: 0.0.7
