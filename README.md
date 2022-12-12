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
    //Sugiyama_phosphosites_updated.txt
        Sugiyama phosphosites, that have had their nomenclature updated so they fit the phosphosite+ Kinase-substrate dataset

//plots
    A folder made for storing and generating plots from data. This folder has some messy code and some plots and is in general more there to segment those away from the rest of the code. I might eventually clean it up, but in the end that's not really the purpose of this folder.

//results
    





Phosphoproteomics
For this project, my goal was to calculate position-weight-matrices for 385 human kinases (354 wt, 21 mutant, 10 lipid) as described by Sugiyama et al. 2019. These updated PWMs are used to predict kinases that act on riboseq phosphosite data from aging mice.

We are still working on recovering the PWMs described above, in the meantime I’m working with PWMs from iKiP-DB instead. The data from iKiP-DB is reported as a .gmt file structured as follow: Kinase name, number of phosphosites as well as their uniprot ID, 13-aa fragments of the phosphosites with the middle AA being the phosphorylated AA.
In some cases, instead of an AA, the sequence had a _ denoted. I assumed this was a missing AA and wrote a program to fetch the missing bits of sequence from uniprot (fetch_from_uniprot.py). After running this program on every single sequence with missing AA’s I quickly realized that these missing AA’s were simply due to the phosphorylation site being near either terminal of the phosphopeptide, and thus this just led to slightly better understanding of the data.
To extract PWMs from the data, I wrote a program (gmt_to_pwm.py) that does counts of each AA for each position along the phosphosite and then calculates the rate of each AA in that position to create a PWM.

Then the goal was to apply these PWMs to the phosphosite data to try and predict what kinase was responsible for each phosphorylation event. First off, we decided to account for background signals in the data by creating PWMs based on the rate of each AA in each position relative to each phosphosite. I did this in the programs (background_pwm_counts.py and bg_pwm_rates.py), one notable decision was that I chose to calculate these separately based on the phosphorylated AA, as I expected there might be slight differences between them. 

Bam bam! Turns out that was (kinda) useless and I’ve been calculating my background wrong. 

Upon revision I just calculated the background signal by counting all amino acids +/- 6 from the phosphosite and taking the chance of any AA to appear in that position.

Prediction To predict i wrote a program (kinase_predition.py). It does a few things:
Takes a input similar to the phosphosites in our DB, converts the index of the phosphorylated AA’s in our sequence data to the corresponding indexes for sequence data from UniProt. 
Calculates the cumulative (multiplicative) probability of a kinase having phosphorylated this phosphosite for each PWM in DB and spits out a sorted list of probabilities.

#As of 3/11 I still need to add the background signal

7/11 - 22
Obtained a new understanding of what I’m actually doing here, so I sort-of started over on writing the whole prediction program in a more “disciplined fashion”.

Some notes: I’m basing the predictions on the paper (https://pubmed.ncbi.nlm.nih.gov/35234914/)
For the likelihood calculations, I’m adding 7 aa’s on either side of the peptide to make sure that each aa IN the phosphosite gets predicted with the entire PWM. This should be fairly simple to remove if it turns out to be dumb.  

Ah shit did some stuff wrong again
I have to search the Kinase-Substrate DB by phosphosite rather than substrate accession ID.
