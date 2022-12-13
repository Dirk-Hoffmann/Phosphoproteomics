Scripts for analysing phosphosite data using positional weight matrices (PWMs)

outdated // again - probably intuitive, but contains scripts that probably don't serve any function anymore

background.py - calculates background for use in kinase prediction

constructFastaDB.py - constructs a dictionary of ACC_ID:FASTA, to reduce overhead by circumventing posting requests to uniprot all the time

evaluate_predictions.py - tests predictions of pwms vs ground truth (eg. phosphosite plus kinase-substrate dataset)

kinase_prediction.py - tries to predict what kinase from a dict of kinase:pwm pairs is most likely to phosphorylate a certain phosphosite. 

pred_data_prep.py - turns output of kinase_prediction.py into pickles which makes evaluate_predictions.py run faster

pwm_cross_entropy.py - calculates cross entropy scores between two pwms for the same kinase

pwms_from_list.py - generates pwm from a list of phosphosites
