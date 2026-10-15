import sys
import re
import numpy as np
import pandas as pd
from itertools import combinations
import torch.nn as nn  
from scipy.optimize import least_squares
import torch
import matplotlib.pyplot as plt
from scipy.stats import *
from thermo_modules import *

### read data
version = sys.argv[1]

ct_conc = np.array([1.140715997,4.202731758,12.55323512,32.15738598,58.73244109,208.5761541])

atac_mean= "/scratch/bclab/yuning/cutTagModel/forTraining/mean_atac_foxa1_neg_dox_pileups_for.new_merged_idr_all6conc_foxa1_1e-5_peaks.txt"
seq_file = "/scratch/bclab/yuning/cutTagModel/forTraining/onehot_noNnoChrM_new_merged_idr_all6conc_foxa1_1e-5_peaks.txt"

cutTag_file = "/scratch/bclab/yuning/cutTagModel/forTraining/cutTag_r1r2_6conc_pileups_for_new_merged_idr_all6conc_foxa1_1e-5_peaks.txt"

atac_mean_df = pd.read_csv(atac_mean, sep="\t", index_col=0, na_values="NaN",header=None)
#seq_df = pd.read_csv(seq_file, sep=",", index_col=0, na_values="NaN",header=None)
cutTag_df = pd.read_csv(cutTag_file, sep="\t", index_col=0, na_values="NaN",header=None)

#normalization_vector =  [0.2048815,0.4008087,1.301086,1.568565,2.495166,1,0.2040202,0.1903435,1.060806,1.594896,2.60861,0.71288]
normalization_vector =  [0.4117973,0.5718221,1.243795,1.46688,2.220708,1,0.4017243,0.4043921,1.05734,1.464847,2.316925,0.7754705]

normalization_vector = np.array(normalization_vector)



valid_coor = list()
seq_forNP = list()

dup_conc = np.concatenate((ct_conc, ct_conc),axis=None)

file1 = open(seq_file, 'r')
for line in file1:
    splited_line = line.strip().split(",")
    i = splited_line[0]
    onehot_seq_here = [int(i) for i in splited_line[1:]]
    curr_ct = cutTag_df.loc[i].to_numpy()/normalization_vector
    rep_cor = pearsonr(curr_ct[0:6], curr_ct[6:12])
    rep_cor2 = spearmanr(curr_ct[0:6], curr_ct[6:12])
    highest_dox_count = curr_ct[5] + curr_ct[11]
    conc_cor = spearmanr(curr_ct, dup_conc)
    if rep_cor.statistic > 0.9 and np.count_nonzero(curr_ct)>9 and conc_cor.statistic > 0.9: ## v9
#    if rep_cor.statistic > 0.9 and np.count_nonzero(curr_ct)>8 and conc_cor.statistic > 0.5: ## v8
    #if rep_cor.statistic > 0.95 and conc_cor.statistic > 0.8: ## v0
    #if rep_cor.statistic > 0.9 and np.count_nonzero(curr_ct)>8 and rep_cor2.statistic>0.9: ## v1
    #if rep_cor.statistic > 0.95 and np.count_nonzero(curr_ct)>8 and rep_cor2.statistic>0.8: ## v2
    #if rep_cor.statistic > 0.85 and np.count_nonzero(curr_ct)>8 and highest_dox_count>100: ## v3
    #if rep_cor.statistic > 0.85 and np.count_nonzero(curr_ct)>8 and highest_dox_count>100 and rep_cor2.statistic>0.85: ## v4
    #if rep_cor.statistic > 0.9 and np.count_nonzero(curr_ct)>8 and highest_dox_count>100 and rep_cor2.statistic>0.8: ## v5
        valid_coor.append(i)
    #print(",".join(line.strip().split(",")[1:]))
 
# Closing files
file1.close()



file2 = open("/scratch/bclab/yuning/cutTagModel/forTraining/new_training_list_"+version+".txt", "x") 

for L in valid_coor:
    file2.write(L+"\n") 

file2.close() 

print("end reading file\n")
