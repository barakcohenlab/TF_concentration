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

#ct_conc = np.array([0.090980598,0.209652969,0.527752875,2.192295078,4.242407153,17.58626505])
ct_conc = np.array([0.096707194,0.212048449,0.533705059,2.215395846,4.290815763,17.77722604])

#dup_conc = np.concatenate((ct_conc, ct_conc),axis=None)

atac_mean= "/scratch/bclab/yuning/cutTagModel/hnf4a_forTraining/mean_atac_hnf4a_neg_dox_pileups_for.new_merged_idr_all6conc_hnf4a_1e-5_peaks.txt"
seq_file = "/scratch/bclab/yuning/cutTagModel/hnf4a_forTraining/onehot_noNnoChrM_new_merged_idr_all6conc_hnf4a_1e-5_peaks.txt"

cutTag_file = "/scratch/bclab/yuning/cutTagModel/hnf4a_forTraining/cutTag_r1r2_6conc_pileups_for_new_merged_idr_all6conc_hnf4a_1e-5_peaks.txt"

atac_mean_df = pd.read_csv(atac_mean, sep="\t", index_col=0, na_values="NaN",header=None)
#seq_df = pd.read_csv(seq_file, sep=",", index_col=0, na_values="NaN",header=None)
cutTag_df = pd.read_csv(cutTag_file, sep="\t", index_col=0, na_values="NaN",header=None)

normalization_vector =  [0.2054157,
                         0.5706981,
                         0.4493818,
                         0.6656557,
                         0.9329646,
                         1,
                         0.2424389,
                         0.5406769,
                         0.7836418,
                         0.2722282,
                         0.7903418,
                         0.8849345]
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
    if rep_cor.statistic > 0.9 and np.count_nonzero(curr_ct)>9 and conc_cor.statistic > 0.8: ## v10
#if rep_cor.statistic > 0.9 and np.count_nonzero(curr_ct)>9 and conc_cor.statistic > 0.9: ## v9
#    if rep_cor.statistic > 0.9 and np.count_nonzero(curr_ct)>8 and conc_cor.statistic > 0.5: ## v8
    #if rep_cor.statistic > 0.95 and conc_cor.statistic > 0.8: ## v0
    #if rep_cor.statistic > 0.9 and np.count_nonzero(curr_ct)>8 and rep_cor2.statistic>0.8: ## v1
    #if rep_cor.statistic > 0.95 and np.count_nonzero(curr_ct)>8 and rep_cor2.statistic>0.8: ## v2
    #if rep_cor.statistic > 0.85 and np.count_nonzero(curr_ct)>8 : ## v3
    #if rep_cor.statistic > 0.85 and np.count_nonzero(curr_ct)>8 and highest_dox_count>100 and rep_cor2.statistic>0.85: ## v4
    #if rep_cor.statistic > 0.9 and np.count_nonzero(curr_ct)>8 and highest_dox_count>100 and rep_cor2.statistic>0.8: ## v5
        valid_coor.append(i)
    #print(",".join(line.strip().split(",")[1:]))
 
# Closing files
file1.close()



file2 = open("/scratch/bclab/yuning/cutTagModel/hnf4a_forTraining/new_training_list_"+version+".txt", "x") 

for L in valid_coor:
    file2.write(L+"\n") 

file2.close() 

print("end reading file\n")
