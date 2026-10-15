import re
import sys
import numpy as np
import pandas as pd
from itertools import combinations
import torch.nn as nn  
from scipy.optimize import least_squares
import torch
import matplotlib.pyplot as plt
from scipy.stats import *
from thermo_modules import *
import random


#version = sys.argv[1]
#re2g_ver = sys.argv[2]
#tf = sys.argv[3]


seq_file = "/scratch/bclab/yuning/cutTagModel/forTraining/tv8_onehot_noNnoChrM_new_merged_idr_all6conc_foxa1_1e-5_peaks.txt"
output_file = "/scratch/bclab/yuning/cutTagModel/bs_count/tv8_foxa1.txt"
#filename="/scratch/bclab/yuning/cutTagModel/new_outputs/pwm_v2v11e4linear.txt" # trained PWM
filename="/scratch/bclab/yuning/cutTagModel/0523_redo_foxa1/foxa1_pwm_nt9_v8v_v8_1e3_0611.txt"


ct_conc = np.array([1.140715997,4.202731758,12.55323512,32.15738598,58.73244109,208.5761541])

pwm=read_pwm_asflatNParray(filename)
motif_len = pwm.shape[-1]   #int(pwm.shape[-1]/4)
motif_length = int(pwm.shape[-1]/4)

def flat_my_one_hot_and_score(seq):
    mapping = dict(zip("ACGTacgt", [0,1,2,3]*2))    
    seq2 = [mapping[i] for i in seq]
    to_return = np.eye(4)[seq2].flatten()
    ka_all = all3Ka_from_kernel_onehot(to_return,pwm,motif_len,motif_length)
    ka_all = np.array(ka_all)
    return "\t".join(ka_all.astype(str))

#seq_df = pd.read_csv(seq_file, sep="\t", index_col=0, na_values="NaN",header=None)

valid_coor = list()
seq_forNP = list()

file1 = open(seq_file, 'r')
for line in file1:
    splited_line = line.strip().split(",")
    i = splited_line[0]
    onehot_seq_here = [int(i) for i in splited_line[1:]]
    valid_coor.append(i)
    seq_forNP.append(onehot_seq_here)
    #print(",".join(line.strip().split(",")[1:]))
 
# Closing files
file1.close()

def simulate_bind(ka_all):
    count = 0
    conc = 208.5761541
    #conc = 58.73244109
    #conc = 58.73244109
    for i in ka_all:
        b_prob = i*conc/(1+i*conc)
        prob = random.random()
        if prob < b_prob:
            count = count + 1
    return count
        

for seq in seq_forNP:
    ka_all = each_kmer_ka_from_kernel_onehot(seq,pwm,motif_len,motif_length)
    count_here = simulate_bind(ka_all)
    #ka_all = np.array(ka_all)
    print(count_here)


#with open(output_file, 'x') as the_file:
    #the_file.write('Hello\n')
#    for index, row in seq_df.iterrows():
        #print( flat_my_one_hot_encode(row[1]))
#        kas_here = flat_my_one_hot_and_score(row[1])
#        the_file.write(index+'\t'+kas_here+'\n')
