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


version = sys.argv[1]
seq_file = "/scratch/bclab/yuning/sequencing/noRep_padded_pft_v1_name_bc.txt"
output_file = "/scratch/bclab/yuning/sequencing/predicted_kas/foxa1_kas_"+version+".txt"
#filename="/scratch/bclab/yuning/cutTagModel/new_outputs/pwm_v2v11e4linear.txt" # trained PWM
filename="/scratch/bclab/yuning/cutTagModel/0523_redo_foxa1/foxa1_pwm_"+version+".txt"


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

seq_df = pd.read_csv(seq_file, sep="\t", index_col=0, na_values="NaN",header=None)

with open(output_file, 'x') as the_file:
    #the_file.write('Hello\n')
    for index, row in seq_df.iterrows():
        #print( flat_my_one_hot_encode(row[1]))
        kas_here = flat_my_one_hot_and_score(row[1])
        the_file.write(index+'\t'+kas_here+'\n')
