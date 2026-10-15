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

### read data


seq_file = "/scratch/bclab/yuning/cutTagModel/forTraining/10p_tv8_onehot_noNnoChrM_new_merged_idr_all6conc_foxa1_1e-5_peaks.txt"
output_prefix = "short_o3_foxa1_thermo_o3_v8_1e3_v8v"

output_file = "/scratch/bclab/yuning/cutTagModel/0523_redo_foxa1/"+output_prefix+"_predicted.txt"
output_fig = "/scratch/bclab/yuning/cutTagModel/0523_redo_foxa1/"+output_prefix+"_allplots.png"

my_fscale=0.1
ct_conc = np.array([1.140715997,4.202731758,12.55323512,32.15738598,58.73244109,208.5761541])


atac_mean= "/scratch/bclab/yuning/cutTagModel/forTraining/mean_atac_foxa1_neg_dox_pileups_for.new_merged_idr_all6conc_foxa1_1e-5_peaks.txt"

cutTag_file = "/scratch/bclab/yuning/cutTagModel/forTraining/cutTag_r1r2_6conc_pileups_for_new_merged_idr_all6conc_foxa1_1e-5_peaks.txt"

atac_mean_df = pd.read_csv(atac_mean, sep="\t", index_col=0, na_values="NaN",header=None)
#atac_max_df = pd.read_csv(atac_max, sep="\t", index_col=0, na_values="NaN",header=None)

#seq_df = pd.read_csv(seq_file, sep=",", index_col=0, na_values="NaN",header=None)
cutTag_df = pd.read_csv(cutTag_file, sep="\t", index_col=0, na_values="NaN",header=None)

valid_coor = list()
seq_forNP = list()

dup_conc = np.concatenate((ct_conc, ct_conc),axis=None)

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

print("end reading file\n")
        
#valid_coor = all_coor


atac_forNP =  atac_mean_df[1][valid_coor].to_numpy() + .1 ### pseudo count
#seq_forNP = seq_df.loc[valid_coor].to_numpy()
#print(atac_forNP)
tmp_atacNF = np.percentile(atac_forNP,90)
print(tmp_atacNF)

#normalization_vector =  [0.2048815,0.4008087,1.301086,1.568565,2.495166,1,0.2040202,0.1903435,1.060806,1.594896,2.60861,0.71288]
normalization_vector =  [0.4117973,0.5718221,1.243795,1.46688,2.220708,1,0.4017243,0.4043921,1.05734,1.464847,2.316925,0.7754705]

normalization_vector = np.array(normalization_vector)

test_output = cutTag_df.loc[valid_coor].to_numpy()/normalization_vector + .1  ### pseudo count
ctNF_r0 = np.max(test_output, axis=1) 
ctNF_r0 = np.full((test_output.shape[0], ), np.percentile(test_output,90)) ### alt
print(np.percentile(test_output,90))

filename="/scratch/bclab/yuning/sequencing/foxa1_pwm_MA0148d5.txt"
pwm=read_pwm_asflatNParray(filename)
pwm
motif_len = pwm.shape[-1]   #int(pwm.shape[-1]/4)
motif_length = int(pwm.shape[-1]/4)

#transfered_pwm = transfer_pwm(pwm/0.25,0.203344,-6.37013/8)
weight_vec = [0.5,0.5,0.5, 1.,1.,1., 0.5,0.5,0.5, 1.,1.,1.]
print(weight_vec)

def loss_transform(x):
    return x * weight_vec

def cutTag_fromPWM_o2_round0(parS, training_seq, atac, ctNF):
    pwm_here = np.array(transfer_pwm(pwm,parS[0], parS[1]))
#   seq_len = len(training_seq)/4
    #print(seq_len)
    ka_all = all3Ka_from_kernel_onehot(training_seq,pwm_here,motif_len,motif_length)
    cs_here = parS[2]* tmp_atacNF / atac - 1
    cs_here = max(cs_here, 0)
    #print(cs_here)
    ### ka1 = parS[0],ka2 = parS[1],atacNF = parS[2],   gamma = parS[3], c2f = parS[4]
    return order3_equation(ct_conc, ka_all ,cs_here ,ctNF*parS[5] , parS[3] , parS[4])
    #return order3_equation(ct_conc, ka_all ,cs_here ,ctNF ,0.391 , 1.)

#np.percentile(atac_forNP,.5)/np.percentile(ctNF_r0,.5)
parStart = [.3, -0.1,.94 , 0.2, .85,1.]
parLower = [ 0, float("-inf"),0, 0, 0, 0.]
parUpper = [float("inf"),float("inf"), float("inf"), 1, float("inf"), float("inf")]

print("r0 starting params:")
print(parStart)
print("==================")

def residFun_r0(parS,seqs,test_output):
    all_resid = np.zeros(test_output.shape)
    for i in range(0,test_output.shape[0]):
        pred = cutTag_fromPWM_o2_round0(parS, seqs[i],atac_forNP[i], ctNF_r0[i])
        #all_resid[i] = test_output[i] - np.concatenate((pred, pred),axis=None)
        all_resid[i] = loss_transform(test_output[i]) - loss_transform(np.concatenate((pred, pred),axis=None) )
    return all_resid.flatten()

#residFun(parStart, seq_forNP, test_output).sum()

res_lsq_r0 = least_squares(residFun_r0, parStart, args=(seq_forNP, test_output),method='trf',bounds=(parLower, parUpper),verbose=2,ftol=1e-3, max_nfev=50, loss='soft_l1', f_scale=my_fscale)
#res_lsq_lm = least_squares(residFun, parStart, args=(seq_forNP, test_output),method='lm')

print("initial round params: ka1, ka2, atacNF, gamma, c2f, ctNF_factor")
print(res_lsq_r0.x)

transfered_pwm_r0 = transfer_pwm(pwm,res_lsq_r0.x[0], res_lsq_r0.x[1])

def cutTag_fromPWM_o2_one(parS,constantpars, training_seq, atac):
    pwm_here = constantpars[0:32]
    ka_all = all3Ka_from_kernel_onehot(training_seq,pwm_here,motif_len,motif_length)
    cs_here = (   constantpars[32]*tmp_atacNF ) / atac - 1
    cs_here = max(cs_here, 0)
    ### loc_factor = parS[0], bg = parS[1]
    ###  atac_NF = constantpars[32], gamma = constantpars[33], c2f = constantpars[34]
    return order3_equation(ct_conc, ka_all ,cs_here , parS[0] , constantpars[33], constantpars[34] )

 
#parStart = [1e6]
#parLower = [0]
#parUpper = [float("inf")]

#parLower <- c(rep(-Inf,32),max(atac_mean)/2,1,0)
#parUpper <- c(rep(Inf,32),Inf,Inf,1)


def residFun_one(parS,constantpars,seq,atac, test_output):
    pred = cutTag_fromPWM_o2_one(parS,constantpars, seq, atac)
    #all_resid = np.log(test_output) - np.log(np.concatenate((pred, pred),axis=None) )
    #all_resid = test_output**0.5 - np.concatenate((pred, pred),axis=None)**0.5
    all_resid = loss_transform(test_output) - loss_transform(np.concatenate((pred, pred),axis=None) )
    return all_resid.flatten()

def fit_and_plot_one(constantpars,i):
    seq = seq_forNP[i]
    atac = atac_forNP[i]
    test_output_local = test_output [i]
    res_lsq_tmp = least_squares(residFun_one, [ctNF_r0[i]*2.], args=(constantpars,seq,atac, test_output_local),method='trf',bounds=([0], [float("inf")]), loss='linear')
    pred = cutTag_fromPWM_o2_one(res_lsq_tmp.x,constantpars, seq, atac)
    return res_lsq_tmp.x
    #return [test_output_local, np.concatenate((pred, pred),axis=None) ]


round0_par = np.append(transfered_pwm_r0,[res_lsq_r0.x[2],res_lsq_r0.x[3] , res_lsq_r0.x[4]])

trained_ctNF = list()

print("start individual round 1")
for index in valid_coor:
    #print(index+",")
    i = valid_coor.index(index)
    tmp_par = fit_and_plot_one(round0_par,i)
    trained_ctNF.append(tmp_par[0])
    #fit_and_plot_one(res_lsq.x,i)

print("end round 1\n=====================================")
    
def cutTag_fromPWM_o2(parS, training_seq, atac, ctNF):
    pwm_here = np.array(parS[0:32])
    ka_all = all3Ka_from_kernel_onehot(training_seq,pwm_here,motif_len,motif_length)
    cs_here = parS[32]*tmp_atacNF  / atac - 1
    #print(cs_here)
    cs_here = max(cs_here, 0)
    ### atac_NF = parS[32],  gamma = parS[33], c2f = parS[34]
    return order3_equation(ct_conc, ka_all ,cs_here ,ctNF , parS[33] , parS[34])
    #return order3_equation(ct_conc, ka_all ,cs_here ,ctNF , parS[33] , round0_par[34]) ### fix


#parStart = np.append(transfered_pwm_r0,[1,1,1])
parStart = round0_par
parLower = np.append(transfered_pwm_r0*0.001,[ 0, 0, 0])
parUpper = np.append(transfered_pwm_r0*1000.,[float("inf"), 1, float("inf")])

#parLower <- c(rep(-Inf,32),max(atac_mean)/2,1,0)
#parUpper <- c(rep(Inf,32),Inf,Inf,1)


def residFun(parS,seqs,test_output):
    all_resid = np.zeros(test_output.shape)
    for i in range(0,test_output.shape[0]):
        pred = cutTag_fromPWM_o2(parS, seqs[i],atac_forNP[i], trained_ctNF[i])
        all_resid[i] = loss_transform(test_output[i]) - loss_transform(np.concatenate((pred, pred),axis=None))
        #all_resid[i] = np.log(test_output[i]) - np.log(np.concatenate((pred, pred),axis=None) )
    return all_resid.flatten()


res_lsq = least_squares(residFun, parStart, args=(seq_forNP, test_output),method='trf',bounds=(parLower, parUpper),verbose=2,ftol=1e-4, max_nfev=50, loss='soft_l1', f_scale=my_fscale)
#res_lsq_lm = least_squares(residFun, parStart, args=(seq_forNP, test_output),method='lm')

print("pwm matrix")
np.savetxt(sys.stdout, res_lsq.x[0:32].reshape(8,4))
print("=======================\nothers: atacNF, gamma, c2f")
print(res_lsq.x[32:])




file2 = open(output_file, "w") 


def predict_after(parS , seqs, write_to):
    all_resid = np.zeros(test_output.shape)
    for i in range(0,test_output.shape[0]):
        pred = cutTag_fromPWM_o2(parS, seqs[i],atac_forNP[i], trained_ctNF[i])
        all_resid[i] = np.concatenate((pred, pred),axis=None) 
        if write_to:
            file2.write(valid_coor[i]+"\t"+"\t".join(map(str,pred))+"\n") 
    return all_resid.flatten()


t_test = ct_conc
y_true = test_output.flatten()
y_lsq = predict_after(res_lsq.x, seq_forNP, True)
y_start = predict_after(round0_par, seq_forNP, False)

file2.close() 



plt.plot(y_true, y_lsq,  'o',label="revised PWM")
#plt.plot(y_true, y_start,'o', label='starting')
#plt.plot(t_test, y_log, label='cauchy loss')
plt.yscale("log")
plt.xscale("log")
plt.xlabel("t")
plt.ylabel("y")
plt.legend()
plt.savefig(output_fig) 

