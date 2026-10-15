import re
import numpy as np
import pandas as pd
from itertools import combinations
import torch.nn as nn  
from scipy.optimize import least_squares
import torch
import matplotlib.pyplot as plt
from scipy.stats import *


    
def rev_comp(seq):
    """Take the reverse compliment of a sequence

    Parameters
    ----------
    seq : str
        The original sequence.

    Returns
    -------
    new_seq : str
        The reverse compliment.
    """
    compliment = {"A": "T", "C": "G", "G": "C", "T": "A","a": "t", "c": "g", "g": "c", "t": "a"}
    new_seq = seq[::-1]
    new_seq = "".join([compliment[i] for i in new_seq])
    return new_seq


def rev_comp_onehot(seq_array):
    """Take the reverse compliment of a one-hot encoded seq

    Parameters
    ----------
    seq : array
        The original sequence in one hot code.

    Returns
    -------
    new_seq : still one hot encoded
        The reverse compliment.
    """
    return np.flip(seq_array)


def transfer_pwm(score_array,ka1,ka2):
    return score_array**(ka1)*np.exp(ka2)

def kasum(score_array):
    return score_array.sum()

def kaprod2(score_array):
    return ( (score_array.sum())**2-(score_array**2).sum() )/2

def kaprod3(score_array):
    return ( (score_array.sum())**3-(score_array**3).sum()-3*((score_array**2).sum()*score_array.sum()-(score_array**3).sum()) )/6

def kaprod4(score_array):
    comb4 = combinations(score_array, 4)
    kaprod4_ini=0
    for i in comb4:
        kaprod4_ini += np.prod(i)
    return kaprod4_ini


def myMaxpool(score_array,motif_length):
    max_pool = nn.MaxPool1d(kernel_size=motif_length, stride=motif_length, ceil_mode=True)
    maxpooled = max_pool(torch.unsqueeze(torch.tensor(score_array),0))
    return maxpooled.numpy()[0]


def myMaxpool2(score_array,motif_length):
    max_pool = nn.MaxPool1d(kernel_size=2*motif_length, stride=2*motif_length, ceil_mode=True)
    maxpooled = max_pool(torch.unsqueeze(torch.tensor(score_array),0))
    return maxpooled.numpy()[0]


def topNscore(score_array,N=15):
    sorted_index_array = np.argsort(score_array) 
    sorted_array = score_array[sorted_index_array] 
    return sorted_array[-N : ] 



def trimScore(score_array,cutoff):
    return score_array[score_array>cutoff]


def ka_from_pwm_perPos_onehot(seq, pwm, motif_len):
    seq = np.array(seq)
    pwm = np.log(pwm)
    ### PWM is nparray shape (4, motif_len), seq is nparray (4,seq_len)
    #print(motif_len)
    # Number of positions where the motif can be scored
    n_scores = int((seq.shape[-1] - motif_len)/4 + 1)
    fscores = np.ones(n_scores)
    # Reverse compliment scores
    rscores = fscores.copy()
    maxOf2OrienScores = fscores.copy()
    # Calculate occ for forward and reverse k-mer at every position
    for pos in range(n_scores):
        pos4 = 4*pos
        f_kmer = seq[pos4:pos4 + motif_len]
        r_kmer = np.flip(f_kmer)

        # Initialize energy score
        fscore = f_kmer * pwm
        rscore = r_kmer * pwm
        
        fscores[pos] = np.exp(fscore.sum())
        rscores[pos] = np.exp(rscore.sum())
        maxOf2OrienScores[pos] = max(fscores[pos],rscores[pos])
        #maxOf2OrienScores[pos] = fscores[pos]+rscores[pos]

    return maxOf2OrienScores
    #return fscores, rscores



def all3Ka_from_kernel_onehot(seq, pwm, motif_len, motif_length):
    normalized_ka = ka_from_pwm_perPos_onehot(seq,pwm,motif_len)
    #maxpooled = myMaxpool(normalized_ka,motif_length)
    maxpooled = topNscore(normalized_ka,10)
    kas = kasum(maxpooled)
    ka_prod2 = kaprod2(maxpooled)
    ka_prod3 = kaprod3(maxpooled)
    return [kas,ka_prod2,ka_prod3] 

def individual3Ka_from_kernel_onehot(seq, pwm, motif_len, motif_length):
    normalized_ka = ka_from_pwm_perPos_onehot(seq,pwm,motif_len)
    #maxpooled = myMaxpool(normalized_ka,motif_length)
    ind_3 = topNscore(normalized_ka,3)
    return ind_3





def all4Ka_from_kernel_onehot(seq, pwm, motif_len, motif_length):
    normalized_ka = ka_from_pwm_perPos_onehot(seq,pwm,motif_len)
    #maxpooled = myMaxpool(normalized_ka,motif_length)
    maxpooled = topNscore(normalized_ka,10)
    kas = kasum(maxpooled)
    ka_prod2 = kaprod2(maxpooled)
    ka_prod3 = kaprod3(maxpooled)
    ka_prod4 = kaprod4(maxpooled)
    return [kas,ka_prod2,ka_prod3,ka_prod4] 



def each_kmer_ka_from_kernel_onehot(seq, pwm, motif_len, motif_length):
    normalized_ka = ka_from_pwm_perPos_onehot(seq,pwm,motif_len)
    #maxpooled = myMaxpool(normalized_ka,motif_length)
    maxpooled = topNscore(normalized_ka,10)
    return maxpooled


def all3Ka_from_kernel_onehot_multiSeq(seqs, pwm, motif_len, motif_length):
    to_return = list()
    for seq in seqs:
        to_return.append(all3Ka_from_kernel_onehot(seq, pwm, motif_len, motif_length))
    return to_return


def preprocess_seq_atac_input(input_array):
    ### input is nparray([['seq', atac_reads]])
    to_return = list()
    for i in input_array:
        #print(np.append(flat_my_one_hot_encode(i[0]),i[1]))
        to_return.append(np.append(flat_my_one_hot_encode(i[0]),i[1]))
        #print(i)
        #print('cycle')
    return np.array(to_return)



def order1_equation(conc, ka_array,cs,ctNF , gamma ):
    kas = ka_array[0]
    p = ctNF * ( conc * kas + gamma * cs * conc * kas) / ( 1 + cs + conc * kas + gamma * cs * conc * kas) +0.1
    return p ### pseudo count



def order2_equation(conc, ka_array,cs,ctNF , gamma,c2f ):
    shortcut = conc * ka_array[0]
    shortcut2 = conc**2 * ka_array[1] * c2f
    shortcut3 = (1+gamma*cs)*shortcut + (1+cs*gamma**2)*shortcut2
    p = ctNF * ( shortcut3) / ( 1 + cs + shortcut3) +0.1
    return p ### pseudo count

def order3_equation(conc, ka_array,cs,ctNF , gamma,c2f ):
    shortcut = conc * ka_array[0]
    shortcut2 = conc**2 * ka_array[1] * c2f
    shortcut3 = (1+gamma*cs)*shortcut + (1+cs*gamma**2)*shortcut2
    shortcut4 = conc**3 * ka_array[2] * c2f**3 * (1+cs*gamma**3)
    p = ctNF * ( shortcut3+shortcut4) / ( 1 + cs + shortcut3+shortcut4) +0.1
    return p ### pseudo count

def order4_equation(conc, ka_array,cs,ctNF , gamma,c2f ):
    shortcut = conc * ka_array[0]
    shortcut2 = conc**2 * ka_array[1] * c2f
    shortcut3 = (1+gamma*cs)*shortcut + (1+cs*gamma**2)*shortcut2
    shortcut4 = conc**3 * ka_array[2] * c2f**3 * (1+cs*gamma**3)
    shortcut5 = conc**4 * ka_array[3] * c2f**6 * (1+cs*gamma**4)
    p = ctNF * ( shortcut3+shortcut4+shortcut5) / ( 1 + cs + shortcut3+shortcut4+shortcut5) +0.1
    return p ### pseudo count



def read_pwm_asflatNParray(filename):
    with open(filename) as fin:
        pwm = []
        for eachline in fin:
            str_list = eachline.split()
            pwm.append([float(ele) for ele in str_list])
    pwm = np.array(pwm)
    #pwm = pwm.astype(np.float)
    return pwm.flatten()
