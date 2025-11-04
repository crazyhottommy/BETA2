# Author: Cliff Meyer, Dana-Farber Cancer Institute, 2008
# Clustering of position weight matrices based on 
# "A Novel Bayesian DNA Motif Comparison Method for Clustering and Retrieval" 

# Habib et al, 2008
# PLOS Computational Biology, Volume 4, Issue 2

import math
import numpy
import copy
import os,sys

SENSE = 0
ANTISENSE = 1
CUTOFF = {0.25:lambda x:0.1127*x+0.4485, 0.05:lambda x:0.2461*x+0.7060, \
          0.01:lambda x:0.4260*x+0.5751, 0.001:lambda x:0.6706*x+0.1165, \
          0.0001:lambda x:0.7301*x+0.2404}

def pwm2ic(pwm):
	"""Judge the information content of bps in a given pwm, return the number
	of ic bigger than 0.85[0.49,0.49,0.01,0.01]."""
	count = 0
	for pos in pwm:
		ic = 2
		for bp in pos:
			ic += bp*math.log(bp+0.0000001,2)
		if ic > 0.85:
			count += 1
	return count

def log_E_p_standard(n):
    """Log of the standard dirichlet prior."""
    E_p = n + 1
    E_p /= numpy.tile(E_p.sum(1), (4,1)).T
    return numpy.log(E_p)

def BLiC_score(M1, M2, cutoff=0.05):
    """
    if flag_r is ANTISENSE the reverse complement is the better match
    returns alignment score shift and orientation of M1 in the match
    matrix A is shifted before reversing
    """
    n1 = M1.shape[0]
    n2 = M2.shape[0]
 
    if n1 >= n2: # A the dame length with B, no changes
        A,B = M1,M2 # make sure A is longer, or the same length.
    else:
        A,B = M2,M1
        n1,n2 = n2,n1
 
    max_score, i_max= -999, -999
    flag_strand, flag_merge = False, False
    Brev = B[::-1,::-1]

    #align B to A, cut the edge of unaligned sequence.
    #i<0:      B:xxxxx
    #          A:    xxxxxxxx
    #i<=n1-n2: B:xxxxx
    #          A:xxxxxxxx
    #i>n1-n2:  B:    xxxxx
    #          A:xxxxxxxx 
    for i in range(1-n2, n1):
        if i<0:
            Bsub = B[-i:, :]
            Brev_sub = Brev[-i:, :]
            Asub = A[:n2+i, :]
           #ii = n1-i
        elif i <= n1-n2:
            Bsub = B
            Brev_sub = Brev
            Asub = A[i:i+n2, :]
           #ii =  n1
        elif n1-i < n2:
            Bsub = B[:n1-i, :]        #B:    xCATCGCxxx
            Brev_sub = Brev[:n1-i, :]
            Asub = A[i:, :]           #A: xxxxxxTCGC 
           #ii = n2+i
        
        score   = BLiC_score_aligned( Asub, Bsub )
        score_r = BLiC_score_aligned( Asub, Brev_sub )
        
        if score_r > score: 
            flag, score = ANTISENSE, score_r
        else:
            flag = SENSE
        
        if score > max_score:
            max_score = score
            flag_strand = flag
            i_max = i
            
    cutoff_len = max(pwm2ic(A), pwm2ic(B))
    if max_score >= CUTOFF[cutoff](cutoff_len):
        flag_merge = True
        
    return max_score, i_max, flag_strand, flag_merge

def BLiC_score_aligned(M1, M2):

    sum_M1_i = M1.sum(axis=1)
    sum_M1_i = 1.0*(sum_M1_i==0) + sum_M1_i
    A1 = M1.transpose()/sum_M1_i
    A1 = A1.transpose()

    sum_M2_i = M2[0:M1.shape[0],].sum(axis=1)
    sum_M2_i = 1.0*(sum_M2_i==0) + sum_M2_i
    A2 = M2[0:M1.shape[0],].transpose()/sum_M2_i
    A2 = A2.transpose()
 
    A12 = A1 + A2
 
    log_p1 = log_E_p_standard(A1)
    log_p2 = log_E_p_standard(A2)
    log_p12 = log_E_p_standard(A1*A2)
    log_pBG = log_E_p_standard(numpy.ones(A1.shape))
    
    s = 2 * (A12 * log_p12).sum() - (A1 * log_p1 + A2 * log_p2 + A12 * log_pBG).sum()
    return s

