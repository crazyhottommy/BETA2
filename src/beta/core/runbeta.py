#!/usr/bin/env python

#########################
####     Su Wang      ###
####    2-27-2015     ###
####   version 1.0.1  ###
#########################

"""Script Description:

BETA-Binding and Expression Targets Analysis: Use the ChIP-x data and Microarray/RNAseq data to predict a factor's targets

Part1:
use binding data to calc the score of each gene to be regulated by factor.
1. For each refseq gene in genome, input a distance (for example 100kb), then I get the peak center within 100kb from gene TSS. 
2. filter the peaks by p-value < 1e-5 from MACS, and only get top 10,000 peaks if it's more than 10,000
3. Then calculate a sum of 'Score' for each gene use this formula: 
  Sg = lambda ldx: sum([math.exp(-0.5-4*t) for t in ldx])
4. output is in bed format. the 5th column is score.

input bed file should 'chr1' instead of 'chrI'

Part2:
use differentail expression data rank each genes by t-score from limma result
1. input expression file need to be the standard limma result format
2. Rank each genes by the t-score calculated by limma, the positive one represent up regualate,and negative one represent downregulate genes.
3. Mutiply the rank of gene score and the rank of gene expression t score, get the the rank product
4. Do the permutation to get the fdr to dicide the confidence of the target genes we found
5. use the regulate potential score to do up and down regulate gene set association analysis compared to the genes with no differential expression as a background

Part3:
Get the FDR of each target genes via permutation


This code is free software; you can redistribute it and/or modify it.

@version: $0.5$
@author: Su Wang
@contact: wangsu0623@gmail.com
"""
import sys
import os
import time
import re
from pkg_resources import resource_filename

from beta.core.pscore import PScore
from beta.core.expr_combine import expr_combine
from beta.core.permp import permutation
from beta.motif.motif_scan import Motif_Scan
from beta.utils import Info, print_elapsed_time
from beta.core.fileformat_check import checkfileformat
from beta.core.opt_validator import opt_validate_basic, opt_validate_super, opt_validate_noexpre

def basicrun(args):
    start = time.time()
    opts = opt_validate_basic(args)
    #os.chdir(opts.output)#change the directory to the output directory
    opts.genomesequence = ''
    check = checkfileformat(opts)
    peakf = check.bed_format()
    opts.peakfile = peakf
    #print peakf
    expre_info = check.check_expr()
        
    g = PScore(opts)
    g.readfile()
    (genepeaks, peaksoutofDistance) = g.ScoreCalc(opts.distance)
    g.Output2File(opts.name)
    e = expr_combine(opts,expre_info)
    bgtotal = e.gene_classify()
    e.rank_genes()
    e.ChGS()
    (selected,counts,geneInfo,nontarget) = e.combine(bgtotal)
    
    f = permutation(opts,selected,counts,geneInfo, bgtotal)
    #counts is thw intersection of binding and expression, total number will be used as a background
    f.fdr()
    f.prepare_ouput_basic()
    m = Motif_Scan(genepeaks,peaksoutofDistance,opts,selected,nontarget)
    m.get_gene_list()
    m.getpeaks()
    print_elapsed_time(start)
    
def plusrun(args):
    start = time.time()
    opts = opt_validate_super(args)
    #os.chdir(opts.output)#change the directory to the output directory
    check = checkfileformat(opts)
    peakf = check.bed_format()
    opts.peakfile = peakf
    #print peakf
    check.check_fasta_dna()
    expre_info = check.check_expr()

    g = PScore(opts)
    g.readfile()
    (genepeaks, peaksoutofDistance) = g.ScoreCalc(opts.distance)
    g.Output2File(opts.name)
    e = expr_combine(opts,expre_info)
    
    bgtotal = e.gene_classify()
    e.rank_genes()
    e.ChGS()
    (selected,counts,geneInfo,nontarget) = e.combine(bgtotal)
    
    f = permutation(opts,selected,counts,geneInfo, bgtotal)
    #counts is the intersection of binding and expression, total number will be used as a background
    f.fdr()
    f.prepare_ouput_plus()
    
    m = Motif_Scan(genepeaks,peaksoutofDistance,opts,selected,nontarget)
    m.get_gene_list()
    m.getpeaks()
    m.getsequence()
    m.run_mis()
    m.statistical_test()
    m.out2html()
    print_elapsed_time(start)

def minusrun(args):
    start = time.time()
    opts = opt_validate_noexpre(args)
    #os.chdir(opts.output)#change the directory to the output directory
    opts.genomesequence = ''
    opts.gname2 = False
    check = checkfileformat(opts)
    peakf = check.bed_format()
    opts.peakfile = peakf
    #print peakf
        
    g = PScore(opts)
    g.readfile()
    (genepeaks, peaksoutofDistance) = g.ScoreCalc(opts.distance)
    g.noexpreoutput(opts.name)
    g.noexprepeaks(opts.name, genepeaks, opts.distance)
    print_elapsed_time(start)
