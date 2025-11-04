import os
import sys
import math
import re
import subprocess
from subprocess import call
from . import motif_parser as MotifParser
from ..core.corelib import *
from pkg_resources import resource_filename

def getImgdir(outdir):
    Imgdir = os.path.join(outdir, 'motiflogos')
    if os.path.isdir(Imgdir):
        cmd = 'rm -rf %s'%Imgdir
        run_cmd(cmd)
        
    Imgdirzip = resource_filename('beta','templates/motiflogos.zip')
    cmd = 'unzip ' + Imgdirzip + ' -d ' + outdir
    run_cmd(cmd)
        
    return Imgdir

def DictToList(root):
    """extract each node information"""
    #print root
    result = []
    if "node" not in list(root.keys()):
        return []
    if not root['node']:
        return result
    else:
        result.append(root['node'])
        for each in root['children']:
            result.extend(DictToList(each))
        return result

def get_mis_oneside(motif, mistable, Imgdir, cutoff):
    if cutoff <= 1:
        pcut = cutoff
        ncut = 1000
    if cutoff > 1:
        pcut = 1
        ncut = cutoff
        
    motif_file = motif       
    inf = open(motif_file)

    count = 0
    outf = open(mistable,'w')
    outf.write('\t'.join(['id','synonym','dbd', 'zscore', 'species', 'pssm','logoImg','hits']) + '\n')
    for line in inf:
        if line.startswith('MotifID') or line.startswith(' '):
            continue
        else:
            line = line.strip()
            line = line.split('\t')
            #ID, species, symbol, DBD, PSSM, tscore, pvalue
            logoImg = os.path.join(Imgdir,line[0] + '.png')
            pvalue = line[-1]
            tscore = float(line[-2])
            if float(pvalue) <= pcut and count < ncut and tscore > 0: 
                data = '\t'.join([line[0],line[2],line[3],line[-1],line[1],line[4],logoImg,line[-2]]) + '\n'
                outf.write(data)
                count += 1
            else:
                continue
            
    outf.close()

def get_mis_twoside(motif, mistable1, mistable2, Imgdir, cutoff):
                
    if cutoff <= 1:
        pcut = cutoff
        ncut = 1000
    if cutoff > 1:
        pcut = 1
        ncut = cutoff
        
    outf1 = open(mistable1,'w')
    outf2 = open(mistable2,'w')
    outf1.write('\t'.join(['id','synonym','dbd', 'zscore', 'species', 'pssm','logoImg','hits']) + '\n')
    outf2.write('\t'.join(['id','synonym','dbd', 'zscore', 'species', 'pssm','logoImg','hits']) + '\n')
    count1 = 0
    count2 = 0

    motif_file = motif       
    inf = open(motif_file)
    for line in inf:
        if line.startswith('MotifID') or line.startswith(' '):
            continue
        else:
            line = line.strip()
            line = line.split('\t')
            logoImg = os.path.join(Imgdir,line[0] + '.png')
            pvalue = line[-1]
            tscore = float(line[-2])
            if float(pvalue) <= pcut and count1 < ncut and tscore > 0:
                data = '\t'.join([line[0],line[2],line[3],line[-1],line[1],line[4],logoImg,line[-2]]) + '\n'
                outf1.write(data)
                count1 += 1
            if float(pvalue) <= pcut and count2 < ncut and tscore < 0:
                data = '\t'.join([line[0],line[2],line[3],line[-1],line[1],line[4],logoImg,line[-2]]) + '\n'
                outf2.write(data)
                count2 += 1
            else:
                continue
            
    outf1.close()
    outf2.close()
    

def motif_info2_oneside(sqposeTable,Imgdir):
       
    p=MotifParser.MotifParser()
    p.ParserTable(sqposeTable)
    s2 = list(p.motifs.values())
    s2 = sorted(s2, key=lambda x:float(x['zscore'][0]))
    
    i = 0
    
    while i< len(s2):
        
        logo = [s2[i]['synonym'][0]]
        dbd = [s2[i]['dbd'][0]]
        for j in range(len(s2)-1,i,-1):
            if p._Similarity(s2[i]['id'][0],s2[j]['id'][0])[-1]:
                logo.append(s2[j]['synonym'][0])
                dbd.append(s2[j]['dbd'][0])
                id = s2[j]['id'][0]
                del p.motifs[id]
                del s2[j]
        s2[i]['synonym'] = logo
        s2[i]['dbd']=dbd
        i = i+1
    output = []
   
    for i in s2:
        logo = i['synonym']
        #cmd = 'cp %s %s'%(str(i['logoImg'][0]),outdir)
        logor = str(i['logoImg'][0])
        species = i['species'][0]
        dbd = i['dbd']
        tempt = ['; '.join(logo),'; '.join(dbd),species,str(i['zscore'][0]),str(i['hits'][0]),logor]
        output.append(tempt)
        
    return output

def motif_info2_twoside(seqposTable1,seqposTable2,Imgdir):
    
    p=MotifParser.MotifParser()
    p.ParserTable(seqposTable1)
    q=MotifParser.MotifParser()
    q.ParserTable(seqposTable2)
    s2 = list(p.motifs.values())
    s2 = sorted(s2, key=lambda x:float(x['zscore'][0]))

    m2 = list(q.motifs.values())
    m2 = sorted(m2, key=lambda x:float(x['zscore'][0]))
    
    i = 0
    
    while i< len(s2):
        
        logo = [s2[i]['synonym'][0]]
        dbd = [s2[i]['dbd'][0]]
        for j in range(len(s2)-1,i,-1):
            if p._Similarity(s2[i]['id'][0],s2[j]['id'][0])[-1]:
                logo.append(s2[j]['synonym'][0])
                dbd.append(s2[j]['dbd'][0])
                id = s2[j]['id'][0]
                del p.motifs[id]
                del s2[j]
        s2[i]['synonym'] = logo
        s2[i]['dbd']=dbd
        i = i+1

    k = 0
    while k< len(m2):
        
        logo = [m2[k]['synonym'][0]]
        dbd = [m2[k]['dbd'][0]]
        for j in range(len(m2)-1,k,-1):
            if q._Similarity(m2[k]['id'][0],m2[j]['id'][0])[-1]:
                logo.append(m2[j]['synonym'][0])
                dbd.append(m2[j]['dbd'][0])
                id = m2[j]['id'][0]
                del q.motifs[id]
                del m2[j]
        m2[k]['synonym'] = logo
        m2[k]['dbd']=dbd
        k = k+1
    
    output = []
   
    for i in s2:
        logo = i['synonym']
        #cmd = 'cp %s %s'%(str(i['logoImg'][0]),outdir)
        logor = str(i['logoImg'][0])
        species = i['species'][0]
        dbd = i['dbd']
        tempt = ['; '.join(logo),'; '.join(dbd),species,str(i['zscore'][0]),str(i['hits'][0]),logor]
        output.append(tempt)
    for j in m2:
        logo = j['synonym']
        #cmd = 'cp %s %s'%(str(i['logoImg'][0]),outdir)
        logor = str(j['logoImg'][0])
        species = j['species'][0]
        dbd = j['dbd']
        tempt = ['; '.join(logo),'; '.join(dbd),species,str(j['zscore'][0]),str(j['hits'][0]),logor]
        output.append(tempt)
    return output

def runmotifcluster_side1(outdir, motif, mistable, cutoff):
    Imgdir = getImgdir(outdir)
    get_mis_oneside(motif, mistable, Imgdir, cutoff)
    output = motif_info2_oneside(mistable, Imgdir)
    run_cmd('rm %s'%mistable)
    return output

def runmotifcluster_side2(outdir, motif, mistable1, mistable2, cutoff):
    Imgdir = getImgdir(outdir)
    get_mis_twoside(motif, mistable1, mistable2, Imgdir, cutoff)
    output = motif_info2_twoside(mistable1, mistable2, Imgdir)
    run_cmd('rm %s'%mistable1)
    run_cmd('rm %s'%mistable2)
    return output       
