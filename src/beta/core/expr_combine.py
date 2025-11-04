##########################
##        PartII        ##
##########################
from pkg_resources import resource_filename
from .corelib import *
import re,sys
from .up_down_distance import *
from .up_down_score import *

class expr_combine:
    
    # Differential expressed gene discorvered by Limma for Microarray data and Cuffdiff for RNAseq data ranked by q value;
    # Doing the Up and Down regulate test and then combine the candidate genes(both have a higer score and significant differential expressed q value);
    # Combine the two parts and get a Rank Product of each candidate gene.
    
    def __init__(self,options,expre_info):
        
        self.diff_expr = options.exprefile
        self.cutoff = options.cutoff
        self.gname2 = options.gname2
        self.diff_fdr = options.diff_fdr
        self.diff_amount = options.diff_amount
        self.peakfile = options.peakfile
        self.name = options.name
        self.method = options.method
        self.outdir = options.output
        self.expre_info = expre_info
        
        self.upgenes = []
        self.downgenes = []
        self.nochange = []
        self.selected = []
        self.upranks = {}
        self.downranks = {}

        if options.genome == 'hg18':
            self.genome = resource_filename('beta','references/hg18.refseq')
        elif options.genome == 'hg19':
            self.genome = resource_filename('beta','references/hg19.refseq')          
        elif options.genome == 'hg38':
            self.genome = resource_filename('beta','references/hg38.refseq')
        elif options.genome == 'mm9':
            self.genome = resource_filename('beta','references/mm9.refseq')           
        elif options.genome == 'mm10':
            self.genome = resource_filename('beta','references/mm10.refseq')
        else:
            self.genome = options.reference

    def gene_classify(self):

        uptotal = 0
        downtotal = 0
        refseq_column = int(self.expre_info['refseq'])
        fdr_column = int(self.expre_info['fdr'])
        t_column = int(self.expre_info['tscore'])
        
        with open(self.diff_expr, mode='r') as exprefile:

            for line in exprefile:
                if not line.strip() or line.startswith('AFFX') or line.startswith('#'): #skip 'Affy' markers and empty lines
                    continue
                line = line.strip()
                line = line.split('\t') #line = [number, ]
                refseq = line[refseq_column].strip('"_at"')#up[i][0].strip('"_at"')is in order to remove "" and _at in limma result {refgene:[t_scorerank]}
                m = re.search(r'\.\d', refseq)#some refseqID maybe like NM_001077.2, we need to remove the extra part
                if not m:
                    refseq = refseq
                else:
                    start = m.start()
                    refseq = refseq[:start]
                #print line[1]
                if line[t_column] != 'NA':
                    if float(line[t_column]) > 0:
                        uptotal += 1
                        self.upgenes.append([str(refseq),float(line[fdr_column])]) #put all the significantly upregulate genes into a list        
                    if float(line[t_column]) < 0:
                        downtotal += 1
                        self.downgenes.append([str(refseq),float(line[fdr_column])]) #put all the significantly downregulate genes into a list                       
                    if float(line[t_column]) == 0:
                        self.nochange.append([str(refseq),float(line[fdr_column])])
                else:
                    self.nochange.append([str(refseq),float(1)])

        total = [uptotal, downtotal]
        print(total)
        Info("Genes were seprated to two parts: up regulated and down regulated.")
        
        return total
        
                        
    def rank_genes(self):
        
        up = self.upgenes
        up.sort(key=lambda x: x[1])#upgene are which t score > 0, the smaller the qvaule, the more significant it is
        #[gene, qvalue]
        down = self.downgenes
        down.sort(key=lambda x: x[1])#downgenes are which t score < 0, the smaller qvaule, the more significant it is
    
        nochange = self.nochange
        nochange.sort(key=lambda x: x[1], reverse = True)
        
        upout = open('upgene.txt','w')
        downout = open('downgene.txt','w')

        
        if self.diff_amount <= 1: #required is a percentage
            c = int(len(up) * self.diff_amount)   #the number of up regulate genes required
            d = int(len(down) * self.diff_amount)  # the number of down regulate genes required
    
        if self.diff_amount > 1:
            c = d = int(self.diff_amount)
      
        newc = 0
        minqvalue = up[0][1]
        rank = 1
        same = 0
        for i in range(c):
            if float(up[i][1]) <= self.diff_fdr:
                qvalue = up[i][1]
                if qvalue == minqvalue:
                    rank = rank
                    same += 1
                else:
                    rank += same
                    same = 1
                minqvalue = qvalue
                self.upranks[up[i][0]] = [up[i][1], rank]
                upout.write(up[i][0] + '\t' + str(up[i][1]) + '\t' + str(rank) + '\n')
                newc += 1
            else:
                self.upranks[up[i][0]] = [up[i][1], 'NA']
                upout.write(up[i][0] + '\t' + str(up[i][1]) + '\t' + 'NA' + '\n')
        for i in range(c,len(up)):
            self.upranks[up[i][0]] = [up[i][1], 'NA']
            upout.write(up[i][0] + '\t' + str(up[i][1]) + '\t' + 'NA' + '\n')
        newd = 0
        minqvalue = down[0][1]
        rank = 1
        same = 0
        for j in range(d):
            if float(down[j][1]) <= self.diff_fdr:
                qvalue = down[j][1]
                if qvalue == minqvalue:
                    rank = rank
                    same += 1
                else:
                    rank += same
                    same = 1
                minqvalue = qvalue
                self.downranks[down[j][0]] = [down[j][1], rank]
                downout.write(down[j][0] + '\t' + str(down[j][1]) + '\t' + str(rank) + '\n')
                newd += 1
            else:
                self.downranks[down[j][0]] = [down[j][1], 'NA']
                downout.write(down[j][0] + '\t' + str(down[j][1]) + '\t' + 'NA' + '\n')               
 
        for j in range(d,len(down)):
            self.downranks[down[j][0]] = [down[j][1], 'NA']
            downout.write(down[j][0] + '\t' + str(down[j][1]) + '\t' + 'NA' + '\n')
        
        upout.close()
        downout.close()
        run_cmd("cut -f1 upgene.txt | head -" + str(newc) + " > up_list.txt")
        run_cmd("cut -f1 downgene.txt | head -" + str(newd) + " > down_list.txt")
        
        notdiffamount = int((newc + newd)/2)
        
        if len(nochange) != 0 and len(nochange) > notdiffamount:
            notdiff = open('notdiff.txt','w')
            for k in range(notdiffamount):
                notdiff.write(nochange[k][0] + '\n')
            notdiff.close()
        elif len(nochange) != 0 and len(nochange) < notdiffamount:
            notdiff = open('notdiff.txt','w')
            for k in range(len(nochange)):
                notdiff.write(nochange[k][0] + '\n')
            with open('upgene.txt') as upgenef:
                upgenelist = upgenef.readlines()
                for m in range((notdiffamount-len(nochange))//2):
                    g = upgenelist[-m-1].split('\t')
                    notdiff.write(g[0] + '\n')
            with open('downgene.txt') as downgenef:
                downgenelist = downgenef.readlines()
                for m in range((notdiffamount-len(nochange))//2):
                    g = downgenelist[-m-1].split('\t')
                    notdiff.write(g[0] + '\n')
            notdiff.close()
        else:
            run_cmd("cut -f1 upgene.txt | tail -" + str(newc//2) + " > notdiff1.txt")
            run_cmd("cut -f1 downgene.txt | tail -" + str(newd//2) + " > notdiff2.txt")
            run_cmd("cat notdiff1.txt notdiff2.txt > notdiff.txt")

        Info("Prepare file for the Up/Down Test")

    def ChGS(self):
        #use gene's regulate potential score to select a significant gene list to show the rank, from this , it will output 2 files, one is a picture in pdf file /
        #the other is a pvalue file store the pvalue information of the up and down regulate gene list with ks-test in a txt file
        
        command_line = "cp " + self.peakfile + " " + self.name + ".bed" # change the input peak file's name
        run_cmd(command_line)

        if self.method == "score":
            scorerun(self.name + '.txt', 'notdiff.txt', 'up_list.txt', 'down_list.txt', 'upregulate', 'downregulate', self.name, self.gname2)
            #scorefile, bglist, uplist, downlist, label1, label2, name, name2
        if self.method == "distance":
            distrun(self.peakfile, self.genome, 'notdiff.txt', 'up_list.txt', 'down_list.txt', 'upregulate', 'downregulate', self.name, self.gname2)
            #bedfile, reference, bglist, upgenes, downgenes, label1, label2, name, name2
            
        ## remove some redundent files
        run_cmd("rm up_list.txt")
        run_cmd("rm down_list.txt")
        run_cmd("rm upgene.txt")
        run_cmd("rm downgene.txt")
        run_cmd("rm notdiff*.txt")
        run_cmd("rm " + self.name + ".bed")
        run_cmd('mv %s %s'%(self.name+'_function_prediction.pdf', self.outdir))
        run_cmd('mv %s %s'%(self.name+'_function_prediction.R', self.outdir))
        Info("Finished, Find the result in %s_function_prediction.pdf"%self.name)
        
    def combine(self,bgtotal):
        # combine the peak rank file and expression rank file, then got a summary rank (peakrank * exprank)
        nontarget = []
        
        ChGS_Pvalues = open(self.name + "_pvalues.txt")
        for line in ChGS_Pvalues:
            line = line.strip()
            line = line.split('\t')
            if len(line) != 2:
                continue
            else:
                p = float(line[1])
                if p <= self.cutoff:
                    self.selected.append(line[0])
        # get the up or down regulate gene list's pavalue in a list, like ["upregualte"]
        if len(self.selected) == 0:
            Info("Both upregualte and downregualte gene list are not closer than the background, please cheak your data or looser the cutoff!")
            run_cmd('rm ' + self.name + "_pvalues.txt")
            run_cmd('rm ' + self.name + ".txt")
            sys.exit(0)
        else:
            upcounts = 0
            downcounts = 0
            p = 0
            for gene_type in self.selected:
                Info('Get the Rank Product of the %s genes'%gene_type)
                if self.gname2 == False:
                    geneID_col = 3
                else:
                    geneID_col = 6
                    
                bgenes = []#the up regulated genes whose binding rank and expression rank are not 'NA'
                egenes = []
                brank = {}
                erank = {}
                out = "%s_%s_target.txt"%(self.name,gene_type[1:-1]) # remove the ("")
                prank = {}#peak rank
                pinf = open(self.name + '.txt') #peak rank's file result from qianzi's method
                outf = open(out,'w')
                chrom = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13',
                         'chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
                
                for line in pinf:
                    if not line.strip():
                        continue
                    line = line.strip()
                    line = line.split('\t')
                    
                    if line[0] not in chrom: #skip the random chromosome and chrM
                        continue
                    else:
                        if line[7] != 'NA':
                            try:
                                prank[line[geneID_col]].append([line[0], line[1], line[2],line[3], line[4], line[5], line[6], int(line[7])])
                                #{geneID: chr, tss, tts, refseq, score, strand, symbol, rank}
                            except:
                                prank[line[geneID_col]] = [[line[0], line[1], line[2],line[3],line[4], line[5], line[6], int(line[7])]]
                        else:
                            try:
                                prank[line[geneID_col]].append([line[0], line[1], line[2], line[3],line[4], line[5], line[6], line[7]])
                                #{qeneID: chr, tss, tts, refseq, score, strand, symbol, rank}
                            except:
                                prank[line[geneID_col]] = [[line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7]]]
                
                for gene in list(prank.keys()):

                    info = prank[gene]
                    # Sort with NA values at the end (use float('inf') for NA)
                    info.sort(key=lambda x: float('inf') if x[-1] == 'NA' else x[-1])
                    prank[gene] = info[0]
                
                if gene_type == '"upregulate"':
                    keys = list(self.upranks.keys())
                    
                if gene_type == '"downregulate"':
                    keys = list(self.downranks.keys())
                    
                pkeys = list(prank.keys())
                bindinggenes = []
                
                for gene in keys:
                    if gene in pkeys:
                        if prank[gene][7] != 'NA':
                            bindinggenes.append([gene, prank[gene][3],prank[gene][6], prank[gene][4]])
                            #[generefseq,symbol,score] the first gene can be refseqID or genesymbol relied on the gname2 you set
                        else:
                            continue
                        
                bindinggenes.sort(key = lambda x:x[-1], reverse=True)# rank from high to low based on the regulatory potential score.
    
                bindingrank = {}
                rank = 1
                maxscore = bindinggenes[0][-1]
                for g in bindinggenes:
                    ID = g[0]
                    refseq = g[1]
                    symbol = g[2]
                    score = g[-1]
                    if score == maxscore:
                        rank = rank
                    else:
                        rank = rank + 1
                    bindingrank[ID] = [refseq, symbol, rank]
                    maxscore = score
     
                bkeys = list(bindingrank.keys())
                for gene in keys:
                    if gene in pkeys:
                        if gene_type == '"upregulate"':

                            if self.upranks[gene][1] != 'NA' and prank[gene][7] != 'NA':
                                
                                bgenes.append([gene, bindingrank[gene][0], bindingrank[gene][1], bindingrank[gene][2]])
                                egenes.append([gene, bindingrank[gene][0], bindingrank[gene][1], self.upranks[gene][1]])
                                upcounts += 1
                            elif self.upranks[gene][1] == 'NA' and prank[gene][7] == 'NA':
                                infos = prank[gene]
                                nongene = [infos[0],gene,infos[1],infos[2],infos[3],infos[5],infos[6]]
                                nontarget.append(tuple(nongene))
                            else:
                                continue

                        if gene_type == '"downregulate"':
                            if self.downranks[gene][1] != 'NA' and prank[gene][7] != 'NA':
                                
                                bgenes.append([gene, bindingrank[gene][0], bindingrank[gene][1], bindingrank[gene][2]])
                                egenes.append([gene, bindingrank[gene][0], bindingrank[gene][1], self.downranks[gene][1]])
                                downcounts += 1
                            elif self.downranks[gene][1] == 'NA' and prank[gene][7] == 'NA':
                                infos = prank[gene]
                                nongene = [infos[0],gene, infos[1],infos[2],infos[3],infos[5],infos[6]]
                                nontarget.append(tuple(nongene))
                            else:
                                continue                            
                    else:
                        continue
                                                                 
                for gene in bgenes:
                    brank[gene[0]] = [gene[1], gene[2], gene[3]]
    
                for gene in egenes:
                    erank[gene[0]] = [gene[1], gene[2], gene[3]]
 
                
                genes = list(brank.keys())
                genenumber = bgtotal[p]
                
                data = []
                for gene in genes:
                    RP = (float(brank[gene][2])/float(genenumber)) * (float(erank[gene][2])/float(genenumber))
                    data.append([brank[gene][0], brank[gene][1], RP])
                data.sort(key=lambda x:x[-1])
                k = 1
                for n in data:
                    rank = k
                    outf.write(n[0] + '\t' +n[1] + '\t' + str(n[2]) + '\t' + str(k) + '\n')
                    k += 1
                p+=1
        outf.close()
        run_cmd('rm ' + self.name + "_pvalues.txt")
        run_cmd('rm ' + self.name + ".txt")

        counts = [upcounts,downcounts]
        return (self.selected, counts, prank, nontarget)
