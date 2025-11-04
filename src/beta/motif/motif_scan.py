from ..core.corelib import *
from . import motif_parser as MP
from pkg_resources import resource_filename
from .fastafrombed import *
from .motif_clustering import *
import json #make the string to list
import os
import math

#Update: make the motif detecting part availiable in other genome assembily (not limited to hg19 and mm9).

class Motif_Scan:
    '''this module has 2 parts, one is sequence preparation, the other is the motif scan by mis'''
    def __init__(self, genepeaks, peaksoutofDistance, options, selected, nontarget):
        self.genepeaks = genepeaks
        self.distance = options.distance
        self.curdir = os.getcwd()
        self.outdir = options.output
        self.name = options.name
        self.motifnumber = options.motifnumber
        self.nontarget = nontarget
        self.peaksoutofDistance = peaksoutofDistance
        if options.genome and options.genome.startswith('hg'):
            self.Species = 'Homo sapiens'
        elif options.genome and options.genome.startswith('mm'):
            self.Species = 'Mus musculus'
        else:
            self.Species = "NA"
        self.genefiles = []
        self.selected = selected
        self.associate_peaks = []
        for i in self.selected:
            if i == '"upregulate"':
                self.genefiles.append(os.path.join(self.outdir,self.name + '_uptarget.txt'))
                self.associate_peaks.append(self.name + '_uptarget_associate_peaks.bed')
            if i == '"downregulate"':
                self.genefiles.append(os.path.join(self.outdir,self.name + '_downtarget.txt'))
                self.associate_peaks.append(self.name + '_downtarget_associate_peaks.bed')
        
        self.genomesequence = options.genomesequence
                          
    ############################
    #part1 Sequence preparation
    ############################
                
    def get_gene_list(self):
        files = self.genefiles
        self.genelists = []#uptarget and downtarget gene lists
        self.motifgenelists = []#Only the significant target genes associated peaks join to do the motif analysis.
        for f in files:
            genelist = [] #selected target genes
            with open (f) as inf:
                infile = inf.readlines()
                totalnumber = len(infile)
                if totalnumber <= 500:
                    top = totalnumber
                elif totalnumber > 500 and totalnumber <= 1000:
                    top = 250 + 0.5*totalnumber
                else:
                    top = 0.5 * totalnumber
            
                for line in infile:
                    if line.startswith('Chroms'):
                        pass
                    else:
                        line = line.strip()
                        line = line.split('\t')
                        gene = [line[0],line[3],line[1],line[2],line[5],line[6]]#chromrefseqID, start, end, strands, symbol
                        genelist.append(tuple(gene))
                motifgenelist = genelist[:int(top)]
            self.motifgenelists.append(motifgenelist)
            self.genelists.append(genelist)
        
    def getpeaks(self):
        '''get all peaks within the certain distance around this gene'''
        Info("pick out the peaks %sbp around the selected genes"%str(self.distance))
        score = lambda t: math.exp(-0.5-4*t)
        peakdist = []
        for i in range(len(self.genelists)):
            dist = {}
            associatepeakf = open(self.associate_peaks[i],'w')
            associatepeakf.write('chrom\tpStart\tpEnd\tRefseq\tSymbol\tDistance\tScore\n')
            for gene in self.genelists[i]:
                if gene[4] == '+':
                    TSS = gene[2]
                else:
                    TSS = gene[3]
                try:
                    peaks = self.genepeaks[gene]
                    for peak in peaks:
                        d = (int(peak[1]) + int(peak[2]))/2 - int(TSS)
                        if d <= self.distance:
                            regscore = score(abs(float(d))/float(self.distance))
                            associatepeakf.write(str(peak[0]) + '\t' + str(peak[1]) + '\t' + str(peak[2]) + '\t' + gene[1] + '\t' + gene[5] + '\t' + str(d) + '\t' + str(regscore) + '\n')
                            if gene in self.motifgenelists[i]:
                                key = tuple(peak)
                                dist[key] = [gene,d]
                        else:
                            continue
                except KeyError:
                    pass
            associatepeakf.close()
            run_cmd('mv %s %s'%(self.associate_peaks[i], self.outdir))
            
            peakdist.append(dist)
    
        Info("Finished: Find target gene associated peaks in %s"%self.outdir)
        
        self.peaklists = []
        
        if len(peakdist) == 2:
            up = peakdist[0]
            down = peakdist[1]
            uppeaks = list(up.keys())
            downpeaks = list(down.keys())
            for peak in uppeaks:
                if peak in downpeaks:
                    updist = up[peak][-1]
                    downdist = down[peak][-1]
                    if updist > downdist:
                        del up[peak]
                    elif downdist > updist:
                        del down[peak]
                    else:
                        del up[peak]
                        del down[peak]
                else:
                    continue
                
            self.peaklists.append(list(up.keys()))
            self.peaklists.append(list(down.keys()))
        elif len(peakdist) == 1:
            self.peaklists.append(list(peakdist[0].keys()))
        self.peaklists.append(self.peaksoutofDistance)
    def getsequence(self):
        '''extend the peak summit to 3 regions, left(-300,-101), middle(-100,99), right(100,299), get
           the fasta sequece by fastaFromBED--bedtools need to be installed'''
        if self.Species == 'NA':
            Info("Motif database is not provided for your species, Skip the motif part ;)")
            Info("For motif analysis, set the -g with hg19, mm9, hg or mm an rerun your job ;)")
            Info("BETA Finished, Find your result in %s"%self.outdir)
            sys.exit(1)
        flags = []
        for flag in self.selected:
            #'upregulate' or 'downregulate'
            if flag == '"upregulate"':
                flag = 'up'
                flags.append(flag)
            if flag == '"downregulate"':
                flag = 'down'
                flags.append(flag)
        flags.append('non')
        self.fastas = []
        i = 0
        for peaks in self.peaklists:#different gene types(up and down)
            flag = flags[i]
            sequences = []
            temp = open('select_peak.bed','w')
            for peak in peaks:#each gene has a lot of peaks within 100k
                peak_center = (int(peak[1]) + int(peak[2]))/2 #peak center equals to peak summit if the input peak file is the peak summit file
                temp.write('%s\t%s\t%s\n'%(peak[0],str(peak_center),str(peak_center + 1)))
            temp.close()
            location = 'middle'
            commandline = "awk -v OFS='\t' '{print $1, $2-100, $3+99}' %s > %s"%('select_peak.bed', flag + '_' + location + '.bed')
            run_cmd(commandline)
            sequences.append(flag + '_' + location + '.bed')
            
            location = 'left'
            commandline = "awk -v OFS='\t' '{print $1, $2-300, $2-101}' %s > %s"%('select_peak.bed', flag + '_' + location + '.bed')
            run_cmd(commandline)
            sequences.append(flag + '_' +  location + '.bed')
            
            location = 'right'
            commandline = "awk -v OFS='\t' '{print $1, $3+100, $3+299}' %s > %s"%('select_peak.bed', flag + '_' + location + '.bed')
            run_cmd(commandline)
            sequences.append(flag + '_' + location + '.bed')
            Info("get three regions of every peak around the %s genes"%flag)
            
            print(sequences) 
            for sequence in sequences:
                outputfasta = sequence.replace('.bed','.fa')
                self.fastas.append(outputfasta)
                runfastafrombed(self.genomesequence, sequence, outputfasta)
                run_cmd('rm %s'%sequence)
            Info("get the fasta format sequence data of the three regions of %s"%flag)
            i += 1
        run_cmd('rm select_peak.bed')
        
    ######################################
    #part2 Motif Scan and score caculate
    ######################################
            
    def run_mis(self):
        '''mis algorithm is refered of MOODS: Motif Occurrence Detection Suite, we write it into a c code to
           improve the efficiency, before run mis, the make command is required'''
        
        Info("run mis to do the known motif scan with cistrome motif database")
        #misdir = resource_filename('beta','mis')
        #os.chdir(misdir)
        #Usage: ./mis <in.seq> <in.db> <p-value> <motif-id> <output-prefix>
        self.motifscore = []
        
        for seq in self.fastas:
            inseq = os.path.join(self.curdir,seq)
            db = resource_filename('beta','references/cistrome.db')
            p_value = 0.001
            motif_ID = 'all'
            prefix = seq.replace('.fa','')
            misout = os.path.join(self.curdir,prefix)
            commandline = 'misp %s %s %s %s %s'%(inseq,db,p_value,motif_ID,misout)#run mis
            run_cmd(commandline)
            run_cmd('rm %s'%inseq)
            scoref = prefix + '_all'
            self.motifscore.append(scoref)

    def statistical_test(self):
        #we are intersted in the motif enriched in the up regulated genes
        #down regulated genes and the up versus down genes
        #curdir = os.getcwd()
       
        Info("T statistic test performed here to do the significance testing of every motif")
        groups = []
        pairfns = []
        pairs = []
        for group in self.selected:
            if group == '"upregulate"':
                upgroup = [f for f in self.motifscore if f.startswith('up')]
                groups.append(upgroup)
            if group == '"downregulate"':
                downgroup = [f for f in self.motifscore if f.startswith('down')]
                groups.append(downgroup)
        nongroup = [f for f in self.motifscore if f.startswith('non')]
        groups.append(nongroup)               
        if len(groups) == 3:
            #if this factor has both active and repressive funcation, we will scan the motif of up vs down.
            pairfn = "upvsdown_motif.txt"
            pairfns.append(pairfn)
            pairf = open(pairfn,'w')
            pairs.append(pairf)
            pairf.write('MotifID\tSpecies\tSymbol\tDNA BindDom\\tPSSM\tTscore\tPvalue\n')
        else:
            pairfn = "upvsdown_motif.txt"
            pairf = ''
        middlescores = []
        
        a = MP.MotifParser()
        motifdb = resource_filename('beta','references/cistrome.xml')
        a.Parser(motifdb)
        motifinfo = a.motifs

        self.FNs = []#the motif output file
        upnonf = ''
        downnonf = ''
        upnon = "up_non_motif.txt"
        downnon = "down_non_motif.txt"
        
        for group in groups:
            motifscore = {}
            if group[0].startswith('up'):
                fn = "up_motif.txt"               
                outf = open(fn,'w')
                upnonf = open(upnon,'w')
                upnonf.write('MotifID\tSpecies\tSymbol\tDNA BindDom\tPSSM\tTscore\tPvalue\n')
            if group[0].startswith('down'):
                fn = "down_motif.txt"
                outf = open(fn,'w')
                downnonf = open(downnon,'w')
                downnonf.write('MotifID\tSpecies\tSymbol\tDNA BindDom\tPSSM\tTscore\tPvalue\n')
            if group[0].startswith('non'):
                fn = "non_motif.txt"
                outf = open(fn,'w')
                
            outf.write('MotifID\tSpecies\tSymbol\tDNA BindDom\tPSSM\tTscore\tPvalue\n')
            group.sort() #have the order like: left,middleright
            if len(group) != 3:
                Info('MISP step wronging!')
                sys.exit(1)
            for f in group:
                f = os.path.join(self.curdir, f)
                with open(f) as scoref:
                    lines = scoref.readlines()
                    for line in lines:
                        if line.startswith('#'):
                            continue
                        elif line.startswith('@'):
                            #@ factor:EN0055
                            index = lines.index(line)
                            line = line.strip()
                            line = line.split('\t')
                            key = line[0].split(':')[1]
                            scoreline = lines[index + 2]
                            #score line looks like 0,0.20(102),0.23(29),0,0,3.34(-395),
                            scoreline = scoreline.strip()
                            scoreline = scoreline.split(',')[:-1]
                            value = []
                            for score in scoreline:
                                score = score.split('(')[0]
                                value.append(score)
                            try:
                                motifscore[key].append(value)
                            except:
                                motifscore[key] = [value]
                        else:
                            continue
                run_cmd('rm %s'%f)
                
            motifs = list(motifscore.keys())
            mscore = []
            for motif in motifs:
                species = motifinfo[motif]['species']
                if len(species) == 1 and species[0] != self.Species:
                    continue
                else:
                    species = ', '.join(species) #list to string
                    scores = motifscore[motif]
                    leftscore = 'left <- c('
                    middlescore = 'middle <- c('
                    rightscore = 'right <- c('
                    string = [leftscore, middlescorerightscore]
                    rtestf = open('temptest.r','w')
                    rtestf.write('options(warn=-1)\n')
                    for i in range(3):
                        for s in scores[i]:
                            string[i] += str(s)
                            string[i] += ','
                        string[i] = string[i].rstrip(',')
                        string[i] += ')'
                        rtestf.write(string[i] + '\n')
                    mscore.append(string[1])
                    rtestf.write('pvalue=1\n')
                    rtestf.write('summaryml = t.test(middle,left,alternative="greater")\n')
                    rtestf.write('mlpvalue = summaryml$p.value\n')
                    rtestf.write('mltscore = summaryml$statistic\n')
                    rtestf.write('summarymr= t.test(middleright,alternative="greater")\n')
                    rtestf.write('mrpvalue= summarymr$p.value\n')
                    rtestf.write('mrtscore = summarymr$statistic\n')
                    rtestf.write('pvalue = max(mlpvalue, mrpvalue)\n')
                    rtestf.write('tscore = min(mltscore, mrtscore)\n')
                    rtestf.write('print(pvalue)\n')
                    rtestf.write('print(tscore)\n')
                    rtestf.close()
                    
                    cmd = 'Rscript temptest.r'
                    #testinfo = os.popen(cmd).read().strip().split('\n')
                    
                    if not os.popen(cmd).read():
                        pvalue = 'NaN'
                        tscore = 'NaN'
                    else:
                        testinfo = os.popen(cmd).read().strip().split('\n')
                        pvalue = testinfo[0].split()[1]
                        tscore = testinfo[-1].split()[-1]
                    symbol = ', '.join(motifinfo[motif]['symbol'])
                    #TFs = motifinfo[motif]['symbol']
                    DBD = ', '.join(motifinfo[motif]['dbd'])
                    description = ', '.join(motifinfo[motif]['description'])
                    synonym = ', '.join(motifinfo[motif]['synonym'])
                    pssm = str(motifinfo[motif]['pssm'])
                    if pvalue == 'NaN':
                        pvalue = 1000
                    if tscore == 'NaN':
                        tscore = 0
                    #print motif, species, symbol, DBD, description, synonym, pssm, tscore, pvalue
                    outf.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(motif, species, symbol, DBD, pssm, str("%.2f"%float(tscore)), str("%.2e"%float(pvalue))))
            outf.close()
            FN = fn.split('.')[0].upper() + 'S.' + fn.split('.')[1] #the sorted motif file
            run_cmd('sort -t "\t" -k7,7g %s > %s'%(fn,FN))
            self.FNs.append(FN)
            run_cmd('rm %s'%fn)
            middlescores.append(mscore)
        
        run_cmd('rm temptest.r')
        pairfns = [upnon, downnon, pairfn]
        pairs = [upnonf, downnonf, pairf]
        
        k = 0
        for f in pairs:
            
            if f == '':
                pass
            else:
                if f == upnonf:
                    middlescore1 = middlescores[0]
                    middlescore2 = middlescores[-1]
                if f == downnonf:
                    middlescore1 = middlescores[-2]
                    middlescore2 = middlescores[-1]
                if f == pairf:
                    middlescore1 = middlescores[0]
                    middlescore2 = middlescores[1]
                
                for i in range(len(middlescores[0])):
                    motif = motifs[i]
                    species = motifinfo[motif]['species']
                    if len(species) == 1 and species[0] != self.Species:
                        pass
                    else:
                        species = ', '.join(species) #list to string
                        pairtest = open('pairtest.r','w')
                        pairtest.write('options(warn=-1)\n')                      
                        middlescore = middlescore1[i].replace('middle','middle1')
                        pairtest.write(middlescore + '\n')
                        middlescore = middlescore2[i].replace('middle','middle2')
                        pairtest.write(middlescore + '\n')
                        pairtest.write('summary = t.test(middle1,middle2)\n')
                        pairtest.write('pvalue = summary$p.value\n')
                        pairtest.write('tscore= summary$statistic\n')
                        pairtest.write('print(pvalue)\n')
                        pairtest.write('print(tscore)')
                        pairtest.close()
                        cmd = 'Rscript pairtest.r'
                        testinfo = os.popen(cmd).read().strip().split('\n')
                        
                        if not os.popen(cmd).read():
                            pvalue = 'NaN'
                            tscore = 'NaN'
                        else:
                            testinfo = os.popen(cmd).read().strip().split('\n')
                            pvalue = testinfo[0].split()[1]
                            tscore = testinfo[-1].split()[-1]
            
                        if pvalue == 'NaN':
                            pvalue = 1000
                        if tscore == 'NaN':
                            tscore = 0
                        description = ', '.join(motifinfo[motifs[i]]['description'])
                        symbol = ', '.join(motifinfo[motifs[i]]['symbol'])
                        DBD = ', '.join(motifinfo[motifs[i]]['dbd'])
                        synonym = ', '.join(motifinfo[motifs[i]]['synonym'])
                        pssm = motifinfo[motifs[i]]['pssm']
                        f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(motifs[i], species, symbol, DBD, pssm, str("%.2f"%float(tscore)), str("%.2e"%float(pvalue))))
                f.close()
                fn = pairfns[k]
                PAIRFN = fn.split('.')[0].upper() + 'S.' + fn.split('.')[1] 
                run_cmd('sort -t "\t" -k7,7g %s > %s'%(fn, PAIRFN))
                self.FNs.append(PAIRFN)
                run_cmd('rm %s'%fn)
                run_cmd('rm pairtest.r')
            k += 1
            
    def out2html(self):
        #write the top 10 motif into the html format

        Info('motif result will be in html format')
        FNs = self.FNs
        template = '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\n\
                    <html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">\n\
                        <head>\n\
                            <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />\n\
                            <title>BETA Motif Part</title>\n \
                            <link rel="stylesheet" type="text/css" href="styles.css" />\n\
                        </head>\n\
                        <body>\n\
                            <div class="section" id="page"> <!-- Defining the #page section with the section tag -->\n\
                                <div class="header"> <!-- Defining the header section of the page with the appropriate tag -->\n\
                                    <h1>BETA: Motif Analysis</h1>\n\
                                    <h3>Motif Scan on the TF Target Genes</h3>\n\
                                </div>\n\
                                <div class="section" id="articles"> <!-- A new section with the articles -->\n\
                                    <div class="line"></div>  <!-- Dividing line -->\n\
                                    <div class="article" id="article1"> <!-- The new article tag. The id is supplied so it can be scrolled into view. -->\n\
                                        <h2>Part1: Up Target Genes</h2>\n\
                                        <div class="line"></div>\n\
                                        <div class="articleBody clear">\n\
                                            <table style="text-align: left; width: 100%%" cellpadding="2" cellspacing="2" border="1" style="border: 1px solid #000;" >\n\
                                                %s\n\
                                                %s\n\
                                            </table>\n\
                                       </div>\n\
                                    </div>\n\
                                    <div class="line"></div>\n\
                                    <div class="article" id="article2">\n\
                                        <h2>Part2: Down Target Genes</h2>\n\
                                        <div class="line"></div>\n\
                                        <div class="articleBody clear">\n\
                                           <table style="text-align: left; width: 100%%" cellpadding="2" cellspacing="2" border="1" style="border: 1px solid #000;" >\n\
                                                %s\n\
                                                %s\n\
                                           </table>\n\
                                        </div>\n\
                                    </div>\n\
                                    <div class="line"></div>\n\
                                    <div class="article" id="article3">\n\
                                        <h2>Part3: UP vs DOWN Motif Scan</h2>\n\
                                        <div class="line"></div>\n\
                                        <div class="articleBody clear">\n\
                                           <table style="text-align: left; width: 100%%" cellpadding="2" cellspacing="2" border="1" style="border: 1px solid #000;" >\n\
                                                %s\n\
                                                %s\n\
                                           </table>\n\
                                        </div>\n\
                                    </div>\n\
                                    <div class="line"></div>  <!-- Dividing line -->\n\
                                    <div class="article" id="article1"> <!-- The new article tag. The id is supplied so it can be scrolled into view. -->\n\
                                        <h2>Part4: Up VS Non Target Motif</h2>\n\
                                        <div class="line"></div>\n\
                                        <div class="articleBody clear">\n\
                                            <table style="text-align: left; width: 100%%" cellpadding="2" cellspacing="2" border="1" style="border: 1px solid #000;" >\n\
                                                %s\n\
                                                %s\n\
                                            </table>\n\
                                       </div>\n\
                                    </div>\n\
                                    <div class="line"></div>  <!-- Dividing line -->\n\
                                    <div class="article" id="article1"> <!-- The new article tag. The id is supplied so it can be scrolled into view. -->\n\
                                        <h2>Part5: Down VS Non Target Motif</h2>\n\
                                        <div class="line"></div>\n\
                                        <div class="articleBody clear">\n\
                                            <table style="text-align: left; width: 100%%" cellpadding="2" cellspacing="2" border="1" style="border: 1px solid #000;" >\n\
                                                %s\n\
                                                %s\n\
                                            </table>\n\
                                       </div>\n\
                                    </div>\n\
                                </div>\n\
                                <div class="footer"> <!-- Marking the footer section -->\n\
                                   <div class="line"></div>\n\
                                       <p>BETA: Binding and Expression Target Analysis</p> <!-- Change the copyright notice -->\n\
                                       <a href="#" class="up">Go UP</a>\n\
                               </div>\n\
                            </div> <!-- Closing the #page section -->\n\
                        <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.3.2/jquery.min.js"></script>\n\
                        <script type="text/javascript" src="jquery.scrollTo-1.4.2/jquery.scrollTo-min.js"></script>\n\
                        <script type="text/javascript" src="script.js"></script>\n\
                        </body>\n\
                    </html>'

        outhtml = open('betamotif.html','w')
        resultdir = os.path.join(self.outdir,'motifresult')
        if os.path.isdir(resultdir):
            run_cmd('rm -r %s'%resultdir)
            run_cmd('mkdir %s'%resultdir)
        else:
            run_cmd('mkdir %s'%resultdir)
        imgdir = os.path.join(resultdir,'img')
        if os.path.isdir(imgdir):
            run_cmd('rm -r %s'%imgdir)
            run_cmd('mkdir %s'%imgdir)
        else:
            run_cmd('mkdir %s'%imgdir)
        
        tabletitle = '\n\
                                    <tr>\n\
                                        <th>\n\
                                            <font size="4" color="yellow">Symbol</font>\n\
                                        </th>\n\
                                        <th>\n\
                                            <font size="4" color="yellow">DNA BindDom</font>\n\
                                        </th>\n\
                                        <th>\n\
                                            <font size="4" color="yellow">Species</font>\n\
                                        </th>\n\
                                        <th>\n\
                                            <font size="4" color="yellow">Pvalue (T Test)</font>\n\
                                        </th>\n\
                                        <th>\n\
                                            <font size="4" color="yellow">T Score</font>\n\
                                        </th>\n\
                                        <th>\n\
                                            <font size="4" color="yellow">Logo</font>\n\
                                        </th>\n\
                                     </tr>'
        nontable = '<p></p>'
        subtemplate = '\n\
                                    <tr>\n\
                                        <td>\n\
                                            %s\n\
                                        </td>\n\
                                        <td>\n\
                                            %s\n\
                                        </td>\n\
                                        <td>\n\
                                            <font size="4" >%s</font>\n\
                                        </td>\n\
                                        <td>\n\
                                            <font size="4" >%s</font>\n\
                                        </td>\n\
                                        <td>\n\
                                            <font size="4" >%s</font>\n\
                                        </td>\n\
                                        <td>\n\
                                            <img width="500" height="220" border="0" align="left" alt="some_text" src=%s>\n\
                                        </td>\n\
                                    </tr>'
        cluster = '\n                                        <p><font size="4">%s</font></p>\n'

        if 'UP_MOTIFS.txt' in FNs:  
            motif1 = 'UP_MOTIFS.txt'#upmotif
            motif2 = 'UP_NON_MOTIFS.txt'#upvsnonmotif
            mistable1 = 'upmis.txt'#up motif mis table
            mistable2 = 'upnonmis.txt'#upvsnon mis table
            motifs = [motif1,motif2]
            mistables = [mistable1,mistable2]
            subups = []
            
            for q in range(2):
                motif = motifs[q]
                mistable = mistables[q]
                subup = ''
                output = runmotifcluster_side1(self.outdir, motif, mistable, self.motifnumber)
                output.sort(key=lambda x:float(x[3]))
                for motif_class in output:
                    factors = motif_class[0].split('; ')
                    factorinput = ''
                    for factor in factors:
                        factorinput += cluster%factor
                    dbdinput = ''
                    dbds = motif_class[1].split('; ')
                    for dbd in dbds:
                        dbdinput += cluster%dbd                    
                    species = motif_class[2]
                    pvalue = motif_class[3]
                    tscore = motif_class[4]
                    img = motif_class[-1]
                    run_cmd('cp %s %s'%(img, imgdir))
                    imgID = os.path.split(img)[1]
                    img = os.path.join('img',imgID)                  
                    temp = subtemplate%(factorinput,dbdinput,species,pvalue,tscore,img)
                    subup += temp
                subups.append(subup)
            tabletitleup = tabletitle
            
            run_cmd('mv %s %s'%('UP_MOTIFS.txt', resultdir))
            run_cmd('mv %s %s'%('UP_NON_MOTIFS.txt', resultdir))
            run_cmd('mv %s %s'%(self.name + '_uptarget.txt', self.outdir))
        else:
            subups = [nontable,nontable]
            tabletitleup = '<font size="5" color="yellow">This Factor has no Active Function to its target genes, Skipped this part Motif Scan</font>'
            subcompare = nontable
            tabletitlecompare = '<font size="5" color="yellow">This Factor has no Active Function to its target genes, Cannot do the comparision Motif Scan</font>'

        if 'DOWN_MOTIFS.txt' in FNs:
            motif1 = 'DOWN_MOTIFS.txt'#upmotif
            motif2 = 'DOWN_NON_MOTIFS.txt'#upvsnonmotif
            mistable1 = 'downmis.txt'#up motif mis table
            mistable2 = 'downnonmis.txt'#upvsnon mis table
            motifs = [motif1,motif2]
            mistables = [mistable1,mistable2]
            subdowns = []
            for q in range(2):
                subdown = ''
                motif = motifs[q]
                mistable = mistables[q]
                output = runmotifcluster_side1(self.outdir, motif, mistable, self.motifnumber)
                output.sort(key=lambda x:float(x[3]))
                for motif_class in output:
                    factors = motif_class[0].split('; ')
                    factorinput = ''
                    for factor in factors:
                        factorinput += cluster%factor
                    dbdinput = ''
                    dbds = motif_class[1].split('; ')
                    for dbd in dbds:
                        dbdinput += cluster%dbd
                    species = motif_class[2]
                    pvalue = motif_class[3]
                    tscore = motif_class[4]
                    img = motif_class[-1]
                    run_cmd('cp %s %s'%(img, imgdir))
                    imgID = os.path.split(img)[1]
                    img = os.path.join('img',imgID)
                    temp = subtemplate%(factorinput,dbdinput,species,pvalue,tscore,img)
                    subdown += temp
                subdowns.append(subdown)
            tabletitledown = tabletitle
            run_cmd('mv %s %s'%('DOWN_MOTIFS.txt', resultdir))
            run_cmd('mv %s %s'%('DOWN_NON_MOTIFS.txt', resultdir))
            run_cmd('mv %s %s'%(self.name + '_downtarget.txt', self.outdir))
        else:
            subdowns = [nontable,nontable]
            tabletitledown = '<font size="5" color="yellow">This Factor has no Repressive Function to its target genes, Skipped this part Motif Scan</font>'
            subcompare = nontable
            tabletitlecompare = '<font size="5" color="yellow">This Factor has no Repressive Function to its target genes, Cannot do the comparision Motif Scan</font>'
                
        if 'UPVSDOWN_MOTIFS.txt' in FNs:
            motif = 'UPVSDOWN_MOTIFS.txt'
            mistable1 = 'upvsdownmis1.txt'
            mostable2 = 'upvsdownmis2.txt'
            output = runmotifcluster_side2(self.outdir, motif, mistable1, mistable2, self.motifnumber)
            output.sort(key=lambda x:float(x[3]))
            subcompare = ''
            for motif_class in output:
                factors = motif_class[0].split('; ')
                factorinput = ''
                for factor in factors:
                    factorinput += cluster%factor
                dbdinput = ''
                dbds = motif_class[1].split('; ')
                for dbd in dbds:
                    dbdinput += cluster%dbd
                species = motif_class[2]
                pvalue = motif_class[3]
                tscore = motif_class[4]
                img = motif_class[-1]
                run_cmd('cp %s %s'%(img, imgdir))
                imgID = os.path.split(img)[1]
                img = os.path.join('img',imgID)
                temp = subtemplate%(factorinput,dbdinput,species,pvalue,tscore,img)
                subcompare += temp
            
            tabletitledown = tabletitle
            run_cmd('mv %s %s'%(self.name + '_downtarget.txt', self.outdir))
            tabletitlecompare = tabletitle
            run_cmd('mv %s %s'%('UPVSDOWN_MOTIFS.txt', 'DIFFERENTIAL_MOTIF_UP_DOWN.txt'))
            run_cmd('mv %s %s'%('DIFFERENTIAL_MOTIF_UP_DOWN.txt', resultdir))
                   
        b = template%(tabletitleup,subups[0],tabletitledown,subdowns[0],tabletitlecompare,subcompare,tabletitleup,subups[1],tabletitledown,subdowns[1])

        outhtml.write(b)
        outhtml.close()
 
        jsfile = resource_filename('beta','templates/script.js')
        cssfile = resource_filename('beta','templates/styles.css')
        run_cmd('rm %s'%'NON_MOTIFS.txt')
        run_cmd('cp %s %s'%(jsfileresultdir))
        run_cmd('cp %s %s'%(cssfileresultdir))
        motiflogos = os.path.join(self.outdir,'motiflogos')
        run_cmd('rm -rf %s'%motiflogos)
        run_cmd('mv betamotif.html %s'%resultdir)
        Info("Done: find motif result in beatmotif.html file")
        Info("Done: find all BETA result in %s"%self.outdir)
