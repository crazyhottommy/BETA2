########################
##        PartI         ##
##########################
from pkg_resources import resource_filename
from .corelib import *
import math

chroms = ['chrY', 'chrX', 'chr13', 'chr12', 'chr11', 'chr10', 'chr17',\
          'chr16', 'chr15', 'chr14', 'chr19', 'chr18', 'chr22', 'chr20',\
          'chr21', 'chr7', 'chr6', 'chr5', 'chr4', 'chr3', 'chr2', 'chr1',\
          'chr9', 'chr8']

#Score calc function
Sg = lambda ldx: sum([math.exp(-0.5-4*t) for t in ldx])


class PScore:
    
    #Caculate every gene's regulate potential score, peaks within 100kb and in the same CTCF region will be cosidered.
    
    def __init__(self, options):

        self.peakfile = options.peakfile
        if options.genome == 'hg19':
            if options.boundarylimit == False:
                Info('You do not like filter peak by CFCT boundary, it will be filtered only by the distance')
                self.boundaryfile = ''
            else:
                if not options.boundaryfile:
                    self.boundaryfile = resource_filename('beta', 'references/hg19_CTCF_bound.bed')
                else:
                    self.boundaryfile = options.boundaryfile
            if not options.reference:
                self.genome = resource_filename('beta','references/hg19.refseq')
            else:
                self.genome = options.reference

        if options.genome == 'mm9':
            if options.boundarylimit == False:
                Info('You do not like filter peak by CFCT boundary, it will be filtered only by the distance')
                self.boundaryfile = ''
            else:
                if not options.boundaryfile:
                    self.boundaryfile = resource_filename('beta', 'references/mm9_CTCF_bound.bed')
                else:
                    self.boundaryfile = options.boundaryfile
            if not options.reference:
                self.genome = resource_filename('beta','references/mm9.refseq')
            else:
                self.genome = options.reference
                
        if options.genome == 'hg38':
            if options.boundarylimit == False:
                Info('You do not like filter peak by CFCT boundary, it will be filtered only by the distance')
                self.boundaryfile = ''
            else:
                if not options.boundaryfile:
                    Info("We didn't provide a CTCF boundary file in hg38 version, please give the boundary file via --bf")
                    sys.exit(1)
                else:
                    self.boundaryfile = options.boundaryfile
            if not options.reference:
                self.genome = resource_filename('beta','references/hg38.refseq')
            else:
                self.genome = options.reference

        if options.genome == 'mm10':
            if options.boundarylimit == False:
                Info('You do not like filter peak by CFCT boundary, it will be filtered only by the distance')
                self.boundaryfile = ''
            else:
                if not options.boundaryfile:
                    Info("We didn't provide a CTCF boundary file in mm10 version, please give the boundary file via --bf")
                    sys.exit(1)
                else:
                    self.boundaryfile = options.boundaryfile
            if not options.reference:
                self.genome = resource_filename('beta','references/mm10.refseq')
            else:
                self.genome = options.reference
                
        if options.genome == 'hg18':
            if options.boundarylimit == False:
                Info('You do not like filter peak by CFCT boundary, it will be filtered only by the distance')
                self.boundaryfile = ''
            else:
                if not options.boundaryfile:
                    Info("We didn't provide a CTCF boundary file in hg18 version, please give the boundary file via --bf")
                    sys.exit(1)
                else:
                    self.boundaryfile = options.boundaryfile
            if not options.reference:
                self.genome = resource_filename('beta','references/hg18.refseq')
            else:
                self.genome = options.reference
                
        if options.genome == 'hg' or options.genome == 'mm' or (not options.genome):
            self.genome = options.reference
            if options.boundarylimit != False:
                self.boundaryfile = options.boundaryfile
            else:
                Info('We don not provide a CTCF boundary file for %s, the peak will be filtered only by the distance'%options.genome)
                self.boundaryfile = ''
        self.peaknumber = options.peaknumber
        self.outdir = options.output
        self.opts_string = "# Argument List:\n" +\
                           "# Name = %s\n" %options.name +\
                           "# peak file = %s\n" %options.peakfile +\
                           "# distance = %d bp\n" %options.distance 
 
        self.peaklist = {}
        
    def readfile(self): #reads the file and returns a vector: each element is a bed_row. 
        
        self.boundarylist = {}
        if self.boundaryfile != '':
            with open(self.boundaryfile) as boundaryf:
                boundaryf = boundaryf.readlines()
                for line in boundaryf:
                    if line.startswith('#') or not line.strip():
                        continue
                    line = line.strip()
                    line = line.split('\t')
                    line = [line[0], int(line[1]), int(line[2]), line[3], float(line[4]), line[5]]
                    #[chroms, tss, tts, name, score, strand]
                    try:
                        self.boundarylist[line[0]].append(line)
                    except:
                        self.boundarylist[line[0]] = [line]
        else:
            pass

        peaks = []
        with open(self.peakfile) as peakf:
            peakf = peakf.readlines()
            for peak in peakf:
                if peak.startswith("#") or not peak.strip():
                    continue
                peak = peak.strip()
                peak = peak.split('\t')
                peaks.append(peak)
            sortpeaks = sorted(peaks, key=lambda peaks:float(peaks[4]), reverse=True)
        
        if len(sortpeaks) > self.peaknumber:
            selectpeaks = sortpeaks[0:self.peaknumber]
        else:
            selectpeaks = sortpeaks
        count = 0
        self.peaklist = {}
        
        for line in selectpeaks:
            #.bed-> 0:chrom 1:pStart 2:pEnd 3:peakName 4:-log10(qvalue)
            line = [line[0], int(line[1]), int(line[2]), line[3], float(line[4])]
            try:
                self.peaklist[line[0]].append(line)
            except KeyError:
                self.peaklist[line[0]] = [line]
            count += 1
        #peakf.close()
            
        for i in list(self.peaklist.keys()):
            self.peaklist[i].sort()
        Info("Read file <%s> OK! All <%d> peaks." %(self.peakfile, count))

    def ScoreCalc(self, distance):
        #get at most 10k sites, p-value < 1e-5. each gene's regulatory potential sg = ...

        INFO = {}
        genepeaks = {}
        
        with open(self.genome,mode = 'r') as refgene:
            self.geneInfo = []
            for line in refgene:
                if line.startswith('#') or line == '\n' or not line:
                    continue
                else:
                    line = line.strip()
                    line = line.split('\t')
                    if line[1] in chroms:
                        info = [line[1], line[0], str(line[3]), str(line[4]), line[2], line[5]]
                        #info = ['chrom', 'refseq', 'start', 'end', 'strand', 'symbol']
                        self.geneInfo.append(info)
                    else:
                        continue
        count = 0
    
        for igene in self.geneInfo:
            if igene[4] == "+":
                #print igene[2]
                gTSS = int(igene[2])
                gTTS = int(igene[3])
            if igene[4] == "-":
                gTSS = int(igene[3])
                gTTS = int(igene[2])
            try:
                boundarys = self.boundarylist[igene[0]]
            except KeyError:
                boundarys = []
            upstreams = []
            downstreams = []
            for i in boundarys:
                if int(i[1]) > int(gTTS):
                    downstreams.append(int(i[1]))#boundary tss > gene tts is the upstream boundarys
                if int(i[2]) < int(gTSS):
                    upstreams.append(int(i[2]))#boundary tts < gene tss is the downstream boundarys
                else:
                    continue
            upstreams.sort()
            downstreams.sort()
            if len(upstreams) != 0:
                uplimit = upstreams[-1]
            else:
                uplimit = 0 
            if len(downstreams) != 0:
                downlimit = downstreams[0]
               
            else:
                downlimit = 250000000#bigger than chr1's length
            
            try:
                peaks = self.peaklist[igene[0]]
            except KeyError:
                peaks = []
            peaksInDistance = [abs((t[1]+t[2])/2-gTSS)*1.0/distance for t in peaks if abs((t[1]+t[2])/2-gTSS) < distance  and not t[2] < uplimit and not t[1] > downlimit] #peak pvalue<1e-5, and distance < distance
            peaksInDistance.sort()
            if len(peaksInDistance) > 10000: # extract no more than 10k peaks
                peaksInDistance = peaksInDistance[:10000]
            #score = sum(math.exp(-0.5-4*t[-1]) for t in peaksInDistance)
            igene.append(Sg(peaksInDistance))
            count += 1
            
            peakswithinDistance = [t for t in peaks if abs((t[1]+t[2])/2-gTSS) < distance  and not(t[2] < uplimit or t[1] >  downlimit)]
            #peaksoutofDistance = [t for t in peaks if abs((t[1]+t[2])/2-gTSS) > distance]
            if len(peakswithinDistance) > 10000: # extract no more than 10k peaks
                peakswithinDistance = peakswithinDistance[:10000]
            key = tuple(igene[:-1]) 
            genepeaks[key] = peakswithinDistance
           
        Info('Process <%d> genes'%count)
        self.geneInfo.sort(key=lambda x:x[-1], reverse=True)
        a = list(self.peaklist.values())#a=[[[chr1],[chr1],[chr1]],[[chr2],[chr2]]]
        totalPeaks = []
        for m in a:
            for j in m:
                if j in totalPeaks:
                    continue
                else:
                    totalPeaks.append(j)
        withinPeaks = []
        b = list(genepeaks.values())
        for n in b:
            for k in n:
                if k in withinPeaks:
                    continue
                else:
                    withinPeaks.append(k)
        peaksoutofDistance = []
        for i in totalPeaks:
            if i in withinPeaks or i == []:
                continue
            else:
                peaksoutofDistance.append(i)
        return (genepeaks, peaksoutofDistance)
    
    def Output2File(self, name):
            
        outf = open("%s.txt"%name, "w")#peaks score and rank file
        outf.write(self.opts_string)
        outf.write('#chrom\ttxStart\ttxEnd\trefseq\tscore\tstrand\tsymbol\trank\n')
        r = 1
        for line in self.geneInfo:
        
            if str('%.3f'%line[6]) == '0.000':
                #if one gene's score is zero, this gene will not be used to rank
                outf.write('%s\t%d\t%d\t%s\t%.3f\t%s\t%s\t%s\n'%(line[0], int(line[2]), int(line[3]), line[1], line[6], line[4], line[5], 'NA'))
            else:
                outf.write('%s\t%d\t%d\t%s\t%.3f\t%s\t%s\t%d\n'%(line[0], int(line[2]), int(line[3]), line[1], line[6], line[4], line[5], r))
            r += 1
        outf.close()
        Info("Finished! Preliminary results saved into temporary file: <%s.txt>"%name)
        
        #return r

    def noexpreoutput(self, name):
        #['chrom', 'refseq', 'start', 'end', 'strand', 'symbol']
        outf = open("%s_targets.txt"%name, "w")#peaks score and rank fileh
        outf.write(self.opts_string)
        outf.write('#Chromsome\tTSS\tTTS\tRefseqID\tScore\tStrand\tGeneSymbol\n')
 
        for line in self.geneInfo:               
            if float(line[6]) != 0:
                outf.write('%s\t%s\t%s\t%s\t%.3f\t%s\t%s\n'%(line[0],str(line[2]),str(line[3]), line[1], float(line[6]),line[4], line[5]))
            else:
                pass
        outf.close()
        run_cmd('mv %s_targets.txt %s'%(name, self.outdir))
        Info("Finished! result output to <%s_targets.txt>"%name)
        
    def noexprepeaks(self, name, genepeaks, distance):

        score = lambda t: math.exp(-0.5-4*t)
        outf = open("%s_targets_associated_peaks.txt"%name, "w")
        outf.write('chrome\tstart\tend\tpeak\trefseqID\tGene_Symbol\tdistance\tscore\n')
        keys = list(genepeaks.keys())
        for gene in keys:
            if gene[4] == "+":
                gTSS = int(gene[2])
            if gene[4] == "-":
                gTSS = int(gene[3])
            for peak in genepeaks[gene]:
                d = (int(peak[2]) + int(peak[1]))/2 - gTSS
                ps = score(abs(float(d))/float(distance))
                outf.write(peak[0] + '\t' + str(peak[1]) + '\t' + str(peak[2]) + '\t' + peak[3] + '\t' + gene[1] + '\t' + gene[5] + '\t' + str(d) + '\t' + str(ps) + '\n')
        run_cmd("mv %s_targets_associated_peaks.txt %s"%(name, self.outdir))
        outf.close()
