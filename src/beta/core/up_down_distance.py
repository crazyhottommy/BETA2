#This script do the TF/CR activating or repressive function prediction by the nearest distance to the binding peaks.
import os,sys
import math
import re, subprocess
# from optparse import OptionParser  # Deprecated, use argparse

CHROM_CONVERT = {'chrI':'chr1','chrII':'chr2','chrIII':'chr3','chrIV':'chr4','chrV':'chr5','chrVI':'chr6',
                 'chrVII':'chr7','chrVIII':'chr8','chrIX':'chr9','chrX':'chr10','chrXI':'chr11','chrXII':'chr12',
                 'chrXIII':'chr13','chrXIV':'chr14','chrXV':'chr15','chrXVI':'chr16','chrXVII':'chr17','chrXVIII':'chr18',
                 'chrXIX':'chr19','chrXX':'chr20','chrXXI':'chr21','chrXXII':'chr22','chrX':'chrX','chrY':'chrY','chrM':'chrM'}
chroms = list(CHROM_CONVERT.values())

def prepare_optparser():
    
    """Prepare optparser object. New options will be added in this
    function first.
    
    """
    
    usage = "usage: %prog <-g gdb -b bed -g geneset> [options]"
    description = "ChGS -- Gene Set Association Analysis. Note that if the number of genomic coordinates is below 1000 or the number of a gene set is below 500, the p value of ChGS may not represent biological implications well enough."
    
    optparser = OptionParser(version="%prog 0.5",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-b","--bed",dest="bed",type="string",\
                         help="BED file of genomic coordinates (e.g. ChIP-Seq peaks). The center of each peak will be used to compute the distance from a gene.")
    optparser.add_option("-r","--refseq",dest="refseq",type="string",\
                         help="Gene annotation table. This gene annotation table must be the UCSC refgene format, that is name, chrom, strand, txstart, txend, names2, This is required")
    optparser.add_option("-g", "--gset", dest="gset", action="append", type="string", help="Gene set to see the association with the genomic coordinates given through -b. Multiple gene sets can be given by repeatedly using this option (e.g. -g geneset1 -g geneset2 -g geneset3). Refseq accession numbers or official gene symbols can be used.")
    optparser.add_option("--bg", dest="bg",type = "string", help = "the background genes file you want to set as a control, make it as a gene list file as the gset you input.Default = ALL", default="ALL")
    optparser.add_option("--name",dest="name", help="Experiment name. This will be used to name the output file. If an experiment name is not given, input BED file name will be used instead.")
    optparser.add_option("-l", "--lab", dest="label", action="append", type="string", help="Label for each gene set. Likewise, multiple gene set labels can be given by repeatedly using this option (e.g. -l label1 -l label2 -l label3). If labels are not given, 'gene set' will be used by default.", default=None)
    optparser.add_option("--gname2",dest="name2",action="store_true", \
                         help="If this switch is on, gene or transcript IDs in files given through -g will be considered as official gene symbols.", default=False)                           
    return optparser

def read_gn_group( fn ):
    """
    Read a file of gene group.
    """

    genes = []
    for line in open(fn, 'r'):

        if line == "\n": continue
        if not line: continue

        genes.append( line.strip().split("\t")[0] )

    return genes


def read_gene_sets( fns ):
    """
    Read a series of gene group files.
    
    fns must be a list containing gene group file names.
    """

    gene_sets = []
    for fn in fns:
        gene_set = list( set( read_gn_group( fn ) ) )
    
        gene_sets.append( gene_set )

    return gene_sets


def readfile(peakfile):
    
    peakf = open(peakfile)
    peaklist = {}
    
    for line in peakf:
        
        if line.startswith('#') or not line.strip():
            continue
        line = line.split() #chrom, start, end, name, score
        line = [line[0], int(line[1]), int(line[2]), line[3], float(line[4])]

        try:
            line[0] = CHROM_CONVERT[line[0]]
        except KeyError:
            pass
        
        try:
            peaklist[line[0]].append(line)  #peaklist = {chr1:[[start,end,name,score],[],[]...], chr2:[[],[],[]...]...}
        except KeyError:
            peaklist[line[0]] = [line]
            
    peakf.close()
    
    for i in list(peaklist.keys()):
        peaklist[i].sort() #sort the peak by their start position in each chromsome

    return peaklist



def get_the_nearest_dist(refseq,peaklist):
    
    #outf = open('refinfo.txt','w')
    '''
    db = sqlite3.connect(genome)
    
    c = db.cursor()
    
    sql = "select chrom, name, txStart, txEnd, strand, name2 from GeneTable order by chrom;"
    c.execute(sql)
    geneInfo = c.fetchall()
    geneInfo.sort()
    geneInfo = [list(t) for t in geneInfo] #[u'chr1', u'NM_000016', 76190042, 76229353, u'+', u'ACADM']
    '''
    refgene = open(refseq,'r')
    geneInfo = []
    for line in refgene:
        if line.startswith('#') or line == '\n' or not line:
            continue
        else:
            line = line.strip()
            line = line.split('\t')
            if line[1] in chroms:
                info = [line[1], line[0], str(line[3]), str(line[4]), line[2], line[5]]#['chr1', 'NM_032291', '66999824', '67210768', '+', 'SGIP1']
                geneInfo.append(info)
            else:
                continue
            
    refseqgene_dist = {}
    symbolgene_dist = {}
    
    for igene in geneInfo:
        #print igene
        #outf.write(igene[0] + '\t' + igene[1] + '\t' + str(igene[2]) + '\t' + str(igene[3]) + '\t' + str(igene[4]) + '\t' + igene[5] + '\n')
        gTSS =int(igene[2])
        
        try:
            peaks = peaklist[igene[0]]
        except KeyError:
            peaks = []
        #print peaks[0]           
        peaksInDistance = [math.log(abs((t[1] + t[2])/2 - gTSS) + 0.0001,10) for t in peaks]
        peaksInDistance.sort()
        
        if len(peaksInDistance) == 0 and igene[1] in chroms:
            nearest_dist = 'NA'
            print(igene[1])
            #print igene
        elif len(peaksInDistance) == 0 and igene[1] not in chroms:
            #print igene[1]
            pass
        else:
            if peaksInDistance[0] >= 0:
                nearest_dist = peaksInDistance[0]
            elif peaksInDistance[0] < 0:
                nearest_dist = 0
        
            try:
                refseqgene_dist[igene[1]].append(nearest_dist)
                symbolgene_dist[igene[5]].append(nearest_dist)
            except KeyError:
                refseqgene_dist[igene[1]] = [nearest_dist]
                symbolgene_dist[igene[5]] = [nearest_dist]

    refdist = {}
    symdist = {}
    
    
    for gene in list(refseqgene_dist.keys()):
        refseqgene_dist[gene].sort()
        refdist[gene] = refseqgene_dist[gene][0]
        
    for gene in list(symbolgene_dist.keys()):
        symbolgene_dist[gene].sort()
        symdist[gene] = symbolgene_dist[gene][0]
    #outf.close()
    refgene.close()
    return (refdist, symdist)


def assign_geneset_dist(gene_sets,dist):

 
    geneset_dists = []
    
    for geneset in gene_sets:
        geneset_dist = []
        for gene in geneset:
            if gene in list(dist.keys()):
                distance = dist[gene]
                geneset_dist.append(distance)
            else:
                pass

        geneset_dists.append(geneset_dist)

    return geneset_dists

def write_R_code(alldist, geneset_dists, label, name):
    
    rscript = ''
    rscript += '\n'
    rscript += 'options(warn = -1)\n'
    rscript += 'cr <- colorRampPalette(col=c("#C8524D", "#BDD791", "#447CBE", "#775A9C"), bias=1)\n'
    rscript += 'col <- cr(2)\n'
    rscript += 'dist <- list()\n'
    bg = 'd <- as.numeric(c(%s'%alldist[0]
    n = len(alldist)
    for i in range(n-2):
        bg += ',%s'%alldist[i + 1]
    bg += ',%s))\n'%alldist[-1]
    rscript += bg
    rscript += "bg <- d\n"
    rscript += "bg[bg == 'NA'] <- NA\n"
    
    for l in range(len(label)):
        group = 'd <- as.numeric(c(%s'%geneset_dists[l][0]
        for j in range(len(geneset_dists[l]) - 2):
            group += ',%s'%geneset_dists[l][j + 1]
        group += ',%s))\n'%geneset_dists[l][-1]
        rscript += group
        rscript += 'dist[["%s"]] <- d\n'%label[l]
    rscript += "\n"
    rscript += "vmax <- max(bg)\n"
    rscript += "uplimit <- round(vmax) + 1\n"
    rscript += "\n"
    rscript += "breaks <- seq(0, uplimit, by = 0.05)\n"
    rscript += "hs <- list()\n"

    for i in range(len(label)):
        rscript += 'd <- as.vector(na.omit(dist[["%s"]]))\n'%label[i]
        rscript += 'd <- d[d >= 0 & d <= uplimit]\n'
        rscript += 'h <- hist(d, breaks = breaks, plot = FALSE)\n'
        rscript += 'hs[["%s"]] <- h\n'%label[i]
    
    rscript += "dbg <- as.vector(na.omit(bg))\n"
    rscript += "dbg <- dbg[dbg >= 0 & dbg <= uplimit]\n"
    rscript += "hbg <- hist(dbg, breaks=breaks, plot=FALSE)\n"
    rscript += "\n"
    rscript += "x <- h$mids\n"
    rscript += "chs <- list()\n"
    for i in label:
        rscript += 'c <- 100*cumsum(hs[["%s"]]$count)/sum(hs[["%s"]]$count)\n'%(i,i)
        rscript += 'chs[["%s"]] <- c\n'%i
    rscript += 'cbg <- 100* cumsum(hbg$count)/sum(hbg$count)\n'
    rscript += "\n"
    rscript += "p_values <- list()\n"
    for i in label:
        rscript += 'd <- as.vector(na.omit(dist[["%s"]]))\n'%i
        rscript += 'd <- d[d >= 0 & d <= uplimit]\n'
        rscript += 'ks <- ks.test(jitter(d), jitter(dbg), alternative = "greater")\n'
        rscript += 'p_values[["%s"]] <- ks$p.value\n'%i

    rscript += "\n"
    
    rscript += "ps = as.matrix(p_values)\n"
    rscript += "write.table(ps, file = '%s_pvalues.txt', sep = '\t', col.names = F)\n"%name

    rscript += 'pdf("%s_function_prediction.pdf", height = 6, width = 6)\n'%name
    rscript += 'plot(x, cbg, type="l", lty=2, lwd=2, col="black", xlab="Rank of genes based on log(10) distance to nearest peak", ylab="Cumulative Fraction of Genes%%", main="%s Activating/Repressive Function Prediction", xlim=c(0, uplimit), ylim=c(0, 100))\n'%name
    rscript += "\n"
    rscript += "l <- c()\n"
    for i in range(len(label)):
        rscript += 'lines(x, chs[["%s"]], lty=1, lwd=2, col=col[%d])\n'%(label[i],i+1)
        rscript += 'l <- c(l, paste("%s (", format(p_values[["%s"]], digits=3), ")", sep=""))\n'%(label[i],label[i])
    rscript += "legend('topleft', legend=l, col=col, pch=15)\n"
    rscript += "dev.off()"
        
    return rscript            
    
def main():
    optparser = prepare_optparser ()
    
    (options,args) = optparser.parse_args()
    if not options.bed or not options.gset:
        
        optparser.print_help()
        sys.exit(1)
        
    #genome = options.gdb
    peakfile = options.bed
    if type(options.gset) == str:
        options.gset = [options.gset]
    if options.label:
        if len(options.label) != len(options.gset):
            error("The number of the gene set labels (-l or --lab) must be the same as that of the gene sets (-g or --gset).")
            sys.exit(1)
    else:
        options.label =  ["gene set" + str(x) for x in range(1, len(options.gset)+1)]
    peaklist = readfile(peakfile)
    (refdist, symdist) = get_the_nearest_dist(options.refseq, peaklist)
    
    if options.name2 == False:
        dist = refdist
    else:
        dist = symdist
    #print dist       
    bgdist = []
    if options.bg == 'ALL':
        bgdist = list(dist.values())
    else:
        bggene = read_gn_group(options.bg)
        #print bggene
        for gene in bggene:
            if gene in list(dist.keys()):
                bgdist.append(dist[gene])
            else:
                continue
    
    fns = options.gset
    gene_sets = read_gene_sets(fns)    
    geneset_dists = assign_geneset_dist(gene_sets,dist)
    rscript = write_R_code(bgdist, geneset_dists, options.label, options.name)
    rcode = open('%s_function_prediction.R'%options.name,'w')
    rcode.write(rscript)
    rcode.close()
    try:
        p = subprocess.Popen("Rscript %s"%(options.name+'_function_prediction.R'), shell=True)
        sts = os.waitpid(p.pid,0)
    except:
        sys.exit()

def distrun(bedfile, reference, bglist, upgenes, downgenes, label1, label2, name, name2):
    
    #this function is for BETA Import, if you use this script as a seprate one, please ignore this function.
    peaklist = readfile(bedfile)
    (refdist, symdist) = get_the_nearest_dist(reference, peaklist)
    
    if name2 == False:
        dist = refdist
    else:
        dist = symdist
    #print dist       
    bgdist = []
    bggene = read_gn_group(bglist)
    #print bggene
    for gene in bggene:
        if gene in list(dist.keys()):
            bgdist.append(dist[gene])
        else:
            continue
    
    fns = [upgenes, downgenes]
    gene_sets = read_gene_sets(fns)    
    geneset_dists = assign_geneset_dist(gene_sets,dist)
    label = [label1,label2]
    rscript = write_R_code(bgdist, geneset_dists, label, name)
    rcode = open('%s_function_prediction.R'%name,'w')
    rcode.write(rscript)
    rcode.close()
    try:
        p = subprocess.Popen("Rscript %s"%(name+'_function_prediction.R'), shell=True)
        sts = os.waitpid(p.pid,0)
    except:
        sys.exit()
    
if __name__ == '__main__':
    main()


