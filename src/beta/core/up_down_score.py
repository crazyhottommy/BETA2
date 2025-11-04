import os,sys
import math
import re, subprocess
# from optparse import OptionParser  # Deprecated, use argparse

CHROM_CONVERT = {'chrI':'chr1','chrII':'chr2','chrIII':'chr3','chrIV':'chr4','chrV':'chr5','chrVI':'chr6',
                 'chrVII':'chr7','chrVIII':'chr8','chrIX':'chr9','chrX':'chr10','chrXI':'chr11','chrXII':'chr12',
                 'chrXIII':'chr13','chrXIV':'chr14','chrXV':'chr15','chrXVI':'chr16','chrXVII':'chr17','chrXVIII':'chr18',
                 'chrXIX':'chr19','chrXX':'chr20','chrXXI':'chr21','chrXXII':'chr22','chrX':'chrX','chrY':'chrY','chrM':'chrM'}
chroms = list(CHROM_CONVERT.values())

def prepare_optparser ():
    
    """Prepare optparser object. New options will be added in this
    function first.
    
    """
    
    usage = "usage: %prog <-g gdb -b bed -g geneset> [options]"
    description = "ChGS -- Gene Set Association Analysis. Note that if the number of genomic coordinates is below 1000 or the number of a gene set is below 500, the p value of ChGS may not represent biological implications well enough."
    
    optparser = OptionParser(version="%prog 0.5",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-s","--score",dest="score",type="string",\
                         help="score file witch caculated with RegPotential python script")
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
    #print gene_sets
    return gene_sets

def read_score_file( scorefile ):

    scoref = open(scorefile,'r')
    refscores = {}
    symbolscores = {}#a refseq ID or gene symbols may have several locations, hence they have several scores
    
    for line in scoref:
        if line.startswith('#'):
            continue
        else:
            line = line.strip()
            line = line.split('\t')
            try:
                refscores[line[3]].append(float(line[4]))
                symbolscores[line[6]].append(float(line[4]))
            except:
                refscores[line[3]] = [float(line[4])]
                symbolscores[line[6]] = [float(line[4])]
    refscore = {}
    symbolscore = {}#take the largerst score of each gene
    
    for gene in list(refscores.keys()):
        refscore[gene] = max(refscores[gene])
    for gene in list(symbolscores.keys()):
        symbolscore[gene] = max(symbolscores[gene])
        
    ref_score = sorted(list(refscore.items()), key=lambda refscore:refscore[1], reverse=True)
    symbol_score = sorted(list(symbolscores.items()), key=lambda symbolscores:symbolscores[1], reverse=True)
    #[('ad',2),('ab',2),('cd',3)]

    rank = 1
    same = 0
    minscore = ref_score[0][1]
    refrank = {}
    for g in ref_score:
        gene = g[0]
        score = g[1]
        if score == minscore:
            rank = rank
            same += 1
        else:
            rank += same
            same = 1
        minscore = score
        refrank[gene] = rank
        
    rank = 1
    same = 0
    minscore = symbol_score[0][1]
    symbolrank = {}
    for g in symbol_score:
        gene = g[0]
        score = g[1]
        if score == minscore:
            rank = rank
            same += 1
        else:
            rank += same
            same = 1
        minscore = score
        symbolrank[gene] = rank
        
    return (refrank,symbolrank)

def assign_score_to_gene_sets(gene_sets, score):
    
    genesets_score = []
    for geneset in gene_sets:
        geneset_score = []
        for gene in geneset:
            if gene in list(score.keys()):
                genescore = score[gene]
                geneset_score.append(genescore)
            else:
                pass
        genesets_score.append(geneset_score)
    #print genesets_score
    return genesets_score

def write_R_code(bgscore, genesets_score, label, name):
    
    rscript = ''
    rscript += '\n'
    rscript += 'options(warn = -1)\n'
    rscript += 'cr <- colorRampPalette(col=c("#C8524D", "#BDD791", "#447CBE", "#775A9C"), bias=1)\n'
    rscript += 'col <- cr(2)\n'
    rscript += 'scores <- list()\n'

    tempmaxs = []
    tempmins = []
    for l in range(len(label)):
        
        tempmax = max(genesets_score[l])
        nonzero = [t for t in genesets_score[l] if float(t) > 0]
        tempmin = min(nonzero)
        tempmaxs.append(tempmax)
        tempmins.append(tempmin)

    maxscore = max(max(bgscore), max(tempmaxs))
    bgnonzero = [t for t in bgscore if float(t) > 0]
    minscore = max(min(min(bgnonzero),min(tempmins)),0.03)
    #minscore = min(min(bgnonzero),min(tempmins))
    
    rscript += "\n"
    rscript += "vmax <- %f\n"%maxscore

    bg = 'd <- as.numeric(c(%f'%bgscore[0]
    n = len(bgscore)
    for i in range(n-2):
        bg += ',%f'%bgscore[i + 1]
    bg += ',%f))\n'%bgscore[-1]
    rscript += bg
    #rscript += "bg <- vmax - d\n"
    rscript += "bg <- d\n"
    rscript += "bg[bg == 'NA'] <- NA\n"
    
    for l in range(len(label)):
        group = 'd <- as.numeric(c(%s'%genesets_score[l][0]
        for j in range(len(genesets_score[l]) - 2):
            group += ',%s'%genesets_score[l][j + 1]
        group += ',%s))\n'%genesets_score[l][-1]
        rscript += group
        #rscript += 'scores[["%s"]] <- vmax - d\n'%label[l]
        rscript += 'scores[["%s"]] <- d\n'%label[l]
        
    rscript += "uplimit <- trunc(vmax) + 1\n"
    rscript += "\n"
    rscript += "breaks <- seq(0, uplimit, by = 1)\n"
    rscript += "hs <- list()\n"
    
    for i in range(len(label)):
        rscript += 'd <- as.vector(na.omit(scores[["%s"]]))\n'%label[i]
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
        rscript += 'd <- as.vector(na.omit(scores[["%s"]]))\n'%i
        rscript += 'd <- d[d >= 0 & d <= uplimit]\n'
        rscript += 'ks <- ks.test(d, dbg, alternative = "greater")\n'
        rscript += 'p_values[["%s"]] <- ks$p.value\n'%i

    rscript += "\n"
    
    rscript += "ps = as.matrix(p_values)\n"
    rscript += "write.table(ps, file = '%s_pvalues.txt', sep = '\t', col.names = F)\n"%name

    rscript += 'pdf("%s_function_prediction.pdf", height = 7, width = 7)\n'%name
    rscript += 'downlimit <- %f\n'%minscore
    rscript += 'plot(x, cbg, type="l", lty=2, lwd=2, col="black", xlab="Rank of genes based on Regulatory Potential Score (from high to low)", ylab="Cumulative Fraction of Genes%%", main="%s Activating/Repressive Function Prediction", xlim=c(downlimit, uplimit), ylim=c(0, 100))\n'%name
    rscript += "\n"
    rscript += "l <- c('static (background)')\n"
    for i in range(len(label)):
        rscript += 'lines(x, chs[["%s"]], lty=1, lwd=2, col=col[%d])\n'%(label[i],i+1)
        rscript += 'l <- c(l, paste("%s (", format(p_values[["%s"]], digits=3), ")", sep=""))\n'%(label[i],label[i])
    rscript += "legend('topleft', legend=l, col=c('black',col), pch=15)\n"
    #rscript += "axis(side=1,at=seq(0,uplimit),labels=seq(uplimit,0))\n"
    rscript += "dev.off()"
        
    return rscript

def main():
    
    optparser = prepare_optparser ()
    
    (options,args) = optparser.parse_args()
    if not options.score or not options.gset:
        
        optparser.print_help()
        sys.exit(1)
        
    #genome = options.gdb
    scorefile = options.score
    if type(options.gset) == str:
        options.gset = [options.gset]
    if options.label:
        if len(options.label) != len(options.gset):
            error("The number of the gene set labels (-l or --lab) must be the same as that of the gene sets (-g or --gset).")
            sys.exit(1)
    else:
        options.label =  ["gene set" + str(x) for x in range(1, len(options.gset)+1)]
 
    (refscore, symbolscore) = read_score_file(scorefile)
    
    if options.name2 == False:
        score = refscore
    else:
        score = symbolscore
          
    bgscore = []
    if options.bg == 'ALL':
        bgscore = list(score.values())
    else:
        bggene = read_gn_group(options.bg)
        #print bggene
        for gene in bggene:
            if gene in list(score.keys()):
                bgscore.append(score[gene])
            else:
                continue
    
    fns = options.gset
    gene_sets = read_gene_sets(fns)
    genesets_score = assign_score_to_gene_sets(gene_sets, score)
    rscript = write_R_code(bgscore, genesets_score, options.label, options.name)
    rcode = open('%s_function_prediction.R'%options.name,'w')
    rcode.write(rscript)
    rcode.close()
    try:
        p = subprocess.Popen("Rscript %s"%(options.name+'_function_prediction.R'), shell=True)
        sts = os.waitpid(p.pid,0)
    except:
        sys.exit()

def scorerun(scorefile, bglist, uplist, downlist, label1, label2, name, name2):

    gset = [uplist, downlist]
    label = [label1, label2]
 
    (refscore, symbolscore) = read_score_file(scorefile)
    
    if name2 == False:
        score = refscore
    else:
        score = symbolscore
          
    bgscore = []
    bggene = read_gn_group(bglist)
 
    for gene in bggene:
        if gene in list(score.keys()):
            bgscore.append(score[gene])
        else:
            continue
    
    fns = gset
    gene_sets = read_gene_sets(fns)
    genesets_score = assign_score_to_gene_sets(gene_sets, score)
    rscript = write_R_code(bgscore, genesets_score, label, name)
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




