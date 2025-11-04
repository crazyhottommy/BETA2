import sys
from beta.utils import *
def opt_validate_basic(optparser):
    """Validate options from a OptParser object.

Ret: Validated options object.
"""

    options = optparser.parse_args()
    if not options.peakfile or not options.exprefile or not options.kind:
        optparser.print_help()
        #error("peak file, expression file and expression tech type required fro BETA basic run, check your input file of -p -e and -k")
        sys.exit(1)
    if not options.genome and not options.reference:
        optparser.print_help()
        Info('ERROR: Genome (hg38, hg19, hg18, mm10, mm9) or the refseq geneinfo from UCSC is required.')
        sys.exit(1)
    
    if not os.path.isfile(options.peakfile):
        Info('ERROR: Cannot find peak file, a tab-delimit peak file must be given through -p (--peakfile).')
        sys.exit(1)
    if not os.path.isfile(options.exprefile):
        Info('ERROR: Cannot find differential expression file, a tab-delimit differential expression file must be given through -e (--diff_expr).')
        sys.exit(1)
    if options.reference and not os.path.isfile(options.reference):
        Info('ERROR: Cannot find the refseq gene info, a tab-dilimit 6 column (refseqID, chroms, strand, txstart, txend, genesymbol) file required via -r --reference')
        sys.exit(1)
    if options.boundaryfile and not os.path.isfile(options.boundaryfile):
        Info('ERROR: Cannot find the bounaryfile, a tab-dilimit bed file (at least 3 columns) is required')
        sys.exit(1)
    if options.genome != 'hg19' and options.genome != 'mm9' and options.genome != 'hg38' and options.genome != 'mm10' and options.genome != 'hg18' and  options.genome:
        Info('ERROR: Input a correct genome: hg38, hg19, hg18, mm10, mm9 (or provide reference files via -r parameter)')
        sys.exit(1)
    if not options.genome and options.boundarylimit == True and not options.boundaryfile:
        Info('ERROR: Do you really want to use CTCF boundary to filter peaks? If so, please inpit the boundary file via --bf, else, please don not use --bl')
        sys.exit(1)
    if not options.kind or (options.kind != 'LIM' and options.kind != 'CUF' and options.kind != "BSF" and options.kind != "O"):
        Info('ERROR: Please input the right kind of your expression file,  LIM for LIMMA standard format. CUF for CUFDIFF standard format, BSF for BETA specific format and O for other formats')
        sys.exit(1)
    if options.kind == "O" and not options.expreinfo:
        Info('ERROR: The differntial expression infomation is required when the it is other formats')
    if options.method != 'score' and options.method != 'distance':
        Info('ERROR: please choose the correct method to do the CR/TF function prediction, score for regulatory potential and distance for the distance to the proximal binding peak')
        sys.exit(1)
    if not options.name:
        options.name = "NA"
        Info("BETA will use 'NA' as all the ouput files prefix name")
    if not options.output:
        options.output = os.getcwd()
        options.output = os.path.join(options.output,'BETA_OUTPUT')
        run_cmd('mkdir %s'%options.output)
    else:
        if not os.path.isdir(options.output):
            run_cmd('mkdir %s'%options.output)
    if not options.expreinfo:
        if options.kind == 'LIM':
            options.expreinfo = '1,2,6'
        if options.kind == 'CUF':
            options.expreinfo = '2,10,13' 
        if options.kind == 'BSF':
            options.expreinfo = '1,2,3'
            
    #if options.diff_fdr > 1 or options.diff_fdr < 0:
    #    Info("--df options error, please set a number 0~1")
    #    sys.exit(1)
    if options.diff_amount < 0:
        Info("--da options error, please set a number 0~1 or bigger than 1")
        sys.exit(1)
    
    Info("Argument List: ")
    Info("Name = " + options.name)
    Info("Peak File = " + options.peakfile)
    Info("Top Peaks Number = %d"%options.peaknumber)
    Info("Distance = %d bp" %options.distance)
    Info("Genome = %s" %options.genome)
    Info("Expression File = %s"%options.exprefile)
    if options.kind == 'LIM':
        Info("Expression Type = MicroArray, LIMMA result")
    if options.kind == 'CUF':
        Info("Expression Type = RNA-Seq, Cuffdiff result")
    if options.kind == 'BSF':
        Info("BETA specific Expression Type")
    if options.kind == 'O':
        Info("Other formats")
    Info("Number of differential expressed genes = %s"%str(options.diff_amount))
    Info("Differential expressed gene FDR Threshold = %s"%str(options.diff_fdr))
    Info("Up/Down Prediction Cutoff = %f"%options.cutoff)
    if options.method == "score":
        Info("Function prediction based on regulatory potential")
    if options.method == "distance":
        Info("Function prediction based on the distance to the proximal binding peak")
    options.motifnumber = 10

    return options
def opt_validate_super(optparser):
    """Validate options from a OptParser object.

Ret: Validated options object.
"""
    options = optparser.parse_args()
    if not options.peakfile or not options.exprefile or not options.kind or not options.genomesequence:
        optparser.print_help()
        sys.exit(1)
    if (not options.genome or options.genome=='mm' or options.genome=='hg') and not options.reference:
        optparser.print_help()
        Info('ERROR: Genome(hg19, mm9,hg38,hg18,mm10) or the refseq geneinfo from UCSC is request.')
        sys.exit(1)
    
    if not os.path.isfile(options.peakfile):
        Info('ERROR: Cannot find peak file, a tab-delimit peak file must be given through -p (--peakfile).')
        sys.exit(1)
    if not os.path.isfile(options.exprefile):
        Info('ERROR: Cannot find differential expression file, a tab-delimit differential expression file must be given through -e (--diff_expr).')
        sys.exit(1)
    if not os.path.isfile(options.genomesequence):
        Info('ERROR: Cannot find genome sequece data, a fasta format genome sequence data is requested')
        sys.exit(1)
    if options.reference and not os.path.isfile(options.reference):
        Info('ERROR: Cannot find the refseq gene info, a tab-dilimit 6 column (refseqID, chroms, strand, txstart, txend, genesymbol) file required via -r --reference')
        sys.exit(1)
    if options.boundaryfile and not os.path.isfile(options.boundaryfile):
        Info('ERROR: Cannot find the bounaryfile, a tab-dilimit bed file (at least 3 columns) is required')
        sys.exit(1)
    if options.genome != 'hg19' and options.genome != 'mm9' and options.genome != 'hg38' and options.genome != 'mm10' and options.genome != 'hg18' and options.genome != 'hg' and options.genome != 'mm' and options.genome:
        Info('ERROR: Input a correct genome: hg38, hg19, hg18, mm10, mm9, hg, mm (or provide reference files via -r parameter)')
        sys.exit(1)
    if not options.genome and options.boundarylimit == True and not options.boundaryfile:
        Info('ERROR: Do you really want to use CTCF boundary to filter peaks? If so, please inpit the boundary file via --bf, else, please don not use --bl')
        sys.exit(1)
    if not options.kind or (options.kind != 'LIM' and options.kind != 'CUF' and options.kind != "BSF" and options.kind != "O"):
        Info('ERROR: Please input the right kind of your expression file,  LIM for LIMMA standard format. CUF for CUFDIFF standard format, BSF for BETA specific format and O for other formats')
        sys.exit(1)
    if options.method != 'score' and options.method != 'distance':
        Info('ERROR: please choose the correct method to do the CR/TF function prediction, score for regulatory potential and distance for the distance to the proximal binding peak')
        sys.exit(1)
    if not options.name:
        options.name = "NA"
    if not options.output:
        options.output = os.getcwd()
        options.output = os.path.join(options.output,'BETA_OUTPUT')
        run_cmd('mkdir %s'%options.output)
    else:
        if not os.path.isdir(options.output):
            run_cmd('mkdir %s'%options.output)
            
    if not options.expreinfo:
        if options.kind == 'LIM':
            options.expreinfo = '1,2,6'
        if options.kind == 'CUF':
            options.expreinfo = '2,10,13' 
        if options.kind == 'BSF':
            options.expreinfo = '1,2,3'
    #if options.diff_fdr > 1 or options.diff_fdr < 0:
    #    Info("--df options error, please set a number 0~1")
    #    sys.exit(1)
    if options.diff_amount < 0:
        Info("--da options error, please set a number 0~1 or bigger than 1")
        sys.exit(1)
    if options.motifnumber < 0:
        Info("--mn error, please set a number bigger than 0, 0~1 as the pvalue and >1 as the number to get the motif output in html file")

    Info("Argument List: ")
    Info("Name = " + options.name)
    Info("Peak File = " + options.peakfile)
    Info("Top Peaks Number = %d"%options.peaknumber)
    Info("Distance = %d bp" %options.distance)
    Info("Genome = %s" %options.genome)
    Info("Expression File = %s"%options.exprefile)
    Info("Genome Sequence fasta formated data = %s"%options.genomesequence)
    if options.kind == 'LIM':
        Info("Expression Type = MicroArray, LIMMA result")
    if options.kind == 'CUF':
        Info("Expression Type = RNA-Seq, Cuffdiff result")
    if options.kind == 'BSF':
        Info("BETA specific Expression Type")
    if options.kind == 'O':
        Info("Other formats")
    Info("Number of differential expressed genes = %s"%str(options.diff_amount))
    Info("Differential expressed gene FDR Threshold = %s"%str(options.diff_fdr))
    Info("Up/Down Prediction Cutoff = %f"%options.cutoff)
    if options.method == "score":
        Info("Function prediction based on regulatory potential")
    if options.method == "distance":
        Info("Function prediction based on the distance to the proximal binding peak")

    return options

def opt_validate_noexpre(optparser):
    """Validate options from a OptParser object.

Ret: Validated options object.
"""

    options = optparser.parse_args()
    if not options.peakfile:
        optparser.print_help()
        #error("peak file, expression file and expression tech type required fro BETA basic run, check your input file of -p -e and -k")
        sys.exit(1)
    if not options.genome and not options.reference:
        optparser.print_help()
        Info('ERROR: Genome(hg19, mm9, hg38, mm10, hg18) or the refseq geneinfo from UCSC is request.')
        sys.exit(1)
    
    if not os.path.isfile(options.peakfile):
        Info('ERROR: Cannot find peak file, a tab-delimit peak file must be given through -p (--peakfile).')
        sys.exit(1)
    
    if options.reference and not os.path.isfile(options.reference):
        Info('ERROR: Cannot find the refseq gene info, a tab-dilimit 6 column (refseqID, chroms, strand, txstart, txend, genesymbol) file required via -r --reference')
        sys.exit(1)
    if options.boundaryfile and not os.path.isfile(options.boundaryfile):
        Info('ERROR: Cannot find the bounaryfile, a tab-dilimit bed file (at least 3 columns) is required')
        sys.exit(1)
    if options.genome != 'hg19' and options.genome != 'mm9' and options.genome != 'hg38' and options.genome != 'mm10' and options.genome != 'hg18' and options.genome:
        Info('ERROR: Input a right genome, hg19, mm9, hg38, mm10, hg18 or input reference files via other parameter')
        sys.exit(1)
    if not options.genome and options.boundarylimit == True and not options.boundaryfile:
        Info('ERROR: Do you really want to use CTCF boundary to filter peaks? If so, please inpit the boundary file via --bf, else, please don not use --bl')
        sys.exit(1)
    if not options.name:
        options.name = "NA"
        Info("BETA will use 'NA' as all the ouput files prefix name")
    if not options.output:
        options.output = os.getcwd()
        options.output = os.path.join(options.output,'BETA_OUTPUT')
        run_cmd('mkdir %s'%options.output)
    else:
        if not os.path.isdir(options.output):
            run_cmd('mkdir %s'%options.output)

    Info("Argument List: ")
    Info("Name = " + options.name)
    Info("Peak File = " + options.peakfile)
    Info("Top Peaks Number = %d"%options.peaknumber)
    Info("Distance = %d bp" %options.distance)
    Info("Genome = %s" %options.genome)
    options.exprefile = ''
    options.expreinfo = ''
    return options
