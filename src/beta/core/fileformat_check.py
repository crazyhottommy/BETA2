from .corelib import *
import re

class checkfileformat:
    def __init__(self, options):
        self.peakfile = options.peakfile
        self.exprefile = options.exprefile
        if options.reference:
            self.reference = options.reference
        if options.boundaryfile:
            self.boundary = options.boundaryfile
        self.genomesequence = options.genomesequence
        self.expreinfo = options.expreinfo
        
        self.gene = options.gname2
    
    def bed_format(self,check_3col_bed=True):
        fname = self.peakfile #the peakfile format should be at least 5 column, if it is 3 column, we need to do some conversion
        
        with open(fname) as to_check:
            lines = to_check.readlines()
            first_line,last_line = lines[0],lines[-1]
            is_5col_bed = lambda a_line:re.search("chr\S+\s\d+\s\d+\s\S+\s\d+[.]*[\d]*",a_line)
            is_3col_bed = lambda a_line:re.search("chr\S+\s\d+\s\d+$",a_line)
            is_3tab_delimite = lambda a_line:len(a_line.split('\t'))==3
            is_5tab_delimite = lambda a_line:len(a_line.split('\t'))==5
            suffix_5col = lambda a_file: a_file+".5c.bed"
            def threetofive(fname):
                with open(suffix_5col(fname),'w') as three2five_coled:
                    for line_num,a_line in enumerate(open(fname)):
                        a_line=a_line.strip()+"\tunknown_peak_%d\t1.0\n"%line_num
                        three2five_coled.write(a_line)
                
            def check_a_line(a_line):
                if is_5col_bed(last_line):
                    if not is_5tab_delimite :
                        error("Your input peak file is not tab delimited")
                    return True
                else:
                    if check_3col_bed:
                        if not is_3tab_delimite :
                            error("Your input peak file is not tab delimited")
                        if  is_3col_bed(last_line):
                            warn("You input a 3 column bed file like this:\t\t%s"%last_line[:50])
                            Info("[3 Column to 5 Column] %s ==> %s "%(suffix_5col(fname),fname))                                                
                            threetofive(fname)
                            Info("Use %s instead of %s as the input BED and run the executive Again"%(suffix_5col(fname),fname))
                        else:
                            error("The input bed file %s has a wrong format!(3 column checking active)"%fname)
                            Info("Wrong Format:\t\t\t%s"%last_line[:50])
                            Info("Right Format should look like:\t%s"%('chr1\t567577\t567578\tMACS_peak_1\t119.00'))
                            Info("Or the depreciate 3-column format like this:\t%s"%('chr1\t567577\t567578'))
                        return False
                    else:
                        error("The input bed file %s has a wrong format!"%fname)
                        Info("Wrong Format:\t\t\t%s"%last_line[:50])
                        Info("Right Format should look like:\t%s"%('chr1\t567577\t567578\tMACS_peak_1\t119.00'))
                        
            if check_a_line(first_line) and check_a_line(last_line):
                Info("Check %s successfully!"%fname)
                return fname
            else:
                return suffix_5col(fname)
            
    def check_fasta_dna(self):
        fname = self.genomesequence
        """
        Check if a file has the format of fasta

        @return:  whether the file passed the fasta check
        """                           
        with open(fname) as fasta_f:
            first_line = fasta_f.readline()
            if not first_line[0] == ">":
                error("The input fasta file %s has a wrong format!"%fname)
                Info("Wrong Format:\t\t\t%s"%first_line[:50])
                Info("Right Format should look like:\t%s"%('>chr1:1150372-1150572'))
                return False
            
            second_line = fasta_f.readline()
            fasta_pattern_scd = "[AGCTN]+"        
            if not re.search(fasta_pattern_scd,second_line):
                error("The input fasta file %s has a wrong format!"%fname)
                Info("Wrong Format:\t\t\t%s"%second_line[:50])
                Info("Right Format should look like:\t%s"%('NGGGCCATTCA'))
                return False
        return True

    def check_expr(self):
        "Check the expression file format"
        is_refseq = lambda string:re.search('[A-Z][A-Z]_\d+(_at)?',string)
        is_genesymbol = lambda string:re.search('\w+',string)
        is_value = lambda data:re.search('\"*(\-*)(\d+)\.*(\d*)\"*',data)
        expreinfo = self.expreinfo.split(',')
        exprefile = self.exprefile
        expr_info = {}
        expr_info['refseq']=int(expreinfo[0])-1
        expr_info['tscore']=int(expreinfo[1])-1
        expr_info['fdr']=int(expreinfo[2])-1
        with open(exprefile, 'r') as expref:
            first_line = expref.readline()
            if not first_line.startswith('#'):
                Info("%s is not the header of the expression file"%first_line)
                second_line = first_line
                Info("Checking the differential expression infomation...")
                Info("Take the first line with Differential Information as an example: %s"%second_line)                
                second_line = second_line.split('\t')
                refseq = second_line[int(expreinfo[0])-1]
                tscore = second_line[int(expreinfo[1])-1]
                fdr = second_line[int(expreinfo[2])-1]
                if self.gene == False:
                    if not is_refseq(refseq) or not is_value(tscore) or not is_value(fdr):
                        Info("BETA cannot recognize the refseq gene ID, status value(logFC) or FDR. Please give the exact column numbers of the refseq, logFC, and fdr like: 1,2,6 for LIMMA; 2,10,13 for Cufdiff; and 1,2,3 for BETA specific format.")
                        sys.exit()
                    else:
                        Info("Differential Expression file format successful passed")
                if self.gene == True:
                    if not is_genesymbol(refseq) or not is_value(tscore) or not is_value(fdr):
                        Info("BETA cannot recognize the official gene symbol, status value(logFC) or FDR. Please give the exact column numbers of the genesymbol, logFC, and fdr like: 1,2,6 for LIMMA; 2,10,13 for Cufdiff; and 1,2,3 for BETA specific format.")
                        sys.exit()
                    else:
                        Info("Differential Expression file format successful passed")
            if first_line.startswith('#'):
                Info("%s is the header of the expression file"%first_line)
                second_line = expref.readline()
                Info("Checking the differential expression infomation...")
                Info("Take the first line with Differential Information as an example: %s"%second_line)
                second_line = second_line.split('\t')
                refseq = second_line[int(expreinfo[0])-1]
                tscore = second_line[int(expreinfo[1])-1]
                fdr = second_line[int(expreinfo[2])-1]
                if self.gene == False:
                    if not is_refseq(refseq) or not is_value(tscore) or not is_value(fdr):
                        Info("BETA cannot recognize the refseq gene ID, status value(logFC) or FDR. Please give the exact column numbers of the refseq, logFC, and fdr like: 1,2,6 for LIMMA; 2,10,13 for Cufdiff; and 1,2,3 for BETA specific format.")
                        sys.exit()
                    else:
                        Info("Differential Expression file format successful passed")
                if self.gene == True:
                    if not is_genesymbol(refseq) or not is_value(tscore) or not is_value(fdr):
                        Info("BETA cannot recognize the official gene symbol, status value(logFC) or FDR. Please give the exact column numbers of the genesymbol, logFC, and fdr like: 1,2,6 for LIMMA; 2,10,13 for Cufdiff; and 1,2,3 for BETA specific format.")
                        sys.exit()
                    else:
                        Info("Differential Expression file format successful passed")
                        
        return expr_info
    
