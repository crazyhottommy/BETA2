# from numpy import random as r
# import numpy
import time
from .corelib import *


class permutation:

    # Get the FDR via multi permutations
    def __init__(self, options, selected, counts, geneInfo, bgtotal):

        self.name = options.name
        self.gname2 = options.gname2
        self.outdir = options.output
        self.counts = counts
        self.selected = selected
        self.geneInfo = geneInfo
        self.bgtotal = bgtotal
        print(self.selected)

    def fdr(self):

        inf_uprank = "%s_upregulate_target.txt" % self.name
        inf_downrank = "%s_downregulate_target.txt" % self.name
        upoutput = "%s_uptarget.txt" % self.name
        downoutput = "%s_downtarget.txt" % self.name
        selected = self.selected
        counts = self.counts
        self.output = []

        for lists in selected:

            if lists == '"upregulate"' and counts[0] != 0:
                infile = inf_uprank
                inf = open(infile)
                outf = open(upoutput, "w")
                self.output.append(upoutput)

            if lists == '"downregulate"' and counts[1] != 0:
                infile = inf_downrank
                inf = open(infile)
                outf = open(downoutput, "w")
                self.output.append(downoutput)

            outf.write(
                "Chroms\ttxStart\ttxEnd\trefseqID\trank_product\tStrands\tGeneSymbol\treg_potential\tbinding_rank\texpr_rank\tlog2FC\tpadj\n"
            )

            rank = (
                {}
            )  # {refseqID:[symbol,RP, rank, reg_score, binding_rank, expr_rank, logfc, padj]}

            GeneID = []
            if self.gname2 == False:
                ID_col = 0
            else:
                ID_col = 1

            for line in inf:
                if not line.strip():
                    continue
                else:
                    line = line.strip()

                    line = line.split("\t")
                    # Columns: refseq, symbol, RP, rank, reg_score, binding_rank, expr_rank, logfc, padj
                    rank[line[ID_col]] = [
                        line[0],  # refseq
                        line[1],  # symbol
                        line[2],  # RP
                        line[3],  # rank
                        line[4],  # reg_score
                        line[5],  # binding_rank
                        line[6],  # expr_rank
                        line[7],  # logfc
                        line[8],  # padj
                    ]
                    GeneID.append(line[ID_col])

            for gene in GeneID:
                ob_rp = float(rank[gene][2])
                Infos = self.geneInfo[gene]
                # [chr, tss, tts, refseq, score, strand, symbol, rank]
                chrom = Infos[0]
                tss = Infos[1]
                tts = Infos[2]
                refseq = Infos[3]
                strand = Infos[5]
                symbol = Infos[6]
                reg_score = rank[gene][4]
                binding_rank = rank[gene][5]
                expr_rank = rank[gene][6]
                logfc = rank[gene][7]
                padj = rank[gene][8]
                outf.write(
                    "%s\t%s\t%s\t%s\t%.3e\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                    % (
                        chrom,
                        tss,
                        tts,
                        refseq,
                        float(ob_rp),
                        strand,
                        symbol,
                        reg_score,
                        binding_rank,
                        expr_rank,
                        logfc,
                        padj,
                    )
                )

            outf.close()
            inf.close()
            run_cmd("rm %s" % infile)

    def prepare_ouput_basic(self):
        # this will be the last step for basic BETA, move all the result into output directory
        for f in self.output:
            run_cmd("mv %s %s" % (f, self.outdir))

    def prepare_ouput_plus(self):
        # this will be the last step for basic BETA, move all the result into output directory
        for f in self.output:
            run_cmd("cp %s %s" % (f, self.outdir))
