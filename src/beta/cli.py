#!/usr/bin/env python3
"""
BETA Command Line Interface

BETA: Binding and Expression Target Analysis
Integrative analysis of ChIP-seq and RNA-seq/microarray data

Authors: Su Wang, Tommy Tang
Version: 2.0.0
"""

import sys
import argparse as ap
from beta import __version__


def prepare_argparser():
    """
    Prepare argument parser object with all subcommands

    Returns:
        Configured ArgumentParser instance
    """
    description = "BETA --- Binding and Expression Target Analysis"
    epilog = "For command line options of each command, type: %(prog)s COMMAND -h"

    argparser = ap.ArgumentParser(description=description, epilog=epilog)
    argparser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}"
    )

    subparsers = argparser.add_subparsers(
        help="sub-command help",
        dest='subcommand_name'
    )

    # Add subcommand parsers
    add_basic_parser(subparsers)
    add_plus_parser(subparsers)
    add_minus_parser(subparsers)

    return argparser


def add_basic_parser(subparsers):
    """Add 'basic' subcommand argument parser for target prediction"""
    basicparser = subparsers.add_parser(
        "basic",
        help="Main BETA Function: Transcription factor direct targets prediction",
        description="BETA-basic: Predict direct targets of TF and active/repressive function prediction.\n"
                    "EXAMPLE: beta basic -p peaks.bed -e gene_exp.diff -k LIM -g hg38 -n test -o output"
    )

    # Required arguments
    basicparser.add_argument(
        "-p", "--peakfile",
        dest="peakfile",
        type=str,
        required=True,
        help="BED format file of peak binding sites (3 or 5 columns: CHROM, START, END [NAME, SCORE])"
    )
    basicparser.add_argument(
        "-e", "--diff_expr",
        dest="exprefile",
        type=str,
        required=True,
        help="Differential expression file from LIMMA (microarray) or Cuffdiff (RNA-seq)"
    )
    basicparser.add_argument(
        "-k", "--kind",
        dest="kind",
        choices=("LIM", "CUF", "BSF", "O"),
        required=True,
        help="Expression file format: LIM (LIMMA), CUF (Cuffdiff), BSF (BETA specific), O (other, requires --info)"
    )

    # Genome arguments
    basicparser.add_argument(
        "-g", "--genome",
        dest="genome",
        choices=("hg38", "hg19", "hg18", "mm10", "mm9"),
        help="Genome assembly: hg38, hg19, hg18, mm10, mm9. For other assemblies, use -r option"
    )
    basicparser.add_argument(
        "-r", "--reference",
        dest="reference",
        type=str,
        help="RefSeq annotation file from UCSC (only if genome is not hg38/hg19/hg18/mm10/mm9)"
    )
    basicparser.add_argument(
        "--gname2",
        dest="gname2",
        action="store_true",
        default=False,
        help="Gene IDs in expression file are official gene symbols (default: False)"
    )

    # Expression file parsing
    basicparser.add_argument(
        "--info",
        dest="expreinfo",
        type=str,
        help="Column specification for expression data: 'geneID,logFC,FDR' (e.g., '1,2,6'). "
             "Default: 1,2,6 (LIMMA); 2,10,13 (Cuffdiff); 1,2,3 (BSF)"
    )

    # Output arguments
    basicparser.add_argument(
        "-o", "--output",
        dest="output",
        type=str,
        help="Output directory (default: current directory)"
    )
    basicparser.add_argument(
        "-n", "--name",
        dest="name",
        type=str,
        help="Prefix for output files (default: 'NA')"
    )

    # Analysis parameters
    basicparser.add_argument(
        "-d", "--distance",
        dest="distance",
        type=int,
        default=100000,
        help="Maximum distance from TSS to consider peaks (bp, default: 100000)"
    )
    basicparser.add_argument(
        "--pn",
        dest="peaknumber",
        type=int,
        default=10000,
        help="Maximum number of peaks to consider (default: 10000)"
    )
    basicparser.add_argument(
        "--method",
        dest="method",
        choices=("score", "distance"),
        default="score",
        help="Method for TF function prediction: 'score' (regulatory potential) or 'distance' (default: score)"
    )

    # CTCF boundary filtering
    basicparser.add_argument(
        "--bl",
        dest="boundarylimit",
        action="store_true",
        default=False,
        help="Use CTCF boundaries to filter peaks (default: False)"
    )
    basicparser.add_argument(
        "--bf",
        dest="boundaryfile",
        type=str,
        help="CTCF peaks BED file (only if --bl is set and genome is not hg19/mm9)"
    )

    # Differential expression filters
    basicparser.add_argument(
        "--df",
        dest="diff_fdr",
        type=float,
        default=1.0,
        help="FDR threshold for differential expression (0-1, default: 1.0)"
    )
    basicparser.add_argument(
        "--da",
        dest="diff_amount",
        type=float,
        default=0.5,
        help="Fraction (0-1) or number (>1) of top DE genes to consider (default: 0.5)"
    )
    basicparser.add_argument(
        "-c", "--cutoff",
        dest="cutoff",
        type=float,
        default=0.001,
        help="P-value cutoff for target gene selection (KS test, default: 0.001)"
    )


def add_plus_parser(subparsers):
    """Add 'plus' subcommand parser for target prediction with motif analysis"""
    plusparser = subparsers.add_parser(
        'plus',
        help="Target prediction with motif analysis",
        description="BETA-plus: Predict direct targets of TF, function prediction, and motif analysis.\n"
                    "EXAMPLE: beta plus -p peaks.bed -e gene_exp.diff -k LIM -g hg38 --gs hg38.fa -n test -o output"
    )

    # Required arguments
    plusparser.add_argument(
        "-p", "--peakfile",
        dest="peakfile",
        type=str,
        required=True,
        help="BED format file of peak binding sites (3 or 5 columns: CHROM, START, END [NAME, SCORE])"
    )
    plusparser.add_argument(
        "-e", "--diff_expr",
        dest="exprefile",
        type=str,
        required=True,
        help="Differential expression file from LIMMA (microarray) or Cuffdiff (RNA-seq)"
    )
    plusparser.add_argument(
        "-k", "--kind",
        dest="kind",
        choices=("LIM", "CUF", "BSF", "O"),
        required=True,
        help="Expression file format: LIM (LIMMA), CUF (Cuffdiff), BSF (BETA specific), O (other, requires --info)"
    )
    plusparser.add_argument(
        "--gs",
        dest="genomesequence",
        type=str,
        required=True,
        help="Genome sequence file in FASTA format (required for motif analysis)"
    )

    # Genome arguments
    plusparser.add_argument(
        "-g", "--genome",
        dest="genome",
        choices=("hg38", "hg19", "hg18", "mm10", "mm9", "hg", "mm"),
        help="Genome assembly: hg38, hg19, hg18, mm10, mm9, hg (human), mm (mouse)"
    )
    plusparser.add_argument(
        "-r", "--reference",
        dest="reference",
        type=str,
        help="RefSeq annotation file from UCSC (only if genome is not hg38/hg19/hg18/mm10/mm9)"
    )
    plusparser.add_argument(
        "--gname2",
        dest="gname2",
        action="store_true",
        default=False,
        help="Gene IDs in expression file are official gene symbols (default: False)"
    )

    # Expression file parsing
    plusparser.add_argument(
        "--info",
        dest="expreinfo",
        type=str,
        help="Column specification for expression data: 'geneID,logFC,FDR' (e.g., '1,2,6'). "
             "Default: 1,2,6 (LIMMA); 2,10,13 (Cuffdiff); 1,2,3 (BSF)"
    )

    # Output arguments
    plusparser.add_argument(
        "-o", "--output",
        dest="output",
        type=str,
        help="Output directory (default: current directory)"
    )
    plusparser.add_argument(
        "-n", "--name",
        dest="name",
        type=str,
        help="Prefix for output files (default: 'NA')"
    )

    # Analysis parameters
    plusparser.add_argument(
        "-d", "--distance",
        dest="distance",
        type=int,
        default=100000,
        help="Maximum distance from TSS to consider peaks (bp, default: 100000)"
    )
    plusparser.add_argument(
        "--pn",
        dest="peaknumber",
        type=int,
        default=10000,
        help="Maximum number of peaks to consider (default: 10000)"
    )
    plusparser.add_argument(
        "--method",
        dest="method",
        choices=("score", "distance"),
        default="score",
        help="Method for TF function prediction: 'score' (regulatory potential) or 'distance' (default: score)"
    )

    # CTCF boundary filtering
    plusparser.add_argument(
        "--bl",
        dest="boundarylimit",
        action="store_true",
        default=False,
        help="Use CTCF boundaries to filter peaks (default: False)"
    )
    plusparser.add_argument(
        "--bf",
        dest="boundaryfile",
        type=str,
        help="CTCF peaks BED file (only if --bl is set and genome is not hg19/mm9)"
    )

    # Differential expression filters
    plusparser.add_argument(
        "--df",
        dest="diff_fdr",
        type=float,
        default=1.0,
        help="FDR threshold for differential expression (0-1, default: 1.0)"
    )
    plusparser.add_argument(
        "--da",
        dest="diff_amount",
        type=float,
        default=0.5,
        help="Fraction (0-1) or number (>1) of top DE genes to consider (default: 0.5)"
    )
    plusparser.add_argument(
        "-c", "--cutoff",
        dest="cutoff",
        type=float,
        default=0.001,
        help="P-value cutoff for target gene selection (KS test, default: 0.001)"
    )

    # Motif analysis parameters
    plusparser.add_argument(
        "--mn",
        dest="motifnumber",
        type=float,
        default=10,
        help="P-value cutoff (0-1) or number of motifs (>1) to report (default: 10)"
    )


def add_minus_parser(subparsers):
    """Add 'minus' subcommand parser for binding-only analysis"""
    minusparser = subparsers.add_parser(
        "minus",
        help="Target prediction with binding data only (no expression data)",
        description="BETA-minus: Calculate regulatory potential scores from binding data only.\n"
                    "EXAMPLE: beta minus -p peaks.bed -g hg38 -n test -o output"
    )

    # Required arguments
    minusparser.add_argument(
        "-p", "--peakfile",
        dest="peakfile",
        type=str,
        required=True,
        help="BED format file of peak binding sites (3 or 5 columns: CHROM, START, END [NAME, SCORE])"
    )

    # Genome arguments
    minusparser.add_argument(
        "-g", "--genome",
        dest="genome",
        choices=("hg38", "hg19", "hg18", "mm10", "mm9"),
        help="Genome assembly: hg38, hg19, hg18, mm10, mm9"
    )
    minusparser.add_argument(
        "-r", "--reference",
        dest="reference",
        type=str,
        help="RefSeq annotation file from UCSC (only if genome is not hg38/hg19/hg18/mm10/mm9)"
    )

    # Output arguments
    minusparser.add_argument(
        "-o", "--output",
        dest="output",
        type=str,
        help="Output directory (default: current directory)"
    )
    minusparser.add_argument(
        "-n", "--name",
        dest="name",
        type=str,
        help="Prefix for output files (default: 'NA')"
    )

    # Analysis parameters
    minusparser.add_argument(
        "-d", "--distance",
        dest="distance",
        type=int,
        default=100000,
        help="Maximum distance from TSS to consider peaks (bp, default: 100000)"
    )
    minusparser.add_argument(
        "--pn",
        dest="peaknumber",
        type=int,
        default=10000,
        help="Maximum number of peaks to consider (default: 10000)"
    )

    # CTCF boundary filtering
    minusparser.add_argument(
        "--bl",
        dest="boundarylimit",
        action="store_true",
        default=False,
        help="Use CTCF boundaries to filter peaks (default: False)"
    )
    minusparser.add_argument(
        "--bf",
        dest="boundaryfile",
        type=str,
        help="CTCF peaks BED file (only if --bl is set and genome is not hg19/mm9)"
    )


def main():
    """Main entry point for BETA CLI"""
    # Parse command-line arguments
    argparser = prepare_argparser()
    args = argparser.parse_args()

    # Get subcommand
    subcommand = args.subcommand_name

    if not subcommand:
        argparser.print_help()
        sys.exit(1)

    # Import and run appropriate subcommand
    try:
        if subcommand == "basic":
            from beta.core.runbeta import basicrun
            basicrun(argparser)
        elif subcommand == "plus":
            from beta.core.runbeta import plusrun
            plusrun(argparser)
        elif subcommand == "minus":
            from beta.core.runbeta import minusrun
            minusrun(argparser)
    except KeyboardInterrupt:
        sys.stderr.write("\nUser interrupted. Exiting...\n")
        sys.exit(0)
    except Exception as e:
        sys.stderr.write(f"\nError: {e}\n")
        sys.exit(1)


if __name__ == '__main__':
    main()
