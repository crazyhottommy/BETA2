# BETA2: Binding and Expression Target Analysis

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-Artistic%202.0-green)](LICENSE)

BETA is a computational tool for integrative analysis of ChIP-seq and RNA-seq/microarray data to predict transcription factor (TF) direct target genes and identify whether the TF primarily functions as a transcriptional activator or repressor.

## Overview

BETA integrates ChIP-seq binding data with differential gene expression data to:

1. **Predict direct target genes** by combining:
   - Regulatory potential scores based on distance to transcription start sites (TSS)
   - Binding peak signals
   - Gene expression changes

2. **Infer TF function** as activator or repressor through:
   - Regulatory potential analysis
   - Gene set enrichment analysis

3. **Identify enriched motifs** in target regions (optional)

## Key Features

- **Integrative Analysis**: Combines ChIP-seq peaks with gene expression data
- **Regulatory Potential Scoring**: Distance-weighted scoring system
- **Statistical Assessment**: Kolmogorov-Smirnov test and permutation-based FDR
- **Motif Analysis**: Optional motif scanning and enrichment analysis
- **Multiple Input Formats**: Supports LIMMA, Cuffdiff, and custom formats
- **Genome Support**: Human (hg38, hg19, hg18) and Mouse (mm10, mm9)

## Installation

### Requirements

- Python 3.8 or higher
- C compiler (gcc) for motif scanning module

### From Source

```bash
git clone https://github.com/yourusername/BETA2.git
cd BETA2
pip install -e .
```

### Using pip (when available)

```bash
pip install beta-binding-analysis
```

## Quick Start

### Basic Analysis

Predict TF target genes and function (activator/repressor):

```bash
beta basic \
  -p peaks.bed \
  -e diff_expr.txt \
  -k LIM \
  -g hg38 \
  -n my_experiment \
  -o output_dir
```

### Plus Mode (with Motif Analysis)

Include motif analysis:

```bash
beta plus \
  -p peaks.bed \
  -e diff_expr.txt \
  -k LIM \
  -g hg38 \
  --gs hg38.fa \
  -n my_experiment \
  -o output_dir
```

### Minus Mode (Peaks Only)

Analyze binding data without expression data:

```bash
beta minus \
  -p peaks.bed \
  -g hg38 \
  -n my_experiment \
  -o output_dir
```

## Input Files

### ChIP-seq Peaks (required)

BED format file (3 or 5 columns):
```
chr1    1000    2000
chr1    5000    6000    peak1    100
```

### Differential Expression (required for basic/plus modes)

Supported formats:

1. **LIMMA** (`-k LIM`): Standard LIMMA output
2. **Cuffdiff** (`-k CUF`): Cuffdiff gene_exp.diff format
3. **BETA Standard Format** (`-k BSF`):
   ```
   GeneSymbol    log2FoldChange    FDR
   TP53          2.5               0.001
   MYC           -1.8              0.01
   ```
4. **Other** (`-k O`): Specify columns with `--info`

### Genome Sequence (for plus mode)

FASTA format genome sequence file (required for motif analysis)

## Output Files

### Basic Mode

- `{name}_targets.txt`: Predicted target genes with statistics
- `{name}_uptarget.txt`: Up-regulated targets
- `{name}_downtarget.txt`: Down-regulated targets
- `{name}_function.pdf`: TF function prediction plot

### Plus Mode

Additional files:
- `{name}_motif.html`: Motif enrichment results
- `{name}_motif_logo/`: Motif logos
- Motif scanning results

## Algorithm

### Regulatory Potential Score

For each gene, BETA calculates a regulatory potential score based on nearby binding peaks:

```
Score = Σ exp(-0.5 - 4 × distance/max_distance)
```

Where distance is from peak center to TSS.

### Target Prediction

1. Rank genes by regulatory potential
2. Rank genes by differential expression
3. Combine rankings using Kolmogorov-Smirnov test
4. Calculate FDR through permutation testing

### Function Prediction

Assess enrichment of up-regulated vs down-regulated genes among predicted targets using one-sided KS test.

## Command-line Options

### Common Options

| Option | Description | Default |
|--------|-------------|---------|
| `-p, --peakfile` | ChIP-seq peaks (BED format) | Required |
| `-g, --genome` | Genome assembly (hg38/hg19/hg18/mm10/mm9) | Required |
| `-n, --name` | Output prefix | "NA" |
| `-o, --output` | Output directory | Current directory |
| `-d, --distance` | Distance from TSS (bp) | 100000 |
| `--pn` | Number of peaks to consider | 10000 |

### Expression Options

| Option | Description |
|--------|-------------|
| `-e, --diff_expr` | Differential expression file |
| `-k, --kind` | Expression file format (LIM/CUF/BSF/O) |
| `--info` | Column specification (for -k O) |
| `--df` | FDR threshold |
| `--da` | Top genes to consider (fraction or number) |

### Advanced Options

| Option | Description |
|--------|-------------|
| `--method` | Scoring method (score/distance) |
| `-c, --cutoff` | P-value cutoff for targets |
| `--bl` | Use CTCF boundary filtering |
| `--gname2` | Gene IDs are gene symbols |

## Examples

### Example 1: Basic TF Analysis (hg38)

```bash
beta basic \
  -p ERalpha_peaks.bed \
  -e ERalpha_treatment_vs_control.txt \
  -k LIM \
  -g hg38 \
  -n ERalpha \
  -d 100000 \
  -c 0.001
```

### Example 2: With Custom Expression Format

```bash
beta basic \
  -p TF_peaks.bed \
  -e expression.txt \
  -k O \
  --info 1,3,7 \
  -g hg38 \
  -n TF_experiment
```
*(Column 1: gene ID, Column 3: log2FC, Column 7: FDR)*

### Example 3: Mouse Analysis with Motif Scanning

```bash
beta plus \
  -p mm10_peaks.bed \
  -e mm10_expression.txt \
  -k CUF \
  -g mm10 \
  --gs mm10.fa \
  -n mouse_TF \
  --mn 20
```

## Migration from BETA 1.x

This is a modernized Python 3 version of BETA. Key changes:

- **Python 3.8+** required (was Python 2.6/2.7)
- **Default genome**: hg38 (was hg19)
- **Improved performance**: Optimized algorithms
- **Better logging**: Structured logging instead of print statements
- **Modern packaging**: Uses pyproject.toml and pip installable

### Command Compatibility

All BETA 1.x commands should work with BETA 2.0 without changes. However, you may need to update:

- Genome references to hg38
- Python environment to 3.8+

## Reference Data

BETA includes reference gene annotations for:

- **Human**: hg38 (default), hg19, hg18
- **Mouse**: mm10, mm9

Default CTCF boundary data included for hg19 and mm9.

### Custom Genomes

For other genome assemblies, provide your own reference:

```bash
beta basic \
  -p peaks.bed \
  -e expression.txt \
  -k LIM \
  -r custom_refseq.txt \
  -n experiment
```

RefSeq format: tab-delimited with columns:
```
bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds score name2 ...
```

## Citation

If you use BETA in your research, please cite:

Wang S, Sun H, Ma J, et al. Target analysis by integration of transcriptome and ChIP-seq data with BETA. *Nature Protocols*. 2013;8(12):2502-2515. doi:10.1038/nprot.2013.150

## License

BETA is distributed under the Artistic License 2.0.

## Support

- **Documentation**: http://cistrome.org/BETA/tutorial.html
- **Issues**: https://github.com/yourusername/BETA2/issues
- **Original paper**: https://doi.org/10.1038/nprot.2013.150

## Authors

- **Original Author**: Su Wang (wangsu0623atgmail.com)
- **Python 3 Port**: Tommy Tang (tangming2005atgmail.com)

## Changelog

### Version 2.0.0 (2025)

- Python 3.8+ support (dropped Python 2)
- Modern project structure and packaging
- Default genome changed to hg38
- Improved code quality and type hints
- Enhanced logging and error handling
- Updated dependencies
- Performance optimizations

### Version 1.0.7 (2015)

- Original Python 2 version
- Basic, plus, and minus modes
- Support for hg18, hg19, mm9, mm10
