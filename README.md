# BETA2: Binding and Expression Target Analysis

[![PyPI version](https://img.shields.io/pypi/v/beta-binding-analysis.svg)](https://pypi.org/project/beta-binding-analysis/)
[![Python Version](https://img.shields.io/pypi/pyversions/beta-binding-analysis.svg)](https://pypi.org/project/beta-binding-analysis/)
[![License](https://img.shields.io/badge/license-Artistic%202.0-green)](LICENSE)
[![PyPI downloads](https://img.shields.io/pypi/dm/beta-binding-analysis)](https://pypi.org/project/beta-binding-analysis/)

BETA is a computational tool for integrative analysis of ChIP-seq and RNA-seq/microarray data to predict transcription factor (TF) direct target genes and identify whether the TF primarily functions as a transcriptional activator or repressor.

## The Biological Problem

When you perform ChIP-seq to find where a transcription factor binds and RNA-seq to see which genes change expression, a critical question arises: **Which genes are direct targets of your factor versus indirect/secondary effects?**

Several challenges complicate this analysis:
- **No 1-to-1 mapping**: A single binding site can regulate multiple genes, and a gene can be regulated by multiple binding sites
- **Not all binding is functional**: Some ChIP-seq peaks may not actually regulate nearby genes due to lack of cofactors or unfavorable chromatin environment
- **Secondary effects**: Binding to direct target genes causes them to change expression, which then affects other genes downstream (indirect targets)

## What BETA Does

BETA addresses these challenges by integrating binding and expression data to answer three key questions:

1. **Is your factor an activator or repressor?**
   - Determines whether the factor primarily activates or represses gene expression by testing if genes with stronger binding potential are enriched among upregulated or downregulated genes

2. **Which genes are direct targets?**
   - Identifies genes that are most likely to be directly regulated by combining two lines of evidence: proximity/strength of binding AND expression changes
   - Genes with both high binding potential and differential expression are prioritized as direct targets

3. **What cofactors modulate the factor's function?** (optional)
   - Identifies DNA-binding motifs enriched near your factor's binding sites to discover collaborating transcription factors

## How BETA Works

**Regulatory Potential Model**: Instead of simply assigning the nearest gene to each peak, BETA calculates a "regulatory potential score" for each gene based on ALL nearby binding sites within a distance window (default 100kb). Binding sites closer to the transcription start site (TSS) contribute more to the score using an exponentially decaying distance function - this reflects the biological reality that closer regulatory elements generally have stronger effects.

**Rank Product Integration**: BETA ranks genes by two criteria:
1. Regulatory potential score (how much binding is nearby)
2. Differential expression significance (how much expression changed)

The rank product identifies genes that score well on BOTH criteria - these are the most confident direct targets. Genes that only show binding OR only show expression changes are deprioritized, reducing false positives from non-functional binding sites and indirect targets.

**Statistical Testing**: BETA uses the Kolmogorov-Smirnov test to determine if upregulated or downregulated genes have significantly higher regulatory potential scores than non-differentially expressed genes, revealing whether your factor functions as an activator, repressor, or both.

> **For Technical Details**: See [METHODOLOGY.md](METHODOLOGY.md) for a comprehensive step-by-step explanation of all calculations, formulas, and statistical tests with worked examples.

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

### From PyPI (Recommended)

```bash
pip install beta-binding-analysis
```

### From Source

```bash
git clone https://github.com/crazyhottommy/BETA2.git
cd BETA2
pip install -e .
```

## Quick Start

### Basic Analysis

Predict TF target genes and function (activator/repressor):

```bash
# Note: diff_expr.txt should be TF activation/stimulation vs control (not knockdown)
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

### Example with Test Data

Using the included test data for estrogen receptor (ER/ESR1):

```bash
beta basic \
  -p BETA_test_data/ER_349_peaks.bed \
  -e BETA_test_data/ESR1_diff_expr.xls \
  -k O \
  --info 1,2,6 \
  -g hg19 \
  -n ESR1
```

**Note**: When using custom format (`-k O`), the `--info` parameter specifies column positions: gene ID, log fold change, and statistical value (e.g., `1,2,6` for columns 1, 2, and 6). See Input Files section below for format details.

## Input Files

### ChIP-seq Peaks (required)

BED format file (3 or 5 columns):
```
chr1    1000    2000
chr1    5000    6000
```

**Peak Number Cutoff - IMPORTANT**:

By default, BETA only uses the **top 10,000 peaks** even if your BED file contains more peaks. This is controlled by the `--pn` parameter.

**Why this default exists**:
- Focus on high-confidence peaks (strongest peaks are most likely functional)
- Reduce computational time
- Reduce noise from weak/low-confidence binding events

**When to adjust**:
- **Use default (10,000)**: For most analyses, or if you have noisy ChIP-seq data
- **Increase the number**: If you have high-quality ChIP-seq with many strong peaks and want comprehensive analysis
  ```bash
  beta basic -p peaks.bed -e expr.txt -k LIM -g hg38 -n my_TF --pn 34000
  ```
- **Use all peaks**: Set to a very high number to ensure all peaks are included
  ```bash
  beta basic -p peaks.bed -e expr.txt -k LIM -g hg38 -n my_TF --pn 100000
  ```
- **Decrease**: If you want to focus only on the very strongest binding events
  ```bash
  beta basic -p peaks.bed -e expr.txt -k LIM -g hg38 -n my_TF --pn 5000
  ```

**Example**: If your BED file contains 34,000 peaks but you don't specify `--pn`, only the top 10,000 peaks (by score, if available) will be used in the analysis.

### Differential Expression (required for basic/plus modes)

**Experimental Design - IMPORTANT**:

The differential expression file should represent gene expression changes from **activating/stimulating your transcription factor**:

- **Preferred**: TF activation/overexpression/stimulation vs control
  - Example: AR activation (androgen treatment) vs vehicle control
  - Example: ESR1 activation (estrogen treatment) vs vehicle control
  - Example: TF overexpression vs empty vector control

- **If you only have knockdown/knockout/inhibition data**: You can use it by **flipping the sign of log2FC values**
  - If you have: TF knockdown vs control
  - Simply multiply all log2FC values by -1
  - Example: Gene X has log2FC = -2.5 in knockdown → use log2FC = +2.5 for BETA
  - Keep p-values and FDR unchanged (only flip the fold change sign)
  - Biological interpretation: If knocking down the TF decreases a gene's expression, that gene is likely activated by the TF (so activating the TF would increase it)

**Why this matters**: BETA determines whether your TF is an activator or repressor by testing if genes with high binding potential are enriched among upregulated or downregulated genes. This logic requires that you provide the "TF activation" comparison. If your factor is an activator, activating it will increase expression of direct targets (positive logFC). If it's a repressor, activating it will decrease expression of direct targets (negative logFC).

**Gene ID Format**:

BETA supports two types of gene identifiers in differential expression files:

1. **RefSeq IDs** (default): Use transcript/gene accessions like `NM_001002231`, `NR_045762`, `XM_012345`
   - This is the default behavior - no additional flag needed
   - Example: `NM_001002231`, `NR_045762_at`

2. **Official Gene Symbols**: Use gene names like `TP53`, `MYC`, `BRCA1`
   - Requires the `--gname2` flag
   - Example command: `beta basic -p peaks.bed -e expression.txt -k LIM --gname2 -g hg38 -n experiment`

**Important**: All genes in your differential expression file must use the SAME identifier type (either all RefSeq IDs or all gene symbols). Do not mix them.

**Supported formats**:

1. **LIMMA** (`-k LIM`): Standard LIMMA output
2. **Cuffdiff** (`-k CUF`): Cuffdiff gene_exp.diff format
3. **BETA Standard Format** (`-k BSF`):
   ```
   GeneSymbol    log2FoldChange    FDR
   TP53          2.5               0.001
   MYC           -1.8              0.01
   ```
4. **Other/Custom** (`-k O`): Custom tab-delimited format with `--info` to specify columns

   **IMPORTANT**: The first line must start with `#` to be treated as a header. Without the `#`, BETA will try to parse the header as data and fail.

   **Example with RefSeq IDs** (columns 1, 2, 6 contain gene ID, logFC, adj.P.Val):
   ```
   #ID              logFC    AveExpr    t           P.Value      adj.P.Val
   NM_001002231     3.21     9.17       35.33       8.07e-11     4.18e-07
   NM_005551        2.14     8.45       28.15       1.23e-09     2.34e-06
   ```
   Then use: `--info 1,2,6`

   **Example with gene symbols** (same columns, add `--gname2` flag):
   ```
   #GeneSymbol      logFC    AveExpr    t           P.Value      adj.P.Val
   TP53             3.21     9.17       35.33       8.07e-11     4.18e-07
   MYC              2.14     8.45       28.15       1.23e-09     2.34e-06
   ```
   Then use: `--info 1,2,6 --gname2`

### Genome Sequence (for plus mode)

FASTA format genome sequence file (required for motif analysis)

## Output Files

### Basic Mode

- `{name}_uptarget.txt`: Up-regulated direct target genes
- `{name}_downtarget.txt`: Down-regulated direct target genes
- `{name}_function_prediction.pdf`: TF function prediction plot (activator/repressor)
- `{name}_function_prediction.R`: R script used to generate the plot

**Target file columns** (`{name}_uptarget.txt` and `{name}_downtarget.txt`):

| Column | Name | Description |
|--------|------|-------------|
| 1 | Chroms | Chromosome |
| 2 | txStart | Transcription start site |
| 3 | txEnd | Transcription end site |
| 4 | refseqID | RefSeq transcript ID |
| 5 | rank_product | Rank product score (lower = more confident target) |
| 6 | Strands | Strand (+/-) |
| 7 | GeneSymbol | Gene symbol |
| 8 | reg_potential | Regulatory potential score (based on nearby binding sites) |
| 9 | binding_rank | Rank by regulatory potential (1 = highest) |
| 10 | expr_rank | Rank by differential expression significance (1 = most significant) |
| 11 | log2FC | Log2 fold change from differential expression |
| 12 | padj | Adjusted p-value (FDR) from differential expression |

**Recommended cutoff**: Rank product < 0.001 for high-confidence targets.

**Usage tip**: The target files are sorted by rank product (column 5). You can filter by rank product, padj, or both depending on your analysis needs. The additional columns allow you to understand why each gene was ranked as a target and to perform downstream analyses.

**Multiple transcripts per gene**:

When using RefSeq IDs as input (default), you may see the same gene symbol appearing in multiple rows with different RefSeq transcript IDs. This is because:
- Regulatory potential is calculated for **each transcript individually** based on its TSS position
- Different transcripts of the same gene can have different TSS positions, resulting in different regulatory potential scores
- Each transcript is ranked separately

Example output showing multiple transcripts for KLK11:
```
chr19  51525486  51531290  NM_001167605  2.019e-06  -  KLK11  0.684  299  1  2.855  1.06e-08
chr19  51525486  51531290  NM_006853     2.019e-06  -  KLK11  0.684  299  1  2.855  1.06e-08
chr19  51525486  51530885  NM_144947     2.046e-06  -  KLK11  0.673  303  1  2.855  1.06e-08
```

**To get one row per gene**, you have two options:

1. **Use gene symbols as input** with the `--gname2` flag:
   ```bash
   beta basic -p peaks.bed -e expr_symbols.txt -k LIM --gname2 -g hg38 -n experiment
   ```
   This consolidates multiple transcripts by gene symbol, keeping the transcript with the best regulatory potential rank (no averaging).

2. **Post-process the output** to keep only the best-ranked transcript per gene:
   ```bash
   # Keep only the first (best) row for each gene symbol
   head -1 output_uptarget.txt > output_uptarget_unique.txt
   tail -n +2 output_uptarget.txt | sort -k5,5g | awk '!seen[$7]++' >> output_uptarget_unique.txt
   ```

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
| `--gname2` | Use official gene symbols instead of RefSeq IDs (default: False) |

## Examples

### Example 1: Basic TF Analysis (hg38)

```bash
# ChIP-seq peaks from ERalpha (ESR1) ChIP-seq
# Differential expression from estrogen treatment (activating ERalpha) vs vehicle control
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

### Example 4: Using Gene Symbols Instead of RefSeq IDs

```bash
beta basic \
  -p TF_peaks.bed \
  -e expression_symbols.txt \
  -k O \
  --info 1,2,6 \
  --gname2 \
  -g hg38 \
  -n TF_genesymbols
```
*(Use `--gname2` when your expression file uses official gene symbols like TP53, MYC, BRCA1 instead of RefSeq IDs)*

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
- **PyPI Package**: https://pypi.org/project/beta-binding-analysis/
- **Issues**: https://github.com/crazyhottommy/BETA2/issues
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
