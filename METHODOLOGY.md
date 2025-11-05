# BETA Methodology: Technical Details

This document provides a step-by-step technical explanation of how BETA calculates regulatory potential scores, predicts transcription factor function, and identifies direct target genes.

## Table of Contents
- [Overview](#overview)
- [Step 1: Regulatory Potential Score](#step-1-regulatory-potential-score)
- [Step 2: Activator/Repressor Function Prediction](#step-2-activatorrepressor-function-prediction)
- [Step 3: Direct Target Gene Prediction](#step-3-direct-target-gene-prediction)
- [Step 4: Motif Analysis](#step-4-motif-analysis-optional)
- [Implementation Notes](#implementation-notes)

---

## Overview

BETA integrates ChIP-seq binding data with differential gene expression data through three main computational steps:

1. **Regulatory Potential Scoring**: Quantifies the likelihood that a transcription factor regulates each gene based on nearby binding sites
2. **Function Prediction**: Uses statistical testing to determine if the factor acts as an activator, repressor, or both
3. **Target Identification**: Combines binding and expression evidence using rank product to identify direct targets

---

## Step 1: Regulatory Potential Score

### Biological Rationale

Transcription factor binding sites do not have a simple one-to-one relationship with target genes. A single binding site can regulate multiple genes, and the regulatory effect typically decreases with distance from the transcription start site (TSS). BETA models this using a distance-weighted scoring system.

### Mathematical Formula

For each gene `g`, the regulatory potential score `S_g` is calculated as:

```
S_g = Σ exp(-0.5 - 4 × Δ_i)
```

Where:
- `Δ_i` = normalized distance from binding site `i` to the gene TSS
- `Δ` = (absolute distance in bp) / (distance cutoff), default cutoff = 100,000 bp
- The summation is over all binding sites within the distance cutoff

### Implementation Details

**Code reference**: `src/beta/core/pscore.py`, line 36
```python
Sg = lambda ldx: sum([math.exp(-0.5 - 4 * t) for t in ldx])
```

**Distance Normalization**:
- If binding site is 10 kb from TSS and cutoff is 100 kb: Δ = 10,000/100,000 = 0.1
- If binding site is 50 kb from TSS: Δ = 50,000/100,000 = 0.5

**Example Calculation**:

For a gene with three binding sites at 5 kb, 20 kb, and 80 kb from TSS (with 100 kb cutoff):

```
Δ₁ = 5,000/100,000 = 0.05
Δ₂ = 20,000/100,000 = 0.2
Δ₃ = 80,000/100,000 = 0.8

S_g = exp(-0.5 - 4×0.05) + exp(-0.5 - 4×0.2) + exp(-0.5 - 4×0.8)
    = exp(-0.7) + exp(-1.3) + exp(-3.7)
    = 0.4966 + 0.2725 + 0.0247
    = 0.7938
```

**Key Properties**:
- Closer binding sites contribute more to the score (exponential decay)
- Sites very far from TSS contribute minimally
- Multiple binding sites have additive effects
- At TSS (Δ=0): contribution = exp(-0.5) ≈ 0.606
- At cutoff limit (Δ=1): contribution = exp(-4.5) ≈ 0.011

### Optional: CTCF Boundary Filtering

If `--bl` (boundary limit) is enabled, BETA only considers binding sites within the same CTCF-delimited topologically associating domain (TAD) as the gene TSS. This reflects the biological reality that CTCF binding sites can act as insulators, preventing distal regulatory elements from affecting genes in neighboring domains.

**Implementation**: For hg19 and mm9, BETA provides pre-computed CTCF boundary files derived from ENCODE data. Peaks are filtered to only include those in the same boundary block as the gene.

### Ranking Genes by Regulatory Potential

After calculating scores for all genes, BETA ranks them from highest to lowest regulatory potential:
- **Rank 1**: Gene with highest regulatory potential score
- **Rank n**: Gene with lowest score

Genes without any binding sites within the distance cutoff receive a score of 0 and rank of "NA".

---

## Step 2: Activator/Repressor Function Prediction

### Biological Question

Does the transcription factor primarily activate gene expression, repress it, or do both? This is answered by testing whether genes with strong binding potential are enriched among upregulated or downregulated genes.

### Statistical Test: Kolmogorov-Smirnov (KS) Test

BETA uses a one-tailed Kolmogorov-Smirnov test to compare cumulative distributions of regulatory potential scores across three gene groups:

1. **UP**: Upregulated genes (significant positive fold change)
2. **DOWN**: Downregulated genes (significant negative fold change)
3. **NON**: Non-differentially expressed genes (background)

### How the KS Test Works

The KS test calculates the maximum vertical distance (D statistic) between two cumulative distribution functions (CDFs):

**Null Hypothesis**: The regulatory potential scores of UP (or DOWN) genes are drawn from the same distribution as NON genes.

**Alternative Hypothesis (one-tailed)**: UP (or DOWN) genes have systematically higher regulatory potential scores than NON genes.

**Test Statistic**:
```
D = max|CDF_UP(x) - CDF_NON(x)|
```

For a one-tailed test checking if UP genes have higher scores:
```
D⁺ = max(CDF_UP(x) - CDF_NON(x))
```

**P-value Interpretation**:
- Small p-value (e.g., p < 0.001) indicates UP genes have significantly higher regulatory potential
- If p_UP < cutoff (default 0.001): **Factor is an ACTIVATOR**
- If p_DOWN < cutoff: **Factor is a REPRESSOR**
- If both p_UP and p_DOWN < cutoff: **Factor is BOTH**

### Implementation Details

**Code reference**: `src/beta/core/up_down_distance.py` and `up_down_score.py`

The R script generated by BETA performs:

```R
# Genes ranked by regulatory potential (high to low)
d_up <- regulatory_potential_scores[upregulated_genes]
d_down <- regulatory_potential_scores[downregulated_genes]
d_bg <- regulatory_potential_scores[nondiff_genes]

# One-tailed KS test
ks_up <- ks.test(d_up, d_bg, alternative = "greater")
ks_down <- ks.test(d_down, d_bg, alternative = "greater")

p_up <- ks_up$p.value
p_down <- ks_down$p.value
```

**Visualization**: The output PDF shows cumulative distribution curves:
- **X-axis**: Genes ranked by regulatory potential (high to low)
- **Y-axis**: Cumulative fraction of genes (0-100%)
- **Lines**:
  - Dashed black: NON (background)
  - Red: UP genes
  - Blue: DOWN genes

If the red line rises faster than the dashed line, upregulated genes are enriched for high regulatory potential (activator function).

### Gene Classification Parameters

Users control which genes are classified as UP/DOWN vs NON:

- `--df <FDR>`: Use FDR threshold (e.g., 0.05)
- `--da <number or fraction>`:
  - If < 1: proportion (e.g., 0.5 = top 50%)
  - If > 1: absolute count (e.g., 500 = top 500 genes)

---

## Step 3: Direct Target Gene Prediction

### The Rank Product Approach

BETA identifies direct targets by requiring BOTH strong binding evidence AND significant expression changes. This is achieved through the rank product method (Breitling et al., 2004, FEBS Letters).

### Step-by-Step Calculation

**Input**:
- `n` genes with both regulatory potential scores AND differential expression data
- Genes without binding (score = 0) OR without expression data are excluded

**Step 3.1: Calculate Two Rankings**

For each gene `g`:

1. **Binding Rank (R_gb)**:
   - Rank genes by decreasing regulatory potential score
   - R_gb = 1 for highest score
   - R_gb = n for lowest score

2. **Expression Rank (R_ge)**:
   - Rank genes by increasing FDR/p-value (or by decreasing |log fold change|)
   - R_ge = 1 for most significant expression change
   - R_ge = n for least significant change

**Step 3.2: Calculate Rank Product**

```
RP_g = (R_gb / n) × (R_ge / n)
```

**Interpretation**:
- RP ranges from (1/n)² to 1
- Lower RP = more confident direct target
- RP can be interpreted as a probability: the chance that a gene has binding rank ≤ R_gb AND expression rank ≤ R_ge by chance alone

**Example**:

Suppose n = 1000 genes, and for gene X:
- R_gb = 10 (10th highest regulatory potential)
- R_ge = 5 (5th most significant expression change)

```
RP_X = (10/1000) × (5/1000) = 0.01 × 0.005 = 0.00005 = 5×10⁻⁵
```

This gene ranks in the top 1% for binding AND top 0.5% for expression, giving very strong evidence for direct regulation.

### Implementation Details

**Code reference**: `src/beta/core/expr_combine.py`, lines 522-524

```python
genenumber = total_genes_with_binding_and_expression
for gene in genes:
    RP = (float(brank[gene][2]) / float(genenumber)) * \
         (float(erank[gene][2]) / float(genenumber))
```

Where:
- `brank[gene][2]` = binding rank for gene
- `erank[gene][2]` = expression rank for gene
- `genenumber` = total genes being considered

### Output

**Target files** (`{name}_uptarget.txt`, `{name}_downtarget.txt`):

Columns:
1-3: Chromosome, TSS, TTS (BED format)
4: RefSeq ID
5: **Rank Product** (RP value)
6: Strand
7: Gene Symbol

Genes are sorted by increasing RP (best targets first).

**Recommended Cutoff**: RP < 0.001 (genes in top ~3.2% for both binding and expression)

### Biological Interpretation

**Why Rank Product Works**:

1. **Reduces false positives from non-functional binding**: A gene with strong binding (low R_gb) but no expression change (high R_ge) will have moderate RP
2. **Reduces false positives from indirect targets**: A gene with strong expression change (low R_ge) but no nearby binding (high R_gb) will have moderate RP
3. **Prioritizes concordant evidence**: Genes scoring well on BOTH criteria get very low RP values

**Example Scenarios**:

| Gene | R_gb | R_ge | RP | Interpretation |
|------|------|------|-----|----------------|
| A | 5 | 8 | (5/1000)×(8/1000) = 4×10⁻⁵ | **Strong direct target** |
| B | 5 | 500 | (5/1000)×(500/1000) = 2.5×10⁻³ | Strong binding, weak expression - may not be functional |
| C | 500 | 8 | (500/1000)×(8/1000) = 4×10⁻³ | Strong expression, weak binding - likely indirect |
| D | 500 | 500 | (500/1000)×(500/1000) = 0.25 | Not a target |

---

## Step 4: Motif Analysis (Optional)

Available only in `beta plus` mode with `--gs` (genome sequence) parameter.

### Purpose

Identify DNA-binding motifs enriched near predicted target genes to discover:
1. The transcription factor's own binding motif (validation)
2. Motifs of collaborating factors (cofactors)
3. Different motifs associated with activation vs repression

### Method: Position-Specific Scoring Matrix (PSSM) Scanning

**Tool**: MISP (Model-based Interval Scanner with PSSM), which implements the MOODS algorithm (Korhonen et al., 2009)

### Step 4.1: Extract Peak Sequences

For each predicted target gene:
- Extract genomic sequences from peak regions
- Focus on 200 bp centered on peak summit

### Step 4.2: Motif Scanning

Scan sequences with a library of known TF motifs (PSSMs) to calculate motif scores for each position.

### Step 4.3: Summit Enrichment Test

**Key Innovation**: Test whether motif instances are enriched at peak summits vs flanking regions

For each motif:
- **Middle region**: 200 bp centered on peak summit
- **Left region**: 200 bp immediately upstream
- **Right region**: 200 bp immediately downstream

**Test statistic**: One-tailed t-test comparing motif scores in middle vs (left + right)

```
H₀: Mean motif score in middle = mean in flanking regions
H₁: Mean motif score in middle > mean in flanking regions
```

**Why this works**: True TF binding sites should be centered on ChIP-seq peak summits, so their motifs should show summit enrichment.

### Step 4.4: Differential Motif Analysis

Compare motifs between gene sets:

1. **UP vs NON**: Motifs enriched near upregulated targets vs non-targets
2. **DOWN vs NON**: Motifs enriched near downregulated targets vs non-targets
3. **UP vs DOWN**: Motifs differentially enriched (identifies activation vs repression cofactors)

**Output**: HTML report with motif logos, p-values, and t-scores grouped by similarity

---

## Implementation Notes

### Programming Languages

- **Core algorithms**: Python 3.8+
- **Statistical analysis**: R scripts generated dynamically
- **Motif scanning**: C extension module (requires gcc)

### Performance Considerations

**Typical runtime** (single-threaded):
- BETA-basic: ~2 minutes
- BETA-minus: ~1 minute
- BETA-plus: ~20 minutes (motif analysis is expensive)

**Memory usage**: ~1-2 GB RAM for typical datasets

### Preprocessing Requirements

1. **ChIP-seq peaks**:
   - Pre-called with MACS2, SICER, or similar
   - BED format (minimum 3 columns)
   - Default: top 10,000 peaks used (adjustable with `--pn`)

2. **Gene expression**:
   - Pre-analyzed with LIMMA, DESeq2, Cuffdiff, or similar
   - Must include: gene ID, log fold change, statistical significance
   - Gene IDs as RefSeq accessions or official gene symbols

3. **Genome annotation**:
   - RefSeq gene models (provided for hg19, hg38, mm9, mm10, hg18)
   - TSS positions used for distance calculations

### Parameter Tuning

**Distance cutoff (`-d`)**: Default 100,000 bp
- Increase for factors with many distal enhancers
- Decrease for factors binding primarily at promoters

**Differential gene cutoff (`--da`, `--df`)**:
- Too stringent: miss real targets with moderate changes
- Too loose: include indirect targets
- Default (top 50%) works well for most datasets
- For high-quality RNA-seq: consider FDR < 0.05

**KS test cutoff (`-c`)**: Default 0.001
- Controls which gene groups (UP/DOWN/BOTH) are analyzed
- More stringent (e.g., 0.0001): only very clear activators/repressors
- More relaxed (e.g., 0.01): detect weak regulatory effects

---

## References

1. Wang et al. (2013). "Target analysis by integration of transcriptome and ChIP-seq data with BETA." *Nature Protocols* 8(12):2502-2515. [doi:10.1038/nprot.2013.150](https://doi.org/10.1038/nprot.2013.150)

2. Breitling et al. (2004). "Rank products: a simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments." *FEBS Letters* 573:83-92.

3. Subramanian et al. (2005). "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles." *PNAS* 102(43):15545-15550.

4. Korhonen et al. (2009). "MOODS: fast search for position weight matrix matches in DNA sequences." *Bioinformatics* 25(23):3181-3182.

---

## Glossary

**Regulatory Potential Score**: Quantitative measure of how likely a transcription factor regulates a gene based on binding site proximity and number

**Rank Product**: Geometric mean of normalized ranks, used to identify genes scoring highly on multiple independent criteria

**Kolmogorov-Smirnov (KS) Test**: Non-parametric test comparing cumulative distributions of two samples

**Position-Specific Scoring Matrix (PSSM)**: Mathematical representation of a DNA-binding motif showing nucleotide preferences at each position

**TSS (Transcription Start Site)**: Genomic position where transcription begins

**FDR (False Discovery Rate)**: Expected proportion of false positives among called significant results, adjusting for multiple testing

**ChIP-seq Peak**: Genomic region with significant enrichment of ChIP-seq signal indicating protein-DNA binding
