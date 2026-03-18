# Perform differential expression analysis for TRAP-seq data using edgeR

Compares IP vs INPUT fractions for a specified brain region and
treatment condition using edgeR's GLM framework. Accounts for paired
structure via a block variable (e.g., individual animal or tube).
Supports both the likelihood ratio test (`glmLRT`) and the
quasi-likelihood F-test (`glmQLFTest`) for the final hypothesis testing
step. Returns a tibble of differential expression results. Optionally
returns a `kableExtra` HTML table of the top DE genes.

## Usage

``` r
ptrap_de(
  counts_mat,
  sample_df,
  gene_ids,
  region_name,
  treatment_name,
  sample_col = "sample",
  fraction_col = "fraction",
  block_col = "tube",
  region_col = "BrainRegion",
  treatment_col = "Treatment",
  ip_level = "IP",
  input_level = "INPUT",
  lfc_threshold = 1,
  fdr_threshold = 0.05,
  test_method = c("LRT", "QLF"),
  ngenes.out = 20,
  kable.out = FALSE
)
```

## Arguments

- counts_mat:

  A numeric matrix of raw counts with genes as rows and samples as
  columns. Column names must match the values in `sample_col`.

- sample_df:

  A data frame containing sample metadata. Must include columns for
  sample names, IP/INPUT fraction, block variable, brain region, and
  treatment condition.

- gene_ids:

  A character vector of gene identifiers corresponding to the rows of
  `counts_mat`.

- region_name:

  The brain region to subset and analyze (e.g., `"POA"`). Must match a
  value in `region_col`.

- treatment_name:

  The treatment condition to subset and analyze (e.g., `"pb"`). Must
  match a value in `treatment_col`.

- sample_col:

  Name of the column in `sample_df` whose values match the column names
  of `counts_mat`. Default is `"sample"`.

- fraction_col:

  Name of the column in `sample_df` that distinguishes IP from INPUT
  fractions. Default is `"fraction"`.

- block_col:

  Name of the column in `sample_df` used as the blocking / pairing
  variable in the design matrix (e.g., individual animal or tube).
  Default is `"tube"`.

- region_col:

  Name of the column in `sample_df` containing brain region labels.
  Default is `"BrainRegion"`.

- treatment_col:

  Name of the column in `sample_df` containing treatment labels. Default
  is `"Treatment"`.

- ip_level:

  The value in `fraction_col` that identifies the IP fraction. Default
  is `"IP"`.

- input_level:

  The value in `fraction_col` that identifies the INPUT fraction (used
  as the reference level). Default is `"INPUT"`.

- lfc_threshold:

  Minimum absolute log2 fold change required to classify a gene as
  differentially expressed. Default is `1`.

- fdr_threshold:

  Maximum FDR allowed to classify a gene as differentially expressed.
  Default is `0.05`.

- test_method:

  The GLM testing method to use. Either `"LRT"` (likelihood ratio test
  via `glmFit` + `glmLRT`, default) or `"QLF"` (quasi-likelihood F-test
  via `glmQLFit` + `glmQLFTest`). `"QLF"` is generally more conservative
  and recommended when the number of samples per group is small.

- ngenes.out:

  Number of top genes (sorted by p-value) to include in the output when
  `kable.out = TRUE`. Default is `20`.

- kable.out:

  Logical. If `TRUE`, returns a `kableExtra` HTML table of the top
  `ngenes.out` genes instead of the full tibble. Requires the
  `kableExtra` package. Default is `FALSE`.

## Value

When `kable.out = FALSE` (default), a tibble with one row per gene,
sorted by p-value. When `kable.out = TRUE`, an HTML `kableExtra` table
of the top `ngenes.out` genes. The tibble contains:

- Gene:

  Gene identifier from `gene_ids`.

- logFC:

  Log2 fold change (IP vs INPUT).

- logCPM:

  Average log2 counts per million.

- LR / F:

  Test statistic (name depends on `test_method`).

- PValue:

  Raw p-value.

- FDR:

  Benjamini-Hochberg adjusted p-value.

- BrainRegion / region_col:

  Brain region label.

- Treatment / treatment_col:

  Treatment label.

- diffexpressed:

  `"UP"`, `"DOWN"`, or `"NO"` based on `lfc_threshold` and
  `fdr_threshold`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Default call using LRT
res_lrt <- ptrap_de(
  counts_mat     = counts_mat,
  sample_df      = sample_df,
  gene_ids       = gene_ids,
  region_name    = "POA",
  treatment_name = "pb"
)

# Using the quasi-likelihood F-test instead
res_qlf <- ptrap_de(
  counts_mat     = counts_mat,
  sample_df      = sample_df,
  gene_ids       = gene_ids,
  region_name    = "POA",
  treatment_name = "pb",
  test_method    = "QLF"
)

# Custom thresholds
res_custom <- ptrap_de(
  counts_mat     = counts_mat,
  sample_df      = sample_df,
  gene_ids       = gene_ids,
  region_name    = "POA",
  treatment_name = "pb",
  lfc_threshold  = 0.5,
  fdr_threshold  = 0.1
)
} # }
```
