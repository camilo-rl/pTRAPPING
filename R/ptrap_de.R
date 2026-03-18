#' Perform differential expression analysis for TRAP-seq data using edgeR
#'
#' Compares IP vs INPUT fractions for a specified brain region and treatment
#' condition using edgeR's GLM framework. Accounts for paired structure via a
#' block variable (e.g., individual animal or tube). Supports both the
#' likelihood ratio test (`glmLRT`) and the quasi-likelihood F-test
#' (`glmQLFTest`) for the final hypothesis testing step. Returns a tibble of
#' differential expression results. Optionally returns a `kableExtra` HTML table of the top DE genes.
#'
#' @param counts_mat A numeric matrix of raw counts with genes as rows and
#'   samples as columns. Column names must match the values in `sample_col`.
#' @param sample_df A data frame containing sample metadata. Must include
#'   columns for sample names, IP/INPUT fraction, block variable, brain region,
#'   and treatment condition.
#' @param gene_ids A character vector of gene identifiers corresponding to the
#'   rows of `counts_mat`.
#' @param region_name The brain region to subset and analyze (e.g., `"POA"`).
#'   Must match a value in `region_col`.
#' @param treatment_name The treatment condition to subset and analyze (e.g.,
#'   `"pb"`). Must match a value in `treatment_col`.
#' @param sample_col Name of the column in `sample_df` whose values match the
#'   column names of `counts_mat`. Default is `"sample"`.
#' @param fraction_col Name of the column in `sample_df` that distinguishes IP
#'   from INPUT fractions. Default is `"fraction"`.
#' @param block_col Name of the column in `sample_df` used as the blocking /
#'   pairing variable in the design matrix (e.g., individual animal or tube).
#'   Default is `"tube"`.
#' @param region_col Name of the column in `sample_df` containing brain region
#'   labels. Default is `"BrainRegion"`.
#' @param treatment_col Name of the column in `sample_df` containing treatment
#'   labels. Default is `"Treatment"`.
#' @param ip_level The value in `fraction_col` that identifies the IP fraction.
#'   Default is `"IP"`.
#' @param input_level The value in `fraction_col` that identifies the INPUT
#'   fraction (used as the reference level). Default is `"INPUT"`.
#' @param lfc_threshold Minimum absolute log2 fold change required to classify
#'   a gene as differentially expressed. Default is `1`.
#' @param fdr_threshold Maximum FDR allowed to classify a gene as
#'   differentially expressed. Default is `0.05`.
#' @param test_method The GLM testing method to use. Either `"LRT"` (likelihood
#'   ratio test via `glmFit` + `glmLRT`, default) or `"QLF"` (quasi-likelihood
#'   F-test via `glmQLFit` + `glmQLFTest`). `"QLF"` is generally more
#'   conservative and recommended when the number of samples per group is small.
#' @param ngenes.out Number of top genes (sorted by p-value) to include in the
#'   output when `kable.out = TRUE`. Default is `20`.
#' @param kable.out Logical. If `TRUE`, returns a `kableExtra` HTML table of
#'   the top `ngenes.out` genes instead of the full tibble. Requires the
#'   `kableExtra` package. Default is `FALSE`.
#'
#' @return When `kable.out = FALSE` (default), a tibble with one row per gene,
#'   sorted by p-value. When `kable.out = TRUE`, an HTML `kableExtra` table of
#'   the top `ngenes.out` genes. The tibble contains:
#'   \describe{
#'     \item{Gene}{Gene identifier from `gene_ids`.}
#'     \item{logFC}{Log2 fold change (IP vs INPUT).}
#'     \item{logCPM}{Average log2 counts per million.}
#'     \item{LR / F}{Test statistic (name depends on `test_method`).}
#'     \item{PValue}{Raw p-value.}
#'     \item{FDR}{Benjamini-Hochberg adjusted p-value.}
#'     \item{BrainRegion / region_col}{Brain region label.}
#'     \item{Treatment / treatment_col}{Treatment label.}
#'     \item{diffexpressed}{`"UP"`, `"DOWN"`, or `"NO"` based on
#'       `lfc_threshold` and `fdr_threshold`.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Default call using LRT
#' res_lrt <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   region_name    = "POA",
#'   treatment_name = "pb"
#' )
#'
#' # Using the quasi-likelihood F-test instead
#' res_qlf <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   region_name    = "POA",
#'   treatment_name = "pb",
#'   test_method    = "QLF"
#' )
#'
#' # Custom thresholds
#' res_custom <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   region_name    = "POA",
#'   treatment_name = "pb",
#'   lfc_threshold  = 0.5,
#'   fdr_threshold  = 0.1
#' )
#' }
#'
#' @importFrom edgeR DGEList filterByExpr calcNormFactors estimateDisp
#'   glmFit glmLRT glmQLFit glmQLFTest topTags
#' @importFrom dplyr filter mutate case_when relocate slice_head
#' @importFrom rlang .data :=
#' @importFrom tibble as_tibble
#' @importFrom stats as.formula model.matrix
#' @importFrom knitr kable
#' @importFrom kableExtra kable_classic row_spec column_spec

ptrap_de <- function(
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
) {
  test_method <- match.arg(test_method)

  # subset samples for the specified region and treatment
  region_samples <- sample_df |>
    filter(
      .data[[region_col]] == region_name,
      .data[[treatment_col]] == treatment_name
    ) |>
    mutate(
      !!fraction_col := factor(
        .data[[fraction_col]],
        levels = c(input_level, ip_level)
      )
    )

  # validate that the region/treatment combination exists in sample_df
  if (nrow(region_samples) == 0) {
    stop(
      "No samples found for ",
      region_col,
      " = '",
      region_name,
      "' and ",
      treatment_col,
      " = '",
      treatment_name,
      "'. ",
      "Check that these values exist in `sample_df`."
    )
  }

  # validate that sample names from sample_df exist as columns in counts_mat
  missing_samples <- setdiff(region_samples[[sample_col]], colnames(counts_mat))
  if (length(missing_samples) > 0) {
    stop(
      "The following samples are in `sample_df` but not in `counts_mat` columns: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  # subset counts matrix to only include columns for the selected samples
  counts_region <- counts_mat[, region_samples[[sample_col]], drop = FALSE]

  # create DGEList object for edgeR analysis
  dge <- DGEList(counts = counts_region)

  # add gene IDs to DGEList object for later reference
  dge$genes <- data.frame(Gene = gene_ids)

  # filter out lowly expressed genes based on counts and groupings
  keep <- filterByExpr(dge, group = region_samples[[fraction_col]])

  # subset DGEList to only include genes that passed filtering, and reset library sizes
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  # calculate normalization factors to account for differences in library sizes and composition
  dge <- calcNormFactors(dge)

  # create design matrix for the model, including fraction (IP vs INPUT) and block variable as covariates
  design_formula <- as.formula(paste("~", fraction_col, "+", block_col))
  design <- model.matrix(design_formula, data = region_samples)

  # estimate dispersion parameters for the model
  dge <- estimateDisp(dge, design)

  # fit GLM and run the chosen test
  lrt_coef <- paste0(fraction_col, ip_level)

  if (test_method == "LRT") {
    fit <- glmFit(dge, design)
    test <- glmLRT(fit, coef = lrt_coef)
  } else {
    fit <- glmQLFit(dge, design)
    test <- glmQLFTest(fit, coef = lrt_coef)
  }

  # extract results and apply user-defined thresholds for DE classification
  # note: topTags already carries the Gene column from dge$genes, so we do not
  # re-assign it here (doing so would attempt to bind the full unfiltered
  # gene_ids vector, causing a size mismatch error)
  results <- topTags(test, n = Inf)$table |>
    as_tibble() |>
    mutate(
      !!region_col    := region_name,
      !!treatment_col := treatment_name,
      diffexpressed = case_when(
        .data$logFC >  lfc_threshold & .data$FDR < fdr_threshold ~ "UP",
        .data$logFC < -lfc_threshold & .data$FDR < fdr_threshold ~ "DOWN",
        TRUE ~ "NO"
      )
    ) |>
    relocate("Gene")

  if (kable.out) {
    return(
      results |>
        slice_head(n = ngenes.out) |>
        kable(
          digits = 2,
          table.attr = 'data-quarto-disable-processing="true"',
          "html"
        ) |>
        kable_classic(full_width = F, html_font = "Cambria") |>
        row_spec(0, italic = T, bold = T) |>
        column_spec(1, italic = F, bold = T)
    )
  }

  return(results)
}
