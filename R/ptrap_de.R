# --- internal helpers ---------------------------------------------------------
# (not exported; used only by ptrap_de)

# Parse a single sample column name into treatment / block / fraction.
# Strategy: strip separators (_, -, .), split at every alpha<->digit boundary,
# then classify each token as fraction keyword | numeric block | treatment.
#
# Examples (all produce the same output fields):
#   "b1input"   -> treatment="b",  block="1", fraction=input_level
#   "nb1ip"     -> treatment="nb", block="1", fraction=ip_level
#   "Nb_IP_1"   -> treatment="Nb", block="1", fraction=ip_level
#   "B.3.INPUT" -> treatment="B",  block="3", fraction=input_level
.parse_one_sample <- function(nm, ip_level, input_level) {
  # remove separators then split at every letter<->digit transition
  # hyphen must be first in the character class to be treated as a literal
  # character and not as a range operator
  cleaned <- gsub("[-_. ]+", "", nm)
  parts <- unlist(strsplit(
    cleaned,
    "(?<=\\d)(?=\\D)|(?<=\\D)(?=\\d)",
    perl = TRUE
  ))
  parts <- parts[nzchar(parts)]
  parts_l <- tolower(parts)

  # fraction: matches ip_level, input_level, or "in" (short alias for INPUT)
  frac_kws <- unique(tolower(c(ip_level, input_level, "in")))
  frac_idx <- which(parts_l %in% frac_kws)
  # block: purely numeric token
  block_idx <- which(grepl("^\\d+$", parts))
  # treatment: everything that is neither fraction nor block
  treat_idx <- setdiff(seq_along(parts), c(frac_idx, block_idx))

  if (length(frac_idx) == 0L) {
    stop(
      "Cannot identify IP/INPUT fraction in column name '",
      nm,
      "'. ",
      "Column names must contain '",
      ip_level,
      "' or '",
      input_level,
      "' (case-insensitive). See ?ptrap_de for naming conventions."
    )
  }
  if (length(block_idx) == 0L) {
    stop(
      "Cannot identify a replicate number in column name '",
      nm,
      "'. ",
      "Column names must contain a digit identifying the biological replicate. ",
      "See ?ptrap_de for naming conventions."
    )
  }

  fraction <- if (parts_l[frac_idx[1L]] == tolower(ip_level)) {
    ip_level
  } else {
    input_level
  }
  block <- parts[block_idx[1L]]
  treatment <- paste(parts[treat_idx], collapse = "")

  list(sample = nm, treatment = treatment, block = block, fraction = fraction)
}

# Build a tibble of sample metadata by parsing every column name.
.build_sample_df_from_cols <- function(col_names, ip_level, input_level) {
  parsed <- lapply(
    col_names,
    .parse_one_sample,
    ip_level = ip_level,
    input_level = input_level
  )
  tibble::tibble(
    sample = vapply(parsed, `[[`, character(1L), "sample"),
    treatment = vapply(parsed, `[[`, character(1L), "treatment"),
    block = vapply(parsed, `[[`, character(1L), "block"),
    fraction = vapply(parsed, `[[`, character(1L), "fraction")
  )
}


# --- main function ------------------------------------------------------------

#' Perform differential expression analysis for TRAP-seq data using edgeR
#' or a paired t-test
#'
#' Compares IP vs INPUT fractions for a specified treatment condition using
#' one of three statistical approaches (see `test_method`). Both `sample_df`
#' and `gene_ids` are optional: the function can derive sample metadata
#' automatically from the column names of `counts_mat`, and gene identifiers
#' from its first column.
#'
#' * `"LRT"` and `"QLF"` use edgeR's GLM framework, account for the paired
#'   animal structure via a block variable, and are recommended when you have
#'   **4 or more replicates** and raw count data.
#' * `"paired.ttest"` runs a per-gene **paired t-test** between IP and INPUT
#'   counts across experimental repeats, following the method described in
#'   Tan et al. (2016) *Cell* 167, 47–59
#'   \doi{10.1016/j.cell.2016.08.028}. This approach is recommended when you
#'   have **only 3 replicates** (as is the typical minimum in PhosphoTRAP
#'   experiments), where GLM-based dispersion estimates are unreliable. The
#'   enrichment ratio (IP / INPUT) is computed per paired sample and expressed
#'   on the log2 scale, as described in the same paper. **Note:** if you pass
#'   normalised values (e.g., RPKM/TPM) instead of raw counts, set
#'   `pseudocount = 0`.
#'
#' @section Automatic column-name parsing:
#' When `sample_df = NULL`, the column names of `counts_mat` (excluding the
#' first, gene-ID column) are parsed to build sample metadata automatically.
#' Each column name must encode three pieces of information, in any order and
#' with any combination of separators (`_`, `-`, `.`, space) or no separator
#' at all:
#' \describe{
#'   \item{Treatment}{One or more letters identifying the experimental group.}
#'   \item{Replicate number}{A digit identifying the biological replicate.}
#'   \item{Fraction}{`ip_level` or `input_level` (case-insensitive); `"in"`
#'     is also accepted as a short alias for the INPUT fraction.}
#' }
#' Valid column name examples (default `ip_level = "IP"`,
#' `input_level = "INPUT"`):
#' \tabular{ll}{
#'   **Column name** \tab **Parsed as** \cr
#'   `b1input`, `b1ip`    \tab treatment = `b`, block = `1` \cr
#'   `nb2INPUT`           \tab treatment = `nb`, block = `2` \cr
#'   `Nb_IP_1`            \tab treatment = `Nb`, block = `1` \cr
#'   `B_3_INPUT`          \tab treatment = `B`, block = `3` \cr
#'   `PB.2.ip`            \tab treatment = `PB`, block = `2` \cr
#'   `Trim_10-INPUT`      \tab treatment = `Trim`, block = `10` \cr
#'   `SOL1INPUT`          \tab treatment = `SOL`, block = `1` \cr
#' }
#'
#' @param counts_mat A counts matrix in one of two formats:
#'   * **Matrix** — numeric, genes × samples; column names are sample IDs;
#'     gene IDs are in `rownames` or supplied via `gene_ids`.
#'   * **Data frame / tibble** — first column is a character vector of gene
#'     IDs; remaining columns are numeric counts with sample names as column
#'     names. When `sample_df = NULL`, column names must follow the naming
#'     convention described in the *Automatic column-name parsing* section.
#' @param sample_df Optional data frame of sample metadata. When `NULL`
#'   (default), metadata is parsed automatically from the column names of
#'   `counts_mat`. When provided, the arguments `sample_col`, `fraction_col`,
#'   `block_col`, `region_col`, and `treatment_col` specify which columns to
#'   use (falling back to their defaults if the names match).
#' @param gene_ids Optional character vector of gene identifiers corresponding
#'   to the rows of `counts_mat`. When `NULL` (default), gene IDs are
#'   extracted from the first column of `counts_mat` (if it is a data frame)
#'   or from `rownames(counts_mat)` (if it is a matrix).
#' @param region_name The brain region to subset and analyze (e.g., `"POA"`).
#'   Only used when `region_col` is not `NULL`. Can be `NULL` when
#'   `region_col = NULL` (default) or when the data contain a single region.
#' @param treatment_name The treatment condition whose samples will be
#'   **subsetted** for the IP vs INPUT comparison (e.g., `"pb"` to analyse
#'   only the pair-bonded samples). Note that `treatment_name` is **not** the
#'   DE contrast variable — the contrast is always IP vs INPUT. To compare two
#'   treatments you call `ptrap_de()` once per treatment and then pass both
#'   results to [pTRAPPING::ptrap_volcano2()].
#'
#'   When `NULL` (default), the value is derived automatically from
#'   `sample_df` (or from the parsed column names when `sample_df = NULL`):
#'   if a single treatment is found the function proceeds silently; if
#'   multiple treatments are found an error asks the user to specify which
#'   one to analyse.
#' @param sample_col Name of the column in `sample_df` whose values match the
#'   column names of `counts_mat`. Ignored when `sample_df = NULL` (auto-set
#'   internally). Default is `"sample"`.
#' @param fraction_col Name of the column in `sample_df` that distinguishes IP
#'   from INPUT fractions. Default is `"fraction"`.
#' @param block_col Name of the column in `sample_df` used as the blocking /
#'   pairing variable (e.g., individual animal or tube). For `"LRT"` / `"QLF"`
#'   it enters the design matrix; for `"paired.ttest"` it aligns each
#'   animal's IP with its own INPUT. Default is `"tube"`.
#' @param region_col Name of the column in `sample_df` containing brain region
#'   labels. Set to `NULL` (default) to skip region filtering — recommended
#'   when the data already contain only one brain region. Only set this when
#'   `counts_mat` contains multiple regions and you want to analyze one at a
#'   time. Default is `NULL`.
#' @param treatment_col Name of the column in `sample_df` containing treatment
#'   labels. Default is `"Treatment"`.
#' @param ip_level The value in `fraction_col` that identifies the IP fraction.
#'   Default is `"IP"`.
#' @param input_level The value in `fraction_col` that identifies the INPUT
#'   fraction (reference level). Default is `"INPUT"`.
#' @param lfc_threshold Minimum absolute log2 fold change required to classify
#'   a gene as differentially expressed. Default is `1`.
#' @param fdr_threshold Maximum FDR (or adjusted p-value for
#'   `"paired.ttest"`) allowed to classify a gene as differentially expressed.
#'   Default is `0.05`.
#' @param test_method Statistical method to use. One of:
#'   * `"LRT"` — likelihood ratio test via `glmFit` + `glmLRT` (default).
#'     Suitable for ≥ 4 replicates and raw count data.
#'   * `"QLF"` — quasi-likelihood F-test via `glmQLFit` + `glmQLFTest`.
#'     More conservative than LRT; recommended for small sample sizes with
#'     raw count data.
#'   * `"paired.ttest"` — per-gene paired t-test between IP and INPUT across
#'     replicates, following Tan et al. (2016)
#'     \doi{10.1016/j.cell.2016.08.028}. Best suited for n = 3 replicates
#'     typical of PhosphoTRAP experiments. P-values are adjusted with
#'     Benjamini-Hochberg. The enrichment ratio per paired sample (IP/INPUT)
#'     is expressed on a log2 scale.
#' @param pseudocount A small value added to counts before computing
#'   log2(IP/INPUT) in `"paired.ttest"`, to avoid log(0). Default is `1`.
#'   Set to `0` if `counts_mat` already contains normalised values (e.g.,
#'   RPKM/TPM) that are guaranteed to be > 0.
#' @param return_long Logical. Only used when `test_method = "paired.ttest"`.
#'   If `TRUE`, returns a named list with `$results` (the DE tibble) and
#'   `$long_data` (the per-gene, per-animal paired table used for the test).
#'   Default is `FALSE`.
#' @param ngenes.out Number of top genes (sorted by p-value) to include in the
#'   output when `kable.out = TRUE`. Default is `20`.
#' @param kable.out Logical. If `TRUE`, returns a `kableExtra` HTML table of
#'   the top `ngenes.out` genes instead of the full tibble. Default is
#'   `FALSE`.
#'
#' @return When `kable.out = FALSE` (default), a tibble with one row per gene
#'   sorted by p-value, containing `Gene`, `logFC`, test-statistic column(s),
#'   `PValue`, `FDR`, treatment label, optionally a region label (when
#'   `region_col` is set), and `diffexpressed` (`"UP"`, `"DOWN"`, `"NO"`).
#'   When `kable.out = TRUE`, an HTML `kableExtra` table of the top
#'   `ngenes.out` genes. For `"paired.ttest"` with `return_long = TRUE`, a
#'   named list with `$results` and `$long_data`.
#'
#' @references
#' Tan, C.L., Cooke, E.K., Leib, D.E., Lin, Y.C., Daly, G.E., Zimmerman,
#' C.A., and Knight, Z.A. (2016). Warm-Sensitive Neurons that Control Body
#' Temperature. *Cell* 167, 47–59.
#' \doi{10.1016/j.cell.2016.08.028}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## ---- Option A: simplest call — auto-parse from column names ---------------
#' # counts_mat is a data frame where:
#' #   col 1      = gene IDs
#' #   col 2+     = samples named like "b1input", "b1ip", "nb2input", etc.
#' counts <- read.table("counts.txt", header = TRUE)
#'
#' # single treatment in the matrix — treatment_name auto-detected
#' res <- ptrap_de(counts_mat = counts, test_method = "paired.ttest")
#'
#' # multiple treatments — specify which one to analyze
#' res_b <- ptrap_de(counts_mat = counts, treatment_name = "b")
#'
#' ## ---- Option B: provide sample_df explicitly (original workflow) -----------
#' res_lrt <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   region_name    = "POA",
#'   treatment_name = "pb"
#' )
#'
#' # quasi-likelihood F-test
#' res_qlf <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   region_name    = "POA",
#'   treatment_name = "pb",
#'   test_method    = "QLF"
#' )
#'
#' # paired t-test — also return long-format paired table
#' res_pt <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   region_name    = "POA",
#'   treatment_name = "pb",
#'   test_method    = "paired.ttest",
#'   return_long    = TRUE
#' )
#' res_pt$results    # DE tibble
#' res_pt$long_data  # per-gene, per-animal pairs
#' }
#'
#' @importFrom edgeR DGEList filterByExpr calcNormFactors estimateDisp
#'   glmFit glmLRT glmQLFit glmQLFTest topTags
#' @importFrom dplyr filter mutate case_when relocate slice_head arrange
#'   bind_rows
#' @importFrom rlang .data :=
#' @importFrom tibble as_tibble tibble
#' @importFrom stats as.formula model.matrix p.adjust pt sd
#' @importFrom knitr kable
#' @importFrom kableExtra kable_classic row_spec column_spec

ptrap_de <- function(
  counts_mat,
  sample_df = NULL,
  gene_ids = NULL,
  region_name = NULL,
  treatment_name = NULL,
  sample_col = "sample",
  fraction_col = "fraction",
  block_col = "tube",
  region_col = NULL,
  treatment_col = "Treatment",
  ip_level = "IP",
  input_level = "INPUT",
  lfc_threshold = 1,
  fdr_threshold = 0.05,
  test_method = c("LRT", "QLF", "paired.ttest"),
  pseudocount = 1,
  return_long = FALSE,
  ngenes.out = 20,
  kable.out = FALSE
) {
  test_method <- match.arg(test_method)

  # ---- 1. resolve gene IDs and ensure counts_mat is a numeric matrix ---------
  if (is.data.frame(counts_mat)) {
    first_col <- counts_mat[[1L]]
    if (is.character(first_col) || is.factor(first_col)) {
      # first column contains gene identifiers
      if (is.null(gene_ids)) {
        gene_ids <- as.character(first_col)
      }
      counts_mat <- as.matrix(counts_mat[, -1L, drop = FALSE])
      mode(counts_mat) <- "numeric"
    } else {
      counts_mat <- as.matrix(counts_mat)
    }
  }

  if (is.null(gene_ids)) {
    rn <- rownames(counts_mat)
    if (!is.null(rn)) {
      gene_ids <- rn
    } else {
      stop(
        "Cannot determine gene IDs. Please either:\n",
        "  (a) include gene IDs as the first column of `counts_mat`, or\n",
        "  (b) pass them explicitly via `gene_ids`."
      )
    }
  }

  # ---- 2. build sample metadata if not supplied ------------------------------
  if (is.null(sample_df)) {
    sample_df <- .build_sample_df_from_cols(
      colnames(counts_mat),
      ip_level = ip_level,
      input_level = input_level
    )
    # override column name arguments to match the auto-built data frame
    sample_col <- "sample"
    fraction_col <- "fraction"
    block_col <- "block"
    treatment_col <- "treatment"

    message(
      "Auto-parsed sample metadata from column names.\n",
      "  Treatments : ",
      paste(sort(unique(sample_df$treatment)), collapse = ", "),
      "\n",
      "  Blocks     : ",
      paste(sort(unique(sample_df$block)), collapse = ", "),
      "\n",
      "  Fractions  : ",
      paste(sort(unique(sample_df$fraction)), collapse = ", ")
    )
  }

  # ---- 3. resolve treatment_name ---------------------------------------------
  if (is.null(treatment_name)) {
    unique_tx <- unique(sample_df[[treatment_col]])
    if (length(unique_tx) == 1L) {
      treatment_name <- unique_tx
      message("Using the only treatment found: '", treatment_name, "'")
    } else {
      stop(
        "Multiple treatments found (",
        paste(unique_tx, collapse = ", "),
        "). Please specify `treatment_name`."
      )
    }
  }

  # ---- 4. resolve region_name (only when region_col is set) ------------------
  if (!is.null(region_col) && is.null(region_name)) {
    unique_reg <- unique(sample_df[[region_col]])
    if (length(unique_reg) == 1L) {
      region_name <- unique_reg
      message("Using the only region found: '", region_name, "'")
    } else {
      stop(
        "Multiple regions found (",
        paste(unique_reg, collapse = ", "),
        "). Please specify `region_name`."
      )
    }
  }

  # ---- 5. subset samples for the specified treatment (and optionally region) -
  if (!is.null(region_col)) {
    region_samples <- sample_df |>
      filter(
        .data[[region_col]] == region_name,
        .data[[treatment_col]] == treatment_name
      )
  } else {
    region_samples <- sample_df |>
      filter(.data[[treatment_col]] == treatment_name)
  }

  region_samples <- region_samples |>
    mutate(
      !!fraction_col := factor(
        .data[[fraction_col]],
        levels = c(input_level, ip_level)
      )
    )

  if (nrow(region_samples) == 0L) {
    stop(
      "No samples found for treatment '",
      treatment_name,
      "'",
      if (!is.null(region_col)) {
        paste0(" in region '", region_name, "'")
      } else {
        ""
      },
      ". Check that these values exist in `sample_df`."
    )
  }

  # ---- 6. validate sample–column matching ------------------------------------
  missing_samples <- setdiff(region_samples[[sample_col]], colnames(counts_mat))
  if (length(missing_samples) > 0L) {
    stop(
      "The following samples are in `sample_df` but not in `counts_mat` columns: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  # subset counts matrix to selected samples
  counts_region <- counts_mat[, region_samples[[sample_col]], drop = FALSE]

  # ---- paired t-test branch (Tan et al. 2016) --------------------------------
  if (test_method == "paired.ttest") {
    # separate and sort IP / INPUT by block to ensure correct pairing
    ip_samples <- region_samples |>
      filter(.data[[fraction_col]] == ip_level) |>
      arrange(.data[[block_col]])

    input_samples <- region_samples |>
      filter(.data[[fraction_col]] == input_level) |>
      arrange(.data[[block_col]])

    if (
      !identical(
        sort(as.character(ip_samples[[block_col]])),
        sort(as.character(input_samples[[block_col]]))
      )
    ) {
      stop(
        "IP and INPUT samples do not share the same '",
        block_col,
        "' values. ",
        "Each animal/tube must have exactly one IP and one INPUT sample."
      )
    }

    # assign gene IDs as rownames, then drop all-zero genes
    rownames(counts_region) <- gene_ids
    counts_region <- counts_region[rowSums(counts_region) > 0, , drop = FALSE]

    ip_mat <- counts_region[, ip_samples[[sample_col]], drop = FALSE]
    input_mat <- counts_region[, input_samples[[sample_col]], drop = FALSE]
    n_reps <- ncol(ip_mat)

    # ---- long-format paired table --------------------------------------------
    long_data <- bind_rows(
      lapply(seq_len(n_reps), function(j) {
        tibble(
          Gene = rownames(ip_mat),
          !!block_col := ip_samples[[block_col]][j],
          ip_count = as.numeric(ip_mat[, j]),
          input_count = as.numeric(input_mat[, j])
        )
      })
    ) |>
      arrange(.data$Gene)

    # ---- vectorised paired t-test --------------------------------------------
    diff_mat <- ip_mat - input_mat
    mean_diff <- rowMeans(diff_mat)
    sd_diff <- apply(diff_mat, 1L, sd)
    t_stat <- mean_diff / (sd_diff / sqrt(n_reps))
    p_val <- 2 * pt(-abs(t_stat), df = n_reps - 1L)

    # ---- log2 enrichment ratio (IP/INPUT) with pseudocount -------------------
    log2_ratios <- log2((ip_mat + pseudocount) / (input_mat + pseudocount))
    mean_log2fc <- rowMeans(log2_ratios)

    # ---- assemble results tibble ---------------------------------------------
    results <- tibble(
      Gene = rownames(ip_mat),
      logFC = as.numeric(mean_log2fc),
      t_statistic = as.numeric(t_stat),
      PValue = as.numeric(p_val)
    ) |>
      mutate(
        FDR = p.adjust(.data$PValue, method = "BH"),
        !!treatment_col := treatment_name,
        diffexpressed = case_when(
          .data$logFC > lfc_threshold & .data$FDR < fdr_threshold ~ "UP",
          .data$logFC < -lfc_threshold & .data$FDR < fdr_threshold ~ "DOWN",
          TRUE ~ "NO"
        )
      ) |>
      arrange(.data$PValue)

    if (!is.null(region_col)) {
      results <- results |> mutate(!!region_col := region_name)
    }

    # ---- return --------------------------------------------------------------
    if (kable.out) {
      return(
        results |>
          slice_head(n = ngenes.out) |>
          kable(
            digits = 3L,
            table.attr = 'data-quarto-disable-processing="true"',
            "html"
          ) |>
          kable_classic(full_width = FALSE, html_font = "Cambria") |>
          row_spec(0L, italic = TRUE, bold = TRUE) |>
          column_spec(1L, italic = FALSE, bold = TRUE)
      )
    }

    if (return_long) {
      return(list(results = results, long_data = long_data))
    }
    return(results)
  }

  # ---- edgeR branch (LRT / QLF) ----------------------------------------------

  dge <- DGEList(counts = counts_region)
  dge$genes <- data.frame(Gene = gene_ids)

  keep <- filterByExpr(dge, group = region_samples[[fraction_col]])
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  design_formula <- as.formula(paste("~", fraction_col, "+", block_col))
  design <- model.matrix(design_formula, data = region_samples)

  dge <- estimateDisp(dge, design)
  lrt_coef <- paste0(fraction_col, ip_level)

  if (test_method == "LRT") {
    fit <- glmFit(dge, design)
    test <- glmLRT(fit, coef = lrt_coef)
  } else {
    fit <- glmQLFit(dge, design)
    test <- glmQLFTest(fit, coef = lrt_coef)
  }

  # topTags already carries the Gene column from dge$genes — do not re-add it
  results <- topTags(test, n = Inf)$table |>
    as_tibble() |>
    mutate(
      !!treatment_col := treatment_name,
      diffexpressed = case_when(
        .data$logFC > lfc_threshold & .data$FDR < fdr_threshold ~ "UP",
        .data$logFC < -lfc_threshold & .data$FDR < fdr_threshold ~ "DOWN",
        TRUE ~ "NO"
      )
    )

  if (!is.null(region_col)) {
    results <- results |> mutate(!!region_col := region_name)
  }

  results <- results |> relocate("Gene")

  if (kable.out) {
    return(
      results |>
        slice_head(n = ngenes.out) |>
        kable(
          digits = 2L,
          table.attr = 'data-quarto-disable-processing="true"',
          "html"
        ) |>
        kable_classic(full_width = FALSE, html_font = "Cambria") |>
        row_spec(0L, italic = TRUE, bold = TRUE) |>
        column_spec(1L, italic = FALSE, bold = TRUE)
    )
  }

  return(results)
}
