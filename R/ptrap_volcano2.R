#' Dual volcano plot comparing two treatment conditions from TRAP-seq DE results
#'
#' Takes two tibbles produced by [pTRAPPING::ptrap_de()] — one per treatment condition,
#' both from the same brain region — joins them by gene, classifies each gene
#' according to its differential expression status in each condition, and
#' returns a scatter plot of logFC\eqn{_{\text{treatment 1}}} vs
#' logFC\eqn{_{\text{treatment 2}}}. Significant genes are highlighted and
#' labelled; threshold lines are drawn at ±`lfc_threshold` on both axes.
#'
#' @param de_result_1 A tibble returned by [pTRAPPING::ptrap_de()] for the **first**
#'   treatment condition (plotted on the **y-axis**).
#' @param de_result_2 A tibble returned by [pTRAPPING::ptrap_de()] for the **second**
#'   treatment condition (plotted on the **x-axis**).
#' @param lfc_threshold Minimum absolute log2 fold change used to define
#'   significance. Must match (or be stricter than) the threshold used in
#'   [pTRAPPING::ptrap_de()]. Default is `1`.
#' @param fdr_threshold Maximum FDR used to define significance. Default
#'   is `0.05`.
#' @param gene_col Name of the column containing gene identifiers in both
#'   result tibbles. Default is `"Gene"`.
#' @param treatment_col Name of the column containing the treatment label in
#'   both result tibbles. Used to auto-generate axis labels and DE class names.
#'   Default is `"Treatment"`.
#' @param region_col Name of the column containing the brain region label.
#'   Used to set the default plot title. Default is `"BrainRegion"`.
#' @param colors A named character vector mapping each DE class to a colour.
#'   Names must be `"DE in both"`, `"DE only <t1>"`, and `"DE only <t2>"`,
#'   where `<t1>` / `<t2>` are the treatment names found in the data.
#'   If `NULL` (default), a colourblind-friendly palette is used.
#' @param point_size Size of the points. Default is `3.5`.
#' @param point_alpha Opacity of the coloured (significant) points.
#'   Default is `0.7`.
#' @param max_overlaps Passed to [ggrepel::geom_text_repel()]. Controls how
#'   many overlapping labels are allowed before they are hidden. Default
#'   is `20`.
#' @param title Plot title. If `NULL` (default), the brain region name
#'   extracted from `de_result_1` is used.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res_pb  <- ptrap_de(counts_mat, sample_df, gene_ids,
#'                     region_name = "POA", treatment_name = "pb")
#' res_sol <- ptrap_de(counts_mat, sample_df, gene_ids,
#'                     region_name = "POA", treatment_name = "sol")
#'
#' # Default plot — treatment 1 on y-axis, treatment 2 on x-axis
#' ptrap_volcano2(res_pb, res_sol)
#'
#' # Custom thresholds, title and colours
#' ptrap_volcano2(res_pb, res_sol,
#'                lfc_threshold = 0.5,
#'                fdr_threshold = 0.1,
#'                title         = "POA: pb vs sol",
#'                colors        = c("DE in both"  = "#D55E00",
#'                                  "DE only pb"  = "#0072B2",
#'                                  "DE only sol" = "#49a15f"))
#' }
#'
#' @importFrom dplyr select left_join mutate case_when arrange filter all_of
#' @importFrom rlang .data :=
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_vline geom_hline
#'   scale_fill_manual theme_classic labs
#' @importFrom ggrepel geom_text_repel

ptrap_volcano2 <- function(
  de_result_1,
  de_result_2,
  lfc_threshold = 1,
  fdr_threshold = 0.05,
  gene_col = "Gene",
  treatment_col = "Treatment",
  region_col = "BrainRegion",
  colors = NULL,
  point_size = 3.5,
  point_alpha = 0.7,
  max_overlaps = 20,
  title = NULL
) {
  # --- extract treatment and region names from the result tables --------------
  t1 <- unique(de_result_1[[treatment_col]])
  t2 <- unique(de_result_2[[treatment_col]])
  region <- unique(de_result_1[[region_col]])

  # --- DE class labels (derived from treatment names) -------------------------
  class_both <- "DE in both"
  class_only1 <- paste0("DE only ", t1)
  class_only2 <- paste0("DE only ", t2)
  class_none <- "Not DE"

  # --- default colourblind-friendly palette -----------------------------------
  if (is.null(colors)) {
    colors <- c("#D55E00", "#0072B2", "#009E73")
    names(colors) <- c(class_both, class_only1, class_only2)
  }

  # --- join the two results by gene -------------------------------------------
  combined <- de_result_1 |>
    select(all_of(gene_col), logFC_1 = "logFC", padj_1 = "FDR") |>
    left_join(
      de_result_2 |>
        select(all_of(gene_col), logFC_2 = "logFC", padj_2 = "FDR"),
      by = gene_col
    )

  # --- drop genes missing a logFC in either condition ------------------------
  # This happens when filterByExpr removes a gene in one condition but not the
  # other, leaving NAs after the left_join. Those genes cannot be placed on a
  # 2D scatter and are excluded from the plot with an informative message.
  n_missing <- sum(is.na(combined$logFC_1) | is.na(combined$logFC_2))
  if (n_missing > 0) {
    message(
      n_missing,
      " gene(s) were dropped because they were only tested in one ",
      "condition (likely filtered out by filterByExpr in the other). ",
      "They cannot be placed on a 2D scatter plot."
    )
    combined <- combined |>
      filter(!is.na(.data$logFC_1), !is.na(.data$logFC_2))
  }

  # --- classify genes based on significance in each condition -----------------
  combined <- combined |>
    mutate(
      sig_1 = !is.na(.data$padj_1) &
        .data$padj_1 < fdr_threshold &
        abs(.data$logFC_1) > lfc_threshold,
      sig_2 = !is.na(.data$padj_2) &
        .data$padj_2 < fdr_threshold &
        abs(.data$logFC_2) > lfc_threshold
    ) |>
    mutate(
      DE_class = case_when(
        .data$sig_1 &  .data$sig_2  ~ class_both,
        .data$sig_1 & !.data$sig_2  ~ class_only1,
        !.data$sig_1 & .data$sig_2  ~ class_only2,
        TRUE ~ class_none
      )
    ) |>
    arrange(.data$DE_class)

  # --- default title ----------------------------------------------------------
  if (is.null(title)) {
    title <- region
  }

  # --- build plot -------------------------------------------------------------
  ggplot(combined, aes(x = .data$logFC_2, y = .data$logFC_1)) +
    geom_point(
      data = combined[combined$DE_class == class_none, ],
      color = "grey80",
      alpha = 0.5,
      size = point_size
    ) +
    geom_point(
      data = combined[combined$DE_class != class_none, ],
      aes(fill = .data$DE_class),
      alpha = point_alpha,
      size = point_size,
      shape = 21
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    geom_vline(
      xintercept = c(-lfc_threshold, lfc_threshold),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = c(-lfc_threshold, lfc_threshold),
      linetype = "dashed"
    ) +
    scale_fill_manual(values = colors) +
    geom_text_repel(
      data = combined[combined$DE_class != class_none, ],
      aes(label = .data[[gene_col]]),
      size = 3,
      max.overlaps = max_overlaps
    ) +
    theme_classic() +
    labs(
      x = paste0("logFC IP/Input (", t2, ")"),
      y = paste0("logFC IP/Input (", t1, ")"),
      title = title,
      fill = "DE class"
    )
}
