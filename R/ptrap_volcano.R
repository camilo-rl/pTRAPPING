#' Single volcano plot for TRAP-seq or RNA-seq differential expression results
#'
#' Takes a single tibble produced by [pTRAPPING::ptrap_de()], [pTRAPPING::limma_de()]
#' or similar and returns a classic volcano plot (logFC on the x-axis, -log10(FDR) on the y-axis).
#' Non-significant genes are shown in grey; genes classified as `"UP"` or
#' `"DOWN"` in the `diffexpressed` column are highlighted with distinct fill
#' colours and labelled with their gene identifier. Dashed threshold lines are
#' drawn at ±`lfc_threshold` (vertical) and `-log10(fdr_threshold)`
#' (horizontal).
#'
#' @param de_result A tibble returned by [pTRAPPING::ptrap_de()], containing at
#'   minimum columns `Gene` (or `gene_col`), `logFC`, `FDR`, and
#'   `diffexpressed`.
#' @param lfc_threshold Minimum absolute log2 fold change used to draw the
#'   vertical threshold lines. Should match the value used in
#'   [pTRAPPING::ptrap_de()]. Default is `1`.
#' @param fdr_threshold FDR cutoff used to draw the horizontal threshold line
#'   (plotted as `-log10(fdr_threshold)`). Should match the value used in
#'   [pTRAPPING::ptrap_de()]. Default is `0.05`.
#' @param gene_col Name of the column containing gene identifiers. Default is
#'   `"Gene"`.
#' @param treatment_col Name of the column containing the treatment label. Used
#'   to build the default plot title. Default is `"Treatment"`.
#' @param region_col Name of the column containing the brain region label. Used
#'   to build the default plot title. Default is `"BrainRegion"`.
#' @param colors A named character vector mapping `"UP"` and `"DOWN"` to
#'   colours. If `NULL` (default), a colourblind-friendly palette is used
#'   (`"#D55E00"` for `"UP"`, `"#0072B2"` for `"DOWN"`).
#' @param point_size Size of the points. Default is `3.5`.
#' @param point_alpha Opacity of the highlighted (significant) points.
#'   Default is `0.7`.
#' @param max_overlaps Passed to [ggrepel::geom_text_repel()]. Controls how
#'   many overlapping labels are suppressed. Default is `20`.
#' @param title Plot title. If `NULL` (default), a title is auto-generated from
#'   the brain region and treatment labels found in `de_result`.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   region_name    = "POA",
#'   treatment_name = "pb"
#' )
#'
#' # Default plot
#' ptrap_volcano(res)
#'
#' # Custom thresholds and colours
#' ptrap_volcano(res,
#'               lfc_threshold = 0.5,
#'               fdr_threshold = 0.1,
#'               colors        = c("UP" = "#E69F00", "DOWN" = "#56B4E9"))
#' }
#'
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline
#'   scale_fill_manual theme_classic labs
#' @importFrom ggrepel geom_text_repel

ptrap_volcano <- function(
  de_result,
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
  # --- extract labels from the result table ------------------------------------
  treatment <- unique(de_result[[treatment_col]])
  region <- unique(de_result[[region_col]])

  # --- default colourblind-friendly palette ------------------------------------
  if (is.null(colors)) {
    colors <- c("UP" = "#D55E00", "DOWN" = "#0072B2")
  }

  # --- tag non-significant genes -----------------------------------------------
  # diffexpressed already holds "UP", "DOWN", or "NO" from ptrap_de();
  # we rename "NO" to "Not DE" for a cleaner legend label
  plot_data <- de_result |>
    mutate(
      DE = ifelse(.data$diffexpressed == "NO", "Not DE", .data$diffexpressed)
    )

  # --- default title -----------------------------------------------------------
  if (is.null(title)) {
    title <- paste0(region, " (", treatment, ")")
  }

  # --- build plot --------------------------------------------------------------
  ggplot(plot_data, aes(x = .data$logFC, y = -log10(.data$FDR))) +
    geom_point(
      data = plot_data[plot_data$DE == "Not DE", ],
      color = "grey80",
      alpha = 0.5,
      size = point_size
    ) +
    geom_point(
      data = plot_data[plot_data$DE != "Not DE", ],
      aes(fill = .data$DE),
      alpha = point_alpha,
      size = point_size,
      shape = 21
    ) +
    geom_vline(
      xintercept = c(-lfc_threshold, lfc_threshold),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(fdr_threshold),
      linetype = "dashed"
    ) +
    scale_fill_manual(values = colors) +
    geom_text_repel(
      data = plot_data[plot_data$DE != "Not DE", ],
      aes(label = .data[[gene_col]]),
      size = 3,
      max.overlaps = max_overlaps
    ) +
    theme_classic() +
    labs(
      x = expression(log[2] ~ "Fold Change (IP / Input)"),
      y = expression(-log[10](FDR)),
      title = title,
      fill = "DE"
    )
}
