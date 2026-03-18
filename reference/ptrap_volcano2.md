# Dual volcano plot comparing two treatment conditions from TRAP-seq DE results

Takes two tibbles produced by
[`ptrap_de()`](https://camilo-rl.github.io/pTRAPPING/reference/ptrap_de.md)
— one per treatment condition, both from the same brain region — joins
them by gene, classifies each gene according to its differential
expression status in each condition, and returns a scatter plot of
logFC\\\_{\text{treatment 1}}\\ vs logFC\\\_{\text{treatment 2}}\\.
Significant genes are highlighted and labelled; threshold lines are
drawn at ±`lfc_threshold` on both axes.

## Usage

``` r
ptrap_volcano2(
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
)
```

## Arguments

- de_result_1:

  A tibble returned by
  [`ptrap_de()`](https://camilo-rl.github.io/pTRAPPING/reference/ptrap_de.md)
  for the **first** treatment condition (plotted on the **y-axis**).

- de_result_2:

  A tibble returned by
  [`ptrap_de()`](https://camilo-rl.github.io/pTRAPPING/reference/ptrap_de.md)
  for the **second** treatment condition (plotted on the **x-axis**).

- lfc_threshold:

  Minimum absolute log2 fold change used to define significance. Must
  match (or be stricter than) the threshold used in
  [`ptrap_de()`](https://camilo-rl.github.io/pTRAPPING/reference/ptrap_de.md).
  Default is `1`.

- fdr_threshold:

  Maximum FDR used to define significance. Default is `0.05`.

- gene_col:

  Name of the column containing gene identifiers in both result tibbles.
  Default is `"Gene"`.

- treatment_col:

  Name of the column containing the treatment label in both result
  tibbles. Used to auto-generate axis labels and DE class names. Default
  is `"Treatment"`.

- region_col:

  Name of the column containing the brain region label. Used to set the
  default plot title. Default is `"BrainRegion"`.

- colors:

  A named character vector mapping each DE class to a colour. Names must
  be `"DE in both"`, `"DE only <t1>"`, and `"DE only <t2>"`, where
  `<t1>` / `<t2>` are the treatment names found in the data. If `NULL`
  (default), a colourblind-friendly palette is used.

- point_size:

  Size of the points. Default is `3.5`.

- point_alpha:

  Opacity of the coloured (significant) points. Default is `0.7`.

- max_overlaps:

  Passed to
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).
  Controls how many overlapping labels are allowed before they are
  hidden. Default is `20`.

- title:

  Plot title. If `NULL` (default), the brain region name extracted from
  `de_result_1` is used.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
if (FALSE) { # \dontrun{
res_pb  <- ptrap_de(counts_mat, sample_df, gene_ids,
                    region_name = "POA", treatment_name = "pb")
res_sol <- ptrap_de(counts_mat, sample_df, gene_ids,
                    region_name = "POA", treatment_name = "sol")

# Default plot — treatment 1 on y-axis, treatment 2 on x-axis
ptrap_volcano2(res_pb, res_sol)

# Custom thresholds, title and colours
ptrap_volcano2(res_pb, res_sol,
               lfc_threshold = 0.5,
               fdr_threshold = 0.1,
               title         = "POA: pb vs sol",
               colors        = c("DE in both"  = "#D55E00",
                                 "DE only pb"  = "#0072B2",
                                 "DE only sol" = "#49a15f"))
} # }
```
