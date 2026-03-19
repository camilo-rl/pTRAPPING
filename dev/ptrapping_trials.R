counts_beg <- read.csv("/Users/camilorl/Downloads/counts_beg.csv")
counts_beg2 <- counts_beg |>
  select(-X)

ptrap_de(counts_beg2, treatment_name = "nb")

library(readr)
library(dplyr)
counts_pairb <- read_tsv(
  "/Users/camilorl/Quarto_projects/PairBond_PhTRAP/counts_Rimitator_gene.txt",
  comment = "#"
)

# remove extra annotation columns to get counts matrix
counts <- counts_pairb %>%
  select(-Chr, -Start, -End, -Strand, -Length)


ptrap_de(counts)
