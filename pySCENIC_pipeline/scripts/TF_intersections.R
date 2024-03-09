# LINK /home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig5a_e.R
library(dplyr)
library(ggplot2)
library(ComplexUpset)

theme_set(theme_minimal())
theme_replace(axis.text = element_text(color = "black"))

split_and_get_first_element <- function(x) {
    strsplit(x, split = "\\(")[[1]][1]
}

read_and_process_csv <- function(file_path, study) {
    df <- read.csv(file_path)
    if (study == "discovery") {
        threshold <- 80
        df$TF <- df$V1
        df$sum <- df$V2
    } else {
        threshold <- 50
    }
    df <- df[df$sum >= threshold, ]
    sapply(df$TF, split_and_get_first_element)
}

m0 <- intersect(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic0_newTF_TF_Gene_sum_above80.csv", "discovery"),
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/mic0_colonna_newTF_100run_onehot_20231030C.csv", "replication")
)

m1 <- intersect(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic1_newTF_TF_Gene_sum_above80.csv", "discovery"),
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/mic1_colonna_newTF_100run_onehot_20231030C.csv", "replication")
)

m3 <- intersect(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic3_newTF_TF_Gene_sum_above80.csv", "discovery"),
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/mic3_colonna_newTF_100run_onehot_20231031C.csv", "replication")
)

a0 <- intersect(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/astro0_newTF_TF_Gene_sum_above80.csv", "discovery"),
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/astro0_colonna_newTF_100run_onehot_20231023C.csv", "replication")
)

a1 <- intersect(
    read_and_process_csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/astro1_newTF_TF_Gene_sum_above80.csv", "discovery"),
    read_and_process_csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/astro1_colonna_newTF_100run_onehot_20231023C.csv", "replication")
)
# Step 1: Create a list of all the groups
groups <- list("mic-resting" = m0, "mic-activated" = m1, "mic-proinflam" = m3, "astro-resting" = a0, "astro-activated" = a1)

# Step 2: Find the union of all elements in these groups
all_elements <- Reduce(union, groups)

# Step 3: Initialize an empty dataframe with rows as the group names and columns as the union of elements
df <- data.frame(matrix(ncol = length(all_elements), nrow = length(groups)))
colnames(df) <- all_elements
rownames(df) <- names(groups)

# Step 4: Iterate over each group and each element in the union of elements
for (group_name in names(groups)) {
    for (element in all_elements) {
        df[group_name, element] <- element %in% groups[[group_name]]
    }
}

df <- as.data.frame(t(df))

# upset(fromList(TFs), order.by = "degree", main.bar.color = "black", matrix.color = "black", text.scale = c(1, 1, 1, 1, 2, 0.75))
p <- upset(df, colnames(df), name = "Cell state", width_ratio = 0.1, mode = "inclusive_intersection", sort_intersections_by = c("degree", "cardinality"), queries = list(upset_query(intersect = c("astro-resting", "astro-activated"), color = "blue", fill = "blue")))
ggsave("additionalAnalyses/plots/pySCENIC_TF_intersection_upsetplot.pdf")
write.csv(df, "additionalAnalyses/data/pySCENIC_TF_intersection_upsetplot.csv")

# gwas <- read.csv("/home/brasel/GWAS_loci_genes/data/2022.08.19_GWAS_loci_genes_inData.csv")
gwas <- read.csv("/home/brasel/GWAS_loci_genes/finalList/GWAS_locus_gene_list.tsv", sep = "\t")

common <- rownames(df[rowSums(df) == 5, ])
#  [1] "ARID3A"  "BCLAF1"  "BHLHE40" "BPTF"    "CEBPB"   "CEBPD"   "CEBPG"   "ELF1"    "ELK3"    "ELK4"    "EP300"   "FOS"
# [13] "FOSL2"   "FOXO1"   "GABPA"   "GATAD1"  "IKZF2"   "JUNB"    "JUND"    "KLF12"   "KMT2B"   "MZF1"    "NFIA"    "NFIC"
# [25] "NR3C1"   "SP3"     "STAT1"   "YY1"     "MAX"     "NFE2L1"  "SREBF2"  "PURA"    "RXRA"    "MAF"     "TAF1"
common[common %in% gwas$Gene]
gwas[gwas$Gene %in% common[common %in% gwas$Gene], ]
#     Locus   Gene
# 2120  WWOX    MAF
# 2140   MAF    MAF
# 2692 ABCA7 ARID3A
# 3287   APP  GABPA

common <- rownames(df[df$"mic-resting" == T & df$"mic-activated" == T & df$"mic-proinflam" == T, ])
common[common %in% gwas$Gene]
gwas[gwas$Gene %in% common[common %in% gwas$Gene], ]
#      Locus   Gene
# 827   IKZF1  IKZF1
# 1363   SPI1   SPI1
# 1730 SPPL2A GABPB1
# 2120   WWOX    MAF
# 2140    MAF    MAF
# 2692  ABCA7 ARID3A
# 2832  KLF16   TCF3
# 3287    APP  GABPA

common <- rownames(df[df$"astro-resting" == T & df$"astro-activated" == T, ])
common[common %in% gwas$Gene]
gwas[gwas$Gene %in% common[common %in% gwas$Gene], ]
#      Locus    Gene
# 1296   BLNK ZNF518A
# 2120   WWOX     MAF
# 2140    MAF     MAF
# 2402 MYO15A  SREBF1
# 2692  ABCA7  ARID3A
# 3287    APP   GABPA
