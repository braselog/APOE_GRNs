library(Seurat)
library(Matrix)

c <- readRDS(sprintf("../%s/data/%s_%s.rds", snakemake@params[["cellType"]], snakemake@params[["cellType"]], snakemake@params[["study"]]))

if (snakemake@params[["study"]] == "discovery") { # discovery data
    geneFilt <- rownames(c)[rowSums(GetAssayData(c, slot = "counts")) >= ncol(c) * 0.05]
    geneFilt <- unique(c(geneFilt, VariableFeatures(c))) #### This was not in the original code ####
    c$Clusters <- Idents(c)
} else { # study data
    geneFilt <- read.csv(snakemake@input[["geneFilt"]], head = F)$V1
    c$Clusters <- c$newClusters
}

csub <- subset(c, subset = Clusters == snakemake@params[["cellState"]])
csub <- GetAssayData(csub, slot = "counts")
csub <- csub[rownames(csub) %in% geneFilt, ]

write.table(csub, sep = ",", quote = F, file = snakemake@output[["subset"]])

# check if file exists
geneFilt <- sprintf("../%s/data/%s_discovery_genes.csv", snakemake@params[["cellType"]], snakemake@params[["cellType"]])
if (snakemake@params[["study"]] == "discovery" && !file.exists(geneFilt)) {
    write.table(rownames(csub), file = geneFilt, quote = F, row.names = F, col.names = F)
}
