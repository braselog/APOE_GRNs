write.table(rownames(read.csv(snakemake@input[["subset"]])), quote = F, row.names = F, col.names = F, sep = ",", file = snakemake@output[["geneList"]])
