m0 <- intersect(
    unlist(lapply(
        read.csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic0_newTF_TF_Gene_sum_above80.csv")$V1,
        function(x) unlist(strsplit(x, "\\("))[1]
    )),
    unlist(lapply(
        read.csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/mic0_colonna_newTF_TF_Gene_sum_above80.csv")$V1,
        function(x) unlist(strsplit(x, "\\("))[1]
    ))
)
m1 <- intersect(
    unlist(lapply(
        read.csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic1_newTF_TF_Gene_sum_above80.csv")$V1,
        function(x) unlist(strsplit(x, "\\("))[1]
    )),
    unlist(lapply(
        read.csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/mic1_colonna_newTF_TF_Gene_sum_above80.csv")$V1,
        function(x) unlist(strsplit(x, "\\("))[1]
    ))
)
m3 <- intersect(
    unlist(lapply(
        read.csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic3_newTF_TF_Gene_sum_above80.csv")$V1,
        function(x) unlist(strsplit(x, "\\("))[1]
    )),
    unlist(lapply(
        read.csv("/home/yanboyu/1PyScenicProject/duplication/ResultedData/mic3_colonna_newTF_TF_Gene_sum_above80_20231031.csv")$V1,
        function(x) unlist(strsplit(x, "\\("))[1]
    ))
)
TFs = list("mic0" = m0, "mic1" = m1, "mic3" = m3)
# TFs = list( "astro0" = unlist(lapply( read.csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/astro0_newTF_TF_Gene_sum_above80.csv")$V1, function(x) unlist(strsplit(x, "\\("))[1])), "astro1" = unlist(lapply( read.csv("/home/yanboyu/1PyScenicProject/discovery/ResultedData/astro1_newTF_TF_Gene_sum_above80.csv")$V1, function(x) unlist(strsplit(x, "\\("))[1])))

# [ins] r$> intersect(intersect(m0,m1),m3)
#  [1] "CEBPB"  "CEBPD"  "CEBPG"  "ELF1"   "ELK3"   "ELK4"   "EP300"  "ETS1"   "ETV5"
# [10] "FLI1"   "FOS"    "FOXO1"  "FOXO3"  "GABPA"  "IKZF1"  "IKZF2"  "ILF2"   "JUND"
# [19] "KMT2B"  "MZF1"   "NFIA"   "NFIC"   "NR3C1"  "PBX1"   "POU2F1" "SP3"    "SPI1"
# [28] "STAT1"  "TBP"    "TCF7L2" "YY1"    "ZNF148" "ETS2"   "GABPB1" "MAX"    "NFE2L1"
# [37] "SREBF2" "TCF3"   "PURA"   "RXRA"   "ELF2"   "MAF"    "TAF1"   "ZEB1"
#
# [ins] r$> m0[!(m0 %in% intersect(intersect(m0,m1),m3))]
#  [1] "ARID3A" "BACH2"  "BPTF"   "CEBPA"  "E2F3"   "ELF4"   "ETV6"   "FOSL2"  "ING4"
# [10] "IRF2"   "IRF7"   "KLF12"  "NFATC2" "NR2C2"  "RUNX3"  "TAL1"   "XBP1"   "ZFX"
# [19] "PLAG1"  "SP4"    "MBD2"   "FOXJ2"  "ATF1"   "MXD1"   "ZNF333" "STAT6"  "TFE3"
# [28] "IRF5"   "IRF3"
#
# [ins] r$> m1[!(m1 %in% intersect(intersect(m0,m1),m3))]
#  [1] "BACH2"  "CHD1"   "E2F3"   "ING4"   "KLF13"  "MECP2"  "NR2C2"  "RFX1"   "SMAD3"
# [10] "STAT2"  "TAL1"   "XBP1"   "ZBTB20" "ZBTB7A" "SP1"    "CLOCK"  "IRF8"   "HLTF"
# [19] "JUNB"   "USF1"   "ZNF333" "PPARG"  "STAT6"  "BCLAF1" "E2F4"   "VEZF1"  "SREBF1"
# [28] "SP4"    "MAZ"    "IRF5"
#
# [ins] r$> m3[!(m3 %in% intersect(intersect(m0,m1),m3))]
#  [1] "BCLAF1" "BPTF"   "ELF4"   "ETV6"   "IRF7"   "JUN"    "KLF9"   "RFX2"   "SP1"
# [10] "CEBPA"  "PPARG"  "ARID3A" "E2F4"   "GATAD1" "RUNX3"  "SREBF1" "ZBTB7A" "CLOCK"
# [19] "MAFG"   "MAZ"    "ZFY"    "JUNB"   "NFE2L2" "FOSL2"  "VEZF1"  "THRB"   "CUX1"

paths <- list("/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic0_newTF_onehot", "/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic1_newTF_onehot", "/home/yanboyu/1PyScenicProject/discovery/ResultedData/mic3_newTF_onehot")
names(paths) <- c("mic0", "mic1", "mic3")

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5800817/
gene <- c("APOE")
# m <- read.csv("/home/brasel/SingleCellProjects/dataObjects/mic_01.csv")
# gene = m$X

byGene <- lapply(gene, function(g) {
    # clust = 'mic0'
    results <- lapply(names(paths), function(clust) {
        files <- list.files(paths[[clust]], pattern = "csv", full.names = TRUE)
        tfs <- TFs[[clust]]
        tfs <- tfs[order(tfs)]
        files <- files[grep(paste0(tfs, "\\(", collapse = "|"), files)]

        hits <- lapply(files, function(f) {
            df <- read.csv(f, header = TRUE, row.names = 1)
            if (g %in% rownames(df)) {
                return(df[g, "sum"])
            } else {
                return(NA)
            }
        })
        names(hits) <- tfs
        hits <- na.omit(unlist(hits))
        return(hits)
    })
    names(results) <- names(paths)
    # Ensure all lists have the same length
    results <- lapply(results, function(r) r[1:max(lengths(results))])
    return(results)
})
names(byGene) <- gene



# $APOE
# $APOE$mic0
#     CEBPA     CEBPB     CEBPD     FOSL2     GABPA      ILF2      IRF7      JUND    NFE2L1
# 23.000000 87.000000  6.000000  2.000000  4.000000  1.000000  6.000000 28.000000  1.010101
#     STAT6      XBP1
#  1.136364  1.000000
#
# $APOE$mic1
#     CEBPB     CEBPD       FOS      JUNB      JUND      SPI1    ZBTB7A      <NA>      <NA>
# 20.000000 17.000000  1.000000 44.329897  1.010101  1.000000  4.000000        NA        NA
#      <NA>      <NA>
#        NA        NA
#
# $APOE$mic3
#     CEBPB     CEBPD      E2F4      JUNB      MAFG     NR3C1      <NA>      <NA>      <NA>
# 78.000000  6.000000  1.020408  5.263158  1.030928  1.000000        NA        NA        NA
#      <NA>      <NA>
#        NA        NA
