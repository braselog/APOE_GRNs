APOE_GRNs
Abstract (from preprint)
Recent single-cell technologies have revealed diverse transcriptional states in multiple brain cell types. Many of these states are implicated in disease, by either driving or responding to the etiologies. Transcription factors (TFs) play a key role in modulating the expression of related sets of genes, thus regulating these different expression states. APOE, a key risk-modifying gene in Alzheimer disease (AD), is expressed in specific glial transcriptional states associated with AD. However, it is still unknown whether the regulatory programs modulating its expression are shared across brain cell types, or are specific to microglia and astrocytes, the two cell types that display highest APOE expression. To better understand this process, we leveraged single-nucleus RNA-seq to determine gene regulatory networks (GRNs) for critical resting and activated cell states in both microglia and astrocytes. Through our analysis, we identified the CEBP, JUN, FOS, and FOXO TF families as key regulators of APOE in microglia. Some of these TF families were also involved in APOE regulation in resting astrocytes, but the steroid/thyroid hormone receptor families, including the THR TF family, dominated APOE regulation in astrocytes. Furthermore, we report that additional AD GWAS associated genes are co-regulated with APOE, and implicate cytokine/inflammation, glucose metabolism, and lipid metabolism networks. Additionally, many of the APOE-regulating TFs were linked to circadian rhythm. Altogether, these results support the building body of evidence that inflammation-induced cholesterol synthesis and cholesterol transport via APOE is at the core of these AD-related pathways. 


This project is a part of my thesis work. I have modified it so it is not as large, but can still represent my work.
The raw data is protected, so cannot be included. Unfortunately, this prevents the project from being run from scratch.
In the future I will extend it to include automatic downloading of some test data so the full set of analyses can be run by anyone.

Running:
1. enter the pySCENIC_pipeline directory and run the snakemake pipeline using
snakemake --use-conda -c4

2. enter the TF_perturbation directory and run the snakemake pipeline using:
snakemake --use-conda -c4
