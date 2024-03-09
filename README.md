APOE_GRNs

We constructed gene regulatory networks for microglial and astrocytic cell states. We then virtually knocked out the transcription factors regulating APOE to visualize their impact on expression. 

This project is a part of my thesis work. I have modified it so it is not as large, but can still represent my work.
The raw data is protected, so cannot be included. Unfortunately, this prevents the project from being run from scratch.
In the future I will extend it to include automatic downloading of some test data so the full set of analyses can be run by anyone.

Running:

1. enter the pySCENIC_pipeline directory and run the snakemake pipeline using
   
snakemake --use-conda -c4

3. enter the TF_perturbation directory and run the snakemake pipeline using:

snakemake --use-conda -c4
