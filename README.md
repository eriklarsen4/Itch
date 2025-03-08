# Itch
R code for bioinformatics analyses used in the [Tmem184b paper published in Pain](https://pmc.ncbi.nlm.nih.gov/articles/PMC8854445/#S16)

Markdown files encompass the analysis of various Illumina (paired-end and single-end) bulk RNAseq experiments

Source code has been revised since their original publications to facilitate better readability and reproducibility

+ includes standard, high-level bioinformatics/transcriptomics visualizations, including:

  + volcano plots
  + MA plots
  + heatmaps

+ includes higher-level analyses on subsets of the data:

  + gene ontology querying and statistical analysis
    + visualizations as bar plots
   
+ includes wrangling and visualization of time-series fluorescence microscopy

At the time, raw sequencing files were pseudoaligned with Salmon on the UA HPC.
The aligned reads were analyzed for differential expression on UseGalaxy.org, employing the DESeq2 algorithm.

Tutorials for transcriptional analysis pipelines are in development.

Included in these and analyses in other repos are tools I've developed to centralize downstream bioinformatic data regarding gene or protein lists.
These tools will soon be included in/as packages and docker images.
