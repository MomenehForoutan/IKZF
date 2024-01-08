This repo contains Rmd report files (under the report/ folder) to perfrom exploratory data analysis of RNA-seq, ChIP-seq, and Cut&Run data of Ikzf1/3 KO samples in mouse NK cells. Some of the results included in these reports have already been published in [Nature Immunology](https://www.nature.com/articles/s41590-023-01718-4).   

Below are the list of Rmd report files and scripts:

report/**Ikzf1_multi_omics.Rmd**: this is the main script that reads in the results of several in-house RNA-seq, ChIP-seq (Ikzf3) and Cut&Run (Ikzf1) data, as well as several publicly available transcriptomics and epigenomics data sets, and attempts to explore the data using different visualisations and ana analysis aproches.

report/**Ikzf3_CnR_MACS2.Rmd**: We ran MACS2 peak calling on the Cut&Run data, and have analysed it using nextflow, manual shell scripting as well as through the Galaxy. Here, I read in the results of my manual narrow peak calling, perform DB analysis using featureCount and edgeR, and additionally, I perform a quick comparison between the results obtained from nextflow and Galaxy.

report/**Ikzf3_CnR_csaw.Rmd**: This report uses bam files and perform a window-based analysis through the csaw, rtracklayer, and edgeR packages.