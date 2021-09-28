# ibd-blood-reproducibility
scripts to reproduce main analyses of IBD blood transcriptomes

The analysis can be done as following:

1)  Differential expression analysis (DEA) in treatment-naive IBD patients using Illumina's BeadChip array (run:scripts/1_microarray_dea.Rmd);
2)  Differential expression analysis (DEA) in treatment-naive and treatment-receiving IBD patients using small RNA-seq (run:scripts/2_microrna_dea.Rmd);
3)  Gene set enrichment analysis (GSEA) of differentially expressed genes and validated targets of differentially expressed microRNAs (run:scripts/3_gsea_dge.Rmd);
4)  Generation of gene co-expression networks (TOMs) for each trait using weighted gene correlation analysis (WGCNA) and treatment-naive IBD patients (run:scripts/4_microarray_wgcna.Rmd);
5)  Nonnegative tensor decomposition of gene co-expression networks (TOMs) using CP heirarchical alternating least squares method (run:scripts/5_wgcna_tensor_decomp.ipynb). The method includes random state and therefore, the results might slightly differ, however they should not change drastically, i.e. the main co-expression components will be detected although the order might be different.
6)  Functional annotation of decomposed co-expression components using gene set enrichment analysis (GSEA) (run:scripts/6_component_annotation.Rmd)