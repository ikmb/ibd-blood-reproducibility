BeadChip microarray differential expression analysis
================

Load prerequisites that are needed for differential gene expression
analysis (DEA) in R environment.

``` r
# libraries
library(tidyverse)
library(limma)
library(illuminaHumanv4.db)
library(org.Hs.eg.db)
```

Read BeadChip microarray intensity data and respective meta information.
The raw data (.idat and bgx files), which are used in this analysis are
also deployed at GEO under
[GSE169568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169568).

``` r
# set wd to main
setwd('..')

# meta data
meta <- read_tsv(
  paste0(getwd(),"/data/se_rna_meta.tsv")
  )

# intensity data
elist <- readRDS(
  paste0(getwd(),"/data/elist.rds.gz")
  )

# add sample info
elist$samples <- meta
```

Summary of the traits in the dataset:

``` r
# summarize samples per trait
meta %>% group_by(diagnosis) %>% tally() %>%
  bind_rows(., tibble(diagnosis = "TOTAL", n = sum(.$n)))
```

    ## # A tibble: 5 × 2
    ##   diagnosis     n
    ##   <chr>     <int>
    ## 1 CD           52
    ## 2 HC           30
    ## 3 SC           65
    ## 4 UC           58
    ## 5 TOTAL       205

Probe QC step. Remove low expressed genes (detected in less than 10% of
samples), not annotated probes, and duplicated probes targeting the same
gene.

``` r
# background correction and normalization
elist_norm <- neqc(elist)

# map probes to array address
ids <- unlist(
  mget(
    rownames(elist_norm),
    revmap(illuminaHumanv4ARRAYADDRESS),
    ifnotfound = NA)
  )

# annotate probes
elist_norm$genes$Symbol <- unlist(
  mget(
    elist_norm$genes$Probe_Id,
    envir = illuminaHumanv4SYMBOL
    )
  )
elist_norm$genes$ENTREZ <- unlist(
  mget(elist_norm$genes$Probe_Id,
       envir = illuminaHumanv4ENTREZID
       )
  )

# get probe quality
qual <- unlist(
  mget(
    ids,
    illuminaHumanv4PROBEQUALITY,
    ifnotfound = NA
    )
  )

# filter bad quality probes
rem <- qual == "No match" | qual == "Bad"
elist_norm_filt <- elist_norm[!rem, ]

# filter low expressed genes
expressed <- rowSums(
  elist_norm_filt$other$Detection<0.05
  ) >= round(nrow(elist_norm_filt$samples)*0.1)
elist_norm_filt <- elist_norm_filt[expressed, ]

# remove not annotated probes (NAs)
elist_norm_filt <- elist_norm_filt[!is.na(elist_norm_filt$genes$Symbol),]

## remove duplicated probes 
## modified from genefilter::findLargest()
rownames(elist_norm_filt) <- elist_norm_filt$genes$Probe_Id

# probe ids
gN <- intersect(
  ls(illuminaHumanv4SYMBOL),
  elist_norm_filt$genes$Probe_Id
  )

# symbol names
sN <- unlist(
  mget(
    gN,
    illuminaHumanv4SYMBOL
    ),
  use.names = FALSE
  )

# order by probe
elist_norm_filt_ord <- elist_norm_filt[gN,]

# calculate medians
testStat <- rowMedians(elist_norm_filt_ord$E)
names(testStat) <- gN

# split probes by symbol
tSsp <- split.default(testStat, sN)

# get index of highest median
idx <- sapply(tSsp, function(x) names(which.max(x)))

# remove duplicates
elist_qced <- elist_norm_filt_ord[idx,]

# number of probes and samples
dim(elist_qced$E)
```

    ## [1] 11909   205

Prior to differential gene expression analysis, we also remove ribosomal
genes (those beginning with Mrps, Rpl, and Rps).

``` r
# remove ribosomal RPS, MRP and RPL transcripts
elist_qced <- elist_qced[!grepl("^RPS|^RPL|^MRP", elist_qced$genes$Symbol),]


# write values for wgcna
exprs <- elist_qced$E
rownames(exprs) <- elist_qced$genes$Symbol

setwd('..')
write_tsv(
  exprs %>% as_tibble(rownames = "gene"),
  paste0(getwd(), "/data/se_rna_exprs.tsv")
  )

# final number of probes
dim(elist_qced)
```

    ## [1] 11725   205

Differential gene expression analysis using limma’s linear models, which
include age as a covariate.

``` r
### differential analysis

# definition of model design
design <- 
  model.matrix(
    ~0+diagnosis+age+gender,
    data = elist_qced$samples
  )

# weights based on quality
aw <- arrayWeights(
  elist_qced,
  design
  )

# contrast matrix
cnts_matrix <-
  makeContrasts(
    "CD_HC"=diagnosisCD-diagnosisHC,
    "UC_HC"=diagnosisUC-diagnosisHC,
    "CD_UC"=diagnosisCD-diagnosisUC,
    "SC_HC"=diagnosisSC-diagnosisHC,
    "CD_SC"=diagnosisCD-diagnosisSC,
    "UC_SC"=diagnosisUC-diagnosisSC,
    levels=design
    )

# fit to a linear model
fit <- lmFit(
  elist_qced,
  design,
  weights = aw
  )
contr_fit <- eBayes(
  contrasts.fit(
    fit,
    cnts_matrix
    )
  )
```

Results of treatment-naive Swedish dataset (Supplementary Table 3, 1st
sheet).

``` r
# get results
results <- lapply(seq(colnames(cnts_matrix)), function(i){
  topTable(contr_fit,
           coef=colnames(cnts_matrix)[i],
           number=nrow(contr_fit),
           sort.by="logFC", confint=T) %>%
    mutate(
      comparison=gsub("_", " vs ", colnames(cnts_matrix)[i]),
      dea = case_when(
        adj.P.Val<0.05 & abs(logFC)>log2(1.5) ~ "de",
        TRUE ~ "non-de")
      )
}) %>% bind_rows()

# write results
setwd('..')
write_tsv(
  results,
  paste0(getwd(), "/results/se_rna_dea_results.tsv")
  )

# print results
head(results)
```

    ##                      Probe_Id Array_Address_Id  Symbol ENTREZ    logFC     CI.L
    ## ILMN_1651358...1 ILMN_1651358          4830541    HBE1   3046 2.605011 2.186610
    ## ILMN_1694548...2 ILMN_1694548           240594   ANXA3    306 2.063832 1.652707
    ## ILMN_1664706...3 ILMN_1664706          4570725   H3C13 653604 1.972432 1.663019
    ## ILMN_2387385...4 ILMN_2387385          3290368  IGFBP1   3484 1.866047 1.488465
    ## ILMN_1748915...5 ILMN_1748915          1410221 S100A12   6283 1.815149 1.466480
    ## ILMN_1724533...6 ILMN_1724533            70167    LY96  23643 1.783370 1.430133
    ##                      CI.R   AveExpr         t      P.Value    adj.P.Val
    ## ILMN_1651358...1 3.023412 11.607552 12.275777 2.668256e-26 1.042843e-22
    ## ILMN_1694548...2 2.474958  8.946410  9.897642 4.105999e-19 1.719387e-16
    ## ILMN_1664706...3 2.281844  7.247558 12.568873 3.313261e-27 1.942399e-23
    ## ILMN_2387385...4 2.243628  8.392602  9.744145 1.154583e-18 4.081533e-16
    ## ILMN_1748915...5 2.163817 11.098550 10.264355 3.400537e-20 2.491956e-17
    ## ILMN_1724533...6 2.136607  9.143565  9.954222 2.801047e-19 1.263164e-16
    ##                         B comparison dea
    ## ILMN_1651358...1 48.98837   CD vs HC  de
    ## ILMN_1694548...2 32.80451   CD vs HC  de
    ## ILMN_1664706...3 51.02815   CD vs HC  de
    ## ILMN_2387385...4 31.79392   CD vs HC  de
    ## ILMN_1748915...5 35.23997   CD vs HC  de
    ## ILMN_1724533...6 33.17838   CD vs HC  de

``` r
sessionInfo()
```

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] illuminaHumanv4.db_1.26.0 org.Hs.eg.db_3.13.0      
    ##  [3] AnnotationDbi_1.54.1      IRanges_2.26.0           
    ##  [5] S4Vectors_0.30.2          Biobase_2.52.0           
    ##  [7] BiocGenerics_0.38.0       limma_3.48.3             
    ##  [9] forcats_0.5.1             stringr_1.4.0            
    ## [11] dplyr_1.0.7               purrr_0.3.4              
    ## [13] readr_2.0.2               tidyr_1.1.4              
    ## [15] tibble_3.1.5              ggplot2_3.3.5            
    ## [17] tidyverse_1.3.1          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.2             vroom_1.5.5            bit64_4.0.5           
    ##  [4] jsonlite_1.7.2         modelr_0.1.8           assertthat_0.2.1      
    ##  [7] blob_1.2.2             GenomeInfoDbData_1.2.6 cellranger_1.1.0      
    ## [10] yaml_2.2.1             pillar_1.6.3           RSQLite_2.2.8         
    ## [13] backports_1.2.1        glue_1.4.2             digest_0.6.28         
    ## [16] XVector_0.32.0         rvest_1.0.1            colorspace_2.0-2      
    ## [19] htmltools_0.5.2        pkgconfig_2.0.3        broom_0.7.9           
    ## [22] haven_2.4.3            zlibbioc_1.38.0        scales_1.1.1          
    ## [25] tzdb_0.1.2             KEGGREST_1.32.0        generics_0.1.0        
    ## [28] ellipsis_0.3.2         cachem_1.0.6           withr_2.4.2           
    ## [31] cli_3.0.1              magrittr_2.0.1         crayon_1.4.1          
    ## [34] readxl_1.3.1           memoise_2.0.0          evaluate_0.14         
    ## [37] fs_1.5.0               fansi_0.5.0            xml2_1.3.2            
    ## [40] tools_4.1.1            hms_1.1.1              lifecycle_1.0.1       
    ## [43] munsell_0.5.0          reprex_2.0.1           Biostrings_2.60.2     
    ## [46] compiler_4.1.1         GenomeInfoDb_1.28.4    rlang_0.4.11          
    ## [49] grid_4.1.1             RCurl_1.98-1.5         rstudioapi_0.13       
    ## [52] bitops_1.0-7           rmarkdown_2.11         gtable_0.3.0          
    ## [55] DBI_1.1.1              R6_2.5.1               lubridate_1.8.0       
    ## [58] knitr_1.36             fastmap_1.1.0          bit_4.0.4             
    ## [61] utf8_1.2.2             stringi_1.7.5          Rcpp_1.0.7            
    ## [64] vctrs_0.3.8            png_0.1-7              dbplyr_2.1.1          
    ## [67] tidyselect_1.1.1       xfun_0.26
