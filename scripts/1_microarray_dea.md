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

    ## # A tibble: 5 x 2
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

    ## [1] 11911   205

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

    ## [1] 11727   205

Differential gene expression analysis using limma’s linear models, which
include age as a covariate.

``` r
### differential analysis

# definition of model design
design <- 
  model.matrix(
    ~0+diagnosis+age,
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
  paste0(getwd(), "/results/se_mirna_dea_results.tsv")
  )

# print results
head(results)
```

    ##       Probe_Id Array_Address_Id  Symbol ENTREZ    logFC     CI.L     CI.R
    ## 1 ILMN_1651358          4830541    HBE1   3046 2.622609 2.206845 3.038374
    ## 2 ILMN_1694548           240594   ANXA3    306 2.018144 1.604871 2.431417
    ## 3 ILMN_1664706          4570725   H3C13 653604 1.946959 1.638636 2.255282
    ## 4 ILMN_2387385          3290368  IGFBP1   3484 1.892989 1.516345 2.269632
    ## 5 ILMN_1748915          1410221 S100A12   6283 1.807271 1.461023 2.153519
    ## 6 ILMN_1724533            70167    LY96  23643 1.765111 1.411113 2.119110
    ##     AveExpr         t      P.Value    adj.P.Val        B comparison dea
    ## 1 11.607552 12.436725 7.967173e-27 3.114368e-23 50.15212   CD vs HC  de
    ## 2  8.946410  9.627982 2.440829e-18 9.870207e-16 31.05893   CD vs HC  de
    ## 3  7.247558 12.450054 7.244552e-27 3.114368e-23 50.24501   CD vs HC  de
    ## 4  8.392602  9.909174 3.675388e-19 1.959149e-16 32.90813   CD vs HC  de
    ## 5 11.098550 10.290972 2.732240e-20 3.204098e-17 35.44738   CD vs HC  de
    ## 6  9.143565  9.830848 6.239534e-19 2.814270e-16 32.39116   CD vs HC  de

``` r
sessionInfo()
```

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS  10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] illuminaHumanv4.db_1.26.0 org.Hs.eg.db_3.11.4      
    ##  [3] AnnotationDbi_1.50.1      IRanges_2.22.2           
    ##  [5] S4Vectors_0.26.1          Biobase_2.48.0           
    ##  [7] BiocGenerics_0.34.0       limma_3.44.3             
    ##  [9] forcats_0.5.0             stringr_1.4.0            
    ## [11] dplyr_1.0.0               purrr_0.3.4              
    ## [13] readr_1.3.1               tidyr_1.1.3              
    ## [15] tibble_3.0.3              ggplot2_3.3.2            
    ## [17] tidyverse_1.3.0          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.6       lubridate_1.7.9  utf8_1.1.4       assertthat_0.2.1
    ##  [5] digest_0.6.25    R6_2.4.1         cellranger_1.1.0 backports_1.1.8 
    ##  [9] reprex_0.3.0     RSQLite_2.2.0    evaluate_0.14    httr_1.4.1      
    ## [13] pillar_1.4.6     rlang_0.4.11     readxl_1.3.1     rstudioapi_0.11 
    ## [17] blob_1.2.1       rmarkdown_2.3    bit_1.1-15.2     munsell_0.5.0   
    ## [21] broom_0.7.0      compiler_4.0.2   modelr_0.1.8     xfun_0.19       
    ## [25] pkgconfig_2.0.3  htmltools_0.5.0  tidyselect_1.1.0 fansi_0.4.1     
    ## [29] crayon_1.3.4     dbplyr_1.4.4     withr_2.4.2      grid_4.0.2      
    ## [33] jsonlite_1.7.2   gtable_0.3.0     lifecycle_0.2.0  DBI_1.1.0       
    ## [37] magrittr_1.5     scales_1.1.1     cli_2.0.2        stringi_1.4.6   
    ## [41] fs_1.4.2         xml2_1.3.2       ellipsis_0.3.1   generics_0.0.2  
    ## [45] vctrs_0.3.8      tools_4.0.2      bit64_0.9-7      glue_1.4.1      
    ## [49] hms_0.5.3        yaml_2.2.1       colorspace_1.4-1 rvest_0.3.6     
    ## [53] memoise_1.1.0    knitr_1.29       haven_2.3.1
