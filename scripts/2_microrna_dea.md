microRNA differential expression analysis
================

Load prerequisites that are needed for differential expression analysis
(DEA) in R environment.

``` r
# libraries
library(tidyverse)
library(edgeR)
library(limma)
```

Read microRNA count data and respective meta information. The raw data
(fastqs) and the counts, which are used in this analysis are deployed at
GEO under
[GSE169569](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169569)
and
[GSE169570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169570)
accession numbers. The data is QCed in terms of outlier samples, the
ouliers were removed based on IQR of expressed microRNAs in logarithmic
scale (IQR \< 1.5).

``` r
# set wd to main
setwd('..')

# count data
counts <- list.files(
  path = paste0(getwd(), "/data"),
  pattern = "counts.tsv",
  full.names = TRUE
  )
names(counts) <- str_extract(counts, "[a-z]{2}_mirna")

counts <- lapply(counts, function(x){
  read_tsv(x) %>% column_to_rownames("mirna")
})

# meta data
metas <- list.files(
  path = paste0(getwd(), "/data"),
  pattern = "mirna_meta.tsv",
  full.names = TRUE
  )
names(metas) <- str_extract(metas, "[a-z]{2}_mirna")

metas <- lapply(metas, function(x){
  i <- read.table(x, header = TRUE)
  rownames(i) <- i$sample_id
  return(i)
})
```

Summary of the traits in the datasets, where “de” - German
(treatment-receiving) and “se” - Swedish (treatment-naive) microRNA
datasets:

``` r
# summarize samples per trait
metas %>%
  bind_rows(.id = "dataset") %>% 
  group_by(dataset, diagnosis) %>% tally() %>%
  bind_rows(., tibble(diagnosis = "TOTAL", n = sum(.$n)))
```

    ## # A tibble: 8 x 3
    ## # Groups:   dataset [3]
    ##   dataset  diagnosis     n
    ##   <chr>    <chr>     <int>
    ## 1 de_mirna CD          100
    ## 2 de_mirna HC           65
    ## 3 de_mirna UC           77
    ## 4 se_mirna CD           51
    ## 5 se_mirna HC           30
    ## 6 se_mirna SC           53
    ## 7 se_mirna UC           55
    ## 8 <NA>     TOTAL       431

microRNA QC: filtering threshold is average raw counts \> 1 in dataset.

``` r
# summarize mirna expression
qced_mir <- lapply(counts, function(x){
  x %>% as_tibble(rownames = "mir") %>% 
    gather(sample_id, value, -mir) %>% 
    group_by(mir) %>%
    summarise(mean = mean(value)) %>% 
    filter(mean>1) %>% pull(mir)
})

# number of qced mirnas
lapply(qced_mir, length) 
```

    ## $de_mirna
    ## [1] 554
    ## 
    ## $se_mirna
    ## [1] 533

Differential expression analysis of microRNA using limma/voom linear
models, which include age as a covariate.

``` r
# generate edgeR objects
dge_list <- vector("list", length = length(counts))
names(dge_list) <- names(counts)

for(set in names(dge_list)){
  dge_list[[set]] <- DGEList(
    counts = counts[[set]][qced_mir[[set]],],
    samples = metas[[set]],
    group = metas[[set]]$diagnosis
  )
}

# calculate size factors
dge_list <- lapply(
  dge_list,
  calcNormFactors
  )
  
# model design
design_list <- lapply(dge_list, function(x){
  model.matrix(
    ~0+diagnosis+age,
    data = x$samples
  )
})

# normalize using voom
dge_list_norm <- lapply(
  1:length(dge_list),
  function(i){
    voomWithQualityWeights(
      dge_list[[i]],
      design_list[[i]],
      plot=FALSE
      )
  })
names(dge_list_norm) <- names(dge_list)
```

Results of treatment-naive Swedish (se) dataset (Supplementary Table 3,
2nd sheet).

``` r
# contrast matrix
cnts_matrix_se <- makeContrasts(
  "CD_HC"=diagnosisCD-diagnosisHC,
  "UC_HC"=diagnosisUC-diagnosisHC,
  "CD_UC"=diagnosisCD-diagnosisUC,
  "SC_HC"=diagnosisSC-diagnosisHC,
  "CD_SC"=diagnosisCD-diagnosisSC,
  "UC_SC"=diagnosisUC-diagnosisSC,
  levels=colnames(design_list[["se_mirna"]])
)

# fit to a linear model
se_fit <- lmFit(dge_list_norm[["se_mirna"]], design_list[["se_mirna"]])
se_contr_fit <- eBayes(contrasts.fit(se_fit, cnts_matrix_se))

# get results
results_se <- lapply(seq(colnames(cnts_matrix_se)), function(i){
  topTable(se_contr_fit,
           coef=colnames(cnts_matrix_se)[i],
           number=nrow(se_contr_fit),
           sort.by="logFC", confint=T) %>%
    rownames_to_column("Symbol") %>% 
    mutate(
      comparison=gsub("_", " vs ", colnames(cnts_matrix_se)[i]),
      dea = case_when(
        adj.P.Val<0.05 & abs(logFC)>log2(1.5) ~ "de",
        TRUE ~ "non-de")
      )
}) %>% bind_rows()

# write results
setwd('..')
write_tsv(
  results_se,
  paste0(getwd(), "/results/se_mirna_dea_results.tsv")
  )

# print results
head(results_se)
```

    ##            Symbol     logFC      CI.L      CI.R    AveExpr          t
    ## 1    hsa-miR-6131 -5.124518 -5.798150 -4.450885 -1.1958024 -15.003822
    ## 2  hsa-miR-126-5p -2.453199 -3.061189 -1.845209  4.4646526  -7.958083
    ## 3  hsa-miR-126-3p -1.994438 -2.503659 -1.485216  3.5292224  -7.724771
    ## 4  hsa-miR-144-5p -1.972666 -2.527163 -1.418169  6.1935593  -7.016588
    ## 5  hsa-miR-144-3p -1.716904 -2.424244 -1.009565  0.4110588  -4.787294
    ## 6 hsa-miR-3158-3p  1.608468  1.164745  2.052191  1.3639639   7.149444
    ##        P.Value    adj.P.Val        B comparison dea
    ## 1 2.906297e-34 1.549056e-31 66.22967   CD vs HC  de
    ## 2 1.438648e-13 2.555998e-11 20.26474   CD vs HC  de
    ## 3 5.893150e-13 6.282098e-11 18.90598   CD vs HC  de
    ## 4 3.724310e-11 2.205619e-09 14.80551   CD vs HC  de
    ## 5 3.350694e-06 3.882434e-05  4.05445   CD vs HC  de
    ## 6 1.739055e-11 1.158645e-09 15.72592   CD vs HC  de

Results of treatment-receiving German (de) dataset (Supplementary Table
3, 3rd sheet).

``` r
# contrast matrix
cnts_matrix_de <- makeContrasts(
  "CD_HC"=diagnosisCD-diagnosisHC,
  "UC_HC"=diagnosisUC-diagnosisHC,
  "CD_UC"=diagnosisCD-diagnosisUC,
  levels=colnames(design_list[["de_mirna"]])
)

# fit to a linear model
de_fit <- lmFit(dge_list_norm[["de_mirna"]], design_list[["de_mirna"]])
de_contr_fit <- eBayes(contrasts.fit(de_fit, cnts_matrix_de))

# get results
results_de <- lapply(seq(colnames(cnts_matrix_de)), function(i){
  topTable(de_contr_fit,
           coef=colnames(cnts_matrix_de)[i],
           number=nrow(de_contr_fit),
           sort.by="logFC", confint=T) %>%
    rownames_to_column("Symbol") %>% 
    mutate(
      comparison=gsub("_", " vs ", colnames(cnts_matrix_de)[i]),
      dea = case_when(
        adj.P.Val<0.05 & abs(logFC)>log2(1.5) ~ "de",
        TRUE ~ "non-de")
      )
}) %>% bind_rows()

# write results
setwd('..')
write_tsv(
  results_de,
  paste0(getwd(), "/results/de_mirna_dea_results.tsv")
  )

# print results
head(results_de)
```

    ##            Symbol     logFC      CI.L       CI.R    AveExpr         t
    ## 1   hsa-miR-32-5p -1.513977 -1.998711 -1.0292420  3.9606606 -6.151662
    ## 2  hsa-miR-590-3p -1.491169 -2.037450 -0.9448879 -0.7208165 -5.376358
    ## 3  hsa-miR-144-3p -1.471243 -1.900810 -1.0416752  8.6701703 -6.745750
    ## 4 hsa-miR-3613-5p -1.382139 -1.899250 -0.8650282  0.8913162 -5.264360
    ## 5  hsa-miR-624-5p -1.344426 -1.737589 -0.9512627  0.3917598 -6.735060
    ## 6 hsa-miR-374a-3p -1.326797 -1.731754 -0.9218389  3.0702931 -6.453154
    ##        P.Value    adj.P.Val         B comparison dea
    ## 1 3.062225e-09 1.304979e-07 10.764675   CD vs HC  de
    ## 2 1.757273e-07 5.123840e-06  6.713372   CD vs HC  de
    ## 3 1.070012e-10 9.011447e-09 13.957448   CD vs HC  de
    ## 4 3.052882e-07 7.047069e-06  6.398631   CD vs HC  de
    ## 5 1.138630e-10 9.011447e-09 13.550151   CD vs HC  de
    ## 6 5.728405e-10 2.834962e-08 12.351598   CD vs HC  de

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] edgeR_3.30.3    limma_3.44.3    forcats_0.5.0   stringr_1.4.0  
    ##  [5] dplyr_1.0.0     purrr_0.3.4     readr_1.3.1     tidyr_1.1.3    
    ##  [9] tibble_3.0.3    ggplot2_3.3.2   tidyverse_1.3.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] locfit_1.5-9.4   tidyselect_1.1.0 xfun_0.19        lattice_0.20-41 
    ##  [5] haven_2.3.1      colorspace_1.4-1 vctrs_0.3.8      generics_0.0.2  
    ##  [9] htmltools_0.5.0  yaml_2.2.1       utf8_1.1.4       blob_1.2.1      
    ## [13] rlang_0.4.11     pillar_1.4.6     glue_1.4.1       withr_2.4.2     
    ## [17] DBI_1.1.0        dbplyr_1.4.4     modelr_0.1.8     readxl_1.3.1    
    ## [21] lifecycle_0.2.0  munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0
    ## [25] rvest_0.3.6      evaluate_0.14    knitr_1.29       fansi_0.4.1     
    ## [29] broom_0.7.0      Rcpp_1.0.6       scales_1.1.1     backports_1.1.8 
    ## [33] jsonlite_1.7.2   fs_1.4.2         hms_0.5.3        digest_0.6.25   
    ## [37] stringi_1.4.6    grid_4.0.2       cli_2.0.2        tools_4.0.2     
    ## [41] magrittr_1.5     crayon_1.3.4     pkgconfig_2.0.3  ellipsis_0.3.1  
    ## [45] xml2_1.3.2       reprex_0.3.0     lubridate_1.7.9  assertthat_0.2.1
    ## [49] rmarkdown_2.3    httr_1.4.1       rstudioapi_0.11  R6_2.4.1        
    ## [53] compiler_4.0.2
