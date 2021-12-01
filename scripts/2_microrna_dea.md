microRNA differential expression analysis
================

Load prerequisites that are needed for differential expression analysis
(DEA) in R environment.

``` r
# libraries
library(tidyverse)
library(ggrepel)
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

    ## # A tibble: 8 × 3
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
    ~0+diagnosis+age+gender,
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
    ## 1    hsa-miR-6131 -5.195778 -5.844381 -4.547175 -1.1958024 -15.800035
    ## 2  hsa-miR-126-5p -2.452700 -3.061130 -1.844270  4.4646526  -7.950969
    ## 3  hsa-miR-126-3p -1.992680 -2.501992 -1.483369  3.5292224  -7.716855
    ## 4  hsa-miR-144-5p -1.971117 -2.526856 -1.415379  6.1935593  -6.995653
    ## 5  hsa-miR-144-3p -1.717053 -2.425396 -1.008710  0.4110588  -4.781086
    ## 6 hsa-miR-3158-3p  1.609619  1.164105  2.055134  1.3639639   7.126034
    ##        P.Value    adj.P.Val         B comparison dea
    ## 1 1.317319e-36 7.021308e-34 71.505209   CD vs HC  de
    ## 2 1.531353e-13 2.720705e-11 20.190958   CD vs HC  de
    ## 3 6.287989e-13 6.702997e-11 18.830151   CD vs HC  de
    ## 4 4.248968e-11 2.516334e-09 14.661408   CD vs HC  de
    ## 5 3.456036e-06 3.878306e-05  4.014384   CD vs HC  de
    ## 6 2.016622e-11 1.343575e-09 15.579354   CD vs HC  de

``` r
# plot results
ggplot(results_se, aes(x=AveExpr, y= logFC)) + 
  geom_point(aes(color=dea), size=0.2, alpha=0.2) + 
  scale_color_brewer("", palette = "Set1") +
  geom_hline(yintercept = 0, linetype=2, color = "red") +
  ggrepel::geom_text_repel(
    data = results_se %>% filter(dea %in% "de") %>% 
      group_by(comparison) %>% top_n(n=5, wt=logFC) %>% filter(logFC>(log2(1.5))),
    aes(label=Symbol), size = 3, segment.alpha = 0.1, vjust = 0,hjust = 0,
    direction = "y", nudge_x = 15, nudge_y = .8) +
  ggrepel::geom_text_repel(
    data = results_se %>% filter(dea %in% "de") %>% 
      group_by(comparison) %>% top_n(n=-5, wt=logFC) %>% filter(logFC<(-log2(1.5))),
    aes(label=Symbol), size = 3, segment.alpha = 0.1, vjust = 0, hjust = 0,
    direction = "y", nudge_x = 15, nudge_y = -1.4) +
  facet_wrap(~comparison, nrow = 2) +
  labs(x= "average expression", y= "log2 fold change")  + 
  theme(legend.position = "bottom")
```

![](2_microrna_dea_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

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
    ## 1   hsa-miR-32-5p -1.442549 -1.931428 -0.9536694  3.9606606 -5.811860
    ## 2  hsa-miR-144-3p -1.418011 -1.851563 -0.9844586  8.6701703 -6.442054
    ## 3  hsa-miR-590-3p -1.406883 -1.959987 -0.8537792 -0.7208165 -5.009999
    ## 4  hsa-miR-624-5p -1.362361 -1.761826 -0.9628960  0.3917598 -6.717380
    ## 5 hsa-miR-3613-5p -1.350098 -1.874621 -0.8255739  0.8913162 -5.069747
    ## 6 hsa-miR-374a-3p -1.262929 -1.671248 -0.8546087  3.0702931 -6.092066
    ##        P.Value    adj.P.Val         B comparison dea
    ## 1 1.902895e-08 7.530028e-07  8.970035   CD vs HC  de
    ## 2 6.136667e-10 3.777460e-08 12.204694   CD vs HC  de
    ## 3 1.039045e-06 2.213966e-05  5.262042   CD vs HC  de
    ## 4 1.270785e-10 2.813354e-08 13.681560   CD vs HC  de
    ## 5 7.829811e-07 1.873227e-05  5.551058   CD vs HC  de
    ## 6 4.259087e-09 1.966278e-07 10.430113   CD vs HC  de

``` r
# plot results
ggplot(results_de, aes(x=AveExpr, y= logFC)) + 
  geom_point(aes(color=dea), size=0.2, alpha=0.2) + 
  scale_color_brewer("", palette = "Set1") +
  geom_hline(yintercept = 0, linetype=2, color = "red") +
  geom_text_repel(
    data = results_de %>% filter(dea %in% "de") %>% 
      group_by(comparison) %>% top_n(n=5, wt=logFC) %>% filter(logFC>(log2(1.5))),
    aes(label=Symbol), size = 3, segment.alpha = 0.1, vjust = 0,hjust = 0,
    direction = "y", nudge_x = 10, nudge_y = .8) +
  geom_text_repel(
    data = results_de %>% filter(dea %in% "de") %>% 
      group_by(comparison) %>% top_n(n=-5, wt=logFC) %>% filter(logFC<(-log2(1.5))),
    aes(label=Symbol), size = 3, segment.alpha = 0.1, vjust = 0, hjust = 0,
    direction = "y", nudge_x = 20, nudge_y = -1.4) +
  facet_wrap(~comparison, nrow = 2) +
  labs(x= "average expression", y= "log2 fold change")  + 
  theme(legend.position = "bottom")
```

![](2_microrna_dea_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] edgeR_3.34.1    limma_3.48.3    ggrepel_0.9.1   forcats_0.5.1  
    ##  [5] stringr_1.4.0   dplyr_1.0.7     purrr_0.3.4     readr_2.0.2    
    ##  [9] tidyr_1.1.4     tibble_3.1.5    ggplot2_3.3.5   tidyverse_1.3.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.7         locfit_1.5-9.4     lubridate_1.8.0    lattice_0.20-45   
    ##  [5] assertthat_0.2.1   digest_0.6.28      utf8_1.2.2         R6_2.5.1          
    ##  [9] cellranger_1.1.0   backports_1.2.1    reprex_2.0.1       evaluate_0.14     
    ## [13] highr_0.9          httr_1.4.2         pillar_1.6.3       rlang_0.4.11      
    ## [17] readxl_1.3.1       rstudioapi_0.13    rmarkdown_2.11     labeling_0.4.2    
    ## [21] bit_4.0.4          munsell_0.5.0      broom_0.7.9        compiler_4.1.1    
    ## [25] modelr_0.1.8       xfun_0.26          pkgconfig_2.0.3    htmltools_0.5.2   
    ## [29] tidyselect_1.1.1   fansi_0.5.0        crayon_1.4.1       tzdb_0.1.2        
    ## [33] dbplyr_2.1.1       withr_2.4.2        grid_4.1.1         jsonlite_1.7.2    
    ## [37] gtable_0.3.0       lifecycle_1.0.1    DBI_1.1.1          magrittr_2.0.1    
    ## [41] scales_1.1.1       cli_3.0.1          stringi_1.7.5      vroom_1.5.5       
    ## [45] farver_2.1.0       fs_1.5.0           xml2_1.3.2         ellipsis_0.3.2    
    ## [49] generics_0.1.0     vctrs_0.3.8        RColorBrewer_1.1-2 tools_4.1.1       
    ## [53] bit64_4.0.5        glue_1.4.2         hms_1.1.1          parallel_4.1.1    
    ## [57] fastmap_1.1.0      yaml_2.2.1         colorspace_2.0-2   rvest_1.0.1       
    ## [61] knitr_1.36         haven_2.4.3
