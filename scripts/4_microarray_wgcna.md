Weighted gene correlation analysis of BeadChip microarray data
================

Load prerequisites that will be needed for weighted gene co-expression
network analysis (WGCNA) in R environment.

``` r
# libraries
library(tidyverse)
library(genefilter)
library(Biobase)
library(WGCNA)
```

Read normalized gene expression intensity values, which were generated
using microarray\_dea.Rmd script.

``` r
# set wd to main
setwd('..')

# exprs data
exprs_rna <- read_tsv(
  paste0(getwd(), "/data/se_rna_exprs.tsv")
) %>% column_to_rownames("gene")


# meta data
meta_rna <- read_tsv(
  paste0(getwd(),"/data/se_rna_meta.tsv")
  )

meta_mirna <- read_tsv(
  paste0(getwd(),"/data/se_mirna_meta.tsv")
  )

# keep overlapping
meta_rna <- meta_rna %>% 
  filter(sample_id %in% meta_mirna$sample_id)
```

Retain the most variable genes using varFilter function from genefilter.

``` r
# keep variable genes
exprs_rna_filt <- exprs(
  varFilter(
    ExpressionSet(
      assayData = as.matrix(exprs_rna)
      ),
    var.func = "sd",
    var.cutoff=0.5
    )
  )
# print number of genes
dim(exprs_rna_filt)
```

    ## [1] 5863  205

Filter outlying samples based on network connectivity (WGCNA workflow):

``` r
# split expression per trait
exprs_list <- exprs_rna_filt %>% t() %>%
  as_tibble(rownames = "array_id") %>% 
  left_join(
    meta_rna %>%
      select(array_id, diagnosis),., by = "array_id"
    ) %>%
  column_to_rownames("array_id") %>%
  split(., .$diagnosis) %>% 
  lapply(., function(x){
    x %>% select(-diagnosis) %>% t()
    })

# squared Euclidean distance
adj_list <- lapply(
  exprs_list,
  adjacency,
  type = "distance"
  )

# connectivity
k_list <- lapply(adj_list, function(x){
  as.numeric(apply(x, 2, sum)) - 1
  })

# standardized connectivity
Z_k_list <- lapply(k_list, scale)

# threshold
thresholdZ_k <- -2.5

# get outlying samples
remove_list <- 
  lapply(Z_k_list, function(x){
    x < thresholdZ_k | is.na(x)
    })

# remove the outliers from expression data
exprs_list_qced <- 
  lapply(seq(exprs_list), function(x) {
    exprs_list[[x]][,!remove_list[[x]]]
    })
names(exprs_list_qced) <- names(exprs_list)

# print dims
lapply(exprs_list_qced, dim)
```

    ## $CD
    ## [1] 5863   49
    ## 
    ## $HC
    ## [1] 5863   29
    ## 
    ## $SC
    ## [1] 5863   51
    ## 
    ## $UC
    ## [1] 5863   54

Calculate soft-thresholding powers for construction of trait-wise
networks.

``` r
# set of soft thresholding powers
powers = c(1:20)

# power based on SFT criterion
sft_list <- lapply(exprs_list_qced, function(x){
  pickSoftThreshold(
    t(x),
    powerVector = powers,
    networkType = "signed",
    RsquaredCut = 0.8
    )
})
```

    ##    Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1    0.338  3.32          0.557  2980.0    3000.0   3280
    ## 2      2    0.252 -1.95          0.154  1710.0    1630.0   2220
    ## 3      3    0.350 -1.57          0.378  1080.0     993.0   1680
    ## 4      4    0.501 -1.37          0.579   728.0     648.0   1350
    ## 5      5    0.626 -1.26          0.715   516.0     444.0   1120
    ## 6      6    0.694 -1.18          0.789   380.0     315.0    945
    ## 7      7    0.759 -1.15          0.840   289.0     231.0    811
    ## 8      8    0.794 -1.16          0.862   225.0     172.0    704
    ## 9      9    0.831 -1.17          0.892   179.0     131.0    616
    ## 10    10    0.859 -1.18          0.906   145.0     101.0    544
    ## 11    11    0.889 -1.17          0.929   118.0      78.7    482
    ## 12    12    0.908 -1.18          0.945    98.3      61.8    430
    ## 13    13    0.920 -1.19          0.953    82.5      49.3    385
    ## 14    14    0.933 -1.19          0.964    69.8      39.4    346
    ## 15    15    0.935 -1.21          0.969    59.5      32.0    314
    ## 16    16    0.930 -1.23          0.967    51.2      26.0    287
    ## 17    17    0.930 -1.25          0.968    44.2      21.2    263
    ## 18    18    0.932 -1.27          0.972    38.5      17.4    242
    ## 19    19    0.919 -1.29          0.964    33.6      14.4    223
    ## 20    20    0.909 -1.31          0.959    29.5      12.0    206
    ##    Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1    0.115  4.83          0.807 2950.00   2960.00 3170.0
    ## 2      2    0.463 -6.67          0.734 1630.00   1590.00 2000.0
    ## 3      3    0.630 -4.85          0.839  959.00    923.00 1400.0
    ## 4      4    0.686 -3.48          0.902  598.00    569.00 1040.0
    ## 5      5    0.708 -2.74          0.925  391.00    366.00  805.0
    ## 6      6    0.738 -2.36          0.945  266.00    244.00  641.0
    ## 7      7    0.780 -2.13          0.964  187.00    166.00  521.0
    ## 8      8    0.802 -2.00          0.973  135.00    116.00  431.0
    ## 9      9    0.825 -1.91          0.978   99.90     82.80  361.0
    ## 10    10    0.852 -1.84          0.988   75.40     60.10  305.0
    ## 11    11    0.864 -1.81          0.989   57.90     44.10  261.0
    ## 12    12    0.880 -1.80          0.993   45.20     33.00  226.0
    ## 13    13    0.887 -1.78          0.993   35.80     24.90  197.0
    ## 14    14    0.898 -1.77          0.997   28.70     18.90  173.0
    ## 15    15    0.900 -1.77          0.995   23.20     14.50  152.0
    ## 16    16    0.906 -1.76          0.995   19.00     11.20  135.0
    ## 17    17    0.904 -1.78          0.990   15.60      8.72  120.0
    ## 18    18    0.907 -1.78          0.989   13.00      6.83  107.0
    ## 19    19    0.914 -1.77          0.991   10.90      5.42   95.6
    ## 20    20    0.917 -1.77          0.990    9.18      4.30   85.9
    ##    Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1    0.617  3.32          0.602  2980.0    3000.0   3290
    ## 2      2    0.350 -2.03          0.194  1720.0    1640.0   2240
    ## 3      3    0.238 -1.26          0.244  1090.0     990.0   1700
    ## 4      4    0.377 -1.20          0.481   734.0     650.0   1370
    ## 5      5    0.504 -1.18          0.655   521.0     448.0   1140
    ## 6      6    0.607 -1.20          0.754   385.0     317.0    970
    ## 7      7    0.673 -1.22          0.805   293.0     233.0    838
    ## 8      8    0.727 -1.27          0.842   228.0     175.0    732
    ## 9      9    0.762 -1.28          0.867   182.0     133.0    645
    ## 10    10    0.782 -1.32          0.875   147.0     102.0    573
    ## 11    11    0.807 -1.34          0.893   120.0      80.1    511
    ## 12    12    0.829 -1.36          0.908   100.0      63.5    459
    ## 13    13    0.829 -1.38          0.906    83.9      50.7    414
    ## 14    14    0.829 -1.40          0.905    71.0      40.8    374
    ## 15    15    0.830 -1.42          0.906    60.5      33.0    340
    ## 16    16    0.832 -1.43          0.909    52.0      26.8    310
    ## 17    17    0.833 -1.44          0.911    44.9      21.9    283
    ## 18    18    0.804 -1.48          0.891    39.0      18.1    259
    ## 19    19    0.815 -1.49          0.899    34.1      15.0    238
    ## 20    20    0.814 -1.49          0.903    29.9      12.4    219
    ##    Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1    0.374  3.15          0.358  2980.0    2990.0   3290
    ## 2      2    0.246 -1.81          0.044  1710.0    1620.0   2220
    ## 3      3    0.401 -1.57          0.244  1070.0     960.0   1670
    ## 4      4    0.522 -1.31          0.453   717.0     619.0   1330
    ## 5      5    0.637 -1.24          0.639   506.0     420.0   1110
    ## 6      6    0.722 -1.22          0.748   372.0     294.0    947
    ## 7      7    0.763 -1.24          0.807   282.0     213.0    821
    ## 8      8    0.784 -1.27          0.833   219.0     157.0    720
    ## 9      9    0.811 -1.30          0.861   174.0     118.0    637
    ## 10    10    0.821 -1.33          0.870   141.0      89.8    568
    ## 11    11    0.817 -1.38          0.864   115.0      69.7    510
    ## 12    12    0.804 -1.42          0.858    95.6      54.8    460
    ## 13    13    0.801 -1.44          0.859    80.2      43.3    417
    ## 14    14    0.790 -1.47          0.853    67.9      34.4    379
    ## 15    15    0.785 -1.49          0.854    57.9      27.7    346
    ## 16    16    0.786 -1.50          0.860    49.7      22.4    316
    ## 17    17    0.784 -1.52          0.867    43.0      18.3    290
    ## 18    18    0.795 -1.52          0.877    37.4      15.1    267
    ## 19    19    0.789 -1.53          0.879    32.7      12.5    246
    ## 20    20    0.797 -1.54          0.890    28.7      10.4    228

``` r
# add signed value of SFT R^2
sft_list <- lapply(sft_list, function(x) {
  x$fitIndices <- x$fitIndices %>% 
    mutate(SFT.R.sq.signed = -sign(slope)*SFT.R.sq)
  return(x)
  })
```

Calculate topological overlap matrices for each trait

``` r
# get median threshold of networks
sft <- lapply(
  sft_list, function(x) x$powerEstimate
  ) %>% 
  bind_rows() %>%
  gather(diagnosis, soft) %>%
  summarise(soft=round(median(soft))) %>% 
  pull(soft)

# list of topological overlap matrices
toms_list <- lapply(exprs_list_qced, function(x){
  TOMsimilarityFromExpr(
    t(x),
    corType = "bicor",
    networkType = "signed",
    power = sft,
    TOMDenom = "mean",
    nThreads = 4
    )
  })
```

    ## TOM calculation: adjacency..
    ## ..will use 4 parallel threads.
    ##  Fraction of slow calculations: 0.000000
    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.
    ## TOM calculation: adjacency..
    ## ..will use 4 parallel threads.
    ##  Fraction of slow calculations: 0.000000
    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.
    ## TOM calculation: adjacency..
    ## ..will use 4 parallel threads.
    ##  Fraction of slow calculations: 0.000000
    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.
    ## TOM calculation: adjacency..
    ## ..will use 4 parallel threads.
    ##  Fraction of slow calculations: 0.000000
    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.

``` r
# add gene names
toms_list <- lapply(toms_list, function(x){
  colnames(x) <- rownames(exprs_list_qced$HC)
  rownames(x) <- rownames(exprs_list_qced$HC)
  return(x)
})

# write for tensor decomposition
setwd('..')
for(i in seq(toms_list)){
  x <- toms_list[[i]] %>% as.data.frame()
  write_csv(
    x, path = gzfile(
      paste0(getwd(), "/data/tom_mrna_", names(toms_list)[i], ".csv.gz")
      )
  )
}

# print dims
lapply(toms_list, dim)
```

    ## $CD
    ## [1] 5863 5863
    ## 
    ## $HC
    ## [1] 5863 5863
    ## 
    ## $SC
    ## [1] 5863 5863
    ## 
    ## $UC
    ## [1] 5863 5863

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
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] WGCNA_1.69            fastcluster_1.1.25    dynamicTreeCut_1.63-1
    ##  [4] Biobase_2.48.0        BiocGenerics_0.34.0   genefilter_1.70.0    
    ##  [7] forcats_0.5.0         stringr_1.4.0         dplyr_1.0.0          
    ## [10] purrr_0.3.4           readr_1.3.1           tidyr_1.1.3          
    ## [13] tibble_3.0.3          ggplot2_3.3.2         tidyverse_1.3.0      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-7          matrixStats_0.58.0    fs_1.4.2             
    ##  [4] lubridate_1.7.9       bit64_0.9-7           RColorBrewer_1.1-2   
    ##  [7] doParallel_1.0.15     httr_1.4.1            tools_4.0.2          
    ## [10] backports_1.1.8       R6_2.4.1              rpart_4.1-15         
    ## [13] Hmisc_4.4-0           DBI_1.1.0             colorspace_1.4-1     
    ## [16] nnet_7.3-14           withr_2.4.2           gridExtra_2.3        
    ## [19] tidyselect_1.1.0      preprocessCore_1.50.0 bit_1.1-15.2         
    ## [22] compiler_4.0.2        cli_2.0.2             rvest_0.3.6          
    ## [25] htmlTable_2.0.1       xml2_1.3.2            checkmate_2.0.0      
    ## [28] scales_1.1.1          digest_0.6.25         foreign_0.8-80       
    ## [31] rmarkdown_2.3         jpeg_0.1-8.1          base64enc_0.1-3      
    ## [34] pkgconfig_2.0.3       htmltools_0.5.0       dbplyr_1.4.4         
    ## [37] htmlwidgets_1.5.1     rlang_0.4.11          readxl_1.3.1         
    ## [40] impute_1.62.0         rstudioapi_0.11       RSQLite_2.2.0        
    ## [43] generics_0.0.2        jsonlite_1.7.2        acepack_1.4.1        
    ## [46] RCurl_1.98-1.3        magrittr_1.5          GO.db_3.11.4         
    ## [49] Formula_1.2-3         Matrix_1.2-18         Rcpp_1.0.6           
    ## [52] munsell_0.5.0         S4Vectors_0.26.1      fansi_0.4.1          
    ## [55] lifecycle_0.2.0       stringi_1.4.6         yaml_2.2.1           
    ## [58] grid_4.0.2            blob_1.2.1            crayon_1.3.4         
    ## [61] lattice_0.20-41       haven_2.3.1           splines_4.0.2        
    ## [64] annotate_1.66.0       hms_0.5.3             knitr_1.29           
    ## [67] pillar_1.4.6          codetools_0.2-16      stats4_4.0.2         
    ## [70] reprex_0.3.0          XML_3.99-0.4          glue_1.4.1           
    ## [73] evaluate_0.14         latticeExtra_0.6-29   data.table_1.12.8    
    ## [76] modelr_0.1.8          png_0.1-7             vctrs_0.3.8          
    ## [79] foreach_1.5.0         cellranger_1.1.0      gtable_0.3.0         
    ## [82] assertthat_0.2.1      xfun_0.19             xtable_1.8-4         
    ## [85] broom_0.7.0           survival_3.2-3        iterators_1.0.12     
    ## [88] AnnotationDbi_1.50.1  memoise_1.1.0         IRanges_2.22.2       
    ## [91] cluster_2.1.0         ellipsis_0.3.1
