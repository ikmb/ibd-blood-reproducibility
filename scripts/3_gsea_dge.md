Gene set enrichment analysis of differentially expressed genes
================

Load prerequisites that will be needed for gene set enrichment analysis
(GSEA) in R environment.

``` r
# libraries
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(multiMiR)
```

Read results of differential gene expression, which were generated using
microrna\_dea.Rmd and microarray\_dea.Rmd.

``` r
# set wd to main
setwd('..')

# dea results
dea_res_list <- 
  list.files(
    path = paste0(getwd(),"/results"),
    pattern = "dea_results.tsv",
    full.names = TRUE
    ) 

dea_res <- lapply(dea_res_list, read_tsv)
names(dea_res) <- str_extract(
  dea_res_list, ".e_.*_dea"
)
dea_res
```

    ## $de_mirna_dea
    ## # A tibble: 1,662 x 11
    ##    Symbol logFC   CI.L   CI.R AveExpr     t  P.Value adj.P.Val     B comparison
    ##    <chr>  <dbl>  <dbl>  <dbl>   <dbl> <dbl>    <dbl>     <dbl> <dbl> <chr>     
    ##  1 hsa-m… -1.51 -2.00  -1.03    3.96  -6.15 3.06e- 9   1.30e-7 10.8  CD vs HC  
    ##  2 hsa-m… -1.49 -2.04  -0.945  -0.721 -5.38 1.76e- 7   5.12e-6  6.71 CD vs HC  
    ##  3 hsa-m… -1.47 -1.90  -1.04    8.67  -6.75 1.07e-10   9.01e-9 14.0  CD vs HC  
    ##  4 hsa-m… -1.38 -1.90  -0.865   0.891 -5.26 3.05e- 7   7.05e-6  6.40 CD vs HC  
    ##  5 hsa-m… -1.34 -1.74  -0.951   0.392 -6.74 1.14e-10   9.01e-9 13.6  CD vs HC  
    ##  6 hsa-m… -1.33 -1.73  -0.922   3.07  -6.45 5.73e-10   2.83e-8 12.4  CD vs HC  
    ##  7 hsa-m… -1.26 -1.63  -0.896   5.62  -6.78 8.86e-11   9.01e-9 14.1  CD vs HC  
    ##  8 hsa-m… -1.25 -1.72  -0.779   4.89  -5.24 3.38e- 7   7.49e-6  6.28 CD vs HC  
    ##  9 hsa-m…  1.22  0.719  1.73    3.44   4.77 3.21e- 6   5.56e-5  4.19 CD vs HC  
    ## 10 hsa-m… -1.12 -1.59  -0.653  -0.679 -4.70 4.23e- 6   6.72e-5  3.91 CD vs HC  
    ## # … with 1,652 more rows, and 1 more variable: dea <chr>
    ## 
    ## $se_mirna_dea
    ## # A tibble: 3,198 x 11
    ##    Symbol logFC   CI.L   CI.R AveExpr      t  P.Value adj.P.Val     B comparison
    ##    <chr>  <dbl>  <dbl>  <dbl>   <dbl>  <dbl>    <dbl>     <dbl> <dbl> <chr>     
    ##  1 hsa-m… -5.12 -5.80  -4.45   -1.20  -15.0  2.91e-34  1.55e-31 66.2  CD vs HC  
    ##  2 hsa-m… -2.45 -3.06  -1.85    4.46   -7.96 1.44e-13  2.56e-11 20.3  CD vs HC  
    ##  3 hsa-m… -1.99 -2.50  -1.49    3.53   -7.72 5.89e-13  6.28e-11 18.9  CD vs HC  
    ##  4 hsa-m… -1.97 -2.53  -1.42    6.19   -7.02 3.72e-11  2.21e- 9 14.8  CD vs HC  
    ##  5 hsa-m… -1.72 -2.42  -1.01    0.411  -4.79 3.35e- 6  3.88e- 5  4.05 CD vs HC  
    ##  6 hsa-m…  1.61  1.16   2.05    1.36    7.15 1.74e-11  1.16e- 9 15.7  CD vs HC  
    ##  7 hsa-m…  1.46  0.913  2.01   -0.790   5.27 3.69e- 7  5.96e- 6  6.24 CD vs HC  
    ##  8 hsa-m… -1.45 -2.02  -0.873   2.08   -4.98 1.43e- 6  1.85e- 5  4.71 CD vs HC  
    ##  9 hsa-m… -1.41 -1.89  -0.928   7.57   -5.77 3.08e- 8  6.57e- 7  8.25 CD vs HC  
    ## 10 hsa-m… -1.40 -1.89  -0.915  -1.92   -5.68 4.91e- 8  9.35e- 7  8.16 CD vs HC  
    ## # … with 3,188 more rows, and 1 more variable: dea <chr>
    ## 
    ## $se_rna_dea
    ## # A tibble: 70,362 x 14
    ##    Probe_Id Array_Address_Id Symbol ENTREZ logFC  CI.L  CI.R AveExpr     t
    ##    <chr>               <dbl> <chr>   <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>
    ##  1 ILMN_16…          4830541 HBE1   3.05e3  2.62  2.21  3.04   11.6  12.4 
    ##  2 ILMN_16…           240594 ANXA3  3.06e2  2.02  1.60  2.43    8.95  9.63
    ##  3 ILMN_16…          4570725 H3C13  6.54e5  1.95  1.64  2.26    7.25 12.5 
    ##  4 ILMN_23…          3290368 IGFBP1 3.48e3  1.89  1.52  2.27    8.39  9.91
    ##  5 ILMN_17…          1410221 S100A… 6.28e3  1.81  1.46  2.15   11.1  10.3 
    ##  6 ILMN_17…            70167 LY96   2.36e4  1.77  1.41  2.12    9.14  9.83
    ##  7 ILMN_21…          2260576 COX7B  1.35e3  1.74  1.33  2.15    7.04  8.35
    ##  8 ILMN_18…          1230164 MTRNR… 1.00e8  1.66  1.41  1.91    7.83 13.1 
    ##  9 ILMN_17…          4260484 COMMD6 1.71e5  1.64  1.18  2.10    8.45  7.05
    ## 10 ILMN_21…          3310161 TDP2   5.16e4  1.47  1.04  1.90    6.94  6.73
    ## # … with 70,352 more rows, and 5 more variables: P.Value <dbl>,
    ## #   adj.P.Val <dbl>, B <dbl>, comparison <chr>, dea <chr>

GSEA for deferentially expressed genes and microRNAs was performed using
in blood expressed genes as a universe for over-representation analysis
(ORA).

``` r
# all expressed genes in peripheral blood dataset
expressed_genes <- dea_res$se_rna_dea %>% 
  dplyr::select(
    Symbol, ENTREZ
  ) %>% distinct() %>% 
  mutate(
    ENTREZID=as.character(ENTREZ)
    )

# print number of genes
nrow(expressed_genes)
```

    ## [1] 11727

GSEA of deferentially expressed genes using GO terms (treatment-naive
cohort):

``` r
# get deferentially expressed genes
se_dea_gene_list <- dea_res$se_rna_dea %>%
  filter(dea %in% "de") %>% 
  split(., f = .$comparison) %>%
  lapply(., function(x){
    x %>% dplyr::select(Symbol, ENTREZ) %>% 
      mutate(ENTREZID = as.character(ENTREZ))
    })

# remove empty elements (such as CD vs UC)
se_dea_gene_list <- compact(se_dea_gene_list)

# print de genes
lapply(se_dea_gene_list, nrow)
```

    ## $`CD vs HC`
    ## [1] 602
    ## 
    ## $`CD vs SC`
    ## [1] 33
    ## 
    ## $`SC vs HC`
    ## [1] 146
    ## 
    ## $`UC vs HC`
    ## [1] 478
    ## 
    ## $`UC vs SC`
    ## [1] 31

``` r
# perform run over-representation analysis
# using expressed genes as universe
se_dea_gene_go <- lapply(
  se_dea_gene_list, function(x) {
    enrichGO(
      gene = x$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      universe = expressed_genes$ENTREZID
      )
  })

# combine results to table for plotting
se_dea_gene_go_tb <- lapply(
  se_dea_gene_go, function(x){
     x %>% .@result  }
  ) %>% bind_rows(.id = "comparison") %>% 
  filter(p.adjust<0.05)

# visualize top sig terms
se_dea_gene_go_tb %>% 
  group_by(comparison) %>%
  top_n(-15, wt=p.adjust) %>% 
  ggplot(aes(x = Count, y = fct_reorder(Description, Count))) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  facet_wrap(~comparison, nrow = 1) + ylab(NULL) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical")
```

![](3_gsea_dge_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

GSEA of validated target genes of deferentially expressed miRNAs.
Analysis was performed using target genes, that are found to be
expressed in blood (using microarray data to define expression status).
The miRTarbase was was accessed via multiMiR package.

Treatment-naive Swedish (se) cohort:

``` r
# get deferentially expressed mirna
se_dea_mir_list <- dea_res$se_mirna_dea %>%
  filter(dea %in% "de") %>% 
  split(., f = .$comparison) %>%
  lapply(., function(x){x %>% .$Symbol})

# remove empty elements (such as CD vs UC)
se_dea_mir_list <- compact(se_dea_mir_list)

# print de mir
lapply(se_dea_mir_list, length)
```

    ## $`CD vs HC`
    ## [1] 115
    ## 
    ## $`SC vs HC`
    ## [1] 97
    ## 
    ## $`UC vs HC`
    ## [1] 101

``` r
# validated MTIs for all comparisons
se_mirtar_list <- 
  lapply(se_dea_mir_list, function(x){
    y <- get_multimir(
      org = "hsa",
      mirna = x,
      target = expressed_genes$ENTREZID,
      table = "mirtarbase",
      summary = TRUE)
    return(y@data)
})
```

    ## Searching mirtarbase ...
    ## Searching mirtarbase ...
    ## Searching mirtarbase ...

``` r
# filter MTIs based on experiment type
se_mirtar_list_filt <- lapply(se_mirtar_list, function(x){
  x %>% filter(
    grepl("Luciferase|Western|pSILAC|Proteomics", experiment)
    )
  })

# run over-representation analysis
# using expressed genes as universe
se_mirtar_go <- lapply(se_mirtar_list_filt, function(x){
  enrichGO(
    gene = unique(x$target_entrez),
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    universe = expressed_genes$ENTREZID
    )
})

# flatten results to table
se_mirtar_go_tb <- lapply(
  se_mirtar_go, function(x) x %>% .@result
  ) %>% bind_rows(.id = "comparison") %>%
  filter(p.adjust<0.05)

# visualize
se_mirtar_go_tb %>% 
  group_by(comparison) %>% 
  top_n(-15, wt=p.adjust) %>% 
  ggplot(aes(x = Count, y = fct_reorder(Description, Count))) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  facet_wrap(~comparison, nrow = 1) + ylab(NULL) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical")
```

![](3_gsea_dge_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Treatment-receiving German (de) cohort:

``` r
# get deferentially expressed mirna
de_dea_mir_list <- dea_res$de_mirna_dea %>%
  filter(dea %in% "de") %>% 
  split(., f = .$comparison) %>%
  lapply(., function(x){x %>% .$Symbol})

# remove empty elements (such as CD vs UC)
de_dea_mir_list <- compact(de_dea_mir_list)

# print de mir
lapply(de_dea_mir_list, length)
```

    ## $`CD vs HC`
    ## [1] 75
    ## 
    ## $`UC vs HC`
    ## [1] 58

``` r
# validated MTIs for all comparisons
de_mirtar_list <- 
  lapply(de_dea_mir_list, function(x){
    y <- get_multimir(
      org = "hsa",
      mirna = x,
      target = expressed_genes$ENTREZID,
      table = "mirtarbase",
      summary = TRUE)
    return(y@data)
})
```

    ## Searching mirtarbase ...
    ## Searching mirtarbase ...

``` r
# filter MTIs based on experiment type
de_mirtar_list_filt <- lapply(de_mirtar_list, function(x){
  x %>% filter(
    grepl("Luciferase|Western|pSILAC|Proteomics", experiment)
    )
})

# run enrichment analysis
de_mirtar_go <- lapply(de_mirtar_list_filt, function(x){
  enrichGO(
    gene = unique(x$target_entrez),
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    universe = expressed_genes$ENTREZID
    )
})

# flatten results to table
de_mirtar_go_tb <- lapply(
  de_mirtar_go, function(x) x %>% .@result
  ) %>% bind_rows(.id = "comparison") %>%
  filter(p.adjust<0.05)

# visualize
de_mirtar_go_tb %>%
  group_by(comparison) %>%
  top_n(-15, wt=p.adjust) %>% 
  ggplot(aes(x = Count, y = fct_reorder(Description, Count))) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  facet_wrap(~comparison, nrow = 1) + ylab(NULL) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical")
```

![](3_gsea_dge_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Enriched MTI pathways that are overlapping with DEG most (top 20)
enriched pathways:

``` r
# overlapping paths
over_paths <- left_join(
  se_dea_gene_go_tb %>% 
    group_by(comparison) %>% 
    top_n(-20, wt=p.adjust) %>% 
    ungroup(),
  se_mirtar_go_tb,
  by = c("comparison", "ID", "Description")
  ) %>% 
  rename_all(
    ~gsub(".x", "_semrna", .x, fixed=TRUE)
    ) %>% 
  rename_all(
    ~gsub(".y", "_semirna", .x, fixed=TRUE)
    ) %>% 
  left_join(
    ., de_mirtar_go_tb %>% 
      rename_at(
        vars(GeneRatio:Count), funs(paste0(., "_demirna"))
        ),
    by = c("comparison", "ID", "Description")
    ) %>% 
  pivot_longer(
    cols = -c(comparison, ID, Description),
    names_to = c( '.value', 'grp'),
    names_pattern = "^(.*_)(.*)"
    ) %>% 
  rename_all(
    ~gsub("_", "", .x)
    )

# visualize
over_paths %>% rowwise() %>% 
  mutate(GeneRatio=eval(parse(text=GeneRatio[1]))) %>%
  mutate(grp = factor(grp, levels = c("semrna", "semirna", "demirna"))) %>% 
  ggplot(aes(x = comparison, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  facet_grid(.~grp, scales = "free_x", space = "free") +
  ylab(NULL)
```

![](3_gsea_dge_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

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
    ##  [1] multiMiR_1.10.0        org.Hs.eg.db_3.11.4    AnnotationDbi_1.50.1  
    ##  [4] IRanges_2.22.2         S4Vectors_0.26.1       Biobase_2.48.0        
    ##  [7] BiocGenerics_0.34.0    clusterProfiler_3.16.0 forcats_0.5.0         
    ## [10] stringr_1.4.0          dplyr_1.0.0            purrr_0.3.4           
    ## [13] readr_1.3.1            tidyr_1.1.3            tibble_3.0.3          
    ## [16] ggplot2_3.3.2          tidyverse_1.3.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fgsea_1.14.0        colorspace_1.4-1    ellipsis_0.3.1     
    ##   [4] ggridges_0.5.2      qvalue_2.20.0       fs_1.4.2           
    ##   [7] rstudioapi_0.11     farver_2.0.3        urltools_1.7.3     
    ##  [10] graphlayouts_0.7.0  ggrepel_0.9.0.9999  bit64_0.9-7        
    ##  [13] fansi_0.4.1         scatterpie_0.1.4    lubridate_1.7.9    
    ##  [16] xml2_1.3.2          splines_4.0.2       GOSemSim_2.14.0    
    ##  [19] knitr_1.29          polyclip_1.10-0     jsonlite_1.7.2     
    ##  [22] broom_0.7.0         GO.db_3.11.4        dbplyr_1.4.4       
    ##  [25] ggforce_0.3.2       BiocManager_1.30.10 compiler_4.0.2     
    ##  [28] httr_1.4.1          rvcheck_0.1.8       backports_1.1.8    
    ##  [31] assertthat_0.2.1    Matrix_1.2-18       cli_2.0.2          
    ##  [34] tweenr_1.0.1        htmltools_0.5.0     prettyunits_1.1.1  
    ##  [37] tools_4.0.2         igraph_1.2.5        gtable_0.3.0       
    ##  [40] glue_1.4.1          reshape2_1.4.4      DO.db_2.9          
    ##  [43] fastmatch_1.1-0     Rcpp_1.0.6          enrichplot_1.8.1   
    ##  [46] cellranger_1.1.0    vctrs_0.3.8         ggraph_2.0.3       
    ##  [49] xfun_0.19           rvest_0.3.6         lifecycle_0.2.0    
    ##  [52] XML_3.99-0.4        DOSE_3.14.0         europepmc_0.4      
    ##  [55] MASS_7.3-51.6       scales_1.1.1        tidygraph_1.2.0    
    ##  [58] hms_0.5.3           RColorBrewer_1.1-2  yaml_2.2.1         
    ##  [61] memoise_1.1.0       gridExtra_2.3       downloader_0.4     
    ##  [64] triebeard_0.3.0     stringi_1.4.6       RSQLite_2.2.0      
    ##  [67] BiocParallel_1.22.0 bitops_1.0-7        rlang_0.4.11       
    ##  [70] pkgconfig_2.0.3     evaluate_0.14       lattice_0.20-41    
    ##  [73] labeling_0.3        cowplot_1.0.0       bit_1.1-15.2       
    ##  [76] tidyselect_1.1.0    plyr_1.8.6          magrittr_1.5       
    ##  [79] R6_2.4.1            generics_0.0.2      DBI_1.1.0          
    ##  [82] pillar_1.4.6        haven_2.3.1         withr_2.4.2        
    ##  [85] RCurl_1.98-1.3      modelr_0.1.8        crayon_1.3.4       
    ##  [88] utf8_1.1.4          rmarkdown_2.3       viridis_0.5.1      
    ##  [91] progress_1.2.2      grid_4.0.2          readxl_1.3.1       
    ##  [94] data.table_1.12.8   blob_1.2.1          reprex_0.3.0       
    ##  [97] digest_0.6.25       gridGraphics_0.5-0  munsell_0.5.0      
    ## [100] viridisLite_0.3.0   ggplotify_0.0.5
