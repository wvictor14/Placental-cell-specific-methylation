---
title: "16_Ethnicity_ancestry"
author: "Victor Yuan"
date: "29/05/2019"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 4
    toc_float:
      collapsed: false
    theme: spacelab
editor_options: 
  chunk_output_type: console
---

For adding ethnicity/ancestry information to the samples, 


```r
library(planet)
library(minfi)
```

```
## Warning: package 'minfi' was built under R version 3.5.2
```

```
## Warning: package 'GenomeInfoDb' was built under R version 3.5.2
```

```
## Warning: package 'matrixStats' was built under R version 3.5.2
```

```
## Warning: package 'BiocParallel' was built under R version 3.5.2
```

```
## Warning: package 'Biostrings' was built under R version 3.5.2
```

```
## Warning: package 'foreach' was built under R version 3.5.2
```

```
## Warning: package 'iterators' was built under R version 3.5.2
```

```
## Warning: package 'locfit' was built under R version 3.5.2
```

```r
library(wateRmelon)
```

```
## Warning: package 'scales' was built under R version 3.5.2
```

```
## Warning: package 'reshape2' was built under R version 3.5.2
```

```
## Warning: package 'ggplot2' was built under R version 3.5.2
```

```
## Warning: package 'GenomicFeatures' was built under R version 3.5.2
```

```r
library(IlluminaHumanMethylationEPICmanifest)
library(dplyr)
```

```
## Warning: package 'dplyr' was built under R version 3.5.2
```

```r
library(egg)
```

```
## Warning: package 'egg' was built under R version 3.5.3
```

```
## Warning: package 'gridExtra' was built under R version 3.5.3
```

```r
mset_noob <- readRDS('../../data/main/interim/14_mset_noob.rds')
rgset <- readRDS('../../data/main/interim/01_rgset_raw.rds')
pDat <- readRDS('../../data/main/interim/13_pDat.rds')
color_code <- readRDS('../../data/main/interim/11_color_code.rds')
color_code_tissue <- setNames(color_code[[1]]$Colors_Tissue, color_code[[1]]$Tissue)
```


```r
set.seed(1)
betas_bmiq <- BMIQ(mset_noob, nfit = 100000)
```

saveRDS(betas_bmiq, '../../data/main/interim/16_bmiq.rds')


```r
betas_bmiq <- readRDS('../../data/main/interim/16_bmiq.rds')
```


```r
#combine snps
snps <- getSnpBeta(rgset)
data <- rbind(betas_bmiq, snps)
dim(data) # 866150    192
```

```
## [1] 866150    192
```

```r
all(pl_ethnicity_features %in% rownames(data)) #T
```

```
## [1] TRUE
```

```r
colnames(data) <- pDat$Sample_Name

results <- pl_infer_ethnicity(data)
```

```
## [1] "1860 of 1860 predictors present."
```

```r
results$Sample_Name <- pDat$Sample_Name

# add to pData
pDat <- pDat %>% left_join(results %>% select(-Highest_Prob, Predicted_ethnicity))
```

```
## Joining, by = "Sample_Name"
```

# Small analysis


```r
# distribution of calls
pDat %>% filter(Tissue == 'Villi') %>% group_by(Trimester) %>% count(Predicted_ethnicity)
```

```
## # A tibble: 7 x 3
## # Groups:   Trimester [3]
##   Trimester Predicted_ethnicity     n
##   <chr>     <chr>               <int>
## 1 First     Ambiguous               5
## 2 First     Caucasian               3
## 3 Second    Asian                   4
## 4 Second    Caucasian              12
## 5 Third     Ambiguous               3
## 6 Third     Asian                  10
## 7 Third     Caucasian              11
```

```r
# distribution of probabilities, for all samples, including cells
ggplot(pDat, aes(x = Prob_Asian, y = Prob_Caucasian, color = Tissue)) +
  geom_point() + theme_bw() +
  facet_wrap(~Trimester) +
  scale_color_manual(values= color_code_tissue[unique(pDat$Tissue)])
```

![](16_ethnicity_ancestry_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
# just villi
pDat %>% 
  filter(Tissue == 'Villi') %>%
  ggplot(aes(x = Prob_Asian, y = Prob_Caucasian)) +
  geom_point() + theme_bw() + labs(title ='Just Villi') +
  facet_wrap(~Trimester)
```

![](16_ethnicity_ancestry_files/figure-html/unnamed-chunk-5-2.png)<!-- -->

```r
# just first trimester, by case
pDat %>% filter(Trimester == 'First') %>% {
  ggplot(data = ., aes(x = Prob_Asian, y = Prob_Caucasian, color = Tissue)) +
  geom_point() + theme_bw() +
  facet_wrap(~Case_ID) +
  scale_color_manual(values= color_code_tissue[unique(.$Tissue)])
}
```

![](16_ethnicity_ancestry_files/figure-html/unnamed-chunk-5-3.png)<!-- -->

```r
# just third trimester, by case
pDat %>% filter(Trimester == 'Third') %>%
  ggplot(aes(x = Prob_Asian, y = Prob_Caucasian, color = Tissue)) +
  geom_point() + theme_bw() +
  facet_wrap(~Case_ID) +
  scale_color_manual(values= color_code_tissue[unique(pDat$Tissue)])
```

![](16_ethnicity_ancestry_files/figure-html/unnamed-chunk-5-4.png)<!-- -->

```r
# comparing cell predictions and villi
p1 <- pDat %>% filter(Trimester == 'Third') %>%
  ggplot(aes(x = Prob_Asian, y = Prob_Caucasian, color = Tissue)) +
  geom_point() + theme_bw() + labs(color = '') +
  scale_color_manual(values= color_code_tissue[unique(pDat$Tissue)]);p1
```

![](16_ethnicity_ancestry_files/figure-html/unnamed-chunk-5-5.png)<!-- -->

```r
p2 <- pDat %>% filter(Trimester == 'Third', Tissue == 'Villi') %>%
  ggplot(aes(x = Prob_Asian, y = Prob_Caucasian, color = Tissue)) +
  geom_point() + theme_bw() + labs(color = '') +
  scale_color_manual(values= color_code_tissue['Villi']);p2
```

![](16_ethnicity_ancestry_files/figure-html/unnamed-chunk-5-6.png)<!-- -->

```r
ggarrange(p1, p2)
```

![](16_ethnicity_ancestry_files/figure-html/unnamed-chunk-5-7.png)<!-- -->

Saving the results for all tissues, in case I want to investigate further later.

For term analysis, I need to use the predicted ethnicity for villi samples only.


```r
pDat_1 <- readRDS('../../data/main/interim/21_pDat_term_cells.rds')

res <- pDat %>% 
  filter(Trimester == 'Third', Tissue == 'Villi') %>% 
  select(Case_ID, contains('ethnicity'), Prob_African, Prob_Asian, Prob_Caucasian)%>%
  # AVERAGE PROBABILITY BETWEEN REPLICATES
  group_by(Case_ID) %>% 
  mutate(Prob_African = mean(Prob_African),
         Prob_Caucasian = mean(Prob_Caucasian),
         Prob_Asian = mean(Prob_Asian)) %>%
  ungroup() %>%
  distinct() # filter out the replicate entries

n1 <- nrow(pDat_1)
pDat_1 <- pDat_1 %>% 
  left_join(res, by = 'Case_ID') 
n2 <- nrow(pDat_1)
n1 == n2 # make sure no new rows are added, Needs to be TRUE
```

```
## [1] TRUE
```

```r
p3 <- ggplot(pDat_1, aes(x = Prob_Asian, y = Prob_Caucasian)) +
  geom_point() + theme_bw() 

ggarrange(p2, p3) # should look the same
```

![](16_ethnicity_ancestry_files/figure-html/unnamed-chunk-6-1.png)<!-- -->


# Save data


```r
# results for all samples
pDat %>% 
  select(contains('Predicted_ethnicity_nothresh'), Prob_African, Prob_Asian, Prob_Caucasian) %>%
  saveRDS(file = '../../data/main/interim/16_ethnicity_results_all.rds')

# results for term, using villi as final predicted ethnicity
pDat_1 %>%
  saveRDS(file = '../../data/main/interim/16_pDat.rds')
```

# SessionInfo


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows Server x64 (build 14393)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252   
## [3] LC_MONETARY=English_Canada.1252 LC_NUMERIC=C                   
## [5] LC_TIME=English_Canada.1252    
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] egg_0.4.2                                         
##  [2] gridExtra_2.3                                     
##  [3] dplyr_0.8.0.1                                     
##  [4] IlluminaHumanMethylationEPICmanifest_0.3.0        
##  [5] wateRmelon_1.26.0                                 
##  [6] illuminaio_0.24.0                                 
##  [7] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0
##  [8] ROC_1.58.0                                        
##  [9] lumi_2.34.0                                       
## [10] methylumi_2.28.0                                  
## [11] FDb.InfiniumMethylation.hg19_2.2.0                
## [12] org.Hs.eg.db_3.7.0                                
## [13] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2           
## [14] GenomicFeatures_1.34.3                            
## [15] AnnotationDbi_1.44.0                              
## [16] ggplot2_3.1.0                                     
## [17] reshape2_1.4.3                                    
## [18] scales_1.0.0                                      
## [19] limma_3.38.3                                      
## [20] minfi_1.28.3                                      
## [21] bumphunter_1.24.5                                 
## [22] locfit_1.5-9.1                                    
## [23] iterators_1.0.10                                  
## [24] foreach_1.4.4                                     
## [25] Biostrings_2.50.2                                 
## [26] XVector_0.22.0                                    
## [27] SummarizedExperiment_1.12.0                       
## [28] DelayedArray_0.8.0                                
## [29] BiocParallel_1.16.6                               
## [30] matrixStats_0.54.0                                
## [31] Biobase_2.42.0                                    
## [32] GenomicRanges_1.34.0                              
## [33] GenomeInfoDb_1.18.2                               
## [34] IRanges_2.16.0                                    
## [35] S4Vectors_0.20.1                                  
## [36] BiocGenerics_0.28.0                               
## [37] planet_0.0.1                                      
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.4-0         siggenes_1.56.0         
##  [3] mclust_5.4.2             base64_2.0              
##  [5] affyio_1.52.0            bit64_0.9-7             
##  [7] fansi_0.4.0              xml2_1.2.0              
##  [9] codetools_0.2-16         splines_3.5.1           
## [11] knitr_1.21               Rsamtools_1.34.1        
## [13] annotate_1.60.0          HDF5Array_1.10.1        
## [15] BiocManager_1.30.4       readr_1.3.1             
## [17] compiler_3.5.1           httr_1.4.0              
## [19] assertthat_0.2.0         Matrix_1.2-15           
## [21] lazyeval_0.2.1           cli_1.0.1               
## [23] htmltools_0.3.6          prettyunits_1.0.2       
## [25] tools_3.5.1              affy_1.60.0             
## [27] gtable_0.2.0             glue_1.3.0              
## [29] GenomeInfoDbData_1.2.0   doRNG_1.7.1             
## [31] Rcpp_1.0.0               multtest_2.38.0         
## [33] preprocessCore_1.44.0    nlme_3.1-137            
## [35] rtracklayer_1.42.2       DelayedMatrixStats_1.4.0
## [37] xfun_0.5                 stringr_1.4.0           
## [39] rngtools_1.3.1           XML_3.98-1.18           
## [41] beanplot_1.2             nleqslv_3.3.2           
## [43] zlibbioc_1.28.0          MASS_7.3-51.1           
## [45] hms_0.4.2                rhdf5_2.26.2            
## [47] GEOquery_2.50.5          RColorBrewer_1.1-2      
## [49] yaml_2.2.0               memoise_1.1.0           
## [51] pkgmaker_0.27            biomaRt_2.38.0          
## [53] reshape_0.8.8            stringi_1.3.1           
## [55] RSQLite_2.1.1            genefilter_1.64.0       
## [57] bibtex_0.4.2             rlang_0.3.1             
## [59] pkgconfig_2.0.2          bitops_1.0-6            
## [61] nor1mix_1.2-3            evaluate_0.13           
## [63] lattice_0.20-38          purrr_0.3.0             
## [65] Rhdf5lib_1.4.2           labeling_0.3            
## [67] GenomicAlignments_1.18.1 bit_1.1-14              
## [69] tidyselect_0.2.5         plyr_1.8.4              
## [71] magrittr_1.5             R6_2.4.0                
## [73] DBI_1.0.0                pillar_1.3.1            
## [75] withr_2.1.2              mgcv_1.8-27             
## [77] survival_2.43-3          RCurl_1.95-4.12         
## [79] tibble_2.0.1             crayon_1.3.4            
## [81] utf8_1.1.4               KernSmooth_2.23-15      
## [83] rmarkdown_1.11           progress_1.2.0          
## [85] grid_3.5.1               data.table_1.12.0       
## [87] blob_1.1.1               digest_0.6.18           
## [89] xtable_1.8-3             tidyr_0.8.3             
## [91] glmnet_2.0-16            openssl_1.2.2           
## [93] munsell_0.5.0            registry_0.5-1          
## [95] askpass_1.1              quadprog_1.5-5
```
