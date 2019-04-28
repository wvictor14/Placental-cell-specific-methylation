---
title: "01_Load_IDATs"
author: "Victor Yuan"
date: "April 27, 2019"
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

Load week 1 samples (n=96, 12 chips)

# 1.0 Libraries


```r
library(tidyr)
library(plyr)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:plyr':
## 
##     arrange, count, desc, failwith, id, mutate, rename, summarise,
##     summarize
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(minfi)
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind,
##     colMeans, colnames, colSums, dirname, do.call, duplicated,
##     eval, evalq, Filter, Find, get, grep, grepl, intersect,
##     is.unsorted, lapply, Map, mapply, match, mget, order, paste,
##     pmax, pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
##     tapply, union, unique, unsplit, which, which.max, which.min
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: stats4
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
```

```
## The following object is masked from 'package:plyr':
## 
##     rename
```

```
## The following object is masked from 'package:tidyr':
## 
##     expand
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
```

```
## The following object is masked from 'package:plyr':
## 
##     desc
```

```
## The following object is masked from 'package:grDevices':
## 
##     windows
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: DelayedArray
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```
## The following object is masked from 'package:dplyr':
## 
##     count
```

```
## The following object is masked from 'package:plyr':
## 
##     count
```

```
## Loading required package: BiocParallel
```

```
## 
## Attaching package: 'DelayedArray'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
```

```
## The following objects are masked from 'package:base':
## 
##     aperm, apply, rowsum
```

```
## Loading required package: Biostrings
```

```
## Loading required package: XVector
```

```
## 
## Attaching package: 'XVector'
```

```
## The following object is masked from 'package:plyr':
## 
##     compact
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:DelayedArray':
## 
##     type
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## Loading required package: bumphunter
```

```
## Loading required package: foreach
```

```
## Loading required package: iterators
```

```
## Loading required package: locfit
```

```
## locfit 1.5-9.1 	 2013-03-22
```

```
## Registered S3 method overwritten by 'openssl':
##   method      from
##   print.bytes Rcpp
```

```
## Setting options('download.file.method.GEOquery'='auto')
```

```
## Setting options('GEOquery.inmemory.gpl'=FALSE)
```

```r
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
```

# 2.0 Data


```r
ss <- read.csv('../../data/main/raw/NIH EPIC Batch 2 week 1 sample sheet.csv', skip = 7)

rgset <- read.metharray.exp('Z:/ROBLAB6 Infinium450k John/EPIC Raw data/NIH EPIC Batch 2/IDATs - Week 1/', 
                            recursive = T, verbose = T, extended = T)
```

```
## [read.metharray] Reading 203067920143_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203067920143_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203067920143_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203067920143_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203067920143_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203067920143_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203067920143_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203067920143_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203067920144_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203067920144_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203067920144_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203067920144_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203067920144_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203067920144_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203067920144_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203067920144_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203068760084_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203068760084_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203068760084_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203068760084_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203068760084_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203068760084_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203068760084_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203068760084_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203068760095_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203068760095_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203068760095_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203068760095_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203068760095_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203068760095_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203068760095_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203068760095_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203068760191_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203068760191_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203068760191_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203068760191_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203068760191_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203068760191_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203068760191_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203068760191_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203068760192_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203068760192_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203068760192_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203068760192_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203068760192_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203068760192_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203068760192_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203068760192_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203068760202_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203068760202_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203068760202_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203068760202_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203068760202_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203068760202_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203068760202_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203068760202_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203068760204_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203068760204_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203068760204_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203068760204_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203068760204_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203068760204_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203068760204_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203068760204_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203072530027_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203072530027_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203072530027_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203072530027_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203072530027_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203072530027_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203072530027_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203072530027_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203072530046_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203072530046_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203072530046_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203072530046_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203072530046_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203072530046_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203072530046_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203072530046_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203072530098_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203072530098_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203072530098_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203072530098_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203072530098_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203072530098_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203072530098_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203072530098_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203072530157_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203072530157_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203072530157_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203072530157_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203072530157_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203072530157_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203072530157_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203072530157_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203067920143_R01C01_Red.idat
```

```
## [read.metharray] Reading 203067920143_R02C01_Red.idat
```

```
## [read.metharray] Reading 203067920143_R03C01_Red.idat
```

```
## [read.metharray] Reading 203067920143_R04C01_Red.idat
```

```
## [read.metharray] Reading 203067920143_R05C01_Red.idat
```

```
## [read.metharray] Reading 203067920143_R06C01_Red.idat
```

```
## [read.metharray] Reading 203067920143_R07C01_Red.idat
```

```
## [read.metharray] Reading 203067920143_R08C01_Red.idat
```

```
## [read.metharray] Reading 203067920144_R01C01_Red.idat
```

```
## [read.metharray] Reading 203067920144_R02C01_Red.idat
```

```
## [read.metharray] Reading 203067920144_R03C01_Red.idat
```

```
## [read.metharray] Reading 203067920144_R04C01_Red.idat
```

```
## [read.metharray] Reading 203067920144_R05C01_Red.idat
```

```
## [read.metharray] Reading 203067920144_R06C01_Red.idat
```

```
## [read.metharray] Reading 203067920144_R07C01_Red.idat
```

```
## [read.metharray] Reading 203067920144_R08C01_Red.idat
```

```
## [read.metharray] Reading 203068760084_R01C01_Red.idat
```

```
## [read.metharray] Reading 203068760084_R02C01_Red.idat
```

```
## [read.metharray] Reading 203068760084_R03C01_Red.idat
```

```
## [read.metharray] Reading 203068760084_R04C01_Red.idat
```

```
## [read.metharray] Reading 203068760084_R05C01_Red.idat
```

```
## [read.metharray] Reading 203068760084_R06C01_Red.idat
```

```
## [read.metharray] Reading 203068760084_R07C01_Red.idat
```

```
## [read.metharray] Reading 203068760084_R08C01_Red.idat
```

```
## [read.metharray] Reading 203068760095_R01C01_Red.idat
```

```
## [read.metharray] Reading 203068760095_R02C01_Red.idat
```

```
## [read.metharray] Reading 203068760095_R03C01_Red.idat
```

```
## [read.metharray] Reading 203068760095_R04C01_Red.idat
```

```
## [read.metharray] Reading 203068760095_R05C01_Red.idat
```

```
## [read.metharray] Reading 203068760095_R06C01_Red.idat
```

```
## [read.metharray] Reading 203068760095_R07C01_Red.idat
```

```
## [read.metharray] Reading 203068760095_R08C01_Red.idat
```

```
## [read.metharray] Reading 203068760191_R01C01_Red.idat
```

```
## [read.metharray] Reading 203068760191_R02C01_Red.idat
```

```
## [read.metharray] Reading 203068760191_R03C01_Red.idat
```

```
## [read.metharray] Reading 203068760191_R04C01_Red.idat
```

```
## [read.metharray] Reading 203068760191_R05C01_Red.idat
```

```
## [read.metharray] Reading 203068760191_R06C01_Red.idat
```

```
## [read.metharray] Reading 203068760191_R07C01_Red.idat
```

```
## [read.metharray] Reading 203068760191_R08C01_Red.idat
```

```
## [read.metharray] Reading 203068760192_R01C01_Red.idat
```

```
## [read.metharray] Reading 203068760192_R02C01_Red.idat
```

```
## [read.metharray] Reading 203068760192_R03C01_Red.idat
```

```
## [read.metharray] Reading 203068760192_R04C01_Red.idat
```

```
## [read.metharray] Reading 203068760192_R05C01_Red.idat
```

```
## [read.metharray] Reading 203068760192_R06C01_Red.idat
```

```
## [read.metharray] Reading 203068760192_R07C01_Red.idat
```

```
## [read.metharray] Reading 203068760192_R08C01_Red.idat
```

```
## [read.metharray] Reading 203068760202_R01C01_Red.idat
```

```
## [read.metharray] Reading 203068760202_R02C01_Red.idat
```

```
## [read.metharray] Reading 203068760202_R03C01_Red.idat
```

```
## [read.metharray] Reading 203068760202_R04C01_Red.idat
```

```
## [read.metharray] Reading 203068760202_R05C01_Red.idat
```

```
## [read.metharray] Reading 203068760202_R06C01_Red.idat
```

```
## [read.metharray] Reading 203068760202_R07C01_Red.idat
```

```
## [read.metharray] Reading 203068760202_R08C01_Red.idat
```

```
## [read.metharray] Reading 203068760204_R01C01_Red.idat
```

```
## [read.metharray] Reading 203068760204_R02C01_Red.idat
```

```
## [read.metharray] Reading 203068760204_R03C01_Red.idat
```

```
## [read.metharray] Reading 203068760204_R04C01_Red.idat
```

```
## [read.metharray] Reading 203068760204_R05C01_Red.idat
```

```
## [read.metharray] Reading 203068760204_R06C01_Red.idat
```

```
## [read.metharray] Reading 203068760204_R07C01_Red.idat
```

```
## [read.metharray] Reading 203068760204_R08C01_Red.idat
```

```
## [read.metharray] Reading 203072530027_R01C01_Red.idat
```

```
## [read.metharray] Reading 203072530027_R02C01_Red.idat
```

```
## [read.metharray] Reading 203072530027_R03C01_Red.idat
```

```
## [read.metharray] Reading 203072530027_R04C01_Red.idat
```

```
## [read.metharray] Reading 203072530027_R05C01_Red.idat
```

```
## [read.metharray] Reading 203072530027_R06C01_Red.idat
```

```
## [read.metharray] Reading 203072530027_R07C01_Red.idat
```

```
## [read.metharray] Reading 203072530027_R08C01_Red.idat
```

```
## [read.metharray] Reading 203072530046_R01C01_Red.idat
```

```
## [read.metharray] Reading 203072530046_R02C01_Red.idat
```

```
## [read.metharray] Reading 203072530046_R03C01_Red.idat
```

```
## [read.metharray] Reading 203072530046_R04C01_Red.idat
```

```
## [read.metharray] Reading 203072530046_R05C01_Red.idat
```

```
## [read.metharray] Reading 203072530046_R06C01_Red.idat
```

```
## [read.metharray] Reading 203072530046_R07C01_Red.idat
```

```
## [read.metharray] Reading 203072530046_R08C01_Red.idat
```

```
## [read.metharray] Reading 203072530098_R01C01_Red.idat
```

```
## [read.metharray] Reading 203072530098_R02C01_Red.idat
```

```
## [read.metharray] Reading 203072530098_R03C01_Red.idat
```

```
## [read.metharray] Reading 203072530098_R04C01_Red.idat
```

```
## [read.metharray] Reading 203072530098_R05C01_Red.idat
```

```
## [read.metharray] Reading 203072530098_R06C01_Red.idat
```

```
## [read.metharray] Reading 203072530098_R07C01_Red.idat
```

```
## [read.metharray] Reading 203072530098_R08C01_Red.idat
```

```
## [read.metharray] Reading 203072530157_R01C01_Red.idat
```

```
## [read.metharray] Reading 203072530157_R02C01_Red.idat
```

```
## [read.metharray] Reading 203072530157_R03C01_Red.idat
```

```
## [read.metharray] Reading 203072530157_R04C01_Red.idat
```

```
## [read.metharray] Reading 203072530157_R05C01_Red.idat
```

```
## [read.metharray] Reading 203072530157_R06C01_Red.idat
```

```
## [read.metharray] Reading 203072530157_R07C01_Red.idat
```

```
## [read.metharray] Reading 203072530157_R08C01_Red.idat
```

```
## [read.metharray] Read idat files in 53.8 seconds
```

```
## [read.metharray] Creating data matrices ... done in 243.5 seconds
## [read.metharray] Instantiating final object ... done in 0.2 seconds
```

```r
rgset
```

```
## class: RGChannelSetExtended 
## dim: 1051815 96 
## metadata(0):
## assays(5): Green Red GreenSD RedSD NBeads
## rownames(1051815): 1600101 1600111 ... 99810990 99810992
## rowData names(0):
## colnames(96): 203067920143_R01C01 203067920143_R02C01 ...
##   203072530157_R07C01 203072530157_R08C01
## colData names(0):
## Annotation
##   array: IlluminaHumanMethylationEPIC
##   annotation: ilm10b4.hg19
```


```r
saveRDS(rgset, '../../data/main/interim/01_rgset_raw.rds')
```

# Session Info


```r
sessionInfo()
```

```
## R version 3.6.0 (2019-04-26)
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
##  [1] IlluminaHumanMethylationEPICmanifest_0.3.0         
##  [2] IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0
##  [3] minfi_1.29.4                                       
##  [4] bumphunter_1.25.4                                  
##  [5] locfit_1.5-9.1                                     
##  [6] iterators_1.0.10                                   
##  [7] foreach_1.4.4                                      
##  [8] Biostrings_2.51.5                                  
##  [9] XVector_0.23.2                                     
## [10] SummarizedExperiment_1.13.0                        
## [11] DelayedArray_0.9.9                                 
## [12] BiocParallel_1.17.18                               
## [13] matrixStats_0.54.0                                 
## [14] Biobase_2.43.1                                     
## [15] GenomicRanges_1.35.1                               
## [16] GenomeInfoDb_1.19.3                                
## [17] IRanges_2.17.5                                     
## [18] S4Vectors_0.21.24                                  
## [19] BiocGenerics_0.29.2                                
## [20] dplyr_0.8.0.1                                      
## [21] plyr_1.8.4                                         
## [22] tidyr_0.8.3                                        
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-139             bitops_1.0-6            
##  [3] bit64_0.9-7              RColorBrewer_1.1-2      
##  [5] progress_1.2.0           httr_1.4.0              
##  [7] tools_3.6.0              doRNG_1.7.1             
##  [9] nor1mix_1.2-3            R6_2.4.0                
## [11] HDF5Array_1.11.1         DBI_1.0.0               
## [13] withr_2.1.2              tidyselect_0.2.5        
## [15] prettyunits_1.0.2        base64_2.0              
## [17] bit_1.1-14               compiler_3.6.0          
## [19] preprocessCore_1.45.0    xml2_1.2.0              
## [21] pkgmaker_0.27            rtracklayer_1.43.4      
## [23] readr_1.3.1              quadprog_1.5-6          
## [25] genefilter_1.65.0        askpass_1.1             
## [27] stringr_1.4.0            digest_0.6.18           
## [29] Rsamtools_1.99.6         illuminaio_0.25.0       
## [31] rmarkdown_1.12           siggenes_1.57.4         
## [33] GEOquery_2.51.6          pkgconfig_2.0.2         
## [35] htmltools_0.3.6          scrime_1.3.5            
## [37] bibtex_0.4.2             limma_3.39.18           
## [39] rlang_0.3.4              RSQLite_2.1.1           
## [41] DelayedMatrixStats_1.5.2 mclust_5.4.3            
## [43] RCurl_1.95-4.12          magrittr_1.5            
## [45] GenomeInfoDbData_1.2.1   Matrix_1.2-17           
## [47] Rcpp_1.0.1               Rhdf5lib_1.5.4          
## [49] stringi_1.4.3            yaml_2.2.0              
## [51] MASS_7.3-51.4            zlibbioc_1.29.0         
## [53] rhdf5_2.27.19            grid_3.6.0              
## [55] blob_1.1.1               crayon_1.3.4            
## [57] lattice_0.20-38          splines_3.6.0           
## [59] annotate_1.61.1          multtest_2.39.0         
## [61] GenomicFeatures_1.35.11  hms_0.4.2               
## [63] knitr_1.22               beanplot_1.2            
## [65] pillar_1.3.1             rngtools_1.3.1.1        
## [67] codetools_0.2-16         biomaRt_2.39.4          
## [69] XML_3.98-1.19            glue_1.3.1              
## [71] evaluate_0.13            data.table_1.12.2       
## [73] openssl_1.3              purrr_0.3.2             
## [75] reshape_0.8.8            assertthat_0.2.1        
## [77] xfun_0.6                 xtable_1.8-4            
## [79] survival_2.44-1.1        tibble_2.1.1            
## [81] GenomicAlignments_1.19.1 AnnotationDbi_1.45.1    
## [83] registry_0.5-1           memoise_1.1.0
```
