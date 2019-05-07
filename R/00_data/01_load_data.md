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

Load week 1&2 samples (n=192, 24 chips)

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
library(readxl)
```

# 2.0 Data


```r
#sample phenodata 
ss <- read_xlsx('../../data/main/raw/NIH EPIC Batch 2 sample sheet.xlsx', skip = 7)
ss
```

```
## # A tibble: 192 x 19
##    Sample_Name `Chip Number`   Row Well  Case_ID Sex      GA Trimester
##    <chr>               <dbl> <dbl> <chr> <chr>   <chr> <dbl> <chr>    
##  1 PM365_endo~            25     1 A1    PM365   F      NA   Third    
##  2 PL295_disc~            25     2 B1    PL295   F      11.6 First    
##  3 PL293_stro~            25     3 C1    PL293   M       9.6 First    
##  4 PL292_endo~            25     4 D1    PL292   M       7.1 First    
##  5 FT73_v                 25     5 E1    FT73    F      NA   Second   
##  6 PM376_vc               25     6 F1    PM376   F      NA   Third    
##  7 PM370_hofb~            25     7 G1    PM370   F      NA   Third    
##  8 PM368_trop~            25     8 H1    PM368   M      NA   Third    
##  9 PM368_vc               26     1 A2    PM368   M      NA   Third    
## 10 PL293_hofb~            26     2 B2    PL293   M       9.6 First    
## # ... with 182 more rows, and 11 more variables: DNA_QP <chr>, Week <dbl>,
## #   Sample_Plate <chr>, Tissue <chr>, Sentrix_ID <dbl>,
## #   Sentrix_Position <chr>, `Scratches on beadchip` <chr>,
## #   Batch_BSC <chr>, `BSC sample adjustedconcentration` <dbl>, `DNA
## #   concentration after speed vac and dilution (ng/ul)` <dbl>, `Amount of
## #   DNA loaded (ng)` <dbl>
```

```r
# load idats
rgset <- read.metharray.exp('Z:/ROBLAB6 Infinium450k John/EPIC Raw data/NIH EPIC Batch 2/', 
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
## [read.metharray] Reading 203072530204_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203072530204_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203072530204_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203072530204_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203072530204_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203072530204_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203072530204_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203072530204_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203072530227_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203072530227_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203072530227_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203072530227_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203072530227_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203072530227_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203072530227_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203072530227_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203077630004_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203077630004_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203077630004_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203077630004_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203077630004_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203077630004_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203077630004_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203077630004_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203077630012_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203077630012_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203077630012_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203077630012_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203077630012_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203077630012_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203077630012_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203077630012_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203077630052_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203077630052_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203077630052_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203077630052_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203077630052_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203077630052_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203077630052_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203077630052_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203077630080_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203077630080_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203077630080_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203077630080_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203077630080_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203077630080_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203077630080_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203077630080_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203077630103_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203077630103_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203077630103_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203077630103_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203077630103_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203077630103_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203077630103_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203077630103_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203077630104_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203077630104_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203077630104_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203077630104_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203077630104_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203077630104_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203077630104_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203077630104_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203077630172_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203077630172_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203077630172_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203077630172_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203077630172_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203077630172_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203077630172_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203077630172_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203077630181_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203077630181_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203077630181_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203077630181_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203077630181_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203077630181_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203077630181_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203077630181_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203077630199_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203077630199_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203077630199_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203077630199_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203077630199_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203077630199_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203077630199_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203077630199_R08C01_Grn.idat
```

```
## [read.metharray] Reading 203077630206_R01C01_Grn.idat
```

```
## [read.metharray] Reading 203077630206_R02C01_Grn.idat
```

```
## [read.metharray] Reading 203077630206_R03C01_Grn.idat
```

```
## [read.metharray] Reading 203077630206_R04C01_Grn.idat
```

```
## [read.metharray] Reading 203077630206_R05C01_Grn.idat
```

```
## [read.metharray] Reading 203077630206_R06C01_Grn.idat
```

```
## [read.metharray] Reading 203077630206_R07C01_Grn.idat
```

```
## [read.metharray] Reading 203077630206_R08C01_Grn.idat
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
## [read.metharray] Reading 203072530204_R01C01_Red.idat
```

```
## [read.metharray] Reading 203072530204_R02C01_Red.idat
```

```
## [read.metharray] Reading 203072530204_R03C01_Red.idat
```

```
## [read.metharray] Reading 203072530204_R04C01_Red.idat
```

```
## [read.metharray] Reading 203072530204_R05C01_Red.idat
```

```
## [read.metharray] Reading 203072530204_R06C01_Red.idat
```

```
## [read.metharray] Reading 203072530204_R07C01_Red.idat
```

```
## [read.metharray] Reading 203072530204_R08C01_Red.idat
```

```
## [read.metharray] Reading 203072530227_R01C01_Red.idat
```

```
## [read.metharray] Reading 203072530227_R02C01_Red.idat
```

```
## [read.metharray] Reading 203072530227_R03C01_Red.idat
```

```
## [read.metharray] Reading 203072530227_R04C01_Red.idat
```

```
## [read.metharray] Reading 203072530227_R05C01_Red.idat
```

```
## [read.metharray] Reading 203072530227_R06C01_Red.idat
```

```
## [read.metharray] Reading 203072530227_R07C01_Red.idat
```

```
## [read.metharray] Reading 203072530227_R08C01_Red.idat
```

```
## [read.metharray] Reading 203077630004_R01C01_Red.idat
```

```
## [read.metharray] Reading 203077630004_R02C01_Red.idat
```

```
## [read.metharray] Reading 203077630004_R03C01_Red.idat
```

```
## [read.metharray] Reading 203077630004_R04C01_Red.idat
```

```
## [read.metharray] Reading 203077630004_R05C01_Red.idat
```

```
## [read.metharray] Reading 203077630004_R06C01_Red.idat
```

```
## [read.metharray] Reading 203077630004_R07C01_Red.idat
```

```
## [read.metharray] Reading 203077630004_R08C01_Red.idat
```

```
## [read.metharray] Reading 203077630012_R01C01_Red.idat
```

```
## [read.metharray] Reading 203077630012_R02C01_Red.idat
```

```
## [read.metharray] Reading 203077630012_R03C01_Red.idat
```

```
## [read.metharray] Reading 203077630012_R04C01_Red.idat
```

```
## [read.metharray] Reading 203077630012_R05C01_Red.idat
```

```
## [read.metharray] Reading 203077630012_R06C01_Red.idat
```

```
## [read.metharray] Reading 203077630012_R07C01_Red.idat
```

```
## [read.metharray] Reading 203077630012_R08C01_Red.idat
```

```
## [read.metharray] Reading 203077630052_R01C01_Red.idat
```

```
## [read.metharray] Reading 203077630052_R02C01_Red.idat
```

```
## [read.metharray] Reading 203077630052_R03C01_Red.idat
```

```
## [read.metharray] Reading 203077630052_R04C01_Red.idat
```

```
## [read.metharray] Reading 203077630052_R05C01_Red.idat
```

```
## [read.metharray] Reading 203077630052_R06C01_Red.idat
```

```
## [read.metharray] Reading 203077630052_R07C01_Red.idat
```

```
## [read.metharray] Reading 203077630052_R08C01_Red.idat
```

```
## [read.metharray] Reading 203077630080_R01C01_Red.idat
```

```
## [read.metharray] Reading 203077630080_R02C01_Red.idat
```

```
## [read.metharray] Reading 203077630080_R03C01_Red.idat
```

```
## [read.metharray] Reading 203077630080_R04C01_Red.idat
```

```
## [read.metharray] Reading 203077630080_R05C01_Red.idat
```

```
## [read.metharray] Reading 203077630080_R06C01_Red.idat
```

```
## [read.metharray] Reading 203077630080_R07C01_Red.idat
```

```
## [read.metharray] Reading 203077630080_R08C01_Red.idat
```

```
## [read.metharray] Reading 203077630103_R01C01_Red.idat
```

```
## [read.metharray] Reading 203077630103_R02C01_Red.idat
```

```
## [read.metharray] Reading 203077630103_R03C01_Red.idat
```

```
## [read.metharray] Reading 203077630103_R04C01_Red.idat
```

```
## [read.metharray] Reading 203077630103_R05C01_Red.idat
```

```
## [read.metharray] Reading 203077630103_R06C01_Red.idat
```

```
## [read.metharray] Reading 203077630103_R07C01_Red.idat
```

```
## [read.metharray] Reading 203077630103_R08C01_Red.idat
```

```
## [read.metharray] Reading 203077630104_R01C01_Red.idat
```

```
## [read.metharray] Reading 203077630104_R02C01_Red.idat
```

```
## [read.metharray] Reading 203077630104_R03C01_Red.idat
```

```
## [read.metharray] Reading 203077630104_R04C01_Red.idat
```

```
## [read.metharray] Reading 203077630104_R05C01_Red.idat
```

```
## [read.metharray] Reading 203077630104_R06C01_Red.idat
```

```
## [read.metharray] Reading 203077630104_R07C01_Red.idat
```

```
## [read.metharray] Reading 203077630104_R08C01_Red.idat
```

```
## [read.metharray] Reading 203077630172_R01C01_Red.idat
```

```
## [read.metharray] Reading 203077630172_R02C01_Red.idat
```

```
## [read.metharray] Reading 203077630172_R03C01_Red.idat
```

```
## [read.metharray] Reading 203077630172_R04C01_Red.idat
```

```
## [read.metharray] Reading 203077630172_R05C01_Red.idat
```

```
## [read.metharray] Reading 203077630172_R06C01_Red.idat
```

```
## [read.metharray] Reading 203077630172_R07C01_Red.idat
```

```
## [read.metharray] Reading 203077630172_R08C01_Red.idat
```

```
## [read.metharray] Reading 203077630181_R01C01_Red.idat
```

```
## [read.metharray] Reading 203077630181_R02C01_Red.idat
```

```
## [read.metharray] Reading 203077630181_R03C01_Red.idat
```

```
## [read.metharray] Reading 203077630181_R04C01_Red.idat
```

```
## [read.metharray] Reading 203077630181_R05C01_Red.idat
```

```
## [read.metharray] Reading 203077630181_R06C01_Red.idat
```

```
## [read.metharray] Reading 203077630181_R07C01_Red.idat
```

```
## [read.metharray] Reading 203077630181_R08C01_Red.idat
```

```
## [read.metharray] Reading 203077630199_R01C01_Red.idat
```

```
## [read.metharray] Reading 203077630199_R02C01_Red.idat
```

```
## [read.metharray] Reading 203077630199_R03C01_Red.idat
```

```
## [read.metharray] Reading 203077630199_R04C01_Red.idat
```

```
## [read.metharray] Reading 203077630199_R05C01_Red.idat
```

```
## [read.metharray] Reading 203077630199_R06C01_Red.idat
```

```
## [read.metharray] Reading 203077630199_R07C01_Red.idat
```

```
## [read.metharray] Reading 203077630199_R08C01_Red.idat
```

```
## [read.metharray] Reading 203077630206_R01C01_Red.idat
```

```
## [read.metharray] Reading 203077630206_R02C01_Red.idat
```

```
## [read.metharray] Reading 203077630206_R03C01_Red.idat
```

```
## [read.metharray] Reading 203077630206_R04C01_Red.idat
```

```
## [read.metharray] Reading 203077630206_R05C01_Red.idat
```

```
## [read.metharray] Reading 203077630206_R06C01_Red.idat
```

```
## [read.metharray] Reading 203077630206_R07C01_Red.idat
```

```
## [read.metharray] Reading 203077630206_R08C01_Red.idat
```

```
## [read.metharray] Read idat files in 113.2 seconds
```

```
## [read.metharray] Creating data matrices ... done in 509.3 seconds
## [read.metharray] Instantiating final object ... done in 0.2 seconds
```

```r
rgset
```

```
## class: RGChannelSetExtended 
## dim: 1051815 192 
## metadata(0):
## assays(5): Green Red GreenSD RedSD NBeads
## rownames(1051815): 1600101 1600111 ... 99810990 99810992
## rowData names(0):
## colnames(192): 203067920143_R01C01 203067920143_R02C01 ...
##   203077630206_R07C01 203077630206_R08C01
## colData names(0):
## Annotation
##   array: IlluminaHumanMethylationEPIC
##   annotation: ilm10b4.hg19
```

```r
# create sentrix column which matches column names of rgset
ss <- ss %>% 
  mutate(Sentrix = paste0(Sentrix_ID, '_', Sentrix_Position)) %>%
  dplyr::rename(Chip_number = `Chip Number`,
                Scratches = `Scratches on beadchip`,
                DNA_conc_BSC_adjusted = `BSC sample adjustedconcentration`,
                DNA_conc_before_load = `DNA concentration after speed vac and dilution (ng/ul)`,
                DNA_loaded = `Amount of DNA loaded (ng)`) %>%
  arrange(Sentrix)

sum(ss$Sentrix %in% colnames(rgset)) #192
```

```
## [1] 192
```

```r
#rearrange argset based on ss
rgset <- rgset[,ss$Sentrix]

# add pData to rgset
ss <- ss %>% as.data.frame()
rownames(ss) <- ss$Sentrix
pData(rgset) <- DataFrame(ss)
pData(rgset)
```

```
## DataFrame with 192 rows and 20 columns
##                          Sample_Name Chip_number       Row        Well
##                          <character>   <numeric> <numeric> <character>
## 203067920143_R01C01    PM365_endo_cs          25         1          A1
## 203067920143_R02C01 PL295_discard_cs          25         2          B1
## 203067920143_R03C01   PL293_strom_cs          25         3          C1
## 203067920143_R04C01    PL292_endo_cs          25         4          D1
## 203067920143_R05C01           FT73_v          25         5          E1
## ...                              ...         ...       ...         ...
## 203077630206_R04C01         PL149_vc          48         4         D12
## 203077630206_R05C01      PM366_vc_R2          48         5         E12
## 203077630206_R06C01          PL295_v          48         6         F12
## 203077630206_R07C01 PM372_discard_cs          48         7         G12
## 203077630206_R08C01    PM374_hofb_cs          48         8         H12
##                         Case_ID         Sex        GA   Trimester
##                     <character> <character> <numeric> <character>
## 203067920143_R01C01       PM365           F        NA       Third
## 203067920143_R02C01       PL295           F      11.6       First
## 203067920143_R03C01       PL293           M       9.6       First
## 203067920143_R04C01       PL292           M       7.1       First
## 203067920143_R05C01        FT73           F        NA      Second
## ...                         ...         ...       ...         ...
## 203077630206_R04C01       PL149           F        NA      Second
## 203077630206_R05C01       PM366           F        NA       Third
## 203077630206_R06C01       PL295           F      11.6       First
## 203077630206_R07C01       PM372           M        NA       Third
## 203077630206_R08C01       PM374           M        NA       Third
##                          DNA_QP      Week   Sample_Plate
##                     <character> <numeric>    <character>
## 203067920143_R01C01       542.5         1 WG6980707-MSA4
## 203067920143_R02C01       701.5         1 WG6980707-MSA4
## 203067920143_R03C01       412.5         1 WG6980707-MSA4
## 203067920143_R04C01       151.5         1 WG6980707-MSA4
## 203067920143_R05C01          NA         1 WG6980707-MSA4
## ...                         ...       ...            ...
## 203077630206_R04C01          NA         2 WG6980638-MSA4
## 203077630206_R05C01          NA         2 WG6980638-MSA4
## 203077630206_R06C01          NA         2 WG6980638-MSA4
## 203077630206_R07C01        1261         2 WG6980638-MSA4
## 203077630206_R08C01         419         2 WG6980638-MSA4
##                                         Tissue   Sentrix_ID
##                                    <character>    <numeric>
## 203067920143_R01C01                Endothelial 203067920143
## 203067920143_R02C01 Dead Cells and Lymphocytes 203067920143
## 203067920143_R03C01                    Stromal 203067920143
## 203067920143_R04C01                Endothelial 203067920143
## 203067920143_R05C01                      Villi 203067920143
## ...                                        ...          ...
## 203077630206_R04C01                      Villi 203077630206
## 203077630206_R05C01                      Villi 203077630206
## 203077630206_R06C01                      Villi 203077630206
## 203077630206_R07C01 Dead Cells and Lymphocytes 203077630206
## 203077630206_R08C01                   Hofbauer 203077630206
##                     Sentrix_Position   Scratches   Batch_BSC
##                          <character> <character> <character>
## 203067920143_R01C01           R01C01          NA        BSC1
## 203067920143_R02C01           R02C01          NA        BSC1
## 203067920143_R03C01           R03C01          NA        BSC1
## 203067920143_R04C01           R04C01          NA        BSC1
## 203067920143_R05C01           R05C01          NA        BSC1
## ...                              ...         ...         ...
## 203077630206_R04C01           R04C01          NA        BSC3
## 203077630206_R05C01           R05C01          NA        BSC3
## 203077630206_R06C01           R06C01          NA    BSC_TEST
## 203077630206_R07C01           R07C01          NA        BSC3
## 203077630206_R08C01           R08C01          NA        BSC3
##                     DNA_conc_BSC_adjusted DNA_conc_before_load
##                                 <numeric>            <numeric>
## 203067920143_R01C01                 30.19     50.3166666666667
## 203067920143_R02C01                 42.65                42.65
## 203067920143_R03C01                 20.99               52.475
## 203067920143_R04C01                  8.04                 20.1
## 203067920143_R05C01                 49.91                49.91
## ...                                   ...                  ...
## 203077630206_R04C01               53.1225                   50
## 203077630206_R05C01                    47                   47
## 203077630206_R06C01                 54.54                   50
## 203077630206_R07C01                77.925                   50
## 203077630206_R08C01                20.145              50.3625
##                           DNA_loaded             Sentrix
##                            <numeric>         <character>
## 203067920143_R01C01 201.266666666667 203067920143_R01C01
## 203067920143_R02C01            170.6 203067920143_R02C01
## 203067920143_R03C01            209.9 203067920143_R03C01
## 203067920143_R04C01             80.4 203067920143_R04C01
## 203067920143_R05C01           199.64 203067920143_R05C01
## ...                              ...                 ...
## 203077630206_R04C01              200 203077630206_R04C01
## 203077630206_R05C01              188 203077630206_R05C01
## 203077630206_R06C01              200 203077630206_R06C01
## 203077630206_R07C01              200 203077630206_R07C01
## 203077630206_R08C01           201.45 203077630206_R08C01
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
##  [1] readxl_1.3.1                                       
##  [2] IlluminaHumanMethylationEPICmanifest_0.3.0         
##  [3] IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0
##  [4] minfi_1.30.0                                       
##  [5] bumphunter_1.25.4                                  
##  [6] locfit_1.5-9.1                                     
##  [7] iterators_1.0.10                                   
##  [8] foreach_1.4.4                                      
##  [9] Biostrings_2.51.5                                  
## [10] XVector_0.23.2                                     
## [11] SummarizedExperiment_1.13.0                        
## [12] DelayedArray_0.9.9                                 
## [13] BiocParallel_1.17.18                               
## [14] matrixStats_0.54.0                                 
## [15] Biobase_2.43.1                                     
## [16] GenomicRanges_1.36.0                               
## [17] GenomeInfoDb_1.19.3                                
## [18] IRanges_2.17.6                                     
## [19] S4Vectors_0.21.24                                  
## [20] BiocGenerics_0.29.2                                
## [21] dplyr_0.8.0.1                                      
## [22] plyr_1.8.4                                         
## [23] tidyr_0.8.3                                        
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-139             bitops_1.0-6            
##  [3] bit64_0.9-7              RColorBrewer_1.1-2      
##  [5] progress_1.2.0           httr_1.4.0              
##  [7] tools_3.6.0              doRNG_1.7.1             
##  [9] nor1mix_1.2-3            utf8_1.1.4              
## [11] R6_2.4.0                 HDF5Array_1.11.1        
## [13] DBI_1.0.0                withr_2.1.2             
## [15] tidyselect_0.2.5         prettyunits_1.0.2       
## [17] base64_2.0               bit_1.1-14              
## [19] compiler_3.6.0           preprocessCore_1.46.0   
## [21] cli_1.1.0                xml2_1.2.0              
## [23] pkgmaker_0.27            rtracklayer_1.44.0      
## [25] readr_1.3.1              quadprog_1.5-6          
## [27] genefilter_1.66.0        askpass_1.1             
## [29] stringr_1.4.0            digest_0.6.18           
## [31] Rsamtools_2.0.0          illuminaio_0.26.0       
## [33] rmarkdown_1.12           siggenes_1.58.0         
## [35] GEOquery_2.52.0          pkgconfig_2.0.2         
## [37] htmltools_0.3.6          scrime_1.3.5            
## [39] bibtex_0.4.2             limma_3.40.0            
## [41] rlang_0.3.4              RSQLite_2.1.1           
## [43] DelayedMatrixStats_1.6.0 mclust_5.4.3            
## [45] RCurl_1.95-4.12          magrittr_1.5            
## [47] GenomeInfoDbData_1.2.1   Matrix_1.2-17           
## [49] fansi_0.4.0              Rcpp_1.0.1              
## [51] Rhdf5lib_1.6.0           stringi_1.4.3           
## [53] yaml_2.2.0               MASS_7.3-51.4           
## [55] zlibbioc_1.30.0          rhdf5_2.28.0            
## [57] grid_3.6.0               blob_1.1.1              
## [59] crayon_1.3.4             lattice_0.20-38         
## [61] splines_3.6.0            annotate_1.62.0         
## [63] multtest_2.40.0          GenomicFeatures_1.36.0  
## [65] hms_0.4.2                knitr_1.22              
## [67] beanplot_1.2             pillar_1.3.1            
## [69] rngtools_1.3.1.1         codetools_0.2-16        
## [71] biomaRt_2.40.0           XML_3.98-1.19           
## [73] glue_1.3.1               evaluate_0.13           
## [75] data.table_1.12.2        cellranger_1.1.0        
## [77] openssl_1.3              purrr_0.3.2             
## [79] reshape_0.8.8            assertthat_0.2.1        
## [81] xfun_0.6                 xtable_1.8-4            
## [83] survival_2.44-1.1        tibble_2.1.1            
## [85] GenomicAlignments_1.20.0 AnnotationDbi_1.46.0    
## [87] registry_0.5-1           memoise_1.1.0
```
