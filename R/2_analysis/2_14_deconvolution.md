---
title: "2_14_deconvolution"
author: "Victor Yuan"
date: "16/01/2020"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 4
    toc_float:
      collapsed: false
    theme: spacelab
    highlight: tango 
editor_options: 
  chunk_output_type: console
---
# Setup

## Libraries


```r
# libraries and data
library(tidyverse)
library(scales)
library(here)
library(readxl)
library(janitor)
library(minfi)
library(broom)
library(pheatmap)
theme_set(theme_bw())
library(EpiDISH)
library(yardstick)
library(ggbeeswarm)
library(planet)
```

## Data


```r
base_path <- file.path('data', 'main', 'interim')

# pData
pDat <- readRDS(here(base_path, '2_3_pDat_contam.rds'))
pDat <- pDat %>%
  mutate(Tissue = case_when(
    !(Tissue %in% c('Villi', 'Villi maternal', 'Syncytiotrophoblast', 'Mixture')) ~ paste(Tissue, 'cs'),
    Tissue == 'Syncytiotrophoblast' ~ 'Trophoblasts enz',
    TRUE ~ Tissue
  )) 

# update ga
new_ga <- readRDS(here::here('data', 'main', 'interim', '0_3_pDat.rds'))

pDat <- pDat %>%
  # update ga
  select(-GA) %>%
  left_join(new_ga %>% select(Sample_Name, GA))
```

```
## Joining, by = "Sample_Name"
```

```r
# raw methylation data
betas <- readRDS(here(base_path, '1_4_betas_noob_filt.rds'))
mset <- readRDS(here(base_path, '1_4_mset_noob.rds' ))
colnames(mset) <- pData(mset)$Sample_Name

# annotation
anno <- readRDS('Z:/Victor/Repositories/EPIC_annotation/hg19_epic_annotation.rds')
anno <- anno %>%
  as_tibble() %>%
  filter(cpg %in% rownames(betas)) # filter to filtered betas cpgs

# cell proportions for mixtures
mixtures <- read_excel(here::here('data', 'main', 'raw', 'Mixtures.xlsx'))

# color key
color_code <- readRDS(here(base_path, '2_3_color_code.rds'))
color_code_tissue <- setNames(color_code$Colors_Tissue, gsub(' cs', '',color_code$label))

color_code_tissue <- c(color_code_tissue, 'nRBC' = 'grey')
color_code_tissue <- c(color_code_tissue, 'Syncytiotrophoblast' = '#f4702e')

# ancestry information
ancestry <- readRDS(here(base_path, '2_11_ancestry.rds'))
```

## Remove samples


```r
pDat_filt <- pDat %>% 
  filter(maternal_contamination_norm_flip < 0.35,
         !Sample_Name %in% c('PM364_hofb_cs', 'PL293_v_R2', 'PM366_vc_R2', 'P131_hofb_cs', 
                             'PM324_V4', 'PM324_V1', 'PM139_vc', 'PM77_vc'),
         !Tissue %in% c('Villi maternal', 
                        'Trophoblasts enz',
                        'Dead Cells and Lymphocytes cs'),
         Trimester != 'Second')

betas_filt <- betas[,pDat_filt$Sentrix]
colnames(betas_filt) <- pDat_filt$Sample_Name
```

# get nRBC data

Meg's recent publciation provides a curated reference for cord blood data: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6712867/#CR22
"Systematic evaluation and validation of reference and library selection methods for deconvolution of cord blood DNA methylation data"

https://bioconductor.org/packages/release/data/experiment/vignettes/FlowSorted.CordBloodCombined.450k/inst/doc/FlowSorted.CordBloodCombined.450k.html

Note this is a 450k dataset

```r
if (memory.limit()>8000){
library(ExperimentHub)  

hub <- ExperimentHub()  

myfiles <- query(hub, "FlowSorted.CordBloodCombined.450k")

FlowSorted.CordBloodCombined.450k <- myfiles[[1]]

FlowSorted.CordBloodCombined.450k
}

# subset to nrbc
nrbc <- FlowSorted.CordBloodCombined.450k[,pData(FlowSorted.CordBloodCombined.450k)$CellType == 'nRBC']

pData(nrbc)
nrbc_noob <- preprocessNoob(nrbc)


nrbc_pDat <- pData(nrbc_noob) %>%
  as_tibble() %>%
  dplyr::rename(Tissue = CellType) %>%
  mutate(Trimester = 'Third') %>%
  select(Sample_Name, Tissue, Trimester, Sex, Study) 
  
pData(nrbc_noob) <- DataFrame(nrbc_pDat)

saveRDS(nrbc_noob, here::here(base_path, '2_14_nrbc_noob.rds'))
```


# Houseman Deconvolution

I apply to first and third trimester samples separately, each with their own trimester-specific reference.


```r
# filter to project samples
mset_filt <- mset[,pDat_filt$Sample_Name]
pData(mset_filt) <- DataFrame(pDat_filt)

# necessary column for pickcompprobes
mset$CellType <- mset$Tissue

mset_ref <- mset[,pData(mset)$Tissue %in% c('Trophoblasts', 'Stromal', 'Hofbauer', 'Endothelial')]
mset_test <- mset[,pData(mset)$Tissue %in% c('Mixture', 'Villi')]

# pick probes
probes_third <- minfi:::pickCompProbes(
  mset_ref[, pData(mset_ref)$Trimester == 'Third'], 
  cellTypes = c('Trophoblasts', 'Stromal', 'Hofbauer', 'Endothelial'),
  compositeCellType = 'Placenta',
  probeSelect = "both")

probes_first <- minfi:::pickCompProbes(
  mset_ref[, pData(mset_ref)$Trimester == 'First'], 
  cellTypes = c('Trophoblasts', 'Stromal', 'Hofbauer', 'Endothelial'),
  compositeCellType = 'Placenta',
  probeSelect = "both")


# extract coefficients
coefs_third <- probes_third$coefEsts
coefs_first <- probes_first$coefEsts

# project 
counts_third <- minfi:::projectCellType(
  getBeta(mset_test)[rownames(coefs_third),pData(mset_test)$Trimester == 'Third'], 
  coefs_third,
  lessThanOne = FALSE)

counts_first<- minfi:::projectCellType(
  getBeta(mset_test)[rownames(coefs_first),pData(mset_test)$Trimester == 'First'], 
  coefs_first,
  lessThanOne = FALSE)

res <- rbind(counts_third, counts_first) %>%
  as.data.frame() %>%
  bind_cols(Sample_Name = rownames(.), .) %>%
  inner_join(pData(mset_test) %>% as_tibble());res
```

```
## Joining, by = "Sample_Name"
```

```
##    Sample_Name Trophoblasts      Stromal    Hofbauer Endothelial Chip_number
## 1     PM376_vc   0.91898550 2.799734e-02 0.027716038  0.04822980          25
## 2     PM368_vc   0.87179085 6.058860e-02 0.032828154  0.06766962          26
## 3    Mixture 7   0.03657594 6.883150e-01 0.040308815  0.16265603          26
## 4     PM378_vc   0.86628791 6.516424e-02 0.021089717  0.06398130          27
## 5     PM364_vc   0.90875323 3.030651e-02 0.023808334  0.06695470          28
## 6     PM365_vc   0.83679048 7.429442e-02 0.045364086  0.08199404          29
## 7     PM362_vc   0.96394349 3.630820e-02 0.017409475  0.03608503          29
## 8     PM371_vc   0.76304964 1.222862e-01 0.076421048  0.07362329          30
## 9   Mixture 11   0.42639283 3.873348e-01 0.046348015  0.04685752          30
## 10    PM324_V1   0.94895001 5.594233e-03 0.002585087  0.06893579          31
## 11   Mixture 1   0.08179434 6.620907e-01 0.005047137  0.17141782          31
## 12    PM379_vc   0.89081829 6.706790e-02 0.030236147  0.03957125          32
## 13    PM139_vc   0.93214824 8.673617e-19 0.002052590  0.08006455          33
## 14   Mixture 9   0.10905301 6.284505e-01 0.063858459  0.10901726          33
## 15    PM324_V4   0.91299993 2.891876e-02 0.007936173  0.08353584          34
## 16  Mixture 14   0.69674020 7.505516e-03 0.024084730  0.24396500          35
## 17    PM370_vc   0.78104453 1.332050e-01 0.035925203  0.07447703          36
## 18    PM369_vc   0.92655822 5.074537e-02 0.023914850  0.03704933          37
## 19    PM367_vc   0.88947635 4.889253e-02 0.025075074  0.05148355          38
## 20    PM373_vc   0.77969352 1.474665e-01 0.039856974  0.05556435          39
## 21  Mixture 10   0.36197981 4.848873e-01 0.016441254  0.00000000          39
## 22    PM372_vc   0.83025325 1.047411e-01 0.036205204  0.06263540          40
## 23    PM359_vc   0.91183665 5.314697e-02 0.024077557  0.05218598          40
## 24    PM381_vc   0.90193769 7.996318e-02 0.029008880  0.01691608          41
## 25   Mixture 6   0.38356088 2.585059e-01 0.127741683  0.06346041          41
## 26   Mixture 8   0.14932535 5.773169e-01 0.000000000  0.16707690          42
## 27    PM377_vc   0.89982428 6.487529e-02 0.032505446  0.02942294          43
## 28  Mixture 13   0.47404358 5.118478e-01 0.008098884  0.00000000          44
## 29     PM77_vc   1.00519861 0.000000e+00 0.000000000  0.02996973          44
## 30   Mixture 3   0.25713496 1.682099e-01 0.432121177  0.15845754          45
## 31    PM375_vc   0.80717246 1.115836e-01 0.059585471  0.06109076          45
## 32   Mixture 2   0.04570272 5.158763e-01 0.209182255  0.23108708          46
## 33 PM366_vc_R1   0.78371449 1.403134e-01 0.032933491  0.05299281          46
## 34  Mixture 12   0.24247263 4.526419e-01 0.144052136  0.03236086          47
## 35    PM374_vc   0.94052447 2.320682e-02 0.015868804  0.04210504          47
## 36 PM366_vc_R2   0.79751474 1.344906e-01 0.029481323  0.04429916          48
## 37     PL290_v   0.25440124 5.675812e-01 0.082630229  0.10474558          27
## 38     PL296_v   0.46382351 3.696750e-01 0.068020385  0.14101894          28
## 39     PL292_v   0.76925462 1.970351e-01 0.058209999  0.00000000          31
## 40     PL289_v   0.67302085 2.459188e-01 0.065061746  0.08872473          32
## 41  PL293_v_R2   0.52965706 3.449103e-01 0.083817134  0.10422396          33
## 42     PL294_v   0.39685916 4.318164e-01 0.043361373  0.15144440          37
## 43  PL293_v_R1   0.49678158 3.532349e-01 0.076776289  0.11050460          40
## 44     PL295_v   0.48666655 3.695904e-01 0.060180709  0.13334343          48
##    Row Well      Case_ID Sex   GA Trimester DNA_QP Week   Sample_Plate  Tissue
## 1    6   F1        PM376   F   NA     Third   <NA>    1 WG6980707-MSA4   Villi
## 2    1   A2        PM368   M   NA     Third   <NA>    1 WG6980707-MSA4   Villi
## 3    6   F2        PM372   M   NA     Third    400    1 WG6980707-MSA4 Mixture
## 4    1   A3        PM378   M   NA     Third   <NA>    1 WG6980707-MSA4   Villi
## 5    3   C4        PM364   M   NA     Third   <NA>    1 WG6980707-MSA4   Villi
## 6    1   A5        PM365   F   NA     Third   <NA>    1 WG6980707-MSA4   Villi
## 7    8   H5        PM362   M   NA     Third     NA    1 WG6980707-MSA4   Villi
## 8    2   B6        PM371   F   NA     Third     NA    1 WG6980707-MSA4   Villi
## 9    5   E6        PM375   F   NA     Third    400    1 WG6980707-MSA4 Mixture
## 10   5   E7        PM324   M   NA     Third   <NA>    1 WG6980707-MSA4   Villi
## 11   7   G7        PM366   F   NA     Third    400    1 WG6980707-MSA4 Mixture
## 12   2   B8        PM379   F   NA     Third     NA    1 WG6980707-MSA4   Villi
## 13   4   D9        PM139   M   NA     Third   <NA>    1 WG6980707-MSA4   Villi
## 14   7   G9        PM372   M   NA     Third    400    1 WG6980707-MSA4 Mixture
## 15   7  G10        PM324   M   NA     Third   <NA>    1 WG6980707-MSA4   Villi
## 16   1  A11        PM374   M   NA     Third    200    1 WG6980707-MSA4 Mixture
## 17   2  B12        PM370   F   NA     Third     NA    1 WG6980707-MSA4   Villi
## 18   8   H1        PM369   F   NA     Third     NA    2 WG6980638-MSA4   Villi
## 19   4   D2        PM367   F   NA     Third     NA    2 WG6980638-MSA4   Villi
## 20   3   C3        PM373   M   NA     Third     NA    2 WG6980638-MSA4   Villi
## 21   4   D3        PM376   F   NA     Third    400    2 WG6980638-MSA4 Mixture
## 22   6   F4        PM372   M   NA     Third     NA    2 WG6980638-MSA4   Villi
## 23   7   G4        PM359   M   NA     Third     NA    2 WG6980638-MSA4   Villi
## 24   3   C5        PM381   M   NA     Third     NA    2 WG6980638-MSA4   Villi
## 25   6   F5  Tech Ctrl 1   M   NA     Third    400    2 WG6980638-MSA4 Mixture
## 26   8   H6        PM375   F   NA     Third    400    2 WG6980638-MSA4 Mixture
## 27   7   G7        PM377   F   NA     Third     NA    2 WG6980638-MSA4   Villi
## 28   2   B8        PM373   M   NA     Third    200    2 WG6980638-MSA4 Mixture
## 29   8   H8         PM77   M   NA     Third   <NA>    2 WG6980638-MSA4   Villi
## 30   1   A9 PM358 (tech)   M   NA     Third    400    2 WG6980638-MSA4 Mixture
## 31   4   D9        PM375   F   NA     Third     NA    2 WG6980638-MSA4   Villi
## 32   4  D10    Tech Ctrl   M   NA     Third    400    2 WG6980638-MSA4 Mixture
## 33   5  E10        PM366   F   NA     Third     NA    2 WG6980638-MSA4   Villi
## 34   3  C11        PM372   M   NA     Third    400    2 WG6980638-MSA4 Mixture
## 35   6  F11        PM374   M   NA     Third     NA    2 WG6980638-MSA4   Villi
## 36   5  E12        PM366   F   NA     Third   <NA>    2 WG6980638-MSA4   Villi
## 37   7   G3        PL290   F 10.5     First   <NA>    1 WG6980707-MSA4   Villi
## 38   1   A4        PL296   M 11.5     First   <NA>    1 WG6980707-MSA4   Villi
## 39   3   C7        PL292   M  7.1     First     NA    1 WG6980707-MSA4   Villi
## 40   6   F8        PL289   M  6.4     First     NA    1 WG6980707-MSA4   Villi
## 41   8   H9        PL293   M   NA     First   <NA>    1 WG6980707-MSA4   Villi
## 42   4   D1        PL294   F 11.6     First     NA    2 WG6980638-MSA4   Villi
## 43   2   B4        PL293   M  9.6     First     NA    2 WG6980638-MSA4   Villi
## 44   6  F12        PL295   F 11.6     First     NA    2 WG6980638-MSA4   Villi
##      Sentrix_ID Sentrix_Position Scratches Batch_BSC DNA_conc_BSC_adjusted
## 1  203067920143           R06C01      <NA>      BSC1               63.2200
## 2  203067920144           R01C01      <NA>  BSC_TEST               57.7800
## 3  203067920144           R06C01      <NA>      BSC0               12.9800
## 4  203068760084           R01C01      <NA>      BSC1               49.5700
## 5  203068760095           R03C01      <NA>      BSC1               64.7100
## 6  203068760191           R01C01      <NA>      BSC1               61.1700
## 7  203068760191           R08C01      <NA>      BSC2               51.3700
## 8  203068760192           R02C01      <NA>      BSC2               43.1200
## 9  203068760192           R05C01      <NA>      BSC0               17.3900
## 10 203068760202           R05C01      <NA>      BSC1               69.7900
## 11 203068760202           R07C01      <NA>      BSC0               15.0600
## 12 203068760204           R02C01      <NA>      BSC1               54.7000
## 13 203072530027           R04C01      <NA>      BSC1               54.9700
## 14 203072530027           R07C01      <NA>      BSC0               14.0900
## 15 203072530046           R07C01      <NA>      BSC2               78.1600
## 16 203072530098           R01C01      <NA>      BSC2                4.5600
## 17 203072530157           R02C01      <NA>      BSC2               62.0700
## 18 203072530204           R08C01      <NA>      BSC2               56.7750
## 19 203072530227           R04C01      <NA>      BSC2               56.2500
## 20 203077630004           R03C01      <NA>      BSC2               58.2525
## 21 203077630004           R04C01      <NA>      BSC0                9.2400
## 22 203077630012           R06C01      <NA>      BSC2               62.3175
## 23 203077630012           R07C01      <NA>      BSC2               51.2625
## 24 203077630052           R03C01      <NA>      BSC2               69.5475
## 25 203077630052           R06C01      <NA>      BSC0               10.8700
## 26 203077630080           R08C01      <NA>      BSC0               15.2100
## 27 203077630103           R07C01      <NA>      BSC3               45.9375
## 28 203077630104           R02C01      <NA>      BSC3                0.0000
## 29 203077630104           R08C01      <NA>      BSC3               70.2900
## 30 203077630172           R01C01      <NA>      BSC3                0.0000
## 31 203077630172           R04C01      <NA>      BSC3               55.1700
## 32 203077630181           R04C01      <NA>      BSC3                0.0000
## 33 203077630181           R05C01      <NA>      BSC3               47.0000
## 34 203077630199           R03C01      <NA>      BSC0                9.3900
## 35 203077630199           R06C01      <NA>      BSC3               61.2375
## 36 203077630206           R05C01      <NA>      BSC3               47.0000
## 37 203068760084           R07C01       Yes      BSC1               46.6900
## 38 203068760095           R01C01      <NA>  BSC_TEST               60.5500
## 39 203068760202           R03C01      <NA>      BSC2               58.7800
## 40 203068760204           R06C01      <NA>  BSC_TEST               46.0300
## 41 203072530027           R08C01      <NA>      BSC2               34.6350
## 42 203072530204           R04C01      <NA>  BSC_TEST               59.0800
## 43 203077630012           R02C01      <NA>      BSC2               42.0525
## 44 203077630206           R06C01      <NA>  BSC_TEST               54.5400
##    DNA_conc_before_load DNA_loaded             Sentrix CellType
## 1               50.0000     200.00 203067920143_R06C01    Villi
## 2               50.0000     200.00 203067920144_R01C01    Villi
## 3               32.4500     129.80 203067920144_R06C01  Mixture
## 4               49.5700     198.28 203068760084_R01C01    Villi
## 5               50.0000     200.00 203068760095_R03C01    Villi
## 6               50.0000     200.00 203068760191_R01C01    Villi
## 7               50.0000     200.00 203068760191_R08C01    Villi
## 8               43.1200     172.48 203068760192_R02C01    Villi
## 9               43.4750     173.90 203068760192_R05C01  Mixture
## 10              50.0000     200.00 203068760202_R05C01    Villi
## 11              37.6500     150.60 203068760202_R07C01  Mixture
## 12              50.0000     200.00 203068760204_R02C01    Villi
## 13              50.0000     200.00 203072530027_R04C01    Villi
## 14              35.2250     140.90 203072530027_R07C01  Mixture
## 15              50.0000     200.00 203072530046_R07C01    Villi
## 16              11.4000      45.60 203072530098_R01C01  Mixture
## 17              50.0000     200.00 203072530157_R02C01    Villi
## 18              50.0000     200.00 203072530204_R08C01    Villi
## 19              50.0000     200.00 203072530227_R04C01    Villi
## 20              50.0000     200.00 203077630004_R03C01    Villi
## 21              23.1000      92.40 203077630004_R04C01  Mixture
## 22              50.0000     200.00 203077630012_R06C01    Villi
## 23              50.0000     200.00 203077630012_R07C01    Villi
## 24              50.0000     200.00 203077630052_R03C01    Villi
## 25              27.1750     108.70 203077630052_R06C01  Mixture
## 26              38.0250     152.10 203077630080_R08C01  Mixture
## 27              45.9375     183.75 203077630103_R07C01    Villi
## 28               0.0000       0.00 203077630104_R02C01  Mixture
## 29              50.0000     200.00 203077630104_R08C01    Villi
## 30               0.0000       0.00 203077630172_R01C01  Mixture
## 31              50.0000     200.00 203077630172_R04C01    Villi
## 32               0.0000       0.00 203077630181_R04C01  Mixture
## 33              47.0000     188.00 203077630181_R05C01    Villi
## 34              23.4750      93.90 203077630199_R03C01  Mixture
## 35              50.0000     200.00 203077630199_R06C01    Villi
## 36              47.0000     188.00 203077630206_R05C01    Villi
## 37              46.6900     186.76 203068760084_R07C01    Villi
## 38              50.0000     200.00 203068760095_R01C01    Villi
## 39              50.0000     200.00 203068760202_R03C01    Villi
## 40              46.0300     184.12 203068760204_R06C01    Villi
## 41              57.7250     230.90 203072530027_R08C01    Villi
## 42              50.0000     200.00 203072530204_R04C01    Villi
## 43              42.0525     168.21 203077630012_R02C01    Villi
## 44              50.0000     200.00 203077630206_R06C01    Villi
```

```r
# plot
res %>%
  select(Tissue, Sample_Name:Endothelial)  %>%
  pivot_longer(cols = -c(Tissue, Sample_Name),
               names_to = 'Estimates',
               values_to = 'p') %>%
  ggplot(aes(x = Sample_Name, y = p, fill = Estimates)) +
  geom_bar(stat = 'identity', position = 'stack') +
  facet_grid(cols = vars(Tissue), scales = 'free', space = 'free') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0,0))
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


# Compare accuracy

Mixture 13 replaces mixture 4 due to failed BSC conversion

Mixture 4 replaces mixture 8 due to failed BSC conversion
Mixture 14 replaces mixture 7 ^ ^ ^ ^
Mixture 7 replaces mixture 5 due to loss of DNA during mixing


```r
mixtures_results <- mixtures %>% select(Sample, Case_ID, contains('%')) %>%
  janitor::clean_names() %>%
  pivot_longer(cols = -c(sample, case_id),
               names_to = 'component',
               values_to = 'actual_percent') %>%
  mutate(component = gsub('_percent', '', component))  %>%
  dplyr::rename(Sample_Name = sample) %>%
  filter(component !='total') %>%
  mutate(Sample_Name = case_when(
    Sample_Name == 'Mixture 4' ~ 'Mixture 13',
    
    Sample_Name == 'Mixture 5' ~ 'Mixture 14',
    TRUE ~ as.character(Sample_Name)
  ))


mixture_results <- res %>%
  select(Sample_Name:Endothelial) %>%
  pivot_longer(cols = -Sample_Name,
               names_to = 'component',
               values_to = 'predicted_percent') %>%
  mutate(component = str_to_lower(component),
         predicted_percent = predicted_percent*100) %>%
  full_join(mixtures_results, by = c('Sample_Name', 'component'))  %>%
  mutate(component = str_to_sentence(component))

# stats
stats <- mixture_results %>%
  filter(grepl('Mix', Sample_Name)) %>%
  nest(data = -component) %>%
  mutate(lm = map(data, ~lm(predicted_percent ~ actual_percent, .)),
         glanced = map(lm, glance)) %>%
  select(-data, -lm) %>%
  unnest(glanced) %>%
  mutate(label = paste0('R2=', signif(r.squared, digits = 2),
                        '\n', pvalue(p.value, add_p = TRUE)))


# scatter plots
mixture_results %>%
  filter(grepl('Mix', Sample_Name)) %>%
  ggplot() +
  geom_point(aes(x = actual_percent, y = predicted_percent, color = component)) +
  geom_smooth(aes(x = actual_percent, y = predicted_percent, color = component),
              method = 'lm', se = FALSE) +
  geom_label(data = stats, aes(label = label),
             x = 0, y = 80, hjust = 0) +
  scale_y_continuous(limits = c(0,100)) +
  facet_wrap(~component) +
  scale_color_manual(values = color_code_tissue[unique(mixture_results$component)], 
                         na.value = '#636363',
                         guide = 'none') 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
# bar plots
mixture_results %>%
  filter(grepl('Mix', Sample_Name)) %>%
  mutate(Sample_Name = factor(Sample_Name, 
                              levels = paste0('Mixture ', c(1:3, 13, 14, 6:12)))) %>%
  pivot_longer(cols = contains('percent'),
               names_to = 'type',
               values_to = 'percent') %>%
  mutate(type = gsub('_.*', '', type)) %>%
  ggplot() +
  geom_bar(stat=  'identity',
           aes(x = Sample_Name, y = percent, fill = component)) +
  facet_wrap(ncol = 1, vars(type)) +
  labs(x = '') +
  scale_x_discrete(labels = function(x)gsub('Mixture ', '', x))
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-6-2.png)<!-- -->

```r
mixture_results %>%
  filter(grepl('Mix', Sample_Name)) %>%
  mutate(Sample_Name = factor(Sample_Name, 
                              levels = paste0('Mixture ', c(1:3, 13, 14, 6:12)))) %>%
  pivot_longer(cols = contains('percent'),
               names_to = 'type',
               values_to = 'percent') %>%
  mutate(type = gsub('_.*', '', type)) %>%
  ggplot() +
  geom_bar(stat=  'identity',
           aes(x = Sample_Name, y = percent, fill = component)) +
  facet_grid(rows = vars(type),
             cols = vars(case_id), scales = 'free', space = 'free') +
  labs(x = '') +
  theme(strip.background = element_blank()) +
  scale_x_discrete(labels = function(x)gsub('Mixture ', '', x))
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-6-3.png)<!-- -->

# EPIDISH


```r
epidish_res1 <- epidish(
  beta.m = getBeta(mset_test)[rownames(coefs_third),pData(mset_test)$Trimester == 'Third'],
  ref.m = coefs_third,
  method = 'RPC')

epidish_res2 <- epidish(
  beta.m = getBeta(mset_test)[rownames(coefs_third),pData(mset_test)$Trimester == 'Third'],
  ref.m = coefs_third,
  method = 'CBS')
```

```
## 1
```

```
## 2
```

```
## 3
```

```r
epidish_res3 <- epidish(
  beta.m = getBeta(mset_test)[rownames(coefs_third),pData(mset_test)$Trimester == 'Third'],
  ref.m = coefs_third,
  method = 'CP',
  constraint = 'inequality')
```

```
## 1
```

```
## 2
```

```
## 3
```

```
## 4
```

```
## 5
```

```
## 6
```

```
## 7
```

```
## 8
```

```
## 9
```

```
## 10
```

```
## 11
```

```
## 12
```

```
## 13
```

```
## 14
```

```
## 15
```

```
## 16
```

```
## 17
```

```
## 18
```

```
## 19
```

```
## 20
```

```
## 21
```

```
## 22
```

```
## 23
```

```
## 24
```

```
## 25
```

```
## 26
```

```
## 27
```

```
## 28
```

```
## 29
```

```
## 30
```

```
## 31
```

```
## 32
```

```
## 33
```

```
## 34
```

```
## 35
```

```
## 36
```

```r
#first trime
epidish_res11 <- epidish(
  beta.m = getBeta(mset_test)[rownames(coefs_first),pData(mset_test)$Trimester == 'First'],
  ref.m = coefs_first,
  method = 'RPC')

epidish_res22 <- epidish(
  beta.m = getBeta(mset_test)[rownames(coefs_first),pData(mset_test)$Trimester == 'First'],
  ref.m = coefs_first,
  method = 'CBS')
```

```
## 1
```

```
## 2
```

```
## 3
```

```r
epidish_res33 <- epidish(
  beta.m = getBeta(mset_test)[rownames(coefs_first),pData(mset_test)$Trimester == 'First'],
  ref.m = coefs_first,
  method = 'CP',
  constraint = 'inequality')
```

```
## 1
```

```
## 2
```

```
## 3
```

```
## 4
```

```
## 5
```

```
## 6
```

```
## 7
```

```
## 8
```

```r
str(epidish_res1)
```

```
## List of 3
##  $ estF   : num [1:36, 1:4] 0.9278 0.902 0.0681 0.899 0.9421 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:36] "PM376_vc" "PM368_vc" "Mixture 7" "PM378_vc" ...
##   .. ..$ : chr [1:4] "Trophoblasts" "Stromal" "Hofbauer" "Endothelial"
##  $ ref    : num [1:400, 1:4] 0.1014 0.0719 0.1507 0.1321 0.1282 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:400] "cg10590657" "cg00281273" "cg14514160" "cg16156113" ...
##   .. ..$ : chr [1:4] "Trophoblasts" "Stromal" "Hofbauer" "Endothelial"
##  $ dataREF: num [1:400, 1:36] 0.1118 0.0797 0.321 0.2209 0.219 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:400] "cg10590657" "cg00281273" "cg14514160" "cg16156113" ...
##   .. ..$ : chr [1:36] "PM376_vc" "PM368_vc" "Mixture 7" "PM378_vc" ...
```

```r
epidish_res1 <- epidish_res1$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (RPC)')

epidish_res2 <- epidish_res2$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (CBS)')

epidish_res3 <- epidish_res3$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (CP)')

epidish_res11 <- epidish_res11$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (RPC)')

epidish_res22 <- epidish_res22$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (CBS)')

epidish_res33 <- epidish_res33$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (CP)')

epidish_results <- bind_rows(epidish_res1, epidish_res2, epidish_res3, 
                             epidish_res11, epidish_res22, epidish_res33) %>%
   dplyr::rename(Sample_Name = Sample_Names) %>%
  mutate(percent = percent*100) %>%
  left_join(mixture_results %>% select(Sample_Name, case_id))
```

```
## Joining, by = "Sample_Name"
```

```r
# combine epidish results with houseman
mixture_results_all <- mixture_results %>%
  pivot_longer(
    cols = contains('percent'),
    names_to = 'algorithm',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = gsub('_percent', '', algorithm),
         algorithm = ifelse(algorithm == 'predicted', 'minfi (CP)', algorithm)) %>%
  select(-case_id) %>%
  bind_rows(epidish_results %>% select(-case_id)) %>%
  mutate(algorithm = factor(algorithm, 
                            levels = c('actual', 'epidish (RPC)', 'epidish (CBS)', 
                                       'epidish (CP)', 'minfi (CP)'))) %>%
  distinct()

# bar plot
mixture_results_all %>%
  filter(grepl('Mix', Sample_Name)) %>%
  mutate(Sample_Name = factor(Sample_Name, 
                              levels = paste0('Mixture ', c(1:3, 13, 14, 6:12)))) %>%
  ggplot() +
  geom_bar(stat=  'identity',
           aes(x = Sample_Name, y = percent, fill = component)) +
  facet_wrap(ncol = 1, vars(algorithm)) +
  labs(x = '', y= '') +
  theme(strip.background = element_blank())+
  scale_fill_manual(values = color_code_tissue[unique(mixture_results_all$component)],
                    na.value = '#636363') +
  scale_y_continuous(labels = function(x)paste0(x, '%'), expand = c(0,0))+
  scale_x_discrete(labels = function(x)gsub('Mixture ', '', x), expand = c(0,0))
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
# stats
stats <- mixture_results_all %>%
  filter(grepl('Mix', Sample_Name)) %>%
  pivot_wider(id_cols = c(Sample_Name, component),
              names_from = 'algorithm',
              values_from = 'percent') %>%
  pivot_longer(cols = -c(Sample_Name, component, actual),
               names_to = 'algorithm',
               values_to = 'predicted') %>%
  
  nest(data = -c(component, algorithm)) %>%
  mutate(lm = map(data, ~lm(predicted ~ actual, .)),
         glanced = map(lm, glance)) %>%
  select(-data, -lm) %>%
  unnest(glanced) %>%
  mutate(label = paste0('R2=', signif(r.squared, digits = 2),
                        '\n', pvalue(p.value, add_p = TRUE)))

mixture_results_all %>%
  filter(grepl('Mix', Sample_Name)) %>%
  pivot_wider(id_cols = c(Sample_Name, component),
              names_from = 'algorithm',
              values_from = 'percent') %>%
  pivot_longer(cols = -c(Sample_Name, component, actual),
               names_to = 'algorithm',
               values_to = 'predicted') %>%
  group_by(algorithm) %>%
  rmse(truth = actual, estimate = predicted)
```

```
## # A tibble: 4 x 4
##   algorithm     .metric .estimator .estimate
##   <chr>         <chr>   <chr>          <dbl>
## 1 epidish (CBS) rmse    standard        21.1
## 2 epidish (CP)  rmse    standard        20.9
## 3 epidish (RPC) rmse    standard        20.9
## 4 minfi (CP)    rmse    standard        20.9
```

```r
stats %>%
  group_by(algorithm) %>%
  summarize(mean_r2 = mean(r.squared))
```

```
## `summarise()` ungrouping output (override with `.groups` argument)
```

```
## # A tibble: 4 x 2
##   algorithm     mean_r2
##   <chr>           <dbl>
## 1 epidish (CBS) 0.00757
## 2 epidish (CP)  0.00572
## 3 epidish (RPC) 0.00846
## 4 minfi (CP)    0.00561
```

```r
mixture_results_all %>%
  filter(grepl('Mix', Sample_Name)) %>%
  pivot_wider(id_cols = c(Sample_Name, component),
              names_from = 'algorithm',
              values_from = 'percent') %>%
  pivot_longer(cols = -c(Sample_Name, component, actual),
               names_to = 'algorithm',
               values_to = 'predicted')  %>%
  ggplot(aes(x = actual, y = predicted, color = component)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  geom_label(color = 'black', size = 3, data = stats, aes(label = label),
             x = 0, y = 100, hjust = 0, vjust = 1) +
  facet_grid(cols = vars(algorithm), rows = vars(component)) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = color_code_tissue[unique(mixture_results_all$component)], 
                         na.value = '#636363',
                         guide = 'none') +
  scale_x_continuous(labels = function(x)paste0(x, '%'), breaks = c(0,30,60,100))+
  scale_y_continuous(limits = c(0,100),labels = function(x)paste0(x, '%')) +
  labs(x = 'Actual', y = 'Predicted')
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-7-2.png)<!-- -->

# Heatmap of deconvolution probes


```r
probes_third %>% str
```

```
## List of 3
##  $ coefEsts   : num [1:400, 1:4] 0.1014 0.0719 0.1507 0.1321 0.1282 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:400] "cg10590657" "cg00281273" "cg14514160" "cg16156113" ...
##   .. ..$ : chr [1:4] "Trophoblasts" "Stromal" "Hofbauer" "Endothelial"
##  $ compTable  :'data.frame':	866091 obs. of  9 variables:
##   ..$ Fstat       : num [1:866091] 0.061 0.0559 87.0695 0.7957 2.1478 ...
##   ..$ p.value     : num [1:866091] 9.80e-01 9.83e-01 6.76e-24 5.00e-01 1.02e-01 ...
##   ..$ Trophoblasts: num [1:866091] 0.0119 0.0245 0.3432 0.9525 0.8267 ...
##   ..$ Stromal     : num [1:866091] 0.012 0.025 0.0847 0.9478 0.8711 ...
##   ..$ Hofbauer    : num [1:866091] 0.0119 0.0243 0.4373 0.9529 0.8743 ...
##   ..$ Endothelial : num [1:866091] 0.0117 0.0246 0.425 0.9554 0.8537 ...
##   ..$ low         : num [1:866091] 0.00715 0.01229 0.02142 0.90207 0.61026 ...
##   ..$ high        : num [1:866091] 0.0161 0.0362 0.6382 0.9809 0.9374 ...
##   ..$ range       : num [1:866091] 0.00894 0.02393 0.61676 0.07888 0.32718 ...
##  $ sampleMeans: Named num [1:76] 0.492 0.506 0.507 0.502 0.494 ...
##   ..- attr(*, "names")= chr [1:76] "Endothelial" "Hofbauer" "Trophoblasts" "Endothelial" ...
```

```r
pheatmap(probes_third$coefEsts, cluster_cols = FALSE, cluster_rows = TRUE,
         show_rownames = FALSE) 
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

# FACS-determined proportions

So mixture data is potentially garbage. Here I look at the FACS cell proportions.


```r
facs_prop <- read_excel(here::here('data', 'main', 'raw', 'FACS_proportions.xlsx'))

# process
facs_prop <- facs_prop %>%
  as.matrix() %>%
  magrittr::set_rownames(.[,'Cell_type']) %>%
  t() %>%
  as.data.frame %>%
  bind_cols(Sample_Name = rownames(.), .) %>%
  dplyr::slice(-1) %>%
  select(Sample_Name, contains('% of live')) %>%
  pivot_longer(cols = -Sample_Name, names_to = 'component', values_to = 'value') %>%
  mutate(component = gsub(' \\(\\% of live\\)', '', component),
         percent = str_extract(value, '[0-9]?\\.?[0-9]+(?=\\%)') %>% as.numeric(),
         Case_ID = str_extract(Sample_Name, 'PM[0-9]+'),
         component = ifelse(component == 'Fibroblasts', 'Stromal', component)) %>%
  select(-value) %>% 
  
  # REMOVE samples
  filter(!grepl('TB Layer', Sample_Name),
         ) %>% 
  
  # average percentage
  group_by(Case_ID, component) %>%
  summarize(percent = mean(percent)) %>%
  
  # rescale such that sums to 1
  mutate(percent = percent*100 / sum(percent)) %>%
  group_by(Case_ID, component == 'Trophoblasts') %>%
  mutate(percent_constrained = percent*100/ sum(percent)) %>%
  ungroup() %>%
  mutate(percent_constrained = ifelse(percent_constrained == 100, NA,
                                      percent_constrained)) %>%
  select(-`component == "Trophoblasts"`) %>%
  
  pivot_longer(cols = contains('percent'),
               names_to = 'FACS',
               values_to = 'actual') %>%
  mutate(FACS = ifelse(FACS == 'percent', 'FACS',
                       ifelse(FACS == 'percent_constrained', 'FACS - constrained',
                              FACS)))
```

```
## `summarise()` regrouping output by 'Case_ID' (override with `.groups` argument)
```

```r
facs_results <-mixture_results_all %>%
  filter(!grepl('Mix', Sample_Name),
         algorithm != 'actual') %>%
  mutate(Case_ID = str_extract(Sample_Name, 'PM[0-9]+'))  %>%
  
  # constrain to non-troph pops
  filter(component != 'Trophoblasts') %>%
  group_by(Case_ID, algorithm) %>%
  mutate(percent = percent*100/sum(percent)) %>%
  inner_join(facs_prop, by = c('component', 'Case_ID'))
  
plot_scatter <- function(x){
  ggplot(data = x, aes(x = actual, y = percent, color =component)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  facet_grid(cols = vars(algorithm), rows = vars(component), scales = 'fixed') +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = color_code_tissue[unique(facs_results$component)], 
                         na.value = '#636363',
                         guide = 'none') +
  scale_x_continuous(labels = function(x)paste0(x, '%'))+
  scale_y_continuous(labels = function(x)paste0(x, '%')) +
  labs(x = 'Actual', y = 'Predicted') +
    coord_equal()
}

facs_results %>%
  filter(FACS == 'FACS',
         component == 'Endothelial')  %>%
  plot_scatter()
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
facs_results %>%
  filter(FACS == 'FACS',
         component == 'Hofbauer')  %>%
  plot_scatter()+
  scale_y_continuous(limits = c(0, 60), labels = function(x)paste0(x, '%')) 
```

```
## Scale for 'y' is already present. Adding another scale for 'y', which will
## replace the existing scale.
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

```r
facs_results %>%
  filter(FACS == 'FACS',
         component == 'Stromal')  %>%
  plot_scatter()
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-9-3.png)<!-- -->

```r
# stats
facs_stats <- facs_results %>%
  filter(FACS == 'FACS - constrained',
         component != 'Trophoblasts') %>%
  nest(data = -c(algorithm, component)) %>%
  mutate(lm = map(data, ~lm(percent ~ actual, .)),
         glanced = map(lm, glance)) %>%
  select(-data, -lm) %>%
  unnest(glanced) %>%
  mutate(label = paste0('R2=', signif(r.squared, digits = 2),
                        '\n', pvalue(p.value, add_p = TRUE)));facs_stats
```

```
## # A tibble: 12 x 14
## # Groups:   algorithm [4]
##    component algorithm r.squared adj.r.squared sigma statistic p.value    df
##    <chr>     <fct>         <dbl>         <dbl> <dbl>     <dbl>   <dbl> <int>
##  1 Stromal   minfi (C~   0.0660        0.0142  10.7      1.27   0.274      2
##  2 Hofbauer  minfi (C~   0.152         0.105    5.33     3.24   0.0888     2
##  3 Endothel~ minfi (C~   0.0143       -0.0405  13.0      0.261  0.616      2
##  4 Stromal   epidish ~   0.0597        0.00749 18.2      1.14   0.299      2
##  5 Hofbauer  epidish ~   0.00928      -0.0458   7.00     0.169  0.686      2
##  6 Endothel~ epidish ~   0.0144       -0.0404  23.1      0.262  0.615      2
##  7 Stromal   epidish ~   0.0439       -0.00923 15.6      0.826  0.375      2
##  8 Hofbauer  epidish ~   0.0547        0.00217  8.08     1.04   0.321      2
##  9 Endothel~ epidish ~   0.0232       -0.0310  20.8      0.428  0.521      2
## 10 Stromal   epidish ~   0.0945        0.0442  12.1      1.88   0.187      2
## 11 Hofbauer  epidish ~   0.113         0.0636   5.72     2.29   0.147      2
## 12 Endothel~ epidish ~   0.0176       -0.0370  14.6      0.322  0.578      2
## # ... with 6 more variables: logLik <dbl>, AIC <dbl>, BIC <dbl>,
## #   deviance <dbl>, df.residual <int>, label <chr>
```

```r
facs_results %>%
  filter(FACS == 'FACS - constrained',
         component != 'Trophoblasts') %>%
  group_by(algorithm) %>%
  rmse(truth = actual, estimate = percent)
```

```
## # A tibble: 4 x 4
##   algorithm     .metric .estimator .estimate
##   <fct>         <chr>   <chr>          <dbl>
## 1 epidish (RPC) rmse    standard        31.1
## 2 epidish (CBS) rmse    standard        28.0
## 3 epidish (CP)  rmse    standard        23.7
## 4 minfi (CP)    rmse    standard        22.8
```

```r
facs_results %>%
  filter(FACS == 'FACS - constrained',
         component != 'Trophoblasts') %>%
  ggplot(aes(x = actual, y = percent, color =component)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +  
  geom_label(data = facs_stats, 
             color = 'black', size = 3, 
             aes(label = label),
             x = 10, y = 120, hjust = 0, vjust = 1,
             alpha = 0) +

  facet_grid(cols = vars(algorithm), rows = vars(component), scales = 'fixed') +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = color_code_tissue[unique(facs_results$component)], 
                         na.value = '#636363',
                         guide = 'none') +
  scale_x_continuous(labels = function(x)paste0(x, '%'))+
  scale_y_continuous(limits = c(0, 120),labels = function(x)paste0(x, '%')) +
  labs(x = 'Actual', y = 'Predicted') + 
  coord_equal()
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-9-4.png)<!-- -->

```r
facs_stats %>%
  group_by(component, algorithm) %>%
  summarize(r2 = mean(r.squared)) %>%
  group_by(algorithm) %>%
  mutate(overall_r2 = mean(r2)) %>%
  ungroup() %>%
  arrange(algorithm, component)
```

```
## `summarise()` regrouping output by 'component' (override with `.groups` argument)
```

```
## # A tibble: 12 x 4
##    component   algorithm          r2 overall_r2
##    <chr>       <fct>           <dbl>      <dbl>
##  1 Endothelial epidish (RPC) 0.0144      0.0278
##  2 Hofbauer    epidish (RPC) 0.00928     0.0278
##  3 Stromal     epidish (RPC) 0.0597      0.0278
##  4 Endothelial epidish (CBS) 0.0232      0.0406
##  5 Hofbauer    epidish (CBS) 0.0547      0.0406
##  6 Stromal     epidish (CBS) 0.0439      0.0406
##  7 Endothelial epidish (CP)  0.0176      0.0750
##  8 Hofbauer    epidish (CP)  0.113       0.0750
##  9 Stromal     epidish (CP)  0.0945      0.0750
## 10 Endothelial minfi (CP)    0.0143      0.0776
## 11 Hofbauer    minfi (CP)    0.152       0.0776
## 12 Stromal     minfi (CP)    0.0660      0.0776
```

# Villi


```r
# both first and third barplots
mixture_results_all %>%
  left_join(pDat_filt %>% select(Sample_Name, Tissue, Trimester)) %>%
  filter(algorithm != 'actual',
         Tissue == 'Villi')   %>%
  ggplot() +
  geom_bar(stat=  'identity',
           aes(x = Sample_Name, y = percent, fill = component)) +
  facet_grid(rows = vars(algorithm), cols = vars(Trimester), 
             scales = 'free', space = 'free') +
  labs(x = '') +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = color_code_tissue[unique(mixture_results_all$component)], 
                    na.value = '#636363') +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1), limits = c(0,100))
```

```
## Joining, by = "Sample_Name"
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
# averaged
mixture_results_all %>%
  left_join(pDat_filt %>% select(Sample_Name, Tissue, Trimester)) %>%
  filter(algorithm != 'actual',
         Tissue == 'Villi')   %>%
  
  group_by(Trimester, algorithm, component) %>%
  summarize(percent = mean(percent)) %>%
  ggplot() +
  geom_bar(stat=  'identity',
           aes(x = Trimester, y = percent, fill = component)) +
  facet_grid(rows = vars(algorithm), scales = 'free', space = 'free') +
  labs(x = '') +
  theme(strip.background = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = color_code_tissue[unique(mixture_results_all$component)], 
                    na.value = '#636363') +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1), limits = c(0,100))
```

```
## Joining, by = "Sample_Name"
```

```
## `summarise()` regrouping output by 'Trimester', 'algorithm' (override with `.groups` argument)
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-10-2.png)<!-- -->

```r
# between algorithms
mixture_results_all %>%
  left_join(pDat_filt %>% select(Sample_Name, Tissue, Trimester)) %>%
  filter(algorithm != 'actual',
         Tissue == 'Villi')   %>%
  ggplot() +
  geom_boxplot(aes(x = algorithm, y = percent, fill = component)) +
  facet_grid(cols = vars(component), rows = vars(Trimester), space = 'free') +
  labs(x = '', fill = '') +
  theme(strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = color_code_tissue[unique(mixture_results_all$component)], 
                    na.value = '#636363',
                    guide = 'none') +
  scale_x_discrete(labels = function(x)gsub('epidish', 'EpiDISH', x)) +
  scale_y_continuous(labels = function(x)paste0(x, '%'))
```

```
## Joining, by = "Sample_Name"
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-10-3.png)<!-- -->

```r
# first to third beeeswarm
mixture_results_all %>%
  left_join(pDat_filt %>% select(Sample_Name, Tissue, Trimester)) %>%
  filter(algorithm != 'actual',
         Tissue == 'Villi') %>%
  ggplot(aes(x = component, y = percent, color = Trimester)) +
  ggbeeswarm::geom_beeswarm(dodge.width = 0.75) +
  facet_grid(rows = vars(algorithm)) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(labels = function(x)paste0(x, '%'))
```

```
## Joining, by = "Sample_Name"
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-10-4.png)<!-- -->

# Add nRBC + nRBC to deconvolution

and apply to villi


```r
mset_filt
```

```
## class: MethylSet 
## dim: 866091 138 
## metadata(0):
## assays(2): Meth Unmeth
## rownames(866091): cg18478105 cg09835024 ... cg10633746 cg12623625
## rowData names(0):
## colnames: NULL
## colData names(68): Sample_Name Chip_number ...
##   maternal_contamination_norm_flip GA
## Annotation
##   array: IlluminaHumanMethylationEPIC
##   annotation: ilm10b4.hg19
## Preprocessing
##   Method: NA
##   minfi version: NA
##   Manifest version: NA
```

```r
nrbc_noob <- readRDS(here::here(base_path, '2_14_nrbc_noob.rds'))

pData(mset)$Study <- 'Yuan'

# combine nrbc to placenta data
mset_combined <- combineArrays(nrbc_noob, mset)
```

```
## [convertArray] Casting as IlluminaHumanMethylation450k
```

```
## Loading required package: IlluminaHumanMethylation450kmanifest
```

```r
mset_combined
```

```
## class: MethylSet 
## dim: 452567 203 
## metadata(0):
## assays(2): Meth Unmeth
## rownames(452567): cg00212031 cg00213748 ... ch.22.47579720R
##   ch.22.48274842R
## rowData names(0):
## colnames: NULL
## colData names(23): Sample_Name Tissue ... CellType ArrayTypes
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19
## Preprocessing
##   Method: NA
##   minfi version: NA
##   Manifest version: NA
```

```r
pData(mset_combined) <- pData(mset_combined) %>% as.data.frame() %>%
  select(Sample_Name, Sex, Trimester, Tissue, Study) %>% 
  mutate(Tissue = gsub(' cs', '', Tissue),
         CellType = Tissue) %>% DataFrame()
  
pData(mset_combined) 
```

```
## DataFrame with 203 rows and 6 columns
##                  Sample_Name         Sex   Trimester                     Tissue
##                  <character> <character> <character>                <character>
## 1                   BakFW001           M       Third                       nRBC
## 2                   BakFW007           F       Third                       nRBC
## 3                   BakFW002           F       Third                       nRBC
## 4                   BakFW016           F       Third                       nRBC
## 5   nRBC_stringent_FACS_ind6           F       Third                       nRBC
## ...                      ...         ...         ...                        ...
## 199                 PL149_vc           F      Second                      Villi
## 200              PM366_vc_R2           F       Third                      Villi
## 201                  PL295_v           F       First                      Villi
## 202         PM372_discard_cs           M       Third Dead Cells and Lymphocytes
## 203            PM374_hofb_cs           M       Third                   Hofbauer
##           Study                   CellType
##     <character>                <character>
## 1      Bakulski                       nRBC
## 2      Bakulski                       nRBC
## 3      Bakulski                       nRBC
## 4      Bakulski                       nRBC
## 5       deGoede                       nRBC
## ...         ...                        ...
## 199        Yuan                      Villi
## 200        Yuan                      Villi
## 201        Yuan                      Villi
## 202        Yuan Dead Cells and Lymphocytes
## 203        Yuan                   Hofbauer
```

```r
pDat_combined <- pData(mset_combined) %>%
  as_tibble() 
  
# subset to reference
mset_combined_ref <- mset_combined[,pData(mset_combined)$CellType %in% c('Trophoblasts', 'Stromal', 'Hofbauer', 'Endothelial', 'nRBC', 'Syncytiotrophoblast')]

# pick probes
probes_third_nrbc <- minfi:::pickCompProbes(
  mset_combined_ref[, pData(mset_combined_ref)$Trimester == 'Third'], 
  cellTypes = c('Trophoblasts', 'Stromal', 'Hofbauer', 'Endothelial', 'nRBC', 'Syncytiotrophoblast'),
  compositeCellType = 'Placenta',
  probeSelect = "both")

probes_first_nrbc <- minfi:::pickCompProbes(
  mset_combined_ref[, pData(mset_combined_ref)$Trimester == 'First' |
             pData(mset_combined_ref)$CellType %in% c('Syncytiotrophoblast', 'nRBC')], 
  cellTypes = c('Trophoblasts', 'Stromal', 'Hofbauer', 'Endothelial', 'nRBC',  'Syncytiotrophoblast'),
  compositeCellType = 'Placenta',
  probeSelect = "both")

# extract coefficients
coefs_combined_third <- probes_third_nrbc$coefEsts
coefs_combined_first <- probes_first_nrbc$coefEsts
```

## Heatmap


```r
pheatmap(coefs_combined_third, cluster_cols = FALSE, cluster_rows = TRUE,
         show_rownames = FALSE) 
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
pheatmap(coefs_combined_first, cluster_cols = FALSE, cluster_rows = TRUE,
         show_rownames = FALSE) 
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-12-2.png)<!-- -->


## Generate in silico mixtures

For each cell type ("a"),
generate proportions randomly from 0-1.
cell type b proportion is drawn from 0-a
cell type c proportion is drawn from 0-(b+a)
cell type d is 1-(a+b+c)

These are the coefficients for generating mixtures with proportions uniformly distributed between 0-1 for cell type A. 

Apply this to generate mixtures for each cell type.



```r
#### generate in silico mixtures
# filter to coefficient probes
mset_combined_coefs <- 
  mset_combined[unique(c(rownames(coefs_combined_third), rownames(coefs_combined_first))),]

colnames(mset_combined_coefs) <- pData(mset_combined_coefs)$Sample_Name

getBeta(mset_combined_coefs) %>%
  as.data.frame() %>%
  bind_cols(cpg = rownames(.), .) %>%
  pivot_longer(cols = -cpg, names_to = 'Sample_Name', values_to = 'beta') %>%
  left_join(pDat_combined) %>%
  nest(data = c(cpg, beta))
```

```
## Joining, by = "Sample_Name"
```

```
## # A tibble: 203 x 7
##    Sample_Name           Sex   Trimester Tissue Study   CellType data           
##    <chr>                 <chr> <chr>     <chr>  <chr>   <chr>    <list>         
##  1 BakFW001              M     Third     nRBC   Bakuls~ nRBC     <tibble [1,050~
##  2 BakFW007              F     Third     nRBC   Bakuls~ nRBC     <tibble [1,050~
##  3 BakFW002              F     Third     nRBC   Bakuls~ nRBC     <tibble [1,050~
##  4 BakFW016              F     Third     nRBC   Bakuls~ nRBC     <tibble [1,050~
##  5 nRBC_stringent_FACS_~ F     Third     nRBC   deGoede nRBC     <tibble [1,050~
##  6 nRBC_stringent_FACS_~ M     Third     nRBC   deGoede nRBC     <tibble [1,050~
##  7 nRBC_stringent_FACS_~ M     Third     nRBC   deGoede nRBC     <tibble [1,050~
##  8 nRBC_stringent_FACS_~ F     Third     nRBC   deGoede nRBC     <tibble [1,050~
##  9 nRBC_stringent_FACS_~ F     Third     nRBC   deGoede nRBC     <tibble [1,050~
## 10 nRBC_stringent_FACS_~ F     Third     nRBC   deGoede nRBC     <tibble [1,050~
## # ... with 193 more rows
```

```r
# generate beta matrices for each trimester + cell type
b_3_troph <- 
  getBeta(mset_combined_coefs)[,pData(mset_combined_coefs)$Trimester == 'Third' &
                                 pData(mset_combined_coefs)$CellType == 'Trophoblasts']
b_3_endo <- 
  getBeta(mset_combined_coefs)[,pData(mset_combined_coefs)$Trimester == 'Third' &
                                 pData(mset_combined_coefs)$CellType == 'Endothelial']
b_3_strom <- 
  getBeta(mset_combined_coefs)[,pData(mset_combined_coefs)$Trimester == 'Third' &
                                 pData(mset_combined_coefs)$CellType == 'Stromal']
b_3_hofb <- 
  getBeta(mset_combined_coefs)[,pData(mset_combined_coefs)$Trimester == 'Third' &
                                 pData(mset_combined_coefs)$CellType == 'Hofbauer']
b_3_nrbc <- 
  getBeta(mset_combined_coefs)[,pData(mset_combined_coefs)$Trimester == 'Third' &
                                 pData(mset_combined_coefs)$CellType == 'nRBC']
b_3_stb <- 
  getBeta(mset_combined_coefs)[,pData(mset_combined_coefs)$Trimester == 'Third' &
                                 pData(mset_combined_coefs)$CellType == 'Syncytiotrophoblast']

b_1_troph <- 
  getBeta(mset_combined_coefs)[,pData(mset_combined_coefs)$Trimester == 'Third' &
                                 pData(mset_combined_coefs)$CellType == 'Trophoblasts']
b_1_endo <- 
  getBeta(mset_combined_coefs)[,pData(mset_combined_coefs)$Trimester == 'Third' &
                                 pData(mset_combined_coefs)$CellType == 'Endothelial']
b_1_strom <- 
  getBeta(mset_combined_coefs)[,pData(mset_combined_coefs)$Trimester == 'Third' &
                                 pData(mset_combined_coefs)$CellType == 'Stromal']
b_1_hofb <- 
  getBeta(mset_combined_coefs)[,pData(mset_combined_coefs)$Trimester == 'Third' &
                                 pData(mset_combined_coefs)$CellType == 'Hofbauer']

# generate coefficients
set.seed(1993)
n <- 250
a <- runif(n = n)
b <- runif(n = n, min = 0, max = 1-a)
c <- runif(n = n, min = 0, max = 1-(b + a))
d <- runif(n = n, min = 0, max = 1-(c + b + a))
e <- runif(n = n, min = 0, max = 1-(d + c + b + a))
f <- 1-(a+b+c+d+e)

all(a + b + c + d + e + f == 1)
```

```
## [1] TRUE
```

```r
# put coefficients into one table
props <- tibble(a = a, b=b, c=c, d=d, e=e, f=f)
props <- props %>%
  as.matrix() %>%
  rbind(., .[,c(2:6,1)], 
        .[,c(3:6,1:2)], 
        .[,c(4:6,1:3)], 
        .[,c(5:6,1:4)], 
        .[,c(6,1:5)]) %>%
  as_tibble() %>%
  mutate(sample = 1:(n*6))

# generate mixtures
set.seed(1)
silico_mix <- props %>%
  group_by(sample) %>%
  mutate(mix3 = list(
           a*b_3_troph[,sample(1:ncol(b_3_troph),1)] + 
           b*b_3_strom[,sample(1:ncol(b_3_strom),1)] +
           c*b_3_endo[,sample(1:ncol(b_3_endo),1)] +
           d*b_3_hofb[,sample(1:ncol(b_3_hofb),1)] +
           e*b_3_nrbc[,sample(1:ncol(b_3_nrbc),1)] +
           f*b_3_stb[,sample(1:ncol(b_3_stb),1)]),
         mix1 = list(
           a*b_1_troph[,sample(1:ncol(b_1_troph),1)] + 
           b*b_1_strom[,sample(1:ncol(b_1_strom),1)] +
           c*b_1_endo[,sample(1:ncol(b_1_endo),1)] +
           d*b_1_hofb[,sample(1:ncol(b_1_hofb),1)] +
           e*b_3_nrbc[,sample(1:ncol(b_3_nrbc),1)] +
           f*b_3_stb[,sample(1:ncol(b_3_stb),1)]),
         
         cpg= list(rownames(mset_combined_coefs))) %>% ungroup()

# process results
silico_pDat <- silico_mix %>%
  select(a:sample) %>%
  mutate(mix3 = 'mix3', mix1 = 'mix1') %>%
   pivot_longer(cols = c(mix3, mix1),
                names_to = 'sample2',
                values_to = 'remove') %>%
  mutate(sample_name = paste0('sample', sample, '_', sample2),
         trimester = ifelse(grepl('mix3', sample_name), 'Third', 'First')) %>%
  select(-sample, -sample2, -remove) 
  

silico_data <- silico_mix %>%
  select(sample, mix3:cpg) %>%
  unnest(c(mix3, mix1, cpg)) %>%
  
  pivot_longer(cols = c(mix3, mix1),
               names_to = 'sample2',
               values_to = 'beta') %>%
  mutate(sample_name = paste0('sample', sample, '_', sample2)) %>%
  select(-sample, -sample2) %>%
  pivot_wider(id_cols = cpg, names_from = 'sample_name', values_from = 'beta' ) %>%
  as.data.frame() %>%
  magrittr::set_rownames(.$cpg) %>%
  select(-cpg) %>%
  as.matrix()
```

Check composition distribution


```r
silico_pDat %>%
  pivot_longer(cols = a:f,
               names_to = 'component',
               values_to = 'proportion') %>%
  mutate(component = case_when(
    component == 'a' ~ 'Trophoblasts',
    component == 'b' ~ 'Stromal',
    component == 'c' ~ 'Endothelial',
    component == 'd' ~ 'Hofbauer',
    component == 'e' ~ 'nRBC',
    component == 'f' ~ 'Syncytiotrophoblast'
  )) %>% {
    ggplot(data =., aes(x = proportion, color = component)) +
  geom_density() +
  facet_grid(rows = vars(component), cols = vars(trimester))  +
  scale_color_manual(values = color_code_tissue[unique(.$component)], 
                    na.value = '#636363',
                    guide = 'none') 
    }
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


## apply to in silico mixtures



```r
# houseman
houseman_silico <- minfi:::projectCellType(
  silico_data[rownames(coefs_combined_third),], 
  coefs_combined_third,
  lessThanOne = FALSE) %>%
  as.data.frame() %>%
  mutate(Sample_Names = rownames(.), .) %>%
  as_tibble() %>%
  pivot_longer(cols = -Sample_Names,
               names_to = 'component',
               values_to = 'percent') %>%
  mutate(algorithm = 'Houseman (CP)')

epidish_res111 <- epidish(
  beta.m = silico_data,
  ref.m = coefs_combined_third,
  method = 'RPC')

epidish_res111 <- epidish_res111$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (RPC)')

epidish_res222 <- epidish(
  beta.m = silico_data,
  ref.m = coefs_combined_third,
  method = 'CBS')
```

```
## 1
```

```
## 2
```

```
## 3
```

```r
epidish_res222 <- epidish_res222$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (CBS)')

epidish_res333 <- epidish(
  beta.m = silico_data,
  ref.m = coefs_combined_third,
  method = 'CP')
```

```
## 1
```

```
## 2
```

```
## 3
```

```
## 4
```

```
## 5
```

```
## 6
```

```
## 7
```

```
## 8
```

```
## 9
```

```
## 10
```

```
## 11
```

```
## 12
```

```
## 13
```

```
## 14
```

```
## 15
```

```
## 16
```

```
## 17
```

```
## 18
```

```
## 19
```

```
## 20
```

```
## 21
```

```
## 22
```

```
## 23
```

```
## 24
```

```
## 25
```

```
## 26
```

```
## 27
```

```
## 28
```

```
## 29
```

```
## 30
```

```
## 31
```

```
## 32
```

```
## 33
```

```
## 34
```

```
## 35
```

```
## 36
```

```
## 37
```

```
## 38
```

```
## 39
```

```
## 40
```

```
## 41
```

```
## 42
```

```
## 43
```

```
## 44
```

```
## 45
```

```
## 46
```

```
## 47
```

```
## 48
```

```
## 49
```

```
## 50
```

```
## 51
```

```
## 52
```

```
## 53
```

```
## 54
```

```
## 55
```

```
## 56
```

```
## 57
```

```
## 58
```

```
## 59
```

```
## 60
```

```
## 61
```

```
## 62
```

```
## 63
```

```
## 64
```

```
## 65
```

```
## 66
```

```
## 67
```

```
## 68
```

```
## 69
```

```
## 70
```

```
## 71
```

```
## 72
```

```
## 73
```

```
## 74
```

```
## 75
```

```
## 76
```

```
## 77
```

```
## 78
```

```
## 79
```

```
## 80
```

```
## 81
```

```
## 82
```

```
## 83
```

```
## 84
```

```
## 85
```

```
## 86
```

```
## 87
```

```
## 88
```

```
## 89
```

```
## 90
```

```
## 91
```

```
## 92
```

```
## 93
```

```
## 94
```

```
## 95
```

```
## 96
```

```
## 97
```

```
## 98
```

```
## 99
```

```
## 100
```

```
## 101
```

```
## 102
```

```
## 103
```

```
## 104
```

```
## 105
```

```
## 106
```

```
## 107
```

```
## 108
```

```
## 109
```

```
## 110
```

```
## 111
```

```
## 112
```

```
## 113
```

```
## 114
```

```
## 115
```

```
## 116
```

```
## 117
```

```
## 118
```

```
## 119
```

```
## 120
```

```
## 121
```

```
## 122
```

```
## 123
```

```
## 124
```

```
## 125
```

```
## 126
```

```
## 127
```

```
## 128
```

```
## 129
```

```
## 130
```

```
## 131
```

```
## 132
```

```
## 133
```

```
## 134
```

```
## 135
```

```
## 136
```

```
## 137
```

```
## 138
```

```
## 139
```

```
## 140
```

```
## 141
```

```
## 142
```

```
## 143
```

```
## 144
```

```
## 145
```

```
## 146
```

```
## 147
```

```
## 148
```

```
## 149
```

```
## 150
```

```
## 151
```

```
## 152
```

```
## 153
```

```
## 154
```

```
## 155
```

```
## 156
```

```
## 157
```

```
## 158
```

```
## 159
```

```
## 160
```

```
## 161
```

```
## 162
```

```
## 163
```

```
## 164
```

```
## 165
```

```
## 166
```

```
## 167
```

```
## 168
```

```
## 169
```

```
## 170
```

```
## 171
```

```
## 172
```

```
## 173
```

```
## 174
```

```
## 175
```

```
## 176
```

```
## 177
```

```
## 178
```

```
## 179
```

```
## 180
```

```
## 181
```

```
## 182
```

```
## 183
```

```
## 184
```

```
## 185
```

```
## 186
```

```
## 187
```

```
## 188
```

```
## 189
```

```
## 190
```

```
## 191
```

```
## 192
```

```
## 193
```

```
## 194
```

```
## 195
```

```
## 196
```

```
## 197
```

```
## 198
```

```
## 199
```

```
## 200
```

```
## 201
```

```
## 202
```

```
## 203
```

```
## 204
```

```
## 205
```

```
## 206
```

```
## 207
```

```
## 208
```

```
## 209
```

```
## 210
```

```
## 211
```

```
## 212
```

```
## 213
```

```
## 214
```

```
## 215
```

```
## 216
```

```
## 217
```

```
## 218
```

```
## 219
```

```
## 220
```

```
## 221
```

```
## 222
```

```
## 223
```

```
## 224
```

```
## 225
```

```
## 226
```

```
## 227
```

```
## 228
```

```
## 229
```

```
## 230
```

```
## 231
```

```
## 232
```

```
## 233
```

```
## 234
```

```
## 235
```

```
## 236
```

```
## 237
```

```
## 238
```

```
## 239
```

```
## 240
```

```
## 241
```

```
## 242
```

```
## 243
```

```
## 244
```

```
## 245
```

```
## 246
```

```
## 247
```

```
## 248
```

```
## 249
```

```
## 250
```

```
## 251
```

```
## 252
```

```
## 253
```

```
## 254
```

```
## 255
```

```
## 256
```

```
## 257
```

```
## 258
```

```
## 259
```

```
## 260
```

```
## 261
```

```
## 262
```

```
## 263
```

```
## 264
```

```
## 265
```

```
## 266
```

```
## 267
```

```
## 268
```

```
## 269
```

```
## 270
```

```
## 271
```

```
## 272
```

```
## 273
```

```
## 274
```

```
## 275
```

```
## 276
```

```
## 277
```

```
## 278
```

```
## 279
```

```
## 280
```

```
## 281
```

```
## 282
```

```
## 283
```

```
## 284
```

```
## 285
```

```
## 286
```

```
## 287
```

```
## 288
```

```
## 289
```

```
## 290
```

```
## 291
```

```
## 292
```

```
## 293
```

```
## 294
```

```
## 295
```

```
## 296
```

```
## 297
```

```
## 298
```

```
## 299
```

```
## 300
```

```
## 301
```

```
## 302
```

```
## 303
```

```
## 304
```

```
## 305
```

```
## 306
```

```
## 307
```

```
## 308
```

```
## 309
```

```
## 310
```

```
## 311
```

```
## 312
```

```
## 313
```

```
## 314
```

```
## 315
```

```
## 316
```

```
## 317
```

```
## 318
```

```
## 319
```

```
## 320
```

```
## 321
```

```
## 322
```

```
## 323
```

```
## 324
```

```
## 325
```

```
## 326
```

```
## 327
```

```
## 328
```

```
## 329
```

```
## 330
```

```
## 331
```

```
## 332
```

```
## 333
```

```
## 334
```

```
## 335
```

```
## 336
```

```
## 337
```

```
## 338
```

```
## 339
```

```
## 340
```

```
## 341
```

```
## 342
```

```
## 343
```

```
## 344
```

```
## 345
```

```
## 346
```

```
## 347
```

```
## 348
```

```
## 349
```

```
## 350
```

```
## 351
```

```
## 352
```

```
## 353
```

```
## 354
```

```
## 355
```

```
## 356
```

```
## 357
```

```
## 358
```

```
## 359
```

```
## 360
```

```
## 361
```

```
## 362
```

```
## 363
```

```
## 364
```

```
## 365
```

```
## 366
```

```
## 367
```

```
## 368
```

```
## 369
```

```
## 370
```

```
## 371
```

```
## 372
```

```
## 373
```

```
## 374
```

```
## 375
```

```
## 376
```

```
## 377
```

```
## 378
```

```
## 379
```

```
## 380
```

```
## 381
```

```
## 382
```

```
## 383
```

```
## 384
```

```
## 385
```

```
## 386
```

```
## 387
```

```
## 388
```

```
## 389
```

```
## 390
```

```
## 391
```

```
## 392
```

```
## 393
```

```
## 394
```

```
## 395
```

```
## 396
```

```
## 397
```

```
## 398
```

```
## 399
```

```
## 400
```

```
## 401
```

```
## 402
```

```
## 403
```

```
## 404
```

```
## 405
```

```
## 406
```

```
## 407
```

```
## 408
```

```
## 409
```

```
## 410
```

```
## 411
```

```
## 412
```

```
## 413
```

```
## 414
```

```
## 415
```

```
## 416
```

```
## 417
```

```
## 418
```

```
## 419
```

```
## 420
```

```
## 421
```

```
## 422
```

```
## 423
```

```
## 424
```

```
## 425
```

```
## 426
```

```
## 427
```

```
## 428
```

```
## 429
```

```
## 430
```

```
## 431
```

```
## 432
```

```
## 433
```

```
## 434
```

```
## 435
```

```
## 436
```

```
## 437
```

```
## 438
```

```
## 439
```

```
## 440
```

```
## 441
```

```
## 442
```

```
## 443
```

```
## 444
```

```
## 445
```

```
## 446
```

```
## 447
```

```
## 448
```

```
## 449
```

```
## 450
```

```
## 451
```

```
## 452
```

```
## 453
```

```
## 454
```

```
## 455
```

```
## 456
```

```
## 457
```

```
## 458
```

```
## 459
```

```
## 460
```

```
## 461
```

```
## 462
```

```
## 463
```

```
## 464
```

```
## 465
```

```
## 466
```

```
## 467
```

```
## 468
```

```
## 469
```

```
## 470
```

```
## 471
```

```
## 472
```

```
## 473
```

```
## 474
```

```
## 475
```

```
## 476
```

```
## 477
```

```
## 478
```

```
## 479
```

```
## 480
```

```
## 481
```

```
## 482
```

```
## 483
```

```
## 484
```

```
## 485
```

```
## 486
```

```
## 487
```

```
## 488
```

```
## 489
```

```
## 490
```

```
## 491
```

```
## 492
```

```
## 493
```

```
## 494
```

```
## 495
```

```
## 496
```

```
## 497
```

```
## 498
```

```
## 499
```

```
## 500
```

```
## 501
```

```
## 502
```

```
## 503
```

```
## 504
```

```
## 505
```

```
## 506
```

```
## 507
```

```
## 508
```

```
## 509
```

```
## 510
```

```
## 511
```

```
## 512
```

```
## 513
```

```
## 514
```

```
## 515
```

```
## 516
```

```
## 517
```

```
## 518
```

```
## 519
```

```
## 520
```

```
## 521
```

```
## 522
```

```
## 523
```

```
## 524
```

```
## 525
```

```
## 526
```

```
## 527
```

```
## 528
```

```
## 529
```

```
## 530
```

```
## 531
```

```
## 532
```

```
## 533
```

```
## 534
```

```
## 535
```

```
## 536
```

```
## 537
```

```
## 538
```

```
## 539
```

```
## 540
```

```
## 541
```

```
## 542
```

```
## 543
```

```
## 544
```

```
## 545
```

```
## 546
```

```
## 547
```

```
## 548
```

```
## 549
```

```
## 550
```

```
## 551
```

```
## 552
```

```
## 553
```

```
## 554
```

```
## 555
```

```
## 556
```

```
## 557
```

```
## 558
```

```
## 559
```

```
## 560
```

```
## 561
```

```
## 562
```

```
## 563
```

```
## 564
```

```
## 565
```

```
## 566
```

```
## 567
```

```
## 568
```

```
## 569
```

```
## 570
```

```
## 571
```

```
## 572
```

```
## 573
```

```
## 574
```

```
## 575
```

```
## 576
```

```
## 577
```

```
## 578
```

```
## 579
```

```
## 580
```

```
## 581
```

```
## 582
```

```
## 583
```

```
## 584
```

```
## 585
```

```
## 586
```

```
## 587
```

```
## 588
```

```
## 589
```

```
## 590
```

```
## 591
```

```
## 592
```

```
## 593
```

```
## 594
```

```
## 595
```

```
## 596
```

```
## 597
```

```
## 598
```

```
## 599
```

```
## 600
```

```
## 601
```

```
## 602
```

```
## 603
```

```
## 604
```

```
## 605
```

```
## 606
```

```
## 607
```

```
## 608
```

```
## 609
```

```
## 610
```

```
## 611
```

```
## 612
```

```
## 613
```

```
## 614
```

```
## 615
```

```
## 616
```

```
## 617
```

```
## 618
```

```
## 619
```

```
## 620
```

```
## 621
```

```
## 622
```

```
## 623
```

```
## 624
```

```
## 625
```

```
## 626
```

```
## 627
```

```
## 628
```

```
## 629
```

```
## 630
```

```
## 631
```

```
## 632
```

```
## 633
```

```
## 634
```

```
## 635
```

```
## 636
```

```
## 637
```

```
## 638
```

```
## 639
```

```
## 640
```

```
## 641
```

```
## 642
```

```
## 643
```

```
## 644
```

```
## 645
```

```
## 646
```

```
## 647
```

```
## 648
```

```
## 649
```

```
## 650
```

```
## 651
```

```
## 652
```

```
## 653
```

```
## 654
```

```
## 655
```

```
## 656
```

```
## 657
```

```
## 658
```

```
## 659
```

```
## 660
```

```
## 661
```

```
## 662
```

```
## 663
```

```
## 664
```

```
## 665
```

```
## 666
```

```
## 667
```

```
## 668
```

```
## 669
```

```
## 670
```

```
## 671
```

```
## 672
```

```
## 673
```

```
## 674
```

```
## 675
```

```
## 676
```

```
## 677
```

```
## 678
```

```
## 679
```

```
## 680
```

```
## 681
```

```
## 682
```

```
## 683
```

```
## 684
```

```
## 685
```

```
## 686
```

```
## 687
```

```
## 688
```

```
## 689
```

```
## 690
```

```
## 691
```

```
## 692
```

```
## 693
```

```
## 694
```

```
## 695
```

```
## 696
```

```
## 697
```

```
## 698
```

```
## 699
```

```
## 700
```

```
## 701
```

```
## 702
```

```
## 703
```

```
## 704
```

```
## 705
```

```
## 706
```

```
## 707
```

```
## 708
```

```
## 709
```

```
## 710
```

```
## 711
```

```
## 712
```

```
## 713
```

```
## 714
```

```
## 715
```

```
## 716
```

```
## 717
```

```
## 718
```

```
## 719
```

```
## 720
```

```
## 721
```

```
## 722
```

```
## 723
```

```
## 724
```

```
## 725
```

```
## 726
```

```
## 727
```

```
## 728
```

```
## 729
```

```
## 730
```

```
## 731
```

```
## 732
```

```
## 733
```

```
## 734
```

```
## 735
```

```
## 736
```

```
## 737
```

```
## 738
```

```
## 739
```

```
## 740
```

```
## 741
```

```
## 742
```

```
## 743
```

```
## 744
```

```
## 745
```

```
## 746
```

```
## 747
```

```
## 748
```

```
## 749
```

```
## 750
```

```
## 751
```

```
## 752
```

```
## 753
```

```
## 754
```

```
## 755
```

```
## 756
```

```
## 757
```

```
## 758
```

```
## 759
```

```
## 760
```

```
## 761
```

```
## 762
```

```
## 763
```

```
## 764
```

```
## 765
```

```
## 766
```

```
## 767
```

```
## 768
```

```
## 769
```

```
## 770
```

```
## 771
```

```
## 772
```

```
## 773
```

```
## 774
```

```
## 775
```

```
## 776
```

```
## 777
```

```
## 778
```

```
## 779
```

```
## 780
```

```
## 781
```

```
## 782
```

```
## 783
```

```
## 784
```

```
## 785
```

```
## 786
```

```
## 787
```

```
## 788
```

```
## 789
```

```
## 790
```

```
## 791
```

```
## 792
```

```
## 793
```

```
## 794
```

```
## 795
```

```
## 796
```

```
## 797
```

```
## 798
```

```
## 799
```

```
## 800
```

```
## 801
```

```
## 802
```

```
## 803
```

```
## 804
```

```
## 805
```

```
## 806
```

```
## 807
```

```
## 808
```

```
## 809
```

```
## 810
```

```
## 811
```

```
## 812
```

```
## 813
```

```
## 814
```

```
## 815
```

```
## 816
```

```
## 817
```

```
## 818
```

```
## 819
```

```
## 820
```

```
## 821
```

```
## 822
```

```
## 823
```

```
## 824
```

```
## 825
```

```
## 826
```

```
## 827
```

```
## 828
```

```
## 829
```

```
## 830
```

```
## 831
```

```
## 832
```

```
## 833
```

```
## 834
```

```
## 835
```

```
## 836
```

```
## 837
```

```
## 838
```

```
## 839
```

```
## 840
```

```
## 841
```

```
## 842
```

```
## 843
```

```
## 844
```

```
## 845
```

```
## 846
```

```
## 847
```

```
## 848
```

```
## 849
```

```
## 850
```

```
## 851
```

```
## 852
```

```
## 853
```

```
## 854
```

```
## 855
```

```
## 856
```

```
## 857
```

```
## 858
```

```
## 859
```

```
## 860
```

```
## 861
```

```
## 862
```

```
## 863
```

```
## 864
```

```
## 865
```

```
## 866
```

```
## 867
```

```
## 868
```

```
## 869
```

```
## 870
```

```
## 871
```

```
## 872
```

```
## 873
```

```
## 874
```

```
## 875
```

```
## 876
```

```
## 877
```

```
## 878
```

```
## 879
```

```
## 880
```

```
## 881
```

```
## 882
```

```
## 883
```

```
## 884
```

```
## 885
```

```
## 886
```

```
## 887
```

```
## 888
```

```
## 889
```

```
## 890
```

```
## 891
```

```
## 892
```

```
## 893
```

```
## 894
```

```
## 895
```

```
## 896
```

```
## 897
```

```
## 898
```

```
## 899
```

```
## 900
```

```
## 901
```

```
## 902
```

```
## 903
```

```
## 904
```

```
## 905
```

```
## 906
```

```
## 907
```

```
## 908
```

```
## 909
```

```
## 910
```

```
## 911
```

```
## 912
```

```
## 913
```

```
## 914
```

```
## 915
```

```
## 916
```

```
## 917
```

```
## 918
```

```
## 919
```

```
## 920
```

```
## 921
```

```
## 922
```

```
## 923
```

```
## 924
```

```
## 925
```

```
## 926
```

```
## 927
```

```
## 928
```

```
## 929
```

```
## 930
```

```
## 931
```

```
## 932
```

```
## 933
```

```
## 934
```

```
## 935
```

```
## 936
```

```
## 937
```

```
## 938
```

```
## 939
```

```
## 940
```

```
## 941
```

```
## 942
```

```
## 943
```

```
## 944
```

```
## 945
```

```
## 946
```

```
## 947
```

```
## 948
```

```
## 949
```

```
## 950
```

```
## 951
```

```
## 952
```

```
## 953
```

```
## 954
```

```
## 955
```

```
## 956
```

```
## 957
```

```
## 958
```

```
## 959
```

```
## 960
```

```
## 961
```

```
## 962
```

```
## 963
```

```
## 964
```

```
## 965
```

```
## 966
```

```
## 967
```

```
## 968
```

```
## 969
```

```
## 970
```

```
## 971
```

```
## 972
```

```
## 973
```

```
## 974
```

```
## 975
```

```
## 976
```

```
## 977
```

```
## 978
```

```
## 979
```

```
## 980
```

```
## 981
```

```
## 982
```

```
## 983
```

```
## 984
```

```
## 985
```

```
## 986
```

```
## 987
```

```
## 988
```

```
## 989
```

```
## 990
```

```
## 991
```

```
## 992
```

```
## 993
```

```
## 994
```

```
## 995
```

```
## 996
```

```
## 997
```

```
## 998
```

```
## 999
```

```
## 1000
```

```
## 1001
```

```
## 1002
```

```
## 1003
```

```
## 1004
```

```
## 1005
```

```
## 1006
```

```
## 1007
```

```
## 1008
```

```
## 1009
```

```
## 1010
```

```
## 1011
```

```
## 1012
```

```
## 1013
```

```
## 1014
```

```
## 1015
```

```
## 1016
```

```
## 1017
```

```
## 1018
```

```
## 1019
```

```
## 1020
```

```
## 1021
```

```
## 1022
```

```
## 1023
```

```
## 1024
```

```
## 1025
```

```
## 1026
```

```
## 1027
```

```
## 1028
```

```
## 1029
```

```
## 1030
```

```
## 1031
```

```
## 1032
```

```
## 1033
```

```
## 1034
```

```
## 1035
```

```
## 1036
```

```
## 1037
```

```
## 1038
```

```
## 1039
```

```
## 1040
```

```
## 1041
```

```
## 1042
```

```
## 1043
```

```
## 1044
```

```
## 1045
```

```
## 1046
```

```
## 1047
```

```
## 1048
```

```
## 1049
```

```
## 1050
```

```
## 1051
```

```
## 1052
```

```
## 1053
```

```
## 1054
```

```
## 1055
```

```
## 1056
```

```
## 1057
```

```
## 1058
```

```
## 1059
```

```
## 1060
```

```
## 1061
```

```
## 1062
```

```
## 1063
```

```
## 1064
```

```
## 1065
```

```
## 1066
```

```
## 1067
```

```
## 1068
```

```
## 1069
```

```
## 1070
```

```
## 1071
```

```
## 1072
```

```
## 1073
```

```
## 1074
```

```
## 1075
```

```
## 1076
```

```
## 1077
```

```
## 1078
```

```
## 1079
```

```
## 1080
```

```
## 1081
```

```
## 1082
```

```
## 1083
```

```
## 1084
```

```
## 1085
```

```
## 1086
```

```
## 1087
```

```
## 1088
```

```
## 1089
```

```
## 1090
```

```
## 1091
```

```
## 1092
```

```
## 1093
```

```
## 1094
```

```
## 1095
```

```
## 1096
```

```
## 1097
```

```
## 1098
```

```
## 1099
```

```
## 1100
```

```
## 1101
```

```
## 1102
```

```
## 1103
```

```
## 1104
```

```
## 1105
```

```
## 1106
```

```
## 1107
```

```
## 1108
```

```
## 1109
```

```
## 1110
```

```
## 1111
```

```
## 1112
```

```
## 1113
```

```
## 1114
```

```
## 1115
```

```
## 1116
```

```
## 1117
```

```
## 1118
```

```
## 1119
```

```
## 1120
```

```
## 1121
```

```
## 1122
```

```
## 1123
```

```
## 1124
```

```
## 1125
```

```
## 1126
```

```
## 1127
```

```
## 1128
```

```
## 1129
```

```
## 1130
```

```
## 1131
```

```
## 1132
```

```
## 1133
```

```
## 1134
```

```
## 1135
```

```
## 1136
```

```
## 1137
```

```
## 1138
```

```
## 1139
```

```
## 1140
```

```
## 1141
```

```
## 1142
```

```
## 1143
```

```
## 1144
```

```
## 1145
```

```
## 1146
```

```
## 1147
```

```
## 1148
```

```
## 1149
```

```
## 1150
```

```
## 1151
```

```
## 1152
```

```
## 1153
```

```
## 1154
```

```
## 1155
```

```
## 1156
```

```
## 1157
```

```
## 1158
```

```
## 1159
```

```
## 1160
```

```
## 1161
```

```
## 1162
```

```
## 1163
```

```
## 1164
```

```
## 1165
```

```
## 1166
```

```
## 1167
```

```
## 1168
```

```
## 1169
```

```
## 1170
```

```
## 1171
```

```
## 1172
```

```
## 1173
```

```
## 1174
```

```
## 1175
```

```
## 1176
```

```
## 1177
```

```
## 1178
```

```
## 1179
```

```
## 1180
```

```
## 1181
```

```
## 1182
```

```
## 1183
```

```
## 1184
```

```
## 1185
```

```
## 1186
```

```
## 1187
```

```
## 1188
```

```
## 1189
```

```
## 1190
```

```
## 1191
```

```
## 1192
```

```
## 1193
```

```
## 1194
```

```
## 1195
```

```
## 1196
```

```
## 1197
```

```
## 1198
```

```
## 1199
```

```
## 1200
```

```
## 1201
```

```
## 1202
```

```
## 1203
```

```
## 1204
```

```
## 1205
```

```
## 1206
```

```
## 1207
```

```
## 1208
```

```
## 1209
```

```
## 1210
```

```
## 1211
```

```
## 1212
```

```
## 1213
```

```
## 1214
```

```
## 1215
```

```
## 1216
```

```
## 1217
```

```
## 1218
```

```
## 1219
```

```
## 1220
```

```
## 1221
```

```
## 1222
```

```
## 1223
```

```
## 1224
```

```
## 1225
```

```
## 1226
```

```
## 1227
```

```
## 1228
```

```
## 1229
```

```
## 1230
```

```
## 1231
```

```
## 1232
```

```
## 1233
```

```
## 1234
```

```
## 1235
```

```
## 1236
```

```
## 1237
```

```
## 1238
```

```
## 1239
```

```
## 1240
```

```
## 1241
```

```
## 1242
```

```
## 1243
```

```
## 1244
```

```
## 1245
```

```
## 1246
```

```
## 1247
```

```
## 1248
```

```
## 1249
```

```
## 1250
```

```
## 1251
```

```
## 1252
```

```
## 1253
```

```
## 1254
```

```
## 1255
```

```
## 1256
```

```
## 1257
```

```
## 1258
```

```
## 1259
```

```
## 1260
```

```
## 1261
```

```
## 1262
```

```
## 1263
```

```
## 1264
```

```
## 1265
```

```
## 1266
```

```
## 1267
```

```
## 1268
```

```
## 1269
```

```
## 1270
```

```
## 1271
```

```
## 1272
```

```
## 1273
```

```
## 1274
```

```
## 1275
```

```
## 1276
```

```
## 1277
```

```
## 1278
```

```
## 1279
```

```
## 1280
```

```
## 1281
```

```
## 1282
```

```
## 1283
```

```
## 1284
```

```
## 1285
```

```
## 1286
```

```
## 1287
```

```
## 1288
```

```
## 1289
```

```
## 1290
```

```
## 1291
```

```
## 1292
```

```
## 1293
```

```
## 1294
```

```
## 1295
```

```
## 1296
```

```
## 1297
```

```
## 1298
```

```
## 1299
```

```
## 1300
```

```
## 1301
```

```
## 1302
```

```
## 1303
```

```
## 1304
```

```
## 1305
```

```
## 1306
```

```
## 1307
```

```
## 1308
```

```
## 1309
```

```
## 1310
```

```
## 1311
```

```
## 1312
```

```
## 1313
```

```
## 1314
```

```
## 1315
```

```
## 1316
```

```
## 1317
```

```
## 1318
```

```
## 1319
```

```
## 1320
```

```
## 1321
```

```
## 1322
```

```
## 1323
```

```
## 1324
```

```
## 1325
```

```
## 1326
```

```
## 1327
```

```
## 1328
```

```
## 1329
```

```
## 1330
```

```
## 1331
```

```
## 1332
```

```
## 1333
```

```
## 1334
```

```
## 1335
```

```
## 1336
```

```
## 1337
```

```
## 1338
```

```
## 1339
```

```
## 1340
```

```
## 1341
```

```
## 1342
```

```
## 1343
```

```
## 1344
```

```
## 1345
```

```
## 1346
```

```
## 1347
```

```
## 1348
```

```
## 1349
```

```
## 1350
```

```
## 1351
```

```
## 1352
```

```
## 1353
```

```
## 1354
```

```
## 1355
```

```
## 1356
```

```
## 1357
```

```
## 1358
```

```
## 1359
```

```
## 1360
```

```
## 1361
```

```
## 1362
```

```
## 1363
```

```
## 1364
```

```
## 1365
```

```
## 1366
```

```
## 1367
```

```
## 1368
```

```
## 1369
```

```
## 1370
```

```
## 1371
```

```
## 1372
```

```
## 1373
```

```
## 1374
```

```
## 1375
```

```
## 1376
```

```
## 1377
```

```
## 1378
```

```
## 1379
```

```
## 1380
```

```
## 1381
```

```
## 1382
```

```
## 1383
```

```
## 1384
```

```
## 1385
```

```
## 1386
```

```
## 1387
```

```
## 1388
```

```
## 1389
```

```
## 1390
```

```
## 1391
```

```
## 1392
```

```
## 1393
```

```
## 1394
```

```
## 1395
```

```
## 1396
```

```
## 1397
```

```
## 1398
```

```
## 1399
```

```
## 1400
```

```
## 1401
```

```
## 1402
```

```
## 1403
```

```
## 1404
```

```
## 1405
```

```
## 1406
```

```
## 1407
```

```
## 1408
```

```
## 1409
```

```
## 1410
```

```
## 1411
```

```
## 1412
```

```
## 1413
```

```
## 1414
```

```
## 1415
```

```
## 1416
```

```
## 1417
```

```
## 1418
```

```
## 1419
```

```
## 1420
```

```
## 1421
```

```
## 1422
```

```
## 1423
```

```
## 1424
```

```
## 1425
```

```
## 1426
```

```
## 1427
```

```
## 1428
```

```
## 1429
```

```
## 1430
```

```
## 1431
```

```
## 1432
```

```
## 1433
```

```
## 1434
```

```
## 1435
```

```
## 1436
```

```
## 1437
```

```
## 1438
```

```
## 1439
```

```
## 1440
```

```
## 1441
```

```
## 1442
```

```
## 1443
```

```
## 1444
```

```
## 1445
```

```
## 1446
```

```
## 1447
```

```
## 1448
```

```
## 1449
```

```
## 1450
```

```
## 1451
```

```
## 1452
```

```
## 1453
```

```
## 1454
```

```
## 1455
```

```
## 1456
```

```
## 1457
```

```
## 1458
```

```
## 1459
```

```
## 1460
```

```
## 1461
```

```
## 1462
```

```
## 1463
```

```
## 1464
```

```
## 1465
```

```
## 1466
```

```
## 1467
```

```
## 1468
```

```
## 1469
```

```
## 1470
```

```
## 1471
```

```
## 1472
```

```
## 1473
```

```
## 1474
```

```
## 1475
```

```
## 1476
```

```
## 1477
```

```
## 1478
```

```
## 1479
```

```
## 1480
```

```
## 1481
```

```
## 1482
```

```
## 1483
```

```
## 1484
```

```
## 1485
```

```
## 1486
```

```
## 1487
```

```
## 1488
```

```
## 1489
```

```
## 1490
```

```
## 1491
```

```
## 1492
```

```
## 1493
```

```
## 1494
```

```
## 1495
```

```
## 1496
```

```
## 1497
```

```
## 1498
```

```
## 1499
```

```
## 1500
```

```
## 1501
```

```
## 1502
```

```
## 1503
```

```
## 1504
```

```
## 1505
```

```
## 1506
```

```
## 1507
```

```
## 1508
```

```
## 1509
```

```
## 1510
```

```
## 1511
```

```
## 1512
```

```
## 1513
```

```
## 1514
```

```
## 1515
```

```
## 1516
```

```
## 1517
```

```
## 1518
```

```
## 1519
```

```
## 1520
```

```
## 1521
```

```
## 1522
```

```
## 1523
```

```
## 1524
```

```
## 1525
```

```
## 1526
```

```
## 1527
```

```
## 1528
```

```
## 1529
```

```
## 1530
```

```
## 1531
```

```
## 1532
```

```
## 1533
```

```
## 1534
```

```
## 1535
```

```
## 1536
```

```
## 1537
```

```
## 1538
```

```
## 1539
```

```
## 1540
```

```
## 1541
```

```
## 1542
```

```
## 1543
```

```
## 1544
```

```
## 1545
```

```
## 1546
```

```
## 1547
```

```
## 1548
```

```
## 1549
```

```
## 1550
```

```
## 1551
```

```
## 1552
```

```
## 1553
```

```
## 1554
```

```
## 1555
```

```
## 1556
```

```
## 1557
```

```
## 1558
```

```
## 1559
```

```
## 1560
```

```
## 1561
```

```
## 1562
```

```
## 1563
```

```
## 1564
```

```
## 1565
```

```
## 1566
```

```
## 1567
```

```
## 1568
```

```
## 1569
```

```
## 1570
```

```
## 1571
```

```
## 1572
```

```
## 1573
```

```
## 1574
```

```
## 1575
```

```
## 1576
```

```
## 1577
```

```
## 1578
```

```
## 1579
```

```
## 1580
```

```
## 1581
```

```
## 1582
```

```
## 1583
```

```
## 1584
```

```
## 1585
```

```
## 1586
```

```
## 1587
```

```
## 1588
```

```
## 1589
```

```
## 1590
```

```
## 1591
```

```
## 1592
```

```
## 1593
```

```
## 1594
```

```
## 1595
```

```
## 1596
```

```
## 1597
```

```
## 1598
```

```
## 1599
```

```
## 1600
```

```
## 1601
```

```
## 1602
```

```
## 1603
```

```
## 1604
```

```
## 1605
```

```
## 1606
```

```
## 1607
```

```
## 1608
```

```
## 1609
```

```
## 1610
```

```
## 1611
```

```
## 1612
```

```
## 1613
```

```
## 1614
```

```
## 1615
```

```
## 1616
```

```
## 1617
```

```
## 1618
```

```
## 1619
```

```
## 1620
```

```
## 1621
```

```
## 1622
```

```
## 1623
```

```
## 1624
```

```
## 1625
```

```
## 1626
```

```
## 1627
```

```
## 1628
```

```
## 1629
```

```
## 1630
```

```
## 1631
```

```
## 1632
```

```
## 1633
```

```
## 1634
```

```
## 1635
```

```
## 1636
```

```
## 1637
```

```
## 1638
```

```
## 1639
```

```
## 1640
```

```
## 1641
```

```
## 1642
```

```
## 1643
```

```
## 1644
```

```
## 1645
```

```
## 1646
```

```
## 1647
```

```
## 1648
```

```
## 1649
```

```
## 1650
```

```
## 1651
```

```
## 1652
```

```
## 1653
```

```
## 1654
```

```
## 1655
```

```
## 1656
```

```
## 1657
```

```
## 1658
```

```
## 1659
```

```
## 1660
```

```
## 1661
```

```
## 1662
```

```
## 1663
```

```
## 1664
```

```
## 1665
```

```
## 1666
```

```
## 1667
```

```
## 1668
```

```
## 1669
```

```
## 1670
```

```
## 1671
```

```
## 1672
```

```
## 1673
```

```
## 1674
```

```
## 1675
```

```
## 1676
```

```
## 1677
```

```
## 1678
```

```
## 1679
```

```
## 1680
```

```
## 1681
```

```
## 1682
```

```
## 1683
```

```
## 1684
```

```
## 1685
```

```
## 1686
```

```
## 1687
```

```
## 1688
```

```
## 1689
```

```
## 1690
```

```
## 1691
```

```
## 1692
```

```
## 1693
```

```
## 1694
```

```
## 1695
```

```
## 1696
```

```
## 1697
```

```
## 1698
```

```
## 1699
```

```
## 1700
```

```
## 1701
```

```
## 1702
```

```
## 1703
```

```
## 1704
```

```
## 1705
```

```
## 1706
```

```
## 1707
```

```
## 1708
```

```
## 1709
```

```
## 1710
```

```
## 1711
```

```
## 1712
```

```
## 1713
```

```
## 1714
```

```
## 1715
```

```
## 1716
```

```
## 1717
```

```
## 1718
```

```
## 1719
```

```
## 1720
```

```
## 1721
```

```
## 1722
```

```
## 1723
```

```
## 1724
```

```
## 1725
```

```
## 1726
```

```
## 1727
```

```
## 1728
```

```
## 1729
```

```
## 1730
```

```
## 1731
```

```
## 1732
```

```
## 1733
```

```
## 1734
```

```
## 1735
```

```
## 1736
```

```
## 1737
```

```
## 1738
```

```
## 1739
```

```
## 1740
```

```
## 1741
```

```
## 1742
```

```
## 1743
```

```
## 1744
```

```
## 1745
```

```
## 1746
```

```
## 1747
```

```
## 1748
```

```
## 1749
```

```
## 1750
```

```
## 1751
```

```
## 1752
```

```
## 1753
```

```
## 1754
```

```
## 1755
```

```
## 1756
```

```
## 1757
```

```
## 1758
```

```
## 1759
```

```
## 1760
```

```
## 1761
```

```
## 1762
```

```
## 1763
```

```
## 1764
```

```
## 1765
```

```
## 1766
```

```
## 1767
```

```
## 1768
```

```
## 1769
```

```
## 1770
```

```
## 1771
```

```
## 1772
```

```
## 1773
```

```
## 1774
```

```
## 1775
```

```
## 1776
```

```
## 1777
```

```
## 1778
```

```
## 1779
```

```
## 1780
```

```
## 1781
```

```
## 1782
```

```
## 1783
```

```
## 1784
```

```
## 1785
```

```
## 1786
```

```
## 1787
```

```
## 1788
```

```
## 1789
```

```
## 1790
```

```
## 1791
```

```
## 1792
```

```
## 1793
```

```
## 1794
```

```
## 1795
```

```
## 1796
```

```
## 1797
```

```
## 1798
```

```
## 1799
```

```
## 1800
```

```
## 1801
```

```
## 1802
```

```
## 1803
```

```
## 1804
```

```
## 1805
```

```
## 1806
```

```
## 1807
```

```
## 1808
```

```
## 1809
```

```
## 1810
```

```
## 1811
```

```
## 1812
```

```
## 1813
```

```
## 1814
```

```
## 1815
```

```
## 1816
```

```
## 1817
```

```
## 1818
```

```
## 1819
```

```
## 1820
```

```
## 1821
```

```
## 1822
```

```
## 1823
```

```
## 1824
```

```
## 1825
```

```
## 1826
```

```
## 1827
```

```
## 1828
```

```
## 1829
```

```
## 1830
```

```
## 1831
```

```
## 1832
```

```
## 1833
```

```
## 1834
```

```
## 1835
```

```
## 1836
```

```
## 1837
```

```
## 1838
```

```
## 1839
```

```
## 1840
```

```
## 1841
```

```
## 1842
```

```
## 1843
```

```
## 1844
```

```
## 1845
```

```
## 1846
```

```
## 1847
```

```
## 1848
```

```
## 1849
```

```
## 1850
```

```
## 1851
```

```
## 1852
```

```
## 1853
```

```
## 1854
```

```
## 1855
```

```
## 1856
```

```
## 1857
```

```
## 1858
```

```
## 1859
```

```
## 1860
```

```
## 1861
```

```
## 1862
```

```
## 1863
```

```
## 1864
```

```
## 1865
```

```
## 1866
```

```
## 1867
```

```
## 1868
```

```
## 1869
```

```
## 1870
```

```
## 1871
```

```
## 1872
```

```
## 1873
```

```
## 1874
```

```
## 1875
```

```
## 1876
```

```
## 1877
```

```
## 1878
```

```
## 1879
```

```
## 1880
```

```
## 1881
```

```
## 1882
```

```
## 1883
```

```
## 1884
```

```
## 1885
```

```
## 1886
```

```
## 1887
```

```
## 1888
```

```
## 1889
```

```
## 1890
```

```
## 1891
```

```
## 1892
```

```
## 1893
```

```
## 1894
```

```
## 1895
```

```
## 1896
```

```
## 1897
```

```
## 1898
```

```
## 1899
```

```
## 1900
```

```
## 1901
```

```
## 1902
```

```
## 1903
```

```
## 1904
```

```
## 1905
```

```
## 1906
```

```
## 1907
```

```
## 1908
```

```
## 1909
```

```
## 1910
```

```
## 1911
```

```
## 1912
```

```
## 1913
```

```
## 1914
```

```
## 1915
```

```
## 1916
```

```
## 1917
```

```
## 1918
```

```
## 1919
```

```
## 1920
```

```
## 1921
```

```
## 1922
```

```
## 1923
```

```
## 1924
```

```
## 1925
```

```
## 1926
```

```
## 1927
```

```
## 1928
```

```
## 1929
```

```
## 1930
```

```
## 1931
```

```
## 1932
```

```
## 1933
```

```
## 1934
```

```
## 1935
```

```
## 1936
```

```
## 1937
```

```
## 1938
```

```
## 1939
```

```
## 1940
```

```
## 1941
```

```
## 1942
```

```
## 1943
```

```
## 1944
```

```
## 1945
```

```
## 1946
```

```
## 1947
```

```
## 1948
```

```
## 1949
```

```
## 1950
```

```
## 1951
```

```
## 1952
```

```
## 1953
```

```
## 1954
```

```
## 1955
```

```
## 1956
```

```
## 1957
```

```
## 1958
```

```
## 1959
```

```
## 1960
```

```
## 1961
```

```
## 1962
```

```
## 1963
```

```
## 1964
```

```
## 1965
```

```
## 1966
```

```
## 1967
```

```
## 1968
```

```
## 1969
```

```
## 1970
```

```
## 1971
```

```
## 1972
```

```
## 1973
```

```
## 1974
```

```
## 1975
```

```
## 1976
```

```
## 1977
```

```
## 1978
```

```
## 1979
```

```
## 1980
```

```
## 1981
```

```
## 1982
```

```
## 1983
```

```
## 1984
```

```
## 1985
```

```
## 1986
```

```
## 1987
```

```
## 1988
```

```
## 1989
```

```
## 1990
```

```
## 1991
```

```
## 1992
```

```
## 1993
```

```
## 1994
```

```
## 1995
```

```
## 1996
```

```
## 1997
```

```
## 1998
```

```
## 1999
```

```
## 2000
```

```
## 2001
```

```
## 2002
```

```
## 2003
```

```
## 2004
```

```
## 2005
```

```
## 2006
```

```
## 2007
```

```
## 2008
```

```
## 2009
```

```
## 2010
```

```
## 2011
```

```
## 2012
```

```
## 2013
```

```
## 2014
```

```
## 2015
```

```
## 2016
```

```
## 2017
```

```
## 2018
```

```
## 2019
```

```
## 2020
```

```
## 2021
```

```
## 2022
```

```
## 2023
```

```
## 2024
```

```
## 2025
```

```
## 2026
```

```
## 2027
```

```
## 2028
```

```
## 2029
```

```
## 2030
```

```
## 2031
```

```
## 2032
```

```
## 2033
```

```
## 2034
```

```
## 2035
```

```
## 2036
```

```
## 2037
```

```
## 2038
```

```
## 2039
```

```
## 2040
```

```
## 2041
```

```
## 2042
```

```
## 2043
```

```
## 2044
```

```
## 2045
```

```
## 2046
```

```
## 2047
```

```
## 2048
```

```
## 2049
```

```
## 2050
```

```
## 2051
```

```
## 2052
```

```
## 2053
```

```
## 2054
```

```
## 2055
```

```
## 2056
```

```
## 2057
```

```
## 2058
```

```
## 2059
```

```
## 2060
```

```
## 2061
```

```
## 2062
```

```
## 2063
```

```
## 2064
```

```
## 2065
```

```
## 2066
```

```
## 2067
```

```
## 2068
```

```
## 2069
```

```
## 2070
```

```
## 2071
```

```
## 2072
```

```
## 2073
```

```
## 2074
```

```
## 2075
```

```
## 2076
```

```
## 2077
```

```
## 2078
```

```
## 2079
```

```
## 2080
```

```
## 2081
```

```
## 2082
```

```
## 2083
```

```
## 2084
```

```
## 2085
```

```
## 2086
```

```
## 2087
```

```
## 2088
```

```
## 2089
```

```
## 2090
```

```
## 2091
```

```
## 2092
```

```
## 2093
```

```
## 2094
```

```
## 2095
```

```
## 2096
```

```
## 2097
```

```
## 2098
```

```
## 2099
```

```
## 2100
```

```
## 2101
```

```
## 2102
```

```
## 2103
```

```
## 2104
```

```
## 2105
```

```
## 2106
```

```
## 2107
```

```
## 2108
```

```
## 2109
```

```
## 2110
```

```
## 2111
```

```
## 2112
```

```
## 2113
```

```
## 2114
```

```
## 2115
```

```
## 2116
```

```
## 2117
```

```
## 2118
```

```
## 2119
```

```
## 2120
```

```
## 2121
```

```
## 2122
```

```
## 2123
```

```
## 2124
```

```
## 2125
```

```
## 2126
```

```
## 2127
```

```
## 2128
```

```
## 2129
```

```
## 2130
```

```
## 2131
```

```
## 2132
```

```
## 2133
```

```
## 2134
```

```
## 2135
```

```
## 2136
```

```
## 2137
```

```
## 2138
```

```
## 2139
```

```
## 2140
```

```
## 2141
```

```
## 2142
```

```
## 2143
```

```
## 2144
```

```
## 2145
```

```
## 2146
```

```
## 2147
```

```
## 2148
```

```
## 2149
```

```
## 2150
```

```
## 2151
```

```
## 2152
```

```
## 2153
```

```
## 2154
```

```
## 2155
```

```
## 2156
```

```
## 2157
```

```
## 2158
```

```
## 2159
```

```
## 2160
```

```
## 2161
```

```
## 2162
```

```
## 2163
```

```
## 2164
```

```
## 2165
```

```
## 2166
```

```
## 2167
```

```
## 2168
```

```
## 2169
```

```
## 2170
```

```
## 2171
```

```
## 2172
```

```
## 2173
```

```
## 2174
```

```
## 2175
```

```
## 2176
```

```
## 2177
```

```
## 2178
```

```
## 2179
```

```
## 2180
```

```
## 2181
```

```
## 2182
```

```
## 2183
```

```
## 2184
```

```
## 2185
```

```
## 2186
```

```
## 2187
```

```
## 2188
```

```
## 2189
```

```
## 2190
```

```
## 2191
```

```
## 2192
```

```
## 2193
```

```
## 2194
```

```
## 2195
```

```
## 2196
```

```
## 2197
```

```
## 2198
```

```
## 2199
```

```
## 2200
```

```
## 2201
```

```
## 2202
```

```
## 2203
```

```
## 2204
```

```
## 2205
```

```
## 2206
```

```
## 2207
```

```
## 2208
```

```
## 2209
```

```
## 2210
```

```
## 2211
```

```
## 2212
```

```
## 2213
```

```
## 2214
```

```
## 2215
```

```
## 2216
```

```
## 2217
```

```
## 2218
```

```
## 2219
```

```
## 2220
```

```
## 2221
```

```
## 2222
```

```
## 2223
```

```
## 2224
```

```
## 2225
```

```
## 2226
```

```
## 2227
```

```
## 2228
```

```
## 2229
```

```
## 2230
```

```
## 2231
```

```
## 2232
```

```
## 2233
```

```
## 2234
```

```
## 2235
```

```
## 2236
```

```
## 2237
```

```
## 2238
```

```
## 2239
```

```
## 2240
```

```
## 2241
```

```
## 2242
```

```
## 2243
```

```
## 2244
```

```
## 2245
```

```
## 2246
```

```
## 2247
```

```
## 2248
```

```
## 2249
```

```
## 2250
```

```
## 2251
```

```
## 2252
```

```
## 2253
```

```
## 2254
```

```
## 2255
```

```
## 2256
```

```
## 2257
```

```
## 2258
```

```
## 2259
```

```
## 2260
```

```
## 2261
```

```
## 2262
```

```
## 2263
```

```
## 2264
```

```
## 2265
```

```
## 2266
```

```
## 2267
```

```
## 2268
```

```
## 2269
```

```
## 2270
```

```
## 2271
```

```
## 2272
```

```
## 2273
```

```
## 2274
```

```
## 2275
```

```
## 2276
```

```
## 2277
```

```
## 2278
```

```
## 2279
```

```
## 2280
```

```
## 2281
```

```
## 2282
```

```
## 2283
```

```
## 2284
```

```
## 2285
```

```
## 2286
```

```
## 2287
```

```
## 2288
```

```
## 2289
```

```
## 2290
```

```
## 2291
```

```
## 2292
```

```
## 2293
```

```
## 2294
```

```
## 2295
```

```
## 2296
```

```
## 2297
```

```
## 2298
```

```
## 2299
```

```
## 2300
```

```
## 2301
```

```
## 2302
```

```
## 2303
```

```
## 2304
```

```
## 2305
```

```
## 2306
```

```
## 2307
```

```
## 2308
```

```
## 2309
```

```
## 2310
```

```
## 2311
```

```
## 2312
```

```
## 2313
```

```
## 2314
```

```
## 2315
```

```
## 2316
```

```
## 2317
```

```
## 2318
```

```
## 2319
```

```
## 2320
```

```
## 2321
```

```
## 2322
```

```
## 2323
```

```
## 2324
```

```
## 2325
```

```
## 2326
```

```
## 2327
```

```
## 2328
```

```
## 2329
```

```
## 2330
```

```
## 2331
```

```
## 2332
```

```
## 2333
```

```
## 2334
```

```
## 2335
```

```
## 2336
```

```
## 2337
```

```
## 2338
```

```
## 2339
```

```
## 2340
```

```
## 2341
```

```
## 2342
```

```
## 2343
```

```
## 2344
```

```
## 2345
```

```
## 2346
```

```
## 2347
```

```
## 2348
```

```
## 2349
```

```
## 2350
```

```
## 2351
```

```
## 2352
```

```
## 2353
```

```
## 2354
```

```
## 2355
```

```
## 2356
```

```
## 2357
```

```
## 2358
```

```
## 2359
```

```
## 2360
```

```
## 2361
```

```
## 2362
```

```
## 2363
```

```
## 2364
```

```
## 2365
```

```
## 2366
```

```
## 2367
```

```
## 2368
```

```
## 2369
```

```
## 2370
```

```
## 2371
```

```
## 2372
```

```
## 2373
```

```
## 2374
```

```
## 2375
```

```
## 2376
```

```
## 2377
```

```
## 2378
```

```
## 2379
```

```
## 2380
```

```
## 2381
```

```
## 2382
```

```
## 2383
```

```
## 2384
```

```
## 2385
```

```
## 2386
```

```
## 2387
```

```
## 2388
```

```
## 2389
```

```
## 2390
```

```
## 2391
```

```
## 2392
```

```
## 2393
```

```
## 2394
```

```
## 2395
```

```
## 2396
```

```
## 2397
```

```
## 2398
```

```
## 2399
```

```
## 2400
```

```
## 2401
```

```
## 2402
```

```
## 2403
```

```
## 2404
```

```
## 2405
```

```
## 2406
```

```
## 2407
```

```
## 2408
```

```
## 2409
```

```
## 2410
```

```
## 2411
```

```
## 2412
```

```
## 2413
```

```
## 2414
```

```
## 2415
```

```
## 2416
```

```
## 2417
```

```
## 2418
```

```
## 2419
```

```
## 2420
```

```
## 2421
```

```
## 2422
```

```
## 2423
```

```
## 2424
```

```
## 2425
```

```
## 2426
```

```
## 2427
```

```
## 2428
```

```
## 2429
```

```
## 2430
```

```
## 2431
```

```
## 2432
```

```
## 2433
```

```
## 2434
```

```
## 2435
```

```
## 2436
```

```
## 2437
```

```
## 2438
```

```
## 2439
```

```
## 2440
```

```
## 2441
```

```
## 2442
```

```
## 2443
```

```
## 2444
```

```
## 2445
```

```
## 2446
```

```
## 2447
```

```
## 2448
```

```
## 2449
```

```
## 2450
```

```
## 2451
```

```
## 2452
```

```
## 2453
```

```
## 2454
```

```
## 2455
```

```
## 2456
```

```
## 2457
```

```
## 2458
```

```
## 2459
```

```
## 2460
```

```
## 2461
```

```
## 2462
```

```
## 2463
```

```
## 2464
```

```
## 2465
```

```
## 2466
```

```
## 2467
```

```
## 2468
```

```
## 2469
```

```
## 2470
```

```
## 2471
```

```
## 2472
```

```
## 2473
```

```
## 2474
```

```
## 2475
```

```
## 2476
```

```
## 2477
```

```
## 2478
```

```
## 2479
```

```
## 2480
```

```
## 2481
```

```
## 2482
```

```
## 2483
```

```
## 2484
```

```
## 2485
```

```
## 2486
```

```
## 2487
```

```
## 2488
```

```
## 2489
```

```
## 2490
```

```
## 2491
```

```
## 2492
```

```
## 2493
```

```
## 2494
```

```
## 2495
```

```
## 2496
```

```
## 2497
```

```
## 2498
```

```
## 2499
```

```
## 2500
```

```
## 2501
```

```
## 2502
```

```
## 2503
```

```
## 2504
```

```
## 2505
```

```
## 2506
```

```
## 2507
```

```
## 2508
```

```
## 2509
```

```
## 2510
```

```
## 2511
```

```
## 2512
```

```
## 2513
```

```
## 2514
```

```
## 2515
```

```
## 2516
```

```
## 2517
```

```
## 2518
```

```
## 2519
```

```
## 2520
```

```
## 2521
```

```
## 2522
```

```
## 2523
```

```
## 2524
```

```
## 2525
```

```
## 2526
```

```
## 2527
```

```
## 2528
```

```
## 2529
```

```
## 2530
```

```
## 2531
```

```
## 2532
```

```
## 2533
```

```
## 2534
```

```
## 2535
```

```
## 2536
```

```
## 2537
```

```
## 2538
```

```
## 2539
```

```
## 2540
```

```
## 2541
```

```
## 2542
```

```
## 2543
```

```
## 2544
```

```
## 2545
```

```
## 2546
```

```
## 2547
```

```
## 2548
```

```
## 2549
```

```
## 2550
```

```
## 2551
```

```
## 2552
```

```
## 2553
```

```
## 2554
```

```
## 2555
```

```
## 2556
```

```
## 2557
```

```
## 2558
```

```
## 2559
```

```
## 2560
```

```
## 2561
```

```
## 2562
```

```
## 2563
```

```
## 2564
```

```
## 2565
```

```
## 2566
```

```
## 2567
```

```
## 2568
```

```
## 2569
```

```
## 2570
```

```
## 2571
```

```
## 2572
```

```
## 2573
```

```
## 2574
```

```
## 2575
```

```
## 2576
```

```
## 2577
```

```
## 2578
```

```
## 2579
```

```
## 2580
```

```
## 2581
```

```
## 2582
```

```
## 2583
```

```
## 2584
```

```
## 2585
```

```
## 2586
```

```
## 2587
```

```
## 2588
```

```
## 2589
```

```
## 2590
```

```
## 2591
```

```
## 2592
```

```
## 2593
```

```
## 2594
```

```
## 2595
```

```
## 2596
```

```
## 2597
```

```
## 2598
```

```
## 2599
```

```
## 2600
```

```
## 2601
```

```
## 2602
```

```
## 2603
```

```
## 2604
```

```
## 2605
```

```
## 2606
```

```
## 2607
```

```
## 2608
```

```
## 2609
```

```
## 2610
```

```
## 2611
```

```
## 2612
```

```
## 2613
```

```
## 2614
```

```
## 2615
```

```
## 2616
```

```
## 2617
```

```
## 2618
```

```
## 2619
```

```
## 2620
```

```
## 2621
```

```
## 2622
```

```
## 2623
```

```
## 2624
```

```
## 2625
```

```
## 2626
```

```
## 2627
```

```
## 2628
```

```
## 2629
```

```
## 2630
```

```
## 2631
```

```
## 2632
```

```
## 2633
```

```
## 2634
```

```
## 2635
```

```
## 2636
```

```
## 2637
```

```
## 2638
```

```
## 2639
```

```
## 2640
```

```
## 2641
```

```
## 2642
```

```
## 2643
```

```
## 2644
```

```
## 2645
```

```
## 2646
```

```
## 2647
```

```
## 2648
```

```
## 2649
```

```
## 2650
```

```
## 2651
```

```
## 2652
```

```
## 2653
```

```
## 2654
```

```
## 2655
```

```
## 2656
```

```
## 2657
```

```
## 2658
```

```
## 2659
```

```
## 2660
```

```
## 2661
```

```
## 2662
```

```
## 2663
```

```
## 2664
```

```
## 2665
```

```
## 2666
```

```
## 2667
```

```
## 2668
```

```
## 2669
```

```
## 2670
```

```
## 2671
```

```
## 2672
```

```
## 2673
```

```
## 2674
```

```
## 2675
```

```
## 2676
```

```
## 2677
```

```
## 2678
```

```
## 2679
```

```
## 2680
```

```
## 2681
```

```
## 2682
```

```
## 2683
```

```
## 2684
```

```
## 2685
```

```
## 2686
```

```
## 2687
```

```
## 2688
```

```
## 2689
```

```
## 2690
```

```
## 2691
```

```
## 2692
```

```
## 2693
```

```
## 2694
```

```
## 2695
```

```
## 2696
```

```
## 2697
```

```
## 2698
```

```
## 2699
```

```
## 2700
```

```
## 2701
```

```
## 2702
```

```
## 2703
```

```
## 2704
```

```
## 2705
```

```
## 2706
```

```
## 2707
```

```
## 2708
```

```
## 2709
```

```
## 2710
```

```
## 2711
```

```
## 2712
```

```
## 2713
```

```
## 2714
```

```
## 2715
```

```
## 2716
```

```
## 2717
```

```
## 2718
```

```
## 2719
```

```
## 2720
```

```
## 2721
```

```
## 2722
```

```
## 2723
```

```
## 2724
```

```
## 2725
```

```
## 2726
```

```
## 2727
```

```
## 2728
```

```
## 2729
```

```
## 2730
```

```
## 2731
```

```
## 2732
```

```
## 2733
```

```
## 2734
```

```
## 2735
```

```
## 2736
```

```
## 2737
```

```
## 2738
```

```
## 2739
```

```
## 2740
```

```
## 2741
```

```
## 2742
```

```
## 2743
```

```
## 2744
```

```
## 2745
```

```
## 2746
```

```
## 2747
```

```
## 2748
```

```
## 2749
```

```
## 2750
```

```
## 2751
```

```
## 2752
```

```
## 2753
```

```
## 2754
```

```
## 2755
```

```
## 2756
```

```
## 2757
```

```
## 2758
```

```
## 2759
```

```
## 2760
```

```
## 2761
```

```
## 2762
```

```
## 2763
```

```
## 2764
```

```
## 2765
```

```
## 2766
```

```
## 2767
```

```
## 2768
```

```
## 2769
```

```
## 2770
```

```
## 2771
```

```
## 2772
```

```
## 2773
```

```
## 2774
```

```
## 2775
```

```
## 2776
```

```
## 2777
```

```
## 2778
```

```
## 2779
```

```
## 2780
```

```
## 2781
```

```
## 2782
```

```
## 2783
```

```
## 2784
```

```
## 2785
```

```
## 2786
```

```
## 2787
```

```
## 2788
```

```
## 2789
```

```
## 2790
```

```
## 2791
```

```
## 2792
```

```
## 2793
```

```
## 2794
```

```
## 2795
```

```
## 2796
```

```
## 2797
```

```
## 2798
```

```
## 2799
```

```
## 2800
```

```
## 2801
```

```
## 2802
```

```
## 2803
```

```
## 2804
```

```
## 2805
```

```
## 2806
```

```
## 2807
```

```
## 2808
```

```
## 2809
```

```
## 2810
```

```
## 2811
```

```
## 2812
```

```
## 2813
```

```
## 2814
```

```
## 2815
```

```
## 2816
```

```
## 2817
```

```
## 2818
```

```
## 2819
```

```
## 2820
```

```
## 2821
```

```
## 2822
```

```
## 2823
```

```
## 2824
```

```
## 2825
```

```
## 2826
```

```
## 2827
```

```
## 2828
```

```
## 2829
```

```
## 2830
```

```
## 2831
```

```
## 2832
```

```
## 2833
```

```
## 2834
```

```
## 2835
```

```
## 2836
```

```
## 2837
```

```
## 2838
```

```
## 2839
```

```
## 2840
```

```
## 2841
```

```
## 2842
```

```
## 2843
```

```
## 2844
```

```
## 2845
```

```
## 2846
```

```
## 2847
```

```
## 2848
```

```
## 2849
```

```
## 2850
```

```
## 2851
```

```
## 2852
```

```
## 2853
```

```
## 2854
```

```
## 2855
```

```
## 2856
```

```
## 2857
```

```
## 2858
```

```
## 2859
```

```
## 2860
```

```
## 2861
```

```
## 2862
```

```
## 2863
```

```
## 2864
```

```
## 2865
```

```
## 2866
```

```
## 2867
```

```
## 2868
```

```
## 2869
```

```
## 2870
```

```
## 2871
```

```
## 2872
```

```
## 2873
```

```
## 2874
```

```
## 2875
```

```
## 2876
```

```
## 2877
```

```
## 2878
```

```
## 2879
```

```
## 2880
```

```
## 2881
```

```
## 2882
```

```
## 2883
```

```
## 2884
```

```
## 2885
```

```
## 2886
```

```
## 2887
```

```
## 2888
```

```
## 2889
```

```
## 2890
```

```
## 2891
```

```
## 2892
```

```
## 2893
```

```
## 2894
```

```
## 2895
```

```
## 2896
```

```
## 2897
```

```
## 2898
```

```
## 2899
```

```
## 2900
```

```
## 2901
```

```
## 2902
```

```
## 2903
```

```
## 2904
```

```
## 2905
```

```
## 2906
```

```
## 2907
```

```
## 2908
```

```
## 2909
```

```
## 2910
```

```
## 2911
```

```
## 2912
```

```
## 2913
```

```
## 2914
```

```
## 2915
```

```
## 2916
```

```
## 2917
```

```
## 2918
```

```
## 2919
```

```
## 2920
```

```
## 2921
```

```
## 2922
```

```
## 2923
```

```
## 2924
```

```
## 2925
```

```
## 2926
```

```
## 2927
```

```
## 2928
```

```
## 2929
```

```
## 2930
```

```
## 2931
```

```
## 2932
```

```
## 2933
```

```
## 2934
```

```
## 2935
```

```
## 2936
```

```
## 2937
```

```
## 2938
```

```
## 2939
```

```
## 2940
```

```
## 2941
```

```
## 2942
```

```
## 2943
```

```
## 2944
```

```
## 2945
```

```
## 2946
```

```
## 2947
```

```
## 2948
```

```
## 2949
```

```
## 2950
```

```
## 2951
```

```
## 2952
```

```
## 2953
```

```
## 2954
```

```
## 2955
```

```
## 2956
```

```
## 2957
```

```
## 2958
```

```
## 2959
```

```
## 2960
```

```
## 2961
```

```
## 2962
```

```
## 2963
```

```
## 2964
```

```
## 2965
```

```
## 2966
```

```
## 2967
```

```
## 2968
```

```
## 2969
```

```
## 2970
```

```
## 2971
```

```
## 2972
```

```
## 2973
```

```
## 2974
```

```
## 2975
```

```
## 2976
```

```
## 2977
```

```
## 2978
```

```
## 2979
```

```
## 2980
```

```
## 2981
```

```
## 2982
```

```
## 2983
```

```
## 2984
```

```
## 2985
```

```
## 2986
```

```
## 2987
```

```
## 2988
```

```
## 2989
```

```
## 2990
```

```
## 2991
```

```
## 2992
```

```
## 2993
```

```
## 2994
```

```
## 2995
```

```
## 2996
```

```
## 2997
```

```
## 2998
```

```
## 2999
```

```
## 3000
```

```r
epidish_res333 <- epidish_res333$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (CP)')

#combine algorithms
epidish_results <- bind_rows(epidish_res111, epidish_res222, epidish_res333, 
                             houseman_silico) %>%
  dplyr::rename(sample_name = Sample_Names) %>%
  left_join(silico_pDat %>%
              pivot_longer(cols = a:f,
                           names_to ='component',
                           values_to = 'actual') %>%
              mutate(component = case_when(
                    component == 'a' ~ 'Trophoblasts',
                    component == 'b' ~ 'Stromal',
                    component == 'c' ~ 'Endothelial',
                    component == 'd' ~ 'Hofbauer',
                    component == 'e' ~ 'nRBC',
                    component == 'f' ~ 'Syncytiotrophoblast'))) %>%
  dplyr::rename(estimated = percent)
```

```
## Joining, by = c("sample_name", "component")
```

```r
# stats
stats <- epidish_results %>%
  group_by(algorithm, component) %>%
  metrics(estimated, actual) %>%
  mutate(label = number(.estimate, accuracy = 0.01)) %>%
  pivot_wider(id_cols = c(algorithm, component),
              names_from = '.metric',
              values_from = c('label', '.estimate')) %>%
  rowwise() %>%
  mutate(label = paste0('RMSE=', label_rmse,
                        '\nR2=', label_rsq,
                        '\nMAE=', label_mae))

# scatter plot
epidish_results %>%
  filter(trimester == 'Third') %>%
  ggplot(aes(x = actual, y = estimated, color = component)) +
  geom_point(alpha = 0.2) +
  geom_text(data = stats, x= 0.05, y= 1, aes(label = label), color = 'black',
             hjust = 0, vjust = 1) +
  facet_grid(cols = vars(algorithm), rows = vars(component)) +
  scale_color_manual(values = color_code_tissue[unique(epidish_results$component)], 
                    na.value = '#636363',
                    guide = 'none') +
  scale_x_continuous(expand = c(0,0), labels = scales::percent) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  theme(panel.grid = element_blank())
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

```r
stats %>%
  group_by(algorithm) %>%
  summarize(rmse = mean(.estimate_rmse),
            rsq = mean(.estimate_rsq),
            mae = mean(.estimate_mae))
```

```
## `summarise()` ungrouping output (override with `.groups` argument)
```

```
## # A tibble: 4 x 4
##   algorithm       rmse   rsq    mae
##   <chr>          <dbl> <dbl>  <dbl>
## 1 epidish (CBS) 0.0447 0.960 0.0235
## 2 epidish (CP)  0.0459 0.957 0.0246
## 3 epidish (RPC) 0.0442 0.960 0.0232
## 4 Houseman (CP) 0.0459 0.957 0.0245
```

```r
# mean difference
epidish_results <- epidish_results %>%
  mutate(deviation = estimated-actual)

epidish_results %>%
  filter(trimester == 'Third') %>%
  ggplot(aes(x = deviation, color = component)) +
  geom_density() +
  geom_vline(data = epidish_results %>%
               filter(trimester == 'Third') %>%
               group_by(algorithm, component) %>%
               summarize(mean_diff = mean(deviation)),
             aes(xintercept = mean_diff),
             linetype = 'dashed') +
  geom_text(data = epidish_results %>%
               filter(trimester == 'Third') %>%
               group_by(algorithm, component) %>%
               summarize(mean_diff = percent(mean(deviation), accuracy = 0.01)),
             aes(label = mean_diff),
            x = -0.05, y = 39, hjust = 0, vjust =1,color = 'black') +
  geom_vline(xintercept = 0) +
  facet_grid(cols = vars(algorithm), rows = vars(component)) +
  scale_color_manual(values = color_code_tissue[unique(epidish_results$component)], 
                    na.value = '#636363',
                    guide = 'none') +
  scale_x_continuous(limits = c(-0.05, 0.05), labels = percent) +
  scale_y_continuous(expand = c(0,0))+
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = '(estimated percentage) - (actual)')
```

```
## `summarise()` regrouping output by 'algorithm' (override with `.groups` argument)
```

```
## `summarise()` regrouping output by 'algorithm' (override with `.groups` argument)
```

```
## Warning: Removed 4383 rows containing non-finite values (stat_density).
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-15-2.png)<!-- -->

## apply to villi


```r
# apply deconvolution and process results
houseman_estimates_third <- minfi:::projectCellType(
  getBeta(mset_test)[rownames(coefs_combined_third)
                                ,pData(mset_test)$Trimester == 'Third'],
  coefs_combined_third,
  lessThanOne = TRUE)

houseman_estimates_third <- houseman_estimates_third %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'Houseman')

houseman_estimates_first <- minfi:::projectCellType(
  getBeta(mset_test)[rownames(coefs_combined_third)
                                ,pData(mset_test)$Trimester == 'First'],
  coefs_combined_first,
  lessThanOne = TRUE)

houseman_estimates_first <- houseman_estimates_first %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'Houseman')

epidish_res111 <- epidish(
  beta.m = getBeta(mset_test)[rownames(coefs_combined_third)
                                ,pData(mset_test)$Trimester == 'Third'],
  ref.m = coefs_combined_third,
  method = 'RPC')

epidish_res1111 <- epidish(
  beta.m = getBeta(mset_test)[rownames(coefs_combined_first)
                                ,pData(mset_test)$Trimester == 'First'],
  ref.m = coefs_combined_first,
  method = 'RPC')

epidish_res111 <- epidish_res111$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (RPC)')

epidish_res1111 <- epidish_res1111$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (RPC)')

epidish_rpc <- bind_rows(houseman_estimates_first, houseman_estimates_third,
                         epidish_res111, epidish_res1111) %>%
  dplyr::rename(Sample_Name = Sample_Names) %>%
  left_join(pDat_filt %>% select(Sample_Name, Sex, Trimester, Tissue, GA))  %>%
  distinct() 
```

```
## Joining, by = "Sample_Name"
```

### plots 


```r
# analyze
epidish_rpc %>%
  filter(Tissue == 'Villi',
         algorithm == 'epidish (RPC)') %>%
  pivot_wider(id_cols = -c(component, percent),
              names_from = 'component',
              values_from = 'percent') %>%
  mutate(Sample_Name = fct_reorder2(Sample_Name, Trophoblasts, Stromal)) %>%
  pivot_longer(cols = Trophoblasts:Syncytiotrophoblast,
               names_to = 'component',
               values_to = 'percent') %>%
  ggplot(aes(x = Sample_Name, y = percent, fill = component)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_code_tissue[unique(epidish_rpc$component)], 
                    na.value = '#636363') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border= element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust =0)) +
  scale_y_continuous(expand = c(0,0), labels = percent) +
  labs(x = '', y = '', fill = '', title = 'epidish (RPC)') +
  facet_grid(cols = vars(Trimester), scales = 'free', space = 'free') 
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
# houseman
epidish_rpc %>%
  filter(Tissue == 'Villi',
         algorithm == 'Houseman') %>%
  pivot_wider(id_cols = -c(component, percent),
              names_from = 'component',
              values_from = 'percent') %>%
  mutate(Sample_Name = fct_reorder2(Sample_Name, Trophoblasts, Stromal)) %>%
  pivot_longer(cols = Trophoblasts:Syncytiotrophoblast,
               names_to = 'component',
               values_to = 'percent') %>%
  ggplot(aes(x = Sample_Name, y = percent, fill = component)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_code_tissue[unique(epidish_rpc$component)], 
                    na.value = '#636363') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border= element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust =0)) +
  scale_y_continuous(expand = c(0,0), labels = percent) +
  labs(x = '', y = '', fill = '',  title = 'Houseman') +
  facet_grid(cols = vars(Trimester), scales = 'free', space = 'free') 
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-17-2.png)<!-- -->

```r
epidish_rpc %>%
  filter(Tissue == 'Villi') %>%
  ggplot(aes(x = component, y = percent, color = component)) +
  geom_boxplot() +
  geom_jitter() +
  scale_color_manual(values = color_code_tissue[unique(epidish_rpc$component)], 
                    na.value = '#636363', guide = 'none') +
  scale_y_continuous(labels = percent) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  facet_grid(cols = vars(Trimester)) +
  labs(x = '', color = '', y = '') 
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-17-3.png)<!-- -->

```r
epidish_rpc %>%
  filter(Tissue == 'Villi',
         algorithm == 'epidish (RPC)') %>%
  group_by(Trimester, component) %>%
  summarize(mean = mean(percent),
            sd = sd(percent)) %>% 
  pivot_wider(id_cols = component, names_from = c('Trimester'), values_from = c('mean', 'sd')) %>%
  mutate(diff = mean_Third-mean_First) %>%
  arrange(desc(component)) %>%
  select(component, contains('First'), contains('Third'), diff)
```

```
## `summarise()` regrouping output by 'Trimester' (override with `.groups` argument)
```

```
## # A tibble: 6 x 6
##   component           mean_First sd_First mean_Third sd_Third    diff
##   <chr>                    <dbl>    <dbl>      <dbl>    <dbl>   <dbl>
## 1 Trophoblasts            0.156    0.117      0.201    0.0553  0.0444
## 2 Syncytiotrophoblast     0.348    0.0868     0.577    0.0832  0.229 
## 3 Stromal                 0.433    0.132      0.124    0.0344 -0.309 
## 4 nRBC                    0        0          0        0       0     
## 5 Hofbauer                0.0339   0.0206     0.0234   0.0141 -0.0105
## 6 Endothelial             0.0285   0.0208     0.0745   0.0145  0.0460
```

### Stats n stuff


```r
# between trimester
epidish_rpc %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  nest(data = -component) %>%
  mutate(lm = map(data, ~lm(percent ~ Trimester, .)),
         glanced = map(lm, glance)) %>%
  select(-data, -lm)  %>%
  unnest(cols = c(glanced)) %>%
  select(component, r.squared, p.value) %>%
  mutate(adj_p = p.adjust(p.value, method = 'bonferroni'),
         p_value = pvalue(adj_p))
```

```
## # A tibble: 6 x 5
##   component           r.squared   p.value     adj_p p_value
##   <chr>                   <dbl>     <dbl>     <dbl> <chr>  
## 1 Trophoblasts           0.0685   1.96e-1   9.82e-1 0.982  
## 2 Stromal                0.794    1.02e-9   5.08e-9 <0.001 
## 3 Hofbauer               0.0848   1.49e-1   7.45e-1 0.745  
## 4 Endothelial            0.629    1.34e-6   6.71e-6 <0.001 
## 5 nRBC                 NaN      NaN       NaN       <NA>   
## 6 Syncytiotrophoblast    0.613    2.27e-6   1.13e-5 <0.001
```

```r
# between sex
stats_by_sex <- epidish_rpc %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  filter(Trimester == 'Third', component != 'nRBC') %>%
  nest(data = -component) %>%
  mutate(lm = map(data, ~lm(percent ~ Sex, .)),
         glanced = map(lm, glance)) %>%
  select(-data, -lm)  %>%
  unnest(cols = c(glanced)) %>%
  select(component, r.squared, p.value) %>%
  mutate(adj_p = p.adjust(p.value, method = 'bonferroni'),
         p_value = pvalue(adj_p),
         testing_variable = 'Sex');stats_by_sex
```

```
## # A tibble: 5 x 6
##   component           r.squared p.value adj_p p_value testing_variable
##   <chr>                   <dbl>   <dbl> <dbl> <chr>   <chr>           
## 1 Trophoblasts           0.0103  0.679  1     >0.999  Sex             
## 2 Stromal                0.0553  0.333  1     >0.999  Sex             
## 3 Hofbauer               0.199   0.0555 0.277 0.277   Sex             
## 4 Endothelial            0.0549  0.334  1     >0.999  Sex             
## 5 Syncytiotrophoblast    0.0789  0.244  1     >0.999  Sex
```

```r
epidish_rpc %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  filter(Trimester == 'Third', component !='nRBC') %>%
  ggplot(aes(x = component, y = percent, color = Sex)) +
  geom_beeswarm(dodge.width = 0.75, cex = 1.5) +
  #geom_jitter(position = position_jitterdodge(dodge.width = 0.75)) +
  scale_y_continuous(labels = function(x)percent(x,accuracy = 1)) +
  scale_color_manual(values = c('#af8dc3', '#7fbf7b')) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  labs(x = '', color = '', y = '') 
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

```r
epidish_rpc %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  filter(Trimester == 'Third') %>%
  select(Sample_Name, Sex) %>%
  distinct() %>%
  dplyr::count(Sex)
```

```
## # A tibble: 2 x 2
##   Sex       n
##   <chr> <int>
## 1 F        10
## 2 M         9
```

```r
# between ancestry
epidish_rpc  %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  filter(Trimester == 'Third') %>%
  left_join(ancestry) %>%
  ggplot(aes(x = PC1_geno, y = PC2_geno, color = Predicted_ethnicity)) +
  geom_point()
```

```
## Joining, by = "Sample_Name"
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-18-2.png)<!-- -->

```r
epidish_rpc %>%
  left_join(ancestry) %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  filter(Trimester == 'Third', Predicted_ethnicity != 'Ambiguous', component != 'nRBC') %>%
  ggplot(aes(x = component, y = percent, color = Predicted_ethnicity)) +
  geom_beeswarm(dodge.width = 0.75,size = 0.9) +
  #geom_boxplot() +
  #geom_jitter(position = position_jitterdodge(dodge.width = 0.75)) +
  scale_y_continuous(labels = percent) +
  scale_color_manual(values = c('#d8b365', '#5ab4ac')) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  labs(x = '', color = '', y = '')  
```

```
## Joining, by = "Sample_Name"
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-18-3.png)<!-- -->

```r
stats_by_ancestry <- epidish_rpc %>%
  left_join(ancestry) %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  filter(Trimester == 'Third', Predicted_ethnicity != 'Ambiguous', component != 'nRBC') %>%
  nest(data = -component) %>%
  mutate(lm = map(data, ~lm(percent ~ Predicted_ethnicity, .)),
         glanced = map(lm, glance)) %>%
  select(-data, -lm)  %>%
  unnest(cols = c(glanced)) %>%
  select(component, r.squared, p.value) %>%
  mutate(adj_p = p.adjust(p.value, method = 'bonferroni'),
         p_value = pvalue(adj_p),
         testing_variable = 'ancestry'); stats_by_ancestry
```

```
## Joining, by = "Sample_Name"
```

```
## # A tibble: 5 x 6
##   component           r.squared p.value adj_p p_value testing_variable
##   <chr>                   <dbl>   <dbl> <dbl> <chr>   <chr>           
## 1 Trophoblasts         0.0594     0.346     1 >0.999  ancestry        
## 2 Stromal              0.00807    0.732     1 >0.999  ancestry        
## 3 Hofbauer             0.000931   0.907     1 >0.999  ancestry        
## 4 Endothelial          0.0285     0.517     1 >0.999  ancestry        
## 5 Syncytiotrophoblast  0.0531     0.373     1 >0.999  ancestry
```

```r
epidish_rpc %>%
  left_join(ancestry) %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  select(Sample_Name, Predicted_ethnicity) %>%
  distinct() %>%
  dplyr::count(Predicted_ethnicity)
```

```
## Joining, by = "Sample_Name"
```

```
## # A tibble: 4 x 2
##   Predicted_ethnicity     n
##   <chr>               <int>
## 1 Ambiguous               2
## 2 Asian                   6
## 3 Caucasian              11
## 4 <NA>                    7
```

# estimate GA

We want to find any association with gestational age and cell composition.

Here I predict Ga using planet


```r
?pl_infer_age()
```

```
## starting httpd help server ... done
```

```r
ga <- tibble(Sample_Name = colnames(mset_test),
  ga_rpc = pl_infer_age(getBeta(mset_test), 'RPC'),
  ga_cpc = pl_infer_age(getBeta(mset_test), 'CPC'),
  ga_rrpc = pl_infer_age(getBeta(mset_test), 'RRPC')
)
```

```
## [1] "558 of 558 predictors present."
## [1] "546 of 546 predictors present."
## [1] "395 of 395 predictors present."
```

```r
# join ga to pDat
pDat_filt <- pDat_filt %>%
  left_join(ga)
```

```
## Joining, by = "Sample_Name"
```

```r
# join deconv estimates to pDat
pDat_filt <- pDat_filt %>%
  left_join(
    by = c('Sample_Name' = 'sample_name'),
    epidish_rpc %>%
      pivot_wider(id_cols = Sample_Name,
                  names_from = c('algorithm', 'component'),
                  values_from = 'percent') %>%
      clean_names() 
  )

pDat_ga_comp <- pDat_filt %>% 
  filter( Tissue == 'Villi') %>%
  select(Sample_Name, Trimester, Tissue, contains('ga'), contains('epidish')) %>%
  pivot_longer(cols = contains('epidish'),
               names_to = 'component',
               names_prefix = 'epidish_rpc_',
               values_to = 'proportion') %>%
  pivot_longer(cols = contains('ga'),
               names_to = 'ga',
               names_prefix = 'ga_',
               values_to = 'weeks') %>%
  mutate(ga = ifelse(ga == 'GA', 'reported', ga) %>%
           as.factor() %>%
           fct_relevel(levels = c('reported')))

stats_by_ga <- pDat_ga_comp %>%
  filter(component != 'n_rbc',
         ga %in% c('reported', 'rpc')) %>%
  nest(data = -c(Trimester, component, ga)) %>%
  mutate(lm = map(data, ~lm(proportion ~ weeks, .)),
         glanced = map(lm, glance)) %>%
  select(-data, -lm)  %>%
  unnest(cols = c(glanced)) %>%
  select(Trimester, component, ga, r.squared, p.value) %>%
  mutate(adj_p = p.adjust(p.value, method = 'fdr'),
         p_value = pvalue(adj_p, 0.01) %>%
           ifelse(. < 0.05, 'p<0.05', .),
         component = str_to_sentence(component) ); stats_by_ga
```

```
## # A tibble: 20 x 7
##    Trimester component           ga       r.squared p.value  adj_p p_value
##    <chr>     <chr>               <fct>        <dbl>   <dbl>  <dbl> <chr>  
##  1 Third     Trophoblasts        reported    0.0415 0.418   0.446  0.45   
##  2 Third     Trophoblasts        rpc         0.162  0.0873  0.146  0.15   
##  3 Third     Stromal             reported    0.105  0.189   0.252  0.25   
##  4 Third     Stromal             rpc         0.210  0.0486  0.124  0.12   
##  5 Third     Hofbauer            reported    0.326  0.0133  0.0698 0.07   
##  6 Third     Hofbauer            rpc         0.303  0.0145  0.0698 0.07   
##  7 Third     Endothelial         reported    0.180  0.0792  0.144  0.14   
##  8 Third     Endothelial         rpc         0.208  0.0496  0.124  0.12   
##  9 Third     Syncytiotrophoblast reported    0.184  0.0754  0.144  0.14   
## 10 Third     Syncytiotrophoblast rpc         0.396  0.00387 0.0698 0.07   
## 11 First     Trophoblasts        reported    0.709  0.0174  0.0698 0.07   
## 12 First     Trophoblasts        rpc         0.296  0.206   0.258  0.26   
## 13 First     Stromal             reported    0.667  0.0250  0.0834 0.08   
## 14 First     Stromal             rpc         0.458  0.0950  0.146  0.15   
## 15 First     Hofbauer            reported    0.0410 0.663   0.663  0.66   
## 16 First     Hofbauer            rpc         0.278  0.224   0.264  0.26   
## 17 First     Endothelial         reported    0.716  0.0164  0.0698 0.07   
## 18 First     Endothelial         rpc         0.538  0.0608  0.135  0.14   
## 19 First     Syncytiotrophoblast reported    0.132  0.424   0.446  0.45   
## 20 First     Syncytiotrophoblast rpc         0.360  0.154   0.220  0.22
```

```r
pDat_ga_comp %>%
  pivot_wider(id_cols = -ga,
              names_from = 'ga',
              values_from = 'weeks') %>%
  ggplot(aes(x= reported, y = rpc)) +
  geom_point() +
  geom_smooth(method ='lm') +
  facet_wrap(vars(Trimester), scales = 'free')
```

```
## `geom_smooth()` using formula 'y ~ x'
```

```
## Warning: Removed 6 rows containing non-finite values (stat_smooth).
```

```
## Warning: Removed 6 rows containing missing values (geom_point).
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

```r
pDat_ga_comp %>% 
  filter(component != 'n_rbc',
         ga %in% c('reported', 'rpc'),
        Trimester == 'Third') %>%
  mutate(component = str_to_sentence(component)) %>%
  {
    ggplot(data = ., aes(x = weeks, y = proportion, color = component)) +
      geom_point() +
      geom_smooth(method = 'lm', se = FALSE,alpha = 0.5) +
      facet_grid(cols = vars(ga), rows = vars(component), scales = 'free_y') +
      scale_y_continuous(labels = percent) +  
      scale_color_manual(values = color_code_tissue[unique(.$component)], 
                       na.value = '#636363',
                       guide = 'none')  +
      theme(panel.grid = element_blank()) +
      labs(title= 'within-trimester association with gestational age')
  }
```

```
## `geom_smooth()` using formula 'y ~ x'
```

```
## Warning: Removed 5 rows containing non-finite values (stat_smooth).
```

```
## Warning: Removed 5 rows containing missing values (geom_point).
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-19-2.png)<!-- -->

# try villi again

amy reports that she gets nRBCs estimated on her villi samples. Let's try again...


```r
library(planet)
pDat_villi <- pDat_filt %>%
  filter(Tissue == 'Villi')


houseman_estimates <- minfi:::projectCellType(
  getBeta(mset)[rownames(coefs_combined_third),pDat_villi$Sample_Name],
  coefs_combined_third,
  lessThanOne = FALSE)

houseman_estimates %>%
  as.data.frame() %>%
  bind_cols(Sample_Name = rownames(.),.) %>%
  pivot_longer(cols =-Sample_Name,
               names_to = 'component',
               values_to = 'proportion') %>%
  left_join(pDat_villi %>% select(Sample_Name, Trimester)) %>%
  ggplot(aes(x = Sample_Name, y = proportion, fill = component)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Trimester, ncol = 2, scales = 'free') +
  scale_fill_brewer(palette = 'Accent') +
  scale_y_continuous(limits = c(-0.1,1.1), breaks = c(0, 0.5, 1), labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = '', fill = '')
```

```
## Joining, by = "Sample_Name"
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

# no STB,

just troph, endo, strom, hc, and nRBC


```r
# subset to reference
mset_combined_ref_nostb <- mset_combined[,pData(mset_combined)$CellType %in% c('Trophoblasts', 'Stromal', 'Hofbauer', 'Endothelial', 'nRBC')]

# pick probes
probes_third_nostb <- minfi:::pickCompProbes(
  mset_combined_ref_nostb[, pData(mset_combined_ref_nostb)$Trimester == 'Third'], 
  cellTypes = c('Trophoblasts', 'Stromal', 'Hofbauer', 'Endothelial', 'nRBC'),
  compositeCellType = 'Placenta',
  probeSelect = "both")

probes_first_nostb <- minfi:::pickCompProbes(
  mset_combined_ref_nostb[, pData(mset_combined_ref_nostb)$Trimester == 'First' |
             pData(mset_combined_ref_nostb)$CellType %in% c('nRBC')], 
  cellTypes = c('Trophoblasts', 'Stromal', 'Hofbauer', 'Endothelial', 'nRBC'),
  compositeCellType = 'Placenta',
  probeSelect = "both")

# extract coefficients
coefs_combined_third_nostb <- probes_third_nostb$coefEsts
coefs_combined_first_nostb <- probes_first_nostb$coefEsts



epidish_nostb_third <- epidish(
  beta.m = getBeta(mset_test)[rownames(coefs_combined_third_nostb)
                                ,pData(mset_test)$Trimester == 'Third'],
  ref.m = coefs_combined_third_nostb,
  method = 'RPC')

epidish_nostb_first <- epidish(
  beta.m = getBeta(mset_test)[rownames(coefs_combined_first_nostb)
                                ,pData(mset_test)$Trimester == 'First'],
  ref.m = coefs_combined_first_nostb,
  method = 'RPC')

epidish_nostb_third <- epidish_nostb_third$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (RPC)')

epidish_nostb_first <- epidish_nostb_first$estF %>%
  as.data.frame() %>%
  bind_cols(Sample_Names = rownames(.), .) %>%
  pivot_longer(
    cols = -Sample_Names,
    names_to = 'component',
    values_to = 'percent'
  ) %>%
  mutate(algorithm = 'epidish (RPC)')

epidish_nostb <- bind_rows(epidish_nostb_first, epidish_nostb_third) %>%
  dplyr::rename(Sample_Name = Sample_Names) %>%
  left_join(pDat_filt %>% select(Sample_Name, Sex, Trimester, Tissue, GA))  %>%
  distinct() 
```

```
## Joining, by = "Sample_Name"
```

## Stats


```r
# between trimester
epidish_nostb %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  nest(data = -component) %>%
  mutate(lm = map(data, ~lm(percent ~ Trimester, .)),
         glanced = map(lm, glance)) %>%
  select(-data, -lm)  %>%
  unnest(cols = c(glanced)) %>%
  select(component, r.squared, p.value) %>%
  mutate(adj_p = p.adjust(p.value, method = 'bonferroni'),
         p_value = pvalue(adj_p))
```

```
## # A tibble: 5 x 5
##   component    r.squared    p.value      adj_p p_value
##   <chr>            <dbl>      <dbl>      <dbl> <chr>  
## 1 Trophoblasts     0.796   9.36e-10   3.75e- 9 <0.001 
## 2 Stromal          0.845   3.23e-11   1.29e-10 <0.001 
## 3 Hofbauer         0.639   9.71e- 7   3.88e- 6 <0.001 
## 4 Endothelial      0.125   7.61e- 2   3.04e- 1 0.304  
## 5 nRBC           NaN     NaN        NaN        <NA>
```

```r
# between sex
stats_by_sex_nostb <- epidish_nostb %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  filter(Trimester == 'Third', component != 'nRBC') %>%
  nest(data = -component) %>%
  mutate(lm = map(data, ~lm(percent ~ Sex, .)),
         glanced = map(lm, glance)) %>%
  select(-data, -lm)  %>%
  unnest(cols = c(glanced)) %>%
  select(component, r.squared, p.value) %>%
  mutate(adj_p = p.adjust(p.value, method = 'bonferroni'),
         p_value = pvalue(adj_p),
         testing_variable = 'Sex');stats_by_sex_nostb
```

```
## # A tibble: 4 x 6
##   component    r.squared p.value adj_p p_value testing_variable
##   <chr>            <dbl>   <dbl> <dbl> <chr>   <chr>           
## 1 Trophoblasts    0.125   0.137  0.548 0.548   Sex             
## 2 Stromal         0.0595  0.314  1     >0.999  Sex             
## 3 Hofbauer        0.236   0.0349 0.140 0.140   Sex             
## 4 Endothelial     0.0770  0.250  1.00  >0.999  Sex
```

```r
epidish_nostb %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  filter(Trimester == 'Third', component !='nRBC') %>%
  ggplot(aes(x = component, y = percent, color = Sex)) +
  geom_beeswarm(dodge.width = 0.75, cex = 1.5) +
  #geom_jitter(position = position_jitterdodge(dodge.width = 0.75)) +
  scale_y_continuous(labels = function(x)percent(x,accuracy = 1)) +
  scale_color_manual(values = c('#af8dc3', '#7fbf7b')) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  labs(x = '', color = '', y = '') 
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```r
epidish_nostb %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  filter(Trimester == 'Third') %>%
  select(Sample_Name, Sex) %>%
  distinct() %>%
  dplyr::count(Sex)
```

```
## # A tibble: 2 x 2
##   Sex       n
##   <chr> <int>
## 1 F        10
## 2 M         9
```

```r
# between ancestry
epidish_nostb  %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  filter(Trimester == 'Third') %>%
  left_join(ancestry) %>%
  ggplot(aes(x = PC1_geno, y = PC2_geno, color = Predicted_ethnicity)) +
  geom_point()
```

```
## Joining, by = "Sample_Name"
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-22-2.png)<!-- -->

```r
epidish_nostb %>%
  left_join(ancestry) %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  filter(Trimester == 'Third', Predicted_ethnicity != 'Ambiguous', component != 'nRBC') %>%
  ggplot(aes(x = component, y = percent, color = Predicted_ethnicity)) +
  geom_beeswarm(dodge.width = 0.75,size = 0.9) +
  #geom_boxplot() +
  #geom_jitter(position = position_jitterdodge(dodge.width = 0.75)) +
  scale_y_continuous(labels = percent) +
  scale_color_manual(values = c('#d8b365', '#5ab4ac')) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  labs(x = '', color = '', y = '')  
```

```
## Joining, by = "Sample_Name"
```

![](2_14_deconvolution_files/figure-html/unnamed-chunk-22-3.png)<!-- -->

```r
stats_by_ancestry_nostb <- epidish_nostb %>%
  left_join(ancestry) %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  filter(Trimester == 'Third', Predicted_ethnicity != 'Ambiguous', component != 'nRBC') %>%
  nest(data = -component) %>%
  mutate(lm = map(data, ~lm(percent ~ Predicted_ethnicity, .)),
         glanced = map(lm, glance)) %>%
  select(-data, -lm)  %>%
  unnest(cols = c(glanced)) %>%
  select(component, r.squared, p.value) %>%
  mutate(adj_p = p.adjust(p.value, method = 'bonferroni'),
         p_value = pvalue(adj_p),
         testing_variable = 'ancestry'); stats_by_ancestry_nostb
```

```
## Joining, by = "Sample_Name"
```

```
## # A tibble: 4 x 6
##   component    r.squared p.value adj_p p_value testing_variable
##   <chr>            <dbl>   <dbl> <dbl> <chr>   <chr>           
## 1 Trophoblasts   0.0122    0.673     1 >0.999  ancestry        
## 2 Stromal        0.0102    0.699     1 >0.999  ancestry        
## 3 Hofbauer       0.00125   0.893     1 >0.999  ancestry        
## 4 Endothelial    0.0363    0.464     1 >0.999  ancestry
```

```r
epidish_nostb %>%
  left_join(ancestry) %>%
  filter(Tissue == 'Villi', algorithm == 'epidish (RPC)') %>%
  select(Sample_Name, Predicted_ethnicity) %>%
  distinct() %>%
  dplyr::count(Predicted_ethnicity)
```

```
## Joining, by = "Sample_Name"
```

```
## # A tibble: 4 x 2
##   Predicted_ethnicity     n
##   <chr>               <int>
## 1 Ambiguous               2
## 2 Asian                   6
## 3 Caucasian              11
## 4 <NA>                    7
```

# Save

Join to pDat_filt and save


```r
pDat_filt %>%
  saveRDS(here(base_path, '2_14_pDat_filt.rds'))
```


```r
# 3rd and 1st trimester reference cpgs
coefs_combined_third %>%
  write_csv( na = 'NA',
             path = here('outs', '2_14_deconvolution_reference_cpgs_third.csv'))

coefs_combined_first %>%
  write_csv(na = 'NA',
            path = here('outs', '2_14_deconvolution_reference_cpgs_first.csv'))

# 3rd and 1st trimester reference cpgs
coefs_combined_third_nostb %>%
  saveRDS(here::here('data', 'main', 'processed', '2_14_coefs_combined_third_nostb.rds'))

coefs_combined_first_nostb %>%
  saveRDS(here::here('data', 'main', 'processed', '2_14_coefs_combined_first_nostb.rds'))

# save deconvolution (nrbc included) estimates on villi
epidish_rpc %>%
  left_join(ancestry) %>%
  saveRDS(here(base_path, '2_14_deconvolution_results.rds'))
epidish_nostb %>%
  left_join(ancestry) %>%
  saveRDS(here(base_path, '2_14_deconvolution_results_nostb.rds'))

# save statistical testing
stats_by_ga %>%
  mutate(ga = ifelse(ga == 'rpc', 'estimated', as.character(ga))) %>%
  mutate(testing_variable = paste('within ', Trimester, ' trimester-', ga, ' - ga', sep = '')) %>%
  select(-Trimester, -ga)  %>%
  bind_rows(., stats_by_sex, stats_by_ancestry) %>%
  write_csv(na = '', path = here('outs', '2_14_stats_by_sex_ancestry.csv'))
bind_rows(stats_by_sex_nostb, stats_by_ancestry_nostb) %>%
  write_csv(na = '', path = here('outs', '2_14_stats_by_sex_ancestry_nostb.csv'))

# in silico stuff
epidish_results %>%
  saveRDS(here(base_path, '2_14_in_silico_deconvolution_results.rds'))
stats %>%
  select(-label) %>%
  write_csv(here('outs', '2_14_in_silico_deconvolution_stats.csv'))
```

