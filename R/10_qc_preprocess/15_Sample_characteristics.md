---
title: "15_Sample_Characteristics"
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

I realized some variables were had an incorrect class (e.g. factor when should be numeric), and some
observations needed to be replaced with NAs (0 for DNA concentration means we did not measure).

I made edits in my 11_QC script to address these issues. 
I also forgot to test for confounders / relationships between sample characteristics.


# Setup


```r
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(viridis)
library(yahew)
library(egg)

pDat <- readRDS('../../data/main/interim/13_pDat.rds')
# PCA association data for plotting
pca <- readRDS('../../data/main/interim/11_pca_association_plot.rds')
```

# Relationships between variables


```r
glimpse(pDat)
```

```
## Observations: 192
## Variables: 65
## $ Sample_Name                  <chr> "PM365_endo_cs", "PL295_discard_c...
## $ Chip_number                  <fct> 25, 25, 25, 25, 25, 25, 25, 25, 2...
## $ Well                         <chr> "A1", "B1", "C1", "D1", "E1", "F1...
## $ Case_ID                      <chr> "PM365", "PL295", "PL293", "PL292...
## $ Sex                          <chr> "F", "F", "M", "M", "F", "F", "F"...
## $ GA                           <dbl> NA, 11.6, 9.6, 7.1, NA, NA, NA, N...
## $ Trimester                    <chr> "Third", "First", "First", "First...
## $ DNA_QP                       <dbl> 542.5, 701.5, 412.5, 151.5, NA, N...
## $ Week                         <fct> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
## $ Sample_Plate                 <chr> "WG6980707-MSA4", "WG6980707-MSA4...
## $ Tissue                       <chr> "Endothelial", "Dead Cells and Ly...
## $ Sentrix_ID                   <dbl> 203067920143, 203067920143, 20306...
## $ Sentrix_Position             <chr> "R01C01", "R02C01", "R03C01", "R0...
## $ Scratches                    <chr> NA, NA, NA, NA, NA, NA, NA, NA, N...
## $ Batch_BSC                    <chr> "BSC1", "BSC1", "BSC1", "BSC1", "...
## $ DNA_conc_BSC_adjusted        <dbl> 30.19, 42.65, 20.99, 8.04, 49.91,...
## $ DNA_conc_before_load         <dbl> 50.31667, 42.65000, 52.47500, 20....
## $ DNA_loaded                   <dbl> 201.2667, 170.6000, 209.9000, 80....
## $ Sentrix                      <chr> "203067920143_R01C01", "203067920...
## $ Colors_Tissue                <chr> "#6A1B9A", "#424242", "#388E3C", ...
## $ Colors_Sex                   <chr> "#F8BBD0", "#F8BBD0", "#BBDEFB", ...
## $ Colors_Trimester             <chr> "#212121", "#E0E0E0", "#E0E0E0", ...
## $ Row_numeric                  <dbl> 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, ...
## $ Row_factor                   <fct> 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, ...
## $ detP_01_minfi                <dbl> 1161, 770, 1815, 874, 1711, 1254,...
## $ beadcount                    <dbl> 2943, 4829, 2755, 2838, 2743, 400...
## $ failed_probes                <dbl> 4001, 5464, 4474, 3632, 4349, 512...
## $ log2_median_meth             <dbl> 11.67154, 11.98050, 11.57033, 11....
## $ log2_median_unmeth           <dbl> 12.21310, 12.34346, 11.99965, 12....
## $ Average_intensity            <dbl> 9904.426, 10495.819, 9338.025, 10...
## $ median_X_intensity           <dbl> 10309.5, 11025.0, 6131.0, 8722.5,...
## $ median_Y_intensity           <dbl> 424, 328, 5283, 7065, 519, 487, 5...
## $ normalized_X_intensity       <dbl> 1.0155989, 1.0403758, 0.7658799, ...
## $ normalized_Y_intensity       <dbl> 0.2040305, 0.1936692, 0.7447600, ...
## $ Flag_Sex                     <lgl> FALSE, FALSE, FALSE, FALSE, FALSE...
## $ genotype_cluster             <fct> 1, 2, 3, 4, 5, 6, 7, 8, 8, 3, 9, ...
## $ Flag_genotype                <lgl> FALSE, FALSE, FALSE, FALSE, FALSE...
## $ Average_SNP_cor_within_donor <dbl> 0.9954983, 0.9706685, 0.9792740, ...
## $ Prob_SNP_outlier             <dbl> 0.08027823, 0.36363599, 0.2981987...
## $ Prob_SNP_outlier_Logodds     <dbl> -3.9498163, -0.7531178, -1.432664...
## $ Agreement_to_donor_villi     <dbl> 1.0000000, 0.9999924, 0.9999983, ...
## $ Agreement_to_unrelated       <dbl> 0.3678490, 0.4355669, 0.4306673, ...
## $ cor_to_donor_villi           <dbl> 0.9960673, 0.9545717, 0.9923917, ...
## $ cor_to_unrelated             <dbl> 0.06523967, 0.18140332, 0.0983996...
## $ cor_to_reference             <dbl> 0.9960673, 0.9545717, 0.9923917, ...
## $ PC1_raw                      <dbl> 39.733547, 23.431438, 7.710260, 4...
## $ PC2_raw                      <dbl> 63.342799, -5.405307, 10.856889, ...
## $ PC3_raw                      <dbl> 90.7902133, -5.7615785, -16.20630...
## $ PC4_raw                      <dbl> 1.464261, -33.965841, -53.532079,...
## $ PC5_raw                      <dbl> 5.9254346, 18.2149374, -5.2857755...
## $ PC6_raw                      <dbl> 10.6816352, 3.8324827, -12.003744...
## $ PC7_raw                      <dbl> 1.4805490, 13.9683012, -0.6926468...
## $ PC8_raw                      <dbl> -2.1988821, 0.8546499, 3.9389876,...
## $ PC9_raw                      <dbl> 1.8757038, 1.4294838, -14.6485750...
## $ PC10_raw                     <dbl> -2.8203244, -1.0747330, 0.9101450...
## $ PC11_raw                     <dbl> 0.02443068, 4.51216170, -9.110262...
## $ PC12_raw                     <dbl> -5.2314949, -6.7112737, -10.76333...
## $ PC13_raw                     <dbl> 6.4089765, -2.3649324, 15.2202045...
## $ PC14_raw                     <dbl> -3.3441211, -3.6917725, 2.6233668...
## $ PC15_raw                     <dbl> 0.01310273, 1.55051127, -1.627309...
## $ PC16_raw                     <dbl> -2.89677945, -3.08055626, -5.1630...
## $ PC17_raw                     <dbl> 1.6113142, -0.1693311, -8.6920188...
## $ PC18_raw                     <dbl> 8.0795749, 2.0039087, -3.5425923,...
## $ PC19_raw                     <dbl> -1.1727345, -5.1312003, 4.7512132...
## $ PC20_raw                     <dbl> -1.58186005, -2.98064749, -1.5496...
```

```r
cov_tests <- pDat %>%
  select(Row_numeric, Row_factor, Chip_number, Week,Batch_BSC,
         DNA_conc_before_load, DNA_loaded,  
         failed_probes, detP_01_minfi , beadcount, Average_intensity, 
         Case_ID, Sex, Tissue,  Trimester,
         cor_to_reference, cor_to_unrelated, Prob_SNP_outlier) %>%
  as.data.frame() %>% pairtest

# make categories
cov_tests <- cov_tests %>% 
  mutate(pval_cat = if_else(p.value < 0.001, '< 0.001',
                            if_else(p.value < 0.01, '< 0.01',
                                    if_else(p.value < 0.05, '< 0.05', '> 0.05'))))

# plot heatmap of associations
ggplot(cov_tests, aes(x=Column, y = Row, fill = pval_cat)) +
  geom_tile(col = 'grey') + theme_bw() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c('> 0.05' = 'White', '< 0.05' = '#fee8c8', 
                               '< 0.01' = '#fdbb84', '< 0.001' = '#e34a33')) +
  labs(x = '', y = '', fill = 'p') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust= 1),
        panel.grid.major.x = element_blank()) +
  coord_equal()
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-2-1.png)<!-- -->


```r
qplot(data = pDat,x = DNA_loaded, y = failed_probes) + geom_smooth()
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
qplot(data = pDat,x = DNA_conc_before_load, y = failed_probes) + geom_smooth()
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

```r
qplot(data = pDat, x = Row_numeric, y = failed_probes) + geom_smooth()
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-3-3.png)<!-- -->

```r
qplot(data = pDat, x = Row_numeric, y = detP_01_minfi) + geom_smooth()
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-3-4.png)<!-- -->

```r
qplot(data = pDat, x = Week, y = detP_01_minfi, geom = 'boxplot') 
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-3-5.png)<!-- -->

```r
qplot(data = pDat, x = Row_numeric, y = Average_intensity) + geom_smooth()
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-3-6.png)<!-- -->

```r
qplot(data = pDat, x = as.numeric(as.factor(Batch_BSC)), y =Average_intensity) + 
  geom_smooth(method = 'lm') + 
  scale_x_continuous(labels = levels(as.factor(pDat$Batch_BSC))) +
  scale_y_continuous(limits = c(0, 14500))
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-3-7.png)<!-- -->

```r
  labs(x = 'BSC Batch')
```

```
## $x
## [1] "BSC Batch"
## 
## attr(,"class")
## [1] "labels"
```

```r
qplot(data = pDat, x = Average_intensity, y = detP_01_minfi)
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-3-8.png)<!-- -->

```r
qplot(data = pDat, x = Average_intensity, y = beadcount)
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-3-9.png)<!-- -->

```r
qplot(data = pDat, x = Sex, y = detP_01_minfi)
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-3-10.png)<!-- -->

```r
ggplot(pDat, aes(x = Tissue, fill = Batch_BSC)) + geom_bar(position = 'dodge')
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-3-11.png)<!-- -->

# PCA


```r
# create color palette
colpal <- c('white', '#fee8c8', '#fdbb84', '#e34a33')
names(colpal) <- levels(pca$Plot_data$pval_cat)

p1 <- ggplot(pca$Plot_data, aes(x = PC, y = dep, fill = pval_cat)) +
  geom_tile(col = 'lightgrey') + theme_bw() +
  scale_x_discrete(expand = c(0, 0), labels = 1:20) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = colpal)  + 
  labs(y = '', fill = 'P value')

p1b <- ggplot(pca$PC_variance, aes(x = PC, y = Prop_var_raw)) +
  geom_bar(stat = 'identity') + theme_bw() + labs(y = '% variance') +
  scale_x_continuous(breaks = 1:20)

ggarrange(p1, p1b, heights = c(3,1))
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-4-1.png)<!-- -->


```r
qplot(data = pDat, x = Row_numeric, y = PC5_raw)
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
qplot(data = pDat, x = Row_numeric, y = PC6_raw)
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-5-2.png)<!-- -->

```r
qplot(data = pDat, x = Row_numeric, y = PC7_raw)
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-5-3.png)<!-- -->

```r
qplot(data = pDat, x = Row_numeric, y = PC8_raw)
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-5-4.png)<!-- -->

```r
qplot(data = pDat, x = Row_numeric, y = PC10_raw)
```

![](15_Sample_characteristics_files/figure-html/unnamed-chunk-5-5.png)<!-- -->
