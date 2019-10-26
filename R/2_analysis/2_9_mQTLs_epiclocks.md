---
title: "2_9_mQTLs_epiclocks"
author: "Victor Yuan"
date: "24/10/2019"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 4
    toc_float:
      collapsed: false
    theme: spacelab
    highlight: tango 
#tango, pygments, kate, monochrome, espresso, zenburn, haddock, textmate
editor_options: 
  chunk_output_type: console
---

In this script, I analyze the overlap between cell specific sites and those sites in:

1) mQTLs
2) placental epigenetic clocks

# Setup

## Libraries


```r
# libraries and data
library(minfi)
library(tidyverse)
library(scales)
library(here)
library(readxl)
library(planet)
theme_set(theme_bw())
```

## Data


```r
base_path <- file.path('data', 'main', 'interim')

# pData
pDat <- readRDS(here(base_path, '2_3_pDat_contam.rds'))
pDat <- pDat %>%
  mutate(Tissue = case_when(
    !(Tissue %in% c('Villi', 'Villi maternal', 'Syncytiotrophoblast')) ~ paste(Tissue, 'cs'),
    Tissue == 'Syncytiotrophoblast' ~ 'Trophoblasts enz',
    TRUE ~ Tissue
  )) 

# raw methylation data
betas <- readRDS(here(base_path, '1_4_betas_noob_filt.rds'))


# annotation
anno <- readRDS('Z:/Victor/Repositories/EPIC_annotation/hg19_epic_annotation.rds')
anno <- anno %>%
  as_tibble() %>%
  filter(cpg %in% rownames(betas)) # filter to filtered betas cpgs
probe_anno <- readRDS(here(base_path, '1_1_probe_anno.rds'))

# color key
pheatmap_color_code <- readRDS(here(base_path, '1_1_color_code.rds'))

color_code <- readRDS(here(base_path, '2_3_color_code.rds'))
color_code_tissue <- setNames(color_code$Colors_Tissue, color_code$label)

#dmcs
dmcs_cell <- readRDS(here(base_path, '2_4_dmcs.rds'))
dmcs_trim <- readRDS(here(base_path, '2_8_all_third_vs_first_dmcs.rds'))

# enrichment results
func_enrich <- readRDS(here(base_path, '2_8_functional_enrichment_results.rds'))
tests <- readRDS(here(base_path, '2_8_genomic_enrichment_results.rds'))
```

## Remove samples


```r
pDat_filt <- pDat %>% 
  filter(maternal_contamination_norm_flip < 0.35,
         !Sample_Name %in% c('PM364_hofb_cs', 'PL293_v_R2', 'PM366_vc_R2', 'P131_hofb_cs', 
                             'PM324_V4', 'PM324_V1'),
         !Tissue %in% c('Villi maternal', 'Trophoblasts enz', 'Mixture cs', 
                        'Dead Cells and Lymphocytes cs'),)

# filter to first trimester
betas_filt <- betas[,pDat_filt$Sentrix]
colnames(betas_filt) <- pDat_filt$Sample_Name
```

# mQTLs


```r
#load mqtls in
mqtls <- read_excel(here('data', 'external', 'journal.pgen.1007785.s018.xlsx'), skip = 1)
```

```
## Warning in read_fun(path = enc2native(normalizePath(path)), sheet_i =
## sheet, : Expecting numeric in H1586 / R1586C8: got '5.11901415656116e-320'
```

```
## Warning in read_fun(path = enc2native(normalizePath(path)), sheet_i =
## sheet, : Expecting numeric in I1586 / R1586C9: got '2.41092769238639e-314'
```

```r
sum(anno$cpg %in% mqtls$cpgID) # 3313 / 4342 mqtl cpgs are in epic
```

```
## [1] 3313
```

```r
# The number of mQTLs athat are also cell-specific:
dmcs_cell %>%
  filter(bonferroni < 0.01 & abs(delta_b) > 0.25) %>%
  separate(col = Group1, sep = '\\.', into = c('Trimester', 'Celltype')) %>%
  filter(Trimester == 'Third') %>%
  select(gene) %>%
  distinct() %>%
  filter(gene %in% mqtls$cpgID) %>%
  nrow() #1022
```

```
## [1] 1483
```

1022 / 3313 (30.85%) cpgs that are placental-specific mQTLs are also DMCs


```r
# calculate number of DMCs that are mQTLs
dmcs_cell %>%
  mutate(placental_mqtl = gene %in% mqtls$cpgID,
         sig = (bonferroni < 0.01 & abs(delta_b) > 0.25)) %>%
  separate(col = Group1, sep = '\\.', into = c('Trimester', 'Celltype')) %>%
  filter(Trimester == 'Third') %>%
  
  group_by(Celltype) %>%
  summarize(
    # total DMCs / mQTLs
    n_sig = sum(sig == TRUE),
    n_mqtl = sum(placental_mqtl == TRUE),
    
    # DMC that is also mqtl
    dmc_and_mqtl = sum(
      (placental_mqtl == TRUE) & 
      (sig == TRUE)),
    
    # DMC that is mQTL / total DMCs
    p_dmc_in_mqtl = percent(dmc_and_mqtl/n_sig),
    
    # DMC that is mQTL / total mQTLs
    p_mqtl_in_dmc = percent(dmc_and_mqtl/n_mqtl))
```

```
## # A tibble: 4 x 6
##   Celltype  n_sig n_mqtl dmc_and_mqtl p_dmc_in_mqtl p_mqtl_in_dmc
##   <chr>     <int>  <int>        <int> <chr>         <chr>        
## 1 Endo_cs   75525   3313          384 0.508%        11.6%        
## 2 Hofb_cs  130733   3313          705 0.539%        21.3%        
## 3 Strom_cs  80153   3313          431 0.538%        13.0%        
## 4 Troph_cs 135553   3313          710 0.524%        21.4%
```

Let's visualize some of these mQTLs that are also cell specific

It'll be good to have ethnicity information. We don't have this currently, so I will use inferred 
ethnicity from planet for now, which I calculated in `1_6_ethnicity_ancestry.html`.


```r
# get inferred ethnicity info
pDat_eth <- readRDS(here(base_path, '1_6_pDat.rds'))
pDat_filt <- pDat_eth %>% 
  select(Sample_Name, contains('ethnicity'), contains('Prob'), -contains('Logodds')) %>%
  left_join(pDat_filt, .)
```

```
## Joining, by = c("Sample_Name", "failed_probes", "Prob_SNP_outlier")
```

```r
#list of dmcs that are also mqtls
dmc_and_mqtl <- dmcs_cell %>%
  filter(bonferroni < 0.01 & abs(delta_b) > 0.25) %>%
  separate(col = Group1, sep = '\\.', into = c('Trimester', 'Celltype')) %>%
  filter(Trimester == 'Third',
         gene %in% mqtls$cpgID) %>%
  dplyr::rename(dmc_for_class = Celltype) %>%
  select(gene, dmc_for_class, bonferroni, delta_b)

dmcs_and_mqtl_betas <- betas_filt[intersect(rownames(betas_filt), dmc_and_mqtl$gene),] %>%
  
  # transpose, remove rownames
  t() %>%
  as.data.frame() %>%
  bind_cols(Sample_Name = rownames(.), .) %>%
  
  # bind to pData
  left_join(pDat_filt %>% 
              select(Sample_Name,
                     Trimester, 
                     Tissue, 
                     contains('Prob_'), 
                     -contains('SNP'),
                     Predicted_ethnicity_nothresh), 
            .,
            by = 'Sample_Name') %>%
  filter(Trimester == 'Third') %>%
  select(-Trimester) %>%
  
  pivot_longer(cols = -c(Sample_Name:Predicted_ethnicity_nothresh),
               names_to = 'cpg_id',
               values_to = 'beta') %>%
  
  # join to dmc data
  left_join(dmc_and_mqtl, by = c('cpg_id' = 'gene'))

# plot function:
plot_dmcs <- function(data) {
  data %>%
    {
      ggplot(data = ., aes(x = Tissue, y= beta, color = Tissue)) +
        ggbeeswarm::geom_beeswarm(cex = 2,
                                  size = 0.7,
                                  priority = 'density',
                                  dodge.width = 0.8,
                                  aes(shape = Predicted_ethnicity_nothresh),
                                  fill = 'white',
                                  stroke = 1) +
        facet_wrap(~cpg_id, ncol = 2) +
        scale_color_manual(values = color_code_tissue[pDat_filt$Tissue],
                           guide = guide_legend(override.aes = list(size = 2))) +
        scale_y_continuous(limits = c(0,1)) +
        scale_shape_manual(values = c('Asian' = 21, 'Caucasian' = 16),
                           guide = guide_legend(override.aes = list(size = 2))) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank()) +
        labs(x = '', y= '', shape = 'Shape', color = 'Color') 
  }
}

dmcs_and_mqtl_betas %>%
  filter(dmc_for_class == 'Hofb_cs') %>%
  group_by(Sample_Name) %>%
  arrange(cpg_id) %>%
  dplyr::slice(1:8) %>%
  plot_dmcs() +
  labs(title = 'Hofbauer cells DMCs')
```

![](2_9_mQTLs_epiclocks_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
dmcs_and_mqtl_betas %>%
  filter(dmc_for_class == 'Troph_cs') %>%
  group_by(Sample_Name) %>%
  arrange(cpg_id) %>%
  dplyr::slice(1:8) %>%
  plot_dmcs() +
  labs(title = 'Trophoblast DMCs')
```

![](2_9_mQTLs_epiclocks_files/figure-html/unnamed-chunk-6-2.png)<!-- -->

```r
dmcs_and_mqtl_betas %>%
  filter(dmc_for_class == 'Endo_cs') %>%
  group_by(Sample_Name) %>%
  arrange(cpg_id) %>%
  dplyr::slice(1:8) %>%
  plot_dmcs() +
  labs(title = 'Endothelial cells DMCs')
```

![](2_9_mQTLs_epiclocks_files/figure-html/unnamed-chunk-6-3.png)<!-- -->

```r
dmcs_and_mqtl_betas %>%
  filter(dmc_for_class == 'Strom_cs') %>%
  group_by(Sample_Name) %>%
  arrange(cpg_id) %>%
  dplyr::slice(1:8) %>%
  plot_dmcs() +
  labs(title = 'Stromal cells DMCs')
```

![](2_9_mQTLs_epiclocks_files/figure-html/unnamed-chunk-6-4.png)<!-- -->

Now let's investigate some mqtls they showed in Delhaye et al 2019

Only 3 out of the 4 cpgs are in my data. Either do not exist in EPIC, or were filtered out for quality.


```r
dmcs_and_mqtl_betas %>%
  filter(cpg_id %in% c('cg08177731', 'cg21139150', 'cg17519949', 'cg26371521')) %>%
  
  {
    ggplot(data = ., aes(x = Tissue, y= beta, color = Tissue)) +
      ggbeeswarm::geom_beeswarm(cex = 0.6,
                                size = 1, 
                                priority = 'density',
                                dodge.width = 0.5,
                                aes(shape = Predicted_ethnicity_nothresh),
                                fill = 'white',
                                stroke = 1.5) +
      facet_wrap(~cpg_id) +
      scale_color_manual(values = color_code_tissue[pDat_filt$Tissue]) +
      scale_y_continuous(limits = c(0,1)) +
      scale_shape_manual(values = c('Asian' = 21, 'Caucasian' = 16)) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank()) +
      labs(x = '', y= '', shape = 'Shape', color = 'Color')
  }
```

![](2_9_mQTLs_epiclocks_files/figure-html/unnamed-chunk-7-1.png)<!-- -->


