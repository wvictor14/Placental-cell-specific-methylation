---
title: "README"
author: "Victor Yuan"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 4
    toc_float:
      collapsed: false
    theme: spacelab
---

# Placental Cell DNA Methylation

For NIH EPIC cells project

## Directory structure

```
Placental Cell DNA Methylation
|   README.html
│   README.md                   
|   README.Rmd 
│   
└───data
│   └───external                  <- Data from third party sources (e.g. annotations)
│   └───main
|           └───raw               <- The original, immutable data (e.g. DNAme, pData)          
|           └───interim           <- Intermediate data that has been transformed
│           └───processed         <- The final, canonical data sets for analysis/sharing
|
└───R
|   └───00_Data                   <- Downloading & loading data into R
|   └───10_qc_preprocess          <- Quality control and preprocessing data
|   └───20_analysis               <- Analyzing processed data
|   └───functions                 
|
└───outs                          <- Figures, tables, reports
```

## Data

## Analysis

## References

[1] Mayne, B. T., Leemaqz, S. Y., Smith, A. K., Breen, J., Roberts, C. T., & Bianco-Miotto, T. 
(2016). Accelerated placental aging in early onset preeclampsia pregnancies identified by DNA 
methylation. Epigenomics, 9(3), 279-289.


