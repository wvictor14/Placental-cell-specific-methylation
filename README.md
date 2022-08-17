---
title: "Placental Cell DNA Methylation"
author: "Victor Yuan"
output:
    html_document:
        keep_md: yes
        self_contained: yes
        theme: spacelab
        toc: yes
        toc_depth: 3
        toc_float:
          collapsed: no
editor_options:
  chunk_output_type: console
---



# README

This repository contains my R scripts for data processing and analysis of the NIH EPIC cells project, which was published in [BMC genomics](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07186-6).

## Organization of repo

R scripts and folders are named in the following system:

* Directories are numbered by the order their contained scripts in the workflow. For example, scripts
in `0_data` are ran before `1_qc_preprocess`, which is before `2_analysis.`
* Scripts are numbered by their directory and workflow order as well. `0_1_load_data.Rmd` is in the 
`0_data` directory and is ran before `0_2_check_data.Rmd.`
* Output from scripts is stored in either `main/interim` or `main/processed`. Output is prefixed by the
same prefix of the script used to generate them. For example, `0_1_load_data.Rmd` produces 
`0_1_rgset.rds`, which is stored in `main/interim`.

# Directory structure

```
Placental Cell DNA Methylation                  
|   README
│   
└───data
│   └───external                  <- Data from third party sources (e.g. annotations)
│   └───main
|           └───raw               <- The original, immutable data (e.g. DNAme, pData)          
|           └───interim           <- Intermediate data that has been transformed
│           └───processed         <- The final, canonical data sets for analysis/sharing
|
└───R
|   └───0_data                   <- Downloading & loading data into R
|   └───1_qc_preprocess          <- Quality control and preprocessing data
|   └───2_analysis               <- Analyzing processed data
|   └───3_app                    <- Some scripts related to the shiny app
|   └───4_figures                <- Figures for manuscript
|   └───functions                 
|
└───outs                         <- Figures, tables, reports
|
└───Manuscript                   <- text and figures for publication
|
└───Reports                      <- ppt or rmd-generated progress reports
```

# TO USE DATA FROM THESE SAMPLES

for Robinson lab members with access to the shared drive.

for using the data that has been processed and filtered.

This should be all the information about the data files that you need to find. Note that these data files contain all 192 samples ran across the 2 weeks worth of arrays.

**raw betas / rgset**:
```
rgset <- readRDS('data/main/interim/0_1_rgset_raw.rds') # where ../../ is Victor/Projects/NIH - cells/
```

**pDat**:

```
pDat <- readRDS('data/main/interim/2_3_pDat_contam.rds')
```

**noob normalized filtered data**. This one has been filtered for probes.
```
betas <- readRDS('data/main/interim/1_4_betas_noob_filt.rds'))
```

**noob normalized unfiltered data**:
```
mset_noob <- readRDS('data/main/interim/1_4_mset_noob.rds')
```

And lastly, this some code I run every script to remove samples based on QC. You might want to alter this depending on what you're doing (probably you are only interested in a subset of the 192 samples). 

The reasoning on why these samples are being removed are in the QC scripts `R/1_qc_preprocess`, also see `Reports/11_QC_Report.[pdf/html]`, `Reports/21_...html` files, and earliest few `pptx` presentations in `Reports/`.

```
pDat_filt <- pDat %>% 
  filter(maternal_contamination_norm_flip < 0.35,
         !Sample_Name %in% c('PM364_hofb_cs', 'PL293_v_R2', 'PM366_vc_R2', 'P131_hofb_cs')))

# then filter betas 
```
