---
title: "ubi-Brk_3h"
author: "Masha"
date: "2024-04-15"
output: 
  html_document: 
    keep_md: yes
---



## Introduction

This is RNAseq analysis using deseq2. Two conditions: "control" and "ubi-brk", 3 repeats for each condition.

Experimental setup:

**SalE/PV-ShineGal4\>Brk(x)** – dark control or "control" in this analysis;

**ubi-ShineGal4\>Brk(x)** - 3h light induction of UAS-Brk overexpression, called "ubi-brk" here;

SalE/PV-ShineGal4\>Brk(x) - 3h ight induction of UAS-Brk overexpression in Sal compartment - not included in this analysis due to low number of def expressed genes.

All conditions were processed at the same time with several rounds of sample collection (different days, several crosses). For each repeat 50 discs were combined in one tube (could use up to 35 next time).

### Loading libraries


```r
library(DESeq2)
library(apeglm)
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(plotly)
library(magrittr)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(writexl)
```

### Load data

