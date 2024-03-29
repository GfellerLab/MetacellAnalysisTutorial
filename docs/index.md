---
title: "Metacell Analysis Tutorial"
author: "Aurélie Gabriel, Léonard Hérault, Mariia Bilous, David Gfeller"
date: "2024-03-07"
site: bookdown::bookdown_site
documentclass: book
bibliography: [citations.bib, packages.bib]
url: https://GfellerLab.github.io/MetacellAnalysisTutorial
description: |
  This is a tutorial about construction of metacells for single-cell RNA-seq data and analysis of metacells.
link-citations: yes
github-repo: GfellerLab/MetacellAnalysisTutorial
biblio-style: apalike
csl: biochemistry-and-cell-biology.csl
---

# About this tutorial{.unnumbered}

In this tutorial, we describe the different steps to build metacells [@baran_metacell_2019] from single-cell data using three frameworks:
i. SuperCell [@SuperCell] (tutorial in \@ref(SuperCell-construction)),
ii. Metacells version 2 (MC2 [@MC2]) (tutorial in \@ref(MC2-construction)),
iii. SEACells [@SEACells] (tutorial in \@ref(SEACells-construction)).

We provide examples of downstream analyses (section \@ref(downstream-analysis)) performed at the metacell level (e.g clustering, differential analysis, marker identification).

We also show how to obtain metacells by running these methods using a **command line tool** that we provide as part of the MetacellAnalysisToolkit (**MATK**) [github repository](https://github.com/GfellerLab/MetacellAnalysisToolkit).
This repository also contains the **MetacellAnalysisToolkit R package** which provides R functions to compute QC
metrics and visualization functions to evaluate the quality of metacells.

In chapter \@ref(command-line), we use MATK to build metacells on a continuous dataset of CD34+ cells and use the R package to visualize the constructed metacells.
The use of the MetacellAnalysisToolkit R package to visualize and evaluate the quality of the metacells is described in chapter \@ref(QCs) and \@ref(command-line).

Finally, using the MATK command line tool, we show in section \@ref(integration) how we can use metacells to **integrate and analyze around 580,000 cells** from 107 individuals distributed in 166 samples.
