# MEcell

==========
* [Introduction](#Introduction)
* [Installation](#Installation)
* [Function](#Function)
* [Tutorial](#Tutorial)
* [Citation](#Citation)

<a name="Introduction"/>

# Introduction

MEcell is a R package to leverage both intrinsic molecular profiles and external microenvironment features for spatial omics cell clustering.

<a name="Installation"/>

# Installation

```
devtools::install_github("liuqivandy/MEcell")
```

<a name="Function"/>

# Function

Taking a Seurat object as input, MEcell calculates the microenvironment profile for each cell, based on which to adjust the nearest neighbor graph from the intrinsic molecular profiles. MEcell returns a microenvironment-aware SNN graph (named MEcell) stored  in the respective slot.

```R
obj <- MEcell(obj)
obj<-FindClusters(obj,graph.name="MEcell")
```


<a name="Tutorial"/>

# Tutorial

- [Analyze MERFISH data to identify region-specific subtypes](https://htmlpreview.github.io/?https://github.com/liuqivandy/MEcell/blob/main/Tutorial/MERFISH.html)
- [Analyze Xenium mouse pup data](https://htmlpreview.github.io/?https://github.com/liuqivandy/MEcell/blob/main/Tutorial/Xenium.html)



<a name="Citation"/>

# Citation
Qi Liu et al. Microenvironment-Aware Spatial Modeling for Accurate Inference of Cell Identity
