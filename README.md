# Introduction
  
Gene regulatory networks model the underlying gene regulation hierarchies that drive gene expression and cell states. The main function of the `epiregulon` package is to construct gene regulatory networks and infer transcription factor (TF) activity in single cells by integration of scATAC-seq and scRNA-seq data and incorporation of public bulk TF ChIP-seq data.

`epiregulon.archr` is a package maintained as a continuation of the former version of `epiregulon` (before November 2023) with dependence on and compatibility with `ArchR` package. Its development parallels that of `epiregulon` but it also keeps possibility to use `ArchR` project on top o which extended versions of the `calculateP2G`, `addTFMotifInfo` and `addMotifScore` functions perform operations leading to the construction of the gene regulatory network model.

For full documentation, please refer to the epiregulon [book](https://xiaosaiyao.github.io/epiregulon.book/).

# Installation

```
# install devtools
if(!require(devtools)) install.packages("devtools")

# install epiregulon.archr
devtools::install_github(repo='xiaosaiyao/epiregulon.archr')
```

Example data included in the tutorial are available from [scMultiome](https://bioconductor.org/packages/release/data/experiment/html/scMultiome.html) 

```
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scMultiome")
```
# Software Requirements

Users should have R version 4.4 or higher.

# Reference
Tomasz Włodarczyk, Aaron Lun, Diana Wu, Shreya Menon, Shushan Toneyan, Kerstin Seidel, Liang Wang, 
Jenille Tan, Shang-Yang Chen, Timothy Keyes, Aleksander Chlebowski, Yu Guo, Ciara Metcalfe, Marc Hafner, 
Christian W. Siebel, M. Ryan Corces, Robert Yauch, Shiqi Xie, Xiaosai Yao. 2023. "Inference of single-cell transcription factor activity to dissect mechanisms of lineage plasticity and drug response" bioRxiv 2023.11.27.568955; doi: [https://doi.org/10.1101/2023.11.27.568955](https://www.biorxiv.org/content/10.1101/2023.11.27.568955v1)

Contact: [Xiaosai Yao](mailto:yao.xiaosai@gene.com), Genentech Inc.
