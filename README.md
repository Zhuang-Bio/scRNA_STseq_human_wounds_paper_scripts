# Spatiotemporal Dynamics of Skin Wound Healing

This repository contains the data processing and analysis scripts for our single-cell RNA-seq and Spatial transcript (ST-seq) of skin wound healing study **Liu et al. preprint (2024)**.


## Repository Contents

  - [01-scRNAseq Acute Wounds](https://github.com/Zhuang-Bio/scRNA_STseq_human_wounds_paper_scripts/tree/main/01-scRNAseq%20Acute%20Wounds) contains scripts and notebooks used for running CellRanger, Removing Doublets, Seurat, cell-cell crosstalk analysis, subclustering analysis of acute wounds, and Milo cell abundance analysis.
   - [02-STseq Acute Wounds](https://github.com/Zhuang-Bio/scRNA_STseq_human_wounds_paper_scripts/tree/main/02-STseq%20Acute%20Wounds) contains scripts used for running SpaceRanger, standard seurat analysis and plotting, and preparation for deconvolution.
   - [03-Deconvolution](https://github.com/Zhuang-Bio/scRNA_STseq_human_wounds_paper_scripts/tree/main/03-Deconvolution) contains scripts and notebooks used for running the deconvolution analysis of spatial data and bulk RNA-seq data based on the cell type signatures in scRNA-seq.
   - [04-Data Integration](https://github.com/Zhuang-Bio/scRNA_STseq_human_wounds_paper_scripts/tree/main/04-Data%20Integration) contains scripts and notebooks used for integrating scRNA-seq datasets from different sources. e.g. healthy adult skin scRNA-seq, diabetic foot ulcer scRNA-seq, venous ulcer scRNA-seq, and acute wound scRNA-seq.
   - [05-Cross species analysis](https://github.com/Zhuang-Bio/scRNA_STseq_human_wounds_paper_scripts/tree/main/05-Cross%20species%20analysis) contains scripts used for corss-species comparison of human and mouse wound healing.
   - [Rshiny scwoundatlas](https://github.com/Zhuang-Bio/scRNA_STseq_human_wounds_paper_scripts/tree/main/Rshiny%20scwoundatlas) contains scripts used for creating RShiny visualization tool in our website.
   - [custom functions](https://github.com/Zhuang-Bio/scRNA_STseq_human_wounds_paper_scripts/tree/main/custom%20functions) contains scripts used for common functions.
   - [scripts reproduce figures](https://github.com/Zhuang-Bio/scRNA_STseq_human_wounds_paper_scripts/tree/main/scripts%20reproduce%20figures) contains scripts used for reproducing figures in our manuscript.


## Data Available

### Web portal (Rshiny tool)
  Processed data can be accessed on our web portal at https://www.xulandenlab.com/tools [Not publicly available right now]

### Raw sequencing data
  Raw requencing data for this study will be deposited in GEO: https://www.ncbi.nlm.nih.gov/geo/. Please keep an eye on our manuscript.


## Contact and collaboration
Feel free to post an [issue](https://github.com/Zhuang-Bio/scRNA_STseq_human_wounds_paper_scripts/issues) in this repository or contact by email zhuang.liu[at]ki.se or ning.xu[at]ki.se if you have any questions.


