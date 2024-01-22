# Co-Expression Algorithm from Microarray Data

## Overview

This repository contains the code implementation of algriothms and figures in "Comparative analyses of gene regulatory networks mediating cancer metastatic potentials across lineage types. Wang, et al, 2024". 

- **Part I (coexpModule):** R code implementation of an algorithm designed to extract co-expression data from microarray data. 
- **Part II (random_forest_regressor):** Python code implementation of the regularized ramdon forest modeling using the subGO expression levels to predict metastatic potentials. (Fig 3C-G)
- **Part III (met.network):** The raw node and edge information of network 1-4 from the protein-protein interaction analyses.

## coexpModule

- **Co-Expression Analysis:** Utilizes microarray data to identify genes that exhibit similar expression patterns.
- **Scalable:** The algorithm is designed to handle large datasets efficiently.
- **Customizable Parameters:** Users can adjust the following parameters to customize the co-expression analysis:

    1. **Termination Patience:** The number of iterations to wait for the FDR values to stabilize. Stabilization is defined as the difference not being greater than a certain threshold value.

    2. **Stabilization Threshold:** The threshold value used to determine stabilization. If the difference between n iterations is not greater than this threshold, the algorithm considers the values stabilized.

    3. **FDR (False Discovery Rate) Values:** Users can set the FDR values to remove poorly coexpressed modules. Adjusting this parameter allows for a more flexible threshold for module co-expression.

### Requirements

- R (>= 4.0.3)
- python 

### Usage

- **Run the Algorithm:** Rscript coExp.R -d ccle.csv -p 25 -s 0.01 -g CDH3,PTPRU,MPZL2,PAK4,BAIAP2,AGR2,SPINT2,HPSE,RAB10,PKP3,CLDN4,CLDN3,CLDN7,PLEKHA7,CXADR,DLG3,DSC2,DSG2,DSP,EFNA1,EFNB2,EMP2,EPHA1,EPHB3,EPHB4,ERBB2,ERBB3,FOLR1,CD2AP,TRIM29,FUT1,NPNT,LYPD3,TNFRSF21,ANK3,RHOD,FOXA1,HES1,IL18,ITGA6,ITGB4,JAG2,JUP,LAMA5,LAMB3,LAMC2,BCAM,TACSTD2,EPCAM,CD46,NRARP,ASS1,MUC1,CEACAM6,ATP1B1,MINK1,F11R,PKP2,PLXNB1,EPB41L4B,VSIG10,FERMT1,PRKCZ,BAIAP2L1,CCL28,PSEN1,IGSF9,EPB41L5,PTPRF,CEACAM1,PERP,SOX9,STX3,EZR,WNT7B,DDR1,GRHL2,ANXA9,PPFIA1,PKP4,ADAM15,CLDN1,CD9,KLF4,TJP2,SOX13,SLK,DGCR2,CDH1 -f 0.01 -o ./

- **-d --data:** Input microsarry data
- **-g --gene:** Genes of interested
- **-p --patience:** Number of iterations to wait for the FDR values to stabilize. Stabilization is defined as the difference not being greater than a certain threshold value (default 25).
- **-s --stabilization_threshold:** The threshold value used to determine stabilization. If the difference between n iterations is not greater than this threshold, the algorithm considers the values stabilized (default 0.01).
- **-f --fdr:** FDR values to remove poorly coexpressed modules. (default 0.01).
- **-o --output:** Output directory

## random_forest_regressor

python ../random_forest.py x.csv y.csv ./example_output

- **Arg 1** Attibutes table
- **Arg 2** Targets table
- **Arg 3** Output directory

