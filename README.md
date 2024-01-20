# Co-Expression Algorithm from Microarray Data

## Overview

This repository contains the R code implementation of an algorithm designed to extract co-expression data from microarray data. 

## Features

- **Co-Expression Analysis:** Utilizes microarray data to identify genes that exhibit similar expression patterns.
- **Scalable:** The algorithm is designed to handle large datasets efficiently.
- **Customizable Parameters:** Users can adjust the following parameters to customize the co-expression analysis:

    1. **Termination Patience:** The number of iterations to wait for the FDR values to stabilize. Stabilization is defined as the difference not being greater than a certain threshold value.

    2. **Stabilization Threshold:** The threshold value used to determine stabilization. If the difference between n iterations is not greater than this threshold, the algorithm considers the values stabilized.

    3. **FDR (False Discovery Rate) Values:** Users can set the FDR values to remove poorly coexpressed modules. Adjusting this parameter allows for a more flexible threshold for module co-expression.

## Requirements

- R (>= 4.0.3)

## Installation

Clone the repository to your local machine:

```bash
git clone https://github.com/zhejilab/MetNet
```


## Usage

- **Download interested microarray data:** For example, Cancer Cell Line Encyclopedia (CCLE) database (https://portals.broadinstitute.org/ccle/) provides the RNA characterization of more than 1000 cancer cell lines.

- **Run the Algorithm:** Rscript coExp.R -p 25 -s 0.01 -f 0.01 -o ./

