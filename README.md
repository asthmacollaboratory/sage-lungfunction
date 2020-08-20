# sage-lungfunction

Analysis and plotting code for analyses of lung function in the Study of African Americans, Asthma, Genes, and Enviornments (SAGE).
**Nota Bene**: this code is provided as a template, an example, and a receipt, but it is _not_ designed as a batteries-included workflow. 
Interested users will need to configure several file paths to make this code useable.

## Prerequisites

* R version >= 3.5
* [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0) v3.0
* Python 2.7 (no Python 3 code used here)

The R/Bioconductor packages used here are:

* GENESIS
* GWASTools
* RColorBrewer
* SNPRelate
* arsenal
* coda
* corrplot
* data.table
* doParallel
* gdsfmt
* ggcorrplot
* ggrepel
* grid
* gridExtra
* methods
* optparse
* tidyverse (ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats) 

## Code 
Analyses are organized in numerical order under the `./src` directory.

## Results
Results are stored in a Dryad repository: https://doi.org/10.5061/dryad.w3r2280nw

## Preprint
If you wish to refer to these project, then please cite our preprint on [_bioRxiv_](https://doi.org/10.1101/2020.05.01.045468):

Goddard PC, Keys KL, White MJ, et al. (2020) ''Integrative genomic analysis in African American children with asthma finds 3 novel loci associated with lung function''. https://doi.org/10.1101/2020.05.01.045468
