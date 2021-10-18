[![Build Status](https://app.travis-ci.com/gabrielgesteira/QTLpoly.svg?branch=main)](https://app.travis-ci.com/gabrielgesteira/QTLpoly)
 [![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# QTLpoly <img src="hex.png" align="right" width="200" />

The R package `qtlpoly` (v. 0.2.2) is an under development software to map quantitative trait loci (QTL) in full-sib families of outcrossing autopolyploid species based on a random-effect multiple QTL model (Pereira et al. 2019). 

In order to do so, you will need a genetic map from which conditional probabilities of putative QTL can be computed. We recommend [`mappoly`](https://github.com/mmollina/MAPpoly), a hidden Markov model-based R package to construct genetic maps in autopolyploids (Mollinari and Garcia 2019).

Variance components associated with putative QTL are tested using score statistics from the R package `varComp` (v. 0.2-0) (Qu et al. 2013). Final models are fitted using residual maximum likelihood (REML) from the R package `sommer` (v. 4.0 or higher) (Covarrubias-Pazaran 2016). Plots for visualizing the results are based on `ggplot2` (v. 3.1 or higher) (Wickham 2016). 

## Install `qtlpoly` package

`qtlpoly` package is available in its development version here on [GitHub](https://github.com/guilherme-pereira/qtlpoly). You can install all needed packages within R using the R package `devtools`:

```r
install.packages("remotes")
remotes::install_github("gabrielgesteira/qtlpoly") 
```

## Documents 

Tutorials as well as simulated and real data set analyses will be listed here opportunately in order to help users to get familiar with the software and allow them to perform their own analyses:

1. [Tutorial on Multiple QTL Mapping in Autopolyploids with QTLpoly](https://guilherme-pereira.github.io/QTLpoly/1-tutorial)
2. [Tools for Polyploids](https://www.polyploids.org/workshop/2021/january/info) training section: [Multiple QTL Mapping in an Autotetraploid F<sub>1</sub> population with QTLpoly](https://guilherme-pereira.github.io/QTLpoly/2-tetraploid_example.html)

## Acknowledgments

This package has been developed as part of the [Genomic Tools for Sweetpotato Improvement](https://sweetpotatogenomics.cals.ncsu.edu/) (GT4SP) project, funded by [Bill \& Melinda Gates Foundation](https://www.gatesfoundation.org/).

## References

Covarrubias-Pazaran G. 2016. “Genome-assisted prediction of quantitative traits using the R package sommer.” PLoS ONE 11 (6): 1-15. [doi:10.1371/journal.pone.0156744](https://doi.org/10.1371/journal.pone.0156744).

Mollinari M, Garcia AAF. 2019. “Linkage analysis and haplotype phasing in experimental autopolyploid populations with high ploidy level using hidden Markov models.” G3: Genes, Genomes, Genetics 9 (10): 3297-3314. [doi:10.1534/g3.119.400378](https://doi.org/10.1534/g3.119.400378).

Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB. 2020. “Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population.” Genetics 215 (3): 579-595. [doi:10.1534/genetics.120.303080](https://doi.org/10.1534/genetics.120.303080).

Qu L, Guennel T, Marshall SL. 2013. “Linear score tests for variance components in linear mixed models and applications to genetic association studies.” Biometrics 69 (4): 883-892. [doi:10.1111/biom.12095](https://doi.org/10.1111/biom.12095).

Wickham H. 2016. “ggplot2: Elegant Graphics for Data Analysis.” Springer. [doi:10.1007/978-0-387-98141-3](https://www.springer.com/gp/book/9780387981413).
