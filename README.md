[![Travis Build Status](https://app.travis-ci.com/gabrielgesteira/QTLpoly.svg?branch=main)](https://app.travis-ci.com/gabrielgesteira/QTLpoly)
 [![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/qtlpoly)](https://cran.r-project.org/package=qtlpoly)
[![R-universe PolyVerse Status Badge](https://polyploids.r-universe.dev/badges/qtlpoly)](https://polyploids.r-universe.dev/badges/qtlpoly)
[![CRAN_monthly_downloads](https://cranlogs.r-pkg.org/badges/qtlpoly)](https://cranlogs.r-pkg.org/badges/qtlpoly?color=orange)



# QTLpoly <img src="https://raw.githubusercontent.com/gabrielgesteira/QTLpoly/main/hex.png" align="right" width="200" />

The R package `qtlpoly` (v. 0.2.4) is an under development software to map quantitative trait loci (QTL) in full-sib families of outcrossing autopolyploid species based on a random-effect multiple QTL model (Pereira et al. 2020). 

In order to do so, you will need a genetic map from which conditional probabilities of putative QTL can be computed. We recommend [`mappoly`](https://cran.r-project.org/package=mappoly), a hidden Markov model-based R package to construct genetic maps in autopolyploids (Mollinari and Garcia 2019).

Variance components associated with putative QTL are tested using score statistics (Qu et al. 2013), and final models are fitted using residual maximum likelihood (REML, adapted from the R package `sommer`). Plots for visualizing the results are based on `ggplot2` (v. 3.1 or higher) (Wickham 2016). 

# Install `qtlpoly` package

## From CRAN (stable version)

`qtlpoly`'s stable version is available on the Comprehensive R Archive Network (CRAN). You can install it by running:

```r
install.packages("qtlpoly")
```

## From Github (development version)

`qtlpoly`'s development version is available here on Github. You can install it using the R package `devtools`:

```r
install.packages("devtools")
devtools::install_github("gabrielgesteira/qtlpoly") 
```

If you are using Windows, you may need to first install the latest recommended version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

# Documents 

Tutorials as well as simulated and real data set analyses will be listed here in order to help users to get familiar with the software and allow them to perform their own analyses:

<!-- 1. [Tutorial on Multiple QTL Mapping in Autopolyploids with QTLpoly](https://guilherme-pereira.github.io/QTLpoly/1-tutorial) -->
1. [Tutorial on Multiple QTL Mapping in Autopolyploids with QTLpoly](https://gabrielgesteira.github.io/QTLpoly/docs/1-tutorial)
2. [Tools for Polyploids](https://www.polyploids.org/workshop/2021/january/info) training section: [Multiple QTL Mapping in an Autotetraploid F<sub>1</sub> population with QTLpoly](https://guilherme-pereira.github.io/QTLpoly/2-tetraploid_example.html)


# Related software

* [Polyverse](https://polyploids.r-universe.dev/ui#builds) - the polyploid R universe (a Lindsay Clark's initiative)

* Variant Calling
  *  [GBSapp: An automated pipeline for variant calling and filtering.](https://github.com/bodeolukolu/GBSapp)

* Simulations
  * [PedigreeSim: Simulation of genetic marker data in diploid and polyploid pedigreed populations.](https://www.wur.nl/en/show/Software-PedigreeSim.htm)

* Genotype calling
  * [ClusterCall: Automated tetraploid genotype calling by hierarchical clustering](https://potatobreeding.cals.wisc.edu/software/)
  * [fitPoly: Genotype Calling for Bi-Allelic Marker Assays](https://CRAN.R-project.org/package=fitPoly)
  * [polyRAD: Genotype Calling with Uncertainty from Sequencing Data in Polyploids and Diploids](https://CRAN.R-project.org/package=polyRAD)
  * [SuperMASSA: Graphical Bayesian inference tool for genotyping polyploids](https://bitbucket.org/orserang/supermassa)
  * [updog: Flexible Genotyping for Polyploids](https://CRAN.R-project.org/package=updog)
  * [VCF2SM: Python script that integrates VCF files and SuperMASSA](https://github.com/guilherme-pereira/vcf2sm)
 
* Genetic mapping in polyploids
  * [MAPpoly: Genetic maps in complex autopolyploids with even ploidy levels](https://CRAN.R-project.org/package=mappoly)
  * [MDSMap: High Density Genetic Linkage Mapping using Multidimensional Scaling](https://CRAN.R-project.org/package=MDSMap)
  * [polymapR: Linkage Analysis in Outcrossing Polyploids](https://CRAN.R-project.org/package=polymapR)
  * [TetraploidSNPMap: Linkage maps and mapping QTLs for autotetraploid species, using SNP dosage data.](https://www.bioss.ac.uk/knowledge/tetraploidmap/)
  
  
* Haplotype reconstruction
  * [TetraOrigin:haplotype reconstruction in a full-sib tetraploid family](https://github.com/chaozhi/TetraOrigin)
  * [PolyOriginR:haplotype reconstruction in polyploid multiparental populations](https://github.com/chaozhi/PolyOriginR)

* QTL mapping
  * [diaQTL: QTL analysis of diploid and autotetraploid diallel populations](https://github.com/jendelman/diaQTL)
  * [polyqtlR: QTL analysis and exploration of meiotic patterns in autopolyploid bi-parental F1 populations.](https://CRAN.R-project.org/package=polyqtlR)
  
* Data visualization
  * [VIEWpoly: integrate, visualize and explore results from genetic analysis, together with genomic information for autopolyploids](https://CRAN.R-project.org/package=viewpoly)
  

# Acknowledgments

This package has been developed as part of the [Genomic Tools for Sweetpotato Improvement](https://sweetpotatogenomics.cals.ncsu.edu/) (GT4SP) and [SweetGAINS](https://cgspace.cgiar.org/handle/10568/106838) projects, both funded by [Bill \& Melinda Gates Foundation](https://www.gatesfoundation.org/).

# References

Covarrubias-Pazaran G. 2016. “Genome-assisted prediction of quantitative traits using the R package sommer.” PLoS ONE 11 (6): 1-15. [doi:10.1371/journal.pone.0156744](https://doi.org/10.1371/journal.pone.0156744).

Mollinari M, Garcia AAF. 2019. “Linkage analysis and haplotype phasing in experimental autopolyploid populations with high ploidy level using hidden Markov models.” G3: Genes, Genomes, Genetics 9 (10): 3297-3314. [doi:10.1534/g3.119.400378](https://doi.org/10.1534/g3.119.400378).

Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB. 2020. “Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population.” Genetics 215 (3): 579-595. [doi:10.1534/genetics.120.303080](https://doi.org/10.1534/genetics.120.303080).

Qu L, Guennel T, Marshall SL. 2013. “Linear score tests for variance components in linear mixed models and applications to genetic association studies.” Biometrics 69 (4): 883-892.

Wickham H. 2016. “ggplot2: Elegant Graphics for Data Analysis.” Springer. [doi:10.1007/978-0-387-98141-3](https://www.springer.com/gp/book/9780387981413).
