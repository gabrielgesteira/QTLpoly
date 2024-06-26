% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{hexafake}
\alias{hexafake}
\title{Simulated autohexaploid dataset.}
\format{
An object of class \code{mappoly.data} which contains a
    list with the following components:
\describe{
    \item{plody}{ploidy level = 6}
    \item{n.ind}{number individuals = 300}
    \item{n.mrk}{total number of markers = 1500}
    \item{ind.names}{the names of the individuals}
    \item{mrk.names}{the names of the markers}
    \item{dosage.p1}{a vector containing the dosage in
      parent P for all \code{n.mrk} markers}
    \item{dosage.p2}{a vector containing the dosage in
      parent Q for all \code{n.mrk} markers}
    \item{chrom}{a vector indicating the chromosome each marker
      belongs. Zero indicates that the marker was not assigned to any
      chromosome}
    \item{genome.pos}{Physical position of the markers into the
      sequence}
    \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
      for each individual (columns). Missing data are represented by 
      \code{ploidy_level + 1 = 7}}
    \item{n.phen}{There are no phenotypes in this simulation}
    \item{phen}{There are no phenotypes in this simulation}
    \item{chisq.pval}{vector containing p-values for all markers associated to 
                      the chi-square test for the expected segregation patterns 
                      under Mendelian segregation}
}
}
\usage{
hexafake
}
\description{
A dataset of a hypothetical autohexaploid full-sib population 
containing three homology groups
}
\examples{
library(mappoly)
plot(hexafake)
}
\references{
Mollinari M, Garcia AAF (2019) Linkage analysis and haplotype phasing in experimental autopolyploid populations with high ploidy level using hidden Markov models, \emph{G3: Genes|Genomes|Genetics} 9 (10): 3297-3314. \doi{10.1534/g3.119.400378}
}
\author{
Marcelo Mollinari, \email{mmollin@ncsu.edu}
}
\keyword{datasets}
