#' Simulated autohexaploid dataset.
#'
#' A dataset of a hypothetical autohexaploid full-sib population 
#' containing three homology groups
#'
#' @docType data
#'
#' @format An object of class \code{mappoly.data} which contains a
#'     list with the following components:
#' \describe{
#'     \item{plody}{ploidy level = 6}
#'     \item{n.ind}{number individuals = 300}
#'     \item{n.mrk}{total number of markers = 1500}
#'     \item{ind.names}{the names of the individuals}
#'     \item{mrk.names}{the names of the markers}
#'     \item{dosage.p1}{a vector containing the dosage in
#'       parent P for all \code{n.mrk} markers}
#'     \item{dosage.p2}{a vector containing the dosage in
#'       parent Q for all \code{n.mrk} markers}
#'     \item{chrom}{a vector indicating the chromosome each marker
#'       belongs. Zero indicates that the marker was not assigned to any
#'       chromosome}
#'     \item{genome.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1 = 7}}
#'     \item{n.phen}{There are no phenotypes in this simulation}
#'     \item{phen}{There are no phenotypes in this simulation}
#'     \item{chisq.pval}{vector containing p-values for all markers associated to 
#'                       the chi-square test for the expected segregation patterns 
#'                       under Mendelian segregation}
#' }
#'
#' @keywords datasets
#'
#' @author Marcelo Mollinari, \email{mmollin@@ncsu.edu}
#'
#' @references 
#'     Mollinari M, Garcia AAF (2019) Linkage analysis and haplotype phasing in experimental autopolyploid populations with high ploidy level using hidden Markov models, \emph{G3: Genes|Genomes|Genetics} 9 (10): 3297-3314. \doi{10.1534/g3.119.400378}
#'
#' @examples
#' library(mappoly)
#' plot(hexafake)
"hexafake"

#' Simulated autohexaploid map
#'
#' A simulated map containing three homology groups of a hypotetical cross between two autohexaploid individuals.
#'
#' @docType data
#'
#' @format An object of class \code{"mappoly.map"} from the package \pkg{mappoly}, which is a list of three linkage groups (LGs):
#'
#' \describe{
#'   \item{LG 1}{538 markers distributed along 112.2 cM}
#'   \item{LG 2}{329 markers distributed along 54.6 cM}
#'   \item{LG 3}{443 markers distributed along 98.2 cM}
#' }
#'
#' @keywords datasets
#'
#' @seealso \code{\link[mappoly]{hexafake}}, \code{\link[qtlpoly]{pheno6x}}
#'
#' @author Marcelo Mollinari, \email{mmollin@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \doi{10.1534/genetics.120.303080}.
#'
#'     Mollinari M, Garcia AAF (2019) Linkage analysis and haplotype phasing in experimental autopolyploid populations with high ploidy level using hidden Markov models, \emph{G3: Genes|Genomes|Genetics} 9 (10): 3297-3314. \doi{10.1534/g3.119.400378}
#'
#' @examples
#' library(mappoly)
#' plot_map_list(maps6x)
"maps6x"

#' Simulated phenotypes
#'
#' A simulated data set of phenotypes for a hipotetical autohexaploid species map.
#'
#' @docType data
#'
#' @format A data frame of phenotypes with 300 named individuals in rows and three named phenotypes in columns, which are:
#'
#' \describe{
#'   \item{T32}{3 QTLs, with heritabilities of 0.20 (LG 1 at 32.03 cM), 0.15 (LG 1 at 95.02 cM) and 0.30 (LG 2 at 40.01 cM).}
#'   \item{T17}{1 QTL, with heritability of 0.15 (LG 3 at 34.51 cM).}
#'   \item{T45}{no QTLs.}
#' }
#'
#' @keywords datasets
#'
#' @seealso \code{\link[qtlpoly]{simulate_qtl}}, \code{\link[qtlpoly]{pheno4x}}
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \doi{10.1534/genetics.120.303080}.
#'
#' @examples
#' head(pheno6x)
"pheno6x"

#' Autotetraploid potato dataset
#'
#' A dataset of the B2721 population which derived from a cross between 
#' two tetraploid potato varieties: Atlantic × B1829-5. 
#'
#' @docType data
#'
#' @format An object of class \code{mappoly.data} from the package \pkg{mappoly}.
#' 
#' @keywords datasets
#'
#' @author Marcelo Mollinari, \email{mmollin@@ncsu.edu}
#'
#' @references 
#'     Mollinari M, Garcia AAF (2019) Linkage analysis and haplotype phasing in experimental autopolyploid populations with high ploidy level using hidden Markov models, \emph{G3: Genes|Genomes|Genetics} 9 (10): 3297-3314. \doi{10.1534/g3.119.400378}
#'
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \doi{10.1534/genetics.120.303080}.
#'     
#'     Pereira GS, Mollinari M, Schumann MJ, Clough ME, Zeng ZB, Yencho C (2021) The recombination landscape and multiple QTL mapping in a \emph{Solanum tuberosum} cv. ‘Atlantic’-derived F_1 population. Heredity. \doi{10.1038/s41437-021-00416-x}.
#' 
#' @examples
#' library(mappoly)
#' print(B2721)
"B2721"

#' Autotetraploid potato map
#'
#' A real autotetraploid potato map containing 12 homology groups from a tetraploid potato full-sib family (Atlantic x B1829-5).
#'
#' @docType data
#'
#' @format An object of class \code{"mappoly.map"} from the package \pkg{mappoly}, which is a list of 12 linkage groups (LGs)
#'
#' @keywords datasets
#'
#' @seealso \code{\link[mappoly]{hexafake}}, \code{\link[qtlpoly]{pheno6x}}
#'
#' @author Marcelo Mollinari, \email{mmollin@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \doi{10.1534/genetics.120.303080}.
#'
#'     Mollinari M, Garcia AAF (2019) Linkage analysis and haplotype phasing in experimental autopolyploid populations with high ploidy level using hidden Markov models, \emph{G3: Genes|Genomes|Genetics} 9 (10): 3297-3314. \doi{10.1534/g3.119.400378}
#'
#'     Pereira GS, Mollinari M, Schumann MJ, Clough ME, Zeng ZB, Yencho C (2021) The recombination landscape and multiple QTL mapping in a \emph{Solanum tuberosum} cv. ‘Atlantic’-derived F_1 population. Heredity. \doi{10.1038/s41437-021-00416-x}.
#'
#' @examples
#' library(mappoly)
#' plot_map_list(maps4x)
"maps4x"

#' Autotetraploid potato phenotypes
#'
#' A subset of phenotypes from a tetraploid potato full-sib family (Atlantic x B1829-5).
#'
#' @docType data
#'
#' @format A data frame of phenotypes with 156 named individuals in rows and three named phenotypes in columns, which are:
#'
#' \describe{
#'   \item{FM07}{Foliage maturity evaluated in 2007.}
#'   \item{FM08}{Foliage maturity evaluated in 2008.}
#'   \item{FM14}{Foliage maturity evaluated in 2014.}
#' }
#'
#' @keywords datasets
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \doi{10.1534/genetics.120.303080}.
#'     
#'     Pereira GS, Mollinari M, Schumann MJ, Clough ME, Zeng ZB, Yencho C (2021) The recombination landscape and multiple QTL mapping in a \emph{Solanum tuberosum} cv. ‘Atlantic’-derived F_1 population. Heredity. \doi{10.1038/s41437-021-00416-x}.
#'
#' @examples
#' head(pheno4x)
"pheno4x"
