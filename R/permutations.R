#' @title Fixed-effect interval mapping (FEIM) model permutations
#'
#' @description Stores maximum LOD scores for a number of permutations of given phenotypes.
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param offset.data a subset of the data object to be used in permutation calculations.
#'
#' @param pheno.col a numeric vector with the phenotype columns to be analyzed; if \code{NULL} (default), all phenotypes from \code{'data'} will be included.
#'
#' @param n.sim a number of simulations, e.g. 1000 (default).
#'
#' @param n.clusters a number of parallel processes to spawn.
#'
#' @param seed an integer for the \code{set.seed()} function; if \code{NULL}, no reproducible seeds are set.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @param x an object of class \code{qtlpoly.perm} to be printed or plotted.
#'
#' @param probs a vector of probability values in [0, 1] representing the quantiles, e.g. c(0.90, 0.95) for the 90\% and 95\% quantiles.
#'
#' @param ... currently ignored
#'
#' @return An object of class \code{qtlpoly.perm} which contains a list of \code{results} for each trait with the maximum LOD score per permutation.
#'
#' @return LOD score thresholds for given quantiles for each trait.
#'
#' @return A \pkg{ggplot2} histogram with the distribution of ordered maximum LOD scores and thresholds for given quantiles for each trait.
#'
#' @seealso \code{\link[qtlpoly]{feim}}
#'
#' @examples
#'   \donttest{
#'   # Estimate conditional probabilities using mappoly package
#'   library(mappoly)
#'   library(qtlpoly)
#'   genoprob4x = lapply(maps4x[c(5)], calc_genoprob)
#'   data = read_data(ploidy = 4, geno.prob = genoprob4x, pheno = pheno4x, step = 1)
#'
#'   # Perform permutations
#'   perm = permutations(data = data, pheno.col = 1, n.sim = 10, n.clusters = 1)
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Churchill GA, Doerge RW (1994) Empirical threshold values for quantitative trait mapping, \emph{Genetics} 138: 963-971.
#'
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \doi{10.1534/genetics.120.303080}.
#'
#' @export permutations
#' @import parallel

permutations <- function(data, offset.data = NULL, pheno.col = NULL, n.sim = 1000, probs = c(0.90, 0.95), n.clusters = NULL, seed = 123, verbose = TRUE) {
  if(is.null(n.clusters)) n.clusters <- 1
  if(verbose) cat("INFO: Using", n.clusters, "CPUs for calculation\n\n")
  cl <- makeCluster(n.clusters)
  clusterSetRNGStream(cl, seed) # to make this reproducible
  #clusterEvalQ(cl, permuta)
  if(is.null(pheno.col)) pheno.col <- 1:dim(data$pheno)[2]
  results <- vector("list", length(pheno.col))
  names(results) <- colnames(data$pheno)[pheno.col]
  x <- array(1:n.sim, dim=c(1,n.sim))
  for(p in 1:length(results)) {
    start <- proc.time()
    if(verbose) cat("Permutations for trait", pheno.col[p], sQuote(colnames(data$pheno)[pheno.col[p]]), "\n")
    lod <- NULL
    lod <- parCapply(cl, x, FUN = permuta, data = data, offset.data = offset.data, pheno.col = pheno.col, p = p)
    results[[p]] <- lod
    if(verbose) cat("  ", last(probs)*100, "% LOD threshold = ", round(quantile(sort(lod), last(probs)), digits = 2), "\n", sep = "")
    end <- proc.time()
    if(verbose) cat("  Calculation took", round((end - start)[3], digits = 2), "seconds\n\n")
  }
  
  stopCluster(cl)
  
  sig.lod <- vector("list", length(probs))
  names(sig.lod) <- as.character(probs)
  for(i in 1:length(probs)) {
    sig.lod[[i]] <- unlist(lapply(results, FUN = function(x) quantile(sort(x), probs[i])))
  }
  
  structure(list(data=deparse(substitute(data)),
                 pheno.col=pheno.col,
                 n.sim=n.sim,
                 probs=probs,
                 seed=seed,
                 results=results,
                 sig.lod=sig.lod
  ),
  class=c("qtlpoly.perm"))
}

#' @rdname permutations
#' @export

print.qtlpoly.perm <- function(x, pheno.col=NULL, probs=c(0.90, 0.95), ...) {
  if(any(class(x) == "qtlpoly.perm")) cat("This is an object of class 'qtlpoly.perm'\n")
  if(is.null(pheno.col)) {
    pheno.col <- 1:length(x$results)
  } else {
    pheno.col <- which(x$pheno.col %in% pheno.col)
  }
  for(p in pheno.col) {
    cat("\n* Trait", p, sQuote(names(x$results)[[p]]), "\n")
    print(round(quantile(sort(x$results[[p]]), probs), digits=2))
  }
}

#' @rdname permutations
#' @import ggplot2
#' @export

plot.qtlpoly.perm <- function(x, pheno.col=NULL, probs=c(0.90, 0.95), ...) {
  LOD = NULL
  if(is.null(pheno.col)) {
    pheno.col <- 1:length(x$results)
  } else {
    pheno.col <- which(x$pheno.col %in% pheno.col)
  }
  for(p in pheno.col) {
    if(!is.null(x$results[[p]])) {
      data <- data.frame(LOD=sort(x$results[[p]]))
      plot <- ggplot(data, aes(LOD)) +
        geom_histogram(bins = 15) +
        geom_vline(xintercept = quantile(data$LOD, probs), linetype="dashed") +
        annotate("text", x=quantile(data$LOD, probs), y=length(data$LOD)/6, angle=90, vjust=-0.1, hjust=0,
                 label=paste(names(quantile(data$LOD, probs)), round(quantile(data$LOD, probs), digits = 2), sep=" = ")) +
        labs(title=names(x$results)[p], y="Count", x="Maximum LOD scores") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5), title=element_text(face="bold"))
      print(plot)
    }
  }
}

#' @keywords internal
permuta <- function(x, data = data, offset.data = offset.data, pheno.col = pheno.col, p = p) {
  ind <- which(dimnames(data$pheno)[[1]] %in% dimnames(data$X)[[1]])
  LRT <- numeric(data$nmrk)
  Y <- sample(data$pheno[ind,pheno.col[p]])
  if(is.null(offset.data)) {
    offset <- NULL
  } else {
    offset <- offset.data[names(Y),pheno.col[p]]
  }
  for(m in 1:data$nmrk) {
    full.mod <- lm(Y ~ 1 + data$X[ind,-c(1,(data$ploidy+1)),m], offset = offset[names(Y)])
    null.mod <- lm(Y ~ 1, offset = offset[names(Y)])
    LRT[m] <- -2 * (logLik(null.mod) - logLik(full.mod))
  }
  LOD <- LRT/(2*log(10))
  return(max(LOD))
}
