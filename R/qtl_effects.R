
#' QTL allele effect estimation
#'
#' Computes allele specific and allele combination (within-parent) heritable effects from multiple QTL models.
#'
#' @param ploidy a numeric value of ploidy level of the cross (currently, only 2, 4 or 6).
#'
#' @param fitted a fitted multiple QTL model of class \code{qtlpoly.fitted}.
#'
#' @param x an object of class \code{qtlpoly.effects} to be plotted.
#'
#' @param pheno.col a numeric vector with the phenotype column numbers to be plotted; if \code{NULL}, all phenotypes from \code{'fitted'} will be included.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @param p1 a character string with the first parent name, e.g. \code{"P1"} (default).
#'
#' @param p2 a character string with the second parent name, e.g. \code{"P2"} (default).
#'
#' @param ... currently ignored
#'
#' @return An object of class \code{qtlpoly.effects} which is a list of \code{results} for each containing the following components:
#'
#'     \item{pheno.col}{a phenotype column number.}
#'     \item{y.hat}{a vector with the predicted values.}
#'
#' @return A \pkg{ggplot2} barplot with parental allele and allele combination effects.
#'
#' @seealso \code{\link[qtlpoly]{read_data}}, \code{\link[qtlpoly]{remim}}, \code{\link[qtlpoly]{fit_model}}
#'
#' @examples
#'   \donttest{
#'   # Estimate conditional probabilities using mappoly package
#'   library(mappoly)
#'   library(qtlpoly)
#'   genoprob4x = lapply(maps4x[c(5)], calc_genoprob)
#'   data = read_data(ploidy = 4, geno.prob = genoprob4x, pheno = pheno4x, step = 1)
#'
#'   # Search for QTL
#'   remim.mod = remim(data = data, pheno.col = 1, w.size = 15, sig.fwd = 0.0011493379,
#' sig.bwd = 0.0002284465, d.sint = 1.5, n.clusters = 1)
#'
#'   # Fit model
#'   fitted.mod = fit_model(data, model=remim.mod, probs="joint", polygenes="none")
#'
#'   # Estimate effects
#'   est.effects = qtl_effects(ploidy = 4, fitted = fitted.mod, pheno.col = 1)
#'
#'   # Plot results
#'   plot(est.effects)
#'   }
#'   
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}, with modifications by Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \doi{10.1534/genetics.120.303080}.
#'     
#'     Kempthorne O (1955) The correlation between relatives in a simple autotetraploid population, \emph{Genetics} 40: 168-174.
#'
#' @export qtl_effects

qtl_effects <- function(ploidy = 6, fitted, pheno.col = NULL, verbose = TRUE) {
  
  if(is.null(pheno.col)) pheno.col <- fitted$pheno.col
  results <- vector("list", length(pheno.col))
  names(results) <- names(fitted$results)[which(fitted$pheno.col %in% pheno.col)]
  
  for(p in 1:length(results)) {
    
    if(!is.null(fitted$results[[pheno.col[p]]]$qtls)) {
      
      nqtl <- dim(fitted$results[[pheno.col[p]]]$qtls)[1]
      if(nqtl > 1) nqtl <- nqtl - 1
      effects <- vector("list", nqtl)
      qtl.mrk <- unlist(fitted$results[[pheno.col[p]]]$qtls[c(1:nqtl),"Nmrk"])
      
      if(verbose) {
        if(length(qtl.mrk) == 1) cat("There is ", length(qtl.mrk), " QTL in the model for trait ", pheno.col[p], " ", sQuote(names(fitted$results)[pheno.col[p]]), ". Computing effects for QTL ", sep="")
        if(length(qtl.mrk) >= 2) cat("There are ", length(qtl.mrk), " QTL in the model for trait ", pheno.col[p], " ", sQuote(names(fitted$results)[pheno.col[p]]), ". Computing effects for QTL ", sep="")
      }
      
      if(ploidy == 6) {
        
        for(q in 1:nqtl) {
          
          if(verbose) {
            if(q < nqtl) cat("...", qtl.mrk[q], "")
            if(q == nqtl) cat(paste0("... ", qtl.mrk[q]))
          }
          
          blups <- fitted$results[[pheno.col[p]]]$fitted$U[[q]]
          alleles = matrix(unlist(strsplit(fitted$results[[pheno.col[p]]]$fitted$alleles, '')), ncol=7, byrow=TRUE)[,-4]
          ## alleles <- matrix(unlist(strsplit(rownames(blups), '')), ncol=7, byrow=TRUE)[,-4]
          
          A <- t(combn(letters[1:12],1))
          D <- t(combn(letters[1:12],2))
          T <- t(combn(letters[1:12],3))
          F <- t(combn(letters[1:12],4))
          G <- t(combn(letters[1:12],5))
          S <- t(combn(letters[1:12],6))
          
          a <- vector("list", dim(A)[1])
          d <- vector("list", dim(D)[1])
          t <- vector("list", dim(T)[1])
          f <- vector("list", dim(F)[1])
          g <- vector("list", dim(G)[1])
          s <- vector("list", dim(S)[1])
          
          for(i in 1:dim(A)[1]) {
            a[[i]] <- which(alleles == as.character(A[i,1]), arr.ind = TRUE)[,1]
            a[[i]] <- mean(blups[Reduce(intersect, list(a[[i]]))])
          }
          names(a) <- as.character(A)
          for(i in 1:dim(D)[1]) {
            d[[i]] <- which(apply(alleles == as.character(D[i,1]) | alleles == as.character(D[i,2]), 1, sum) == 2)
            d[[i]] <- mean(blups[Reduce(intersect, list(d[[i]]))])
          }
          names(d) <- apply(D, 1, paste, collapse="")
          for(i in 1:dim(T)[1]) {
            t[[i]] <- which(apply(alleles == as.character(T[i,1]) | alleles == as.character(T[i,2]) | alleles == as.character(T[i,3]), 1, sum) == 3)
            t[[i]] <- mean(blups[Reduce(intersect, list(t[[i]]))])
          }
          names(t) <- apply(T, 1, paste, collapse="")
          for(i in 1:dim(F)[1]) {
            f[[i]] <- which(apply(alleles == as.character(F[i,1]) | alleles == as.character(F[i,2]) | alleles == as.character(F[i,3]) | alleles == as.character(F[i,4]), 1, sum) == 4)
            f[[i]] <- mean(blups[Reduce(intersect, list(f[[i]]))])
          }
          names(f) <- apply(F, 1, paste, collapse="")
          for(i in 1:dim(G)[1]) {
            g[[i]] <- which(apply(alleles == as.character(G[i,1]) | alleles == as.character(G[i,2]) | alleles == as.character(G[i,3]) | alleles == as.character(G[i,4]) | alleles == as.character(G[i,5]), 1, sum) == 5)
            g[[i]] <- mean(blups[Reduce(intersect, list(g[[i]]))])
          }
          names(g) <- apply(G, 1, paste, collapse="")
          for(i in 1:dim(S)[1]) {
            s[[i]] <- which(apply(alleles == as.character(S[i,1]) | alleles == as.character(S[i,2]) | alleles == as.character(S[i,3]) | alleles == as.character(S[i,4]) | alleles == as.character(S[i,5]) | alleles == as.character(S[i,6]), 1, sum) == 6)
            s[[i]] <- mean(blups[Reduce(intersect, list(s[[i]]))])
          }
          names(s) <- apply(S, 1, paste, collapse="")
          
          a <- a[!is.nan(unlist(a))]
          d <- d[!is.nan(unlist(d))]
          t <- t[!is.nan(unlist(t))]
          f <- f[!is.nan(unlist(f))]
          g <- g[!is.nan(unlist(g))]
          s <- s[!is.nan(unlist(s))]
          
          for(i in 1:length(d)) {
            d[[i]] <- d[[i]] -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(d), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          for(i in 1:length(t)) {
            t[[i]] <- t[[i]] -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(t), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(t), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          for(i in 1:length(f)) {
            f[[i]] <- f[[i]] -
              sum(unlist(t[which(lapply(lapply(strsplit(names(t), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 3) == TRUE)])) -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          for(i in 1:length(g)) {
            g[[i]] <- g[[i]] -
              sum(unlist(f[which(lapply(lapply(strsplit(names(f), split = ""), function(x) intersect(strsplit(names(g), split = "")[[i]], x)), function(x) length(x) == 4) == TRUE)])) -
              sum(unlist(t[which(lapply(lapply(strsplit(names(t), split = ""), function(x) intersect(strsplit(names(g), split = "")[[i]], x)), function(x) length(x) == 3) == TRUE)])) -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(g), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(g), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          for(i in 1:length(s)) {
            s[[i]] <- s[[i]] -
              sum(unlist(g[which(lapply(lapply(strsplit(names(g), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 5) == TRUE)])) -
              sum(unlist(f[which(lapply(lapply(strsplit(names(f), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 4) == TRUE)])) -
              sum(unlist(t[which(lapply(lapply(strsplit(names(t), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 3) == TRUE)])) -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          effects[[q]] <- list(unlist(a), unlist(d), unlist(t), unlist(f), unlist(g), unlist(s))
          
        }
      }
      
      if(ploidy == 4) {
        
        for(q in 1:nqtl) {
          
          if(verbose) {
            if(q < nqtl) cat("...", qtl.mrk[q], "")
            if(q == nqtl) cat(paste0("... ", qtl.mrk[q]))
          }
          
          blups <- fitted$results[[pheno.col[p]]]$fitted$U[[q]]
          alleles = matrix(unlist(strsplit(fitted$results[[pheno.col[p]]]$fitted$alleles, '')), ncol=5, byrow=TRUE)[,-3]
          ## alleles <- matrix(unlist(strsplit(rownames(blups), '')), ncol=5, byrow=TRUE)[,-3]
          
          A <- t(combn(letters[1:8],1))
          D <- t(combn(letters[1:8],2))
          T <- t(combn(letters[1:8],3))
          F <- t(combn(letters[1:8],4))
          
          a <- vector("list", dim(A)[1])
          d <- vector("list", dim(D)[1])
          t <- vector("list", dim(T)[1])
          f <- vector("list", dim(F)[1])
          
          for(i in 1:dim(A)[1]) {
            a[[i]] <- which(alleles == as.character(A[i,]), arr.ind = TRUE)[,1]
            a[[i]] <- mean(blups[Reduce(intersect, list(a[[i]]))])
          }
          names(a) <- as.character(A)
          for(i in 1:dim(D)[1]) {
            d[[i]] <- which(apply(alleles == as.character(D[i,1]) | alleles == as.character(D[i,2]), 1, sum) == 2)
            d[[i]] <- mean(blups[Reduce(intersect, list(d[[i]]))])
          }
          names(d) <- apply(D, 1, paste, collapse="")
          for(i in 1:dim(T)[1]) {
            t[[i]] <- which(apply(alleles == as.character(T[i,1]) | alleles == as.character(T[i,2]) | alleles == as.character(T[i,3]), 1, sum) == 3)
            t[[i]] <- mean(blups[Reduce(intersect, list(t[[i]]))])
          }
          names(t) <- apply(T, 1, paste, collapse="")
          for(i in 1:dim(F)[1]) {
            f[[i]] <- which(apply(alleles == as.character(F[i,1]) | alleles == as.character(F[i,2]) | alleles == as.character(F[i,3]) | alleles == as.character(F[i,4]), 1, sum) == 4)
            f[[i]] <- mean(blups[Reduce(intersect, list(f[[i]]))])
          }
          names(f) <- apply(F, 1, paste, collapse="")
          
          a <- a[!is.nan(unlist(a))]
          d <- d[!is.nan(unlist(d))]
          t <- t[!is.nan(unlist(t))]
          f <- f[!is.nan(unlist(f))]
          
          for(i in 1:length(d)) {
            d[[i]] <- d[[i]] -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(d), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          for(i in 1:length(t)) {
            t[[i]] <- t[[i]] -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(t), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(t), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          for(i in 1:length(f)) {
            f[[i]] <- f[[i]] -
              sum(unlist(t[which(lapply(lapply(strsplit(names(t), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 3) == TRUE)])) -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          effects[[q]] <- list(unlist(a), unlist(d), unlist(t), unlist(f))
          
        }
        
      }
      if(ploidy == 2) {
        
        for(q in 1:nqtl) {
          
          if(verbose) {
            if(q < nqtl) cat("...", qtl.mrk[q], "")
            if(q == nqtl) cat(paste0("... ", qtl.mrk[q]))
          }
          
          blups <- fitted$results[[pheno.col[p]]]$fitted$U[[q]]
          alleles = matrix(unlist(strsplit(fitted$results[[pheno.col[p]]]$fitted$alleles, '')), ncol=3, byrow=TRUE)[,-2]
          ## alleles <- matrix(unlist(strsplit(rownames(blups), '')), ncol=5, byrow=TRUE)[,-3]
          
          A <- t(combn(letters[1:4],1))
          D <- t(combn(letters[1:4],2))
          
          a <- vector("list", dim(A)[1])
          d <- vector("list", dim(D)[1])
          
          for(i in 1:dim(A)[1]) {
            a[[i]] <- which(alleles == as.character(A[i,]), arr.ind = TRUE)[,1]
            a[[i]] <- mean(blups[Reduce(intersect, list(a[[i]]))])
          }
          names(a) <- as.character(A)
          for(i in 1:dim(D)[1]) {
            d[[i]] <- which(apply(alleles == as.character(D[i,1]) | alleles == as.character(D[i,2]), 1, sum) == 2)
            d[[i]] <- mean(blups[Reduce(intersect, list(d[[i]]))])
          }
          names(d) <- apply(D, 1, paste, collapse="")
          
          a <- a[!is.nan(unlist(a))]
          d <- d[!is.nan(unlist(d))]
          
          for(i in 1:length(d)) {
            d[[i]] <- d[[i]] -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(d), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          effects[[q]] <- list(unlist(a), unlist(d))
          
        }
        
      }
      
      if(verbose) cat(". Done! \n\n", sep="")
      
    } else {
      
      if(verbose) cat("There are no QTL in the model for trait ", pheno.col[p], " ", sQuote(names(fitted$results)[pheno.col[p]]), ". Skipping! \n\n", sep="")
      effects <- NULL
    }
    
    results[[p]] <- list(
      pheno.col=fitted$results[[pheno.col[p]]]$pheno.col,
      effects=effects)
    
  }
  
  structure(list(fitted=deparse(substitute(fitted)),
                 ploidy=ploidy,
                 pheno.col=fitted$pheno.col,
                 results=results),
            class="qtlpoly.effects")
  
}

#' @rdname qtl_effects
#' @import ggplot2
#' @export

plot.qtlpoly.effects <- function(x, pheno.col = NULL, p1 = "P1", p2 = "P2", ...) {
  Alleles = Estimates = NULL
  if(is.null(pheno.col)) {
    pheno.col <- 1:length(x$results)
  } else {
    pheno.col <- which(x$pheno.col %in% pheno.col)
  }
  res = list()
  for(p in pheno.col) {
    nqtl <- length(x$results[[p]]$effects)
    if(nqtl > 0) {
      for(q in 1:nqtl) {
        if(x$ploidy == 2) {
          data <- unlist(x$results[[p]]$effects[[q]])[1:8]
          data <- data.frame(Estimates=as.numeric(data), Alleles=names(data), Parent=c(rep(p1,2),rep(p2,2),rep(p1,2),rep(p2,2)), Effects=c(rep("Additive",4),rep("Digenic",4)))
        }
        if(x$ploidy == 4) {
          data <- unlist(x$results[[p]]$effects[[q]])[1:36]
          data <- data.frame(Estimates=as.numeric(data), Alleles=names(data), Parent=c(rep(p1,4),rep(p2,4),rep(p1,14),rep(p2,14)), Effects=c(rep("Additive",8),rep("Digenic",28)))
          data <- data[-c(12:15,18:21,23:30),]
        }
        if(x$ploidy == 6) {
          data <- unlist(x$results[[p]]$effects[[q]])[-c(18:23,28:33,37:42,45:50,52:63,83:88,92:97,100:105,107:133,137:142,145:150,152:178,181:186,188:214,216:278,299:1763)]
          data <- data.frame(Estimates=as.numeric(data), Alleles=names(data), Parent=c(rep(p1,6),rep(p2,6),rep(p1,15),rep(p2,15),rep(p1,20),rep(p2,20)), Effects=c(rep("Additive",12),rep("Digenic",30),rep("Trigenic",40)))
        }
        data$Parent <- factor(data$Parent, levels=unique(data$Parent))
        plot <- ggplot(data[which(data$Effects == "Additive"),], aes(x = Alleles, y = Estimates, fill = Estimates)) +
          geom_bar(stat="identity") +
          scale_fill_gradient2(low = "red", high = "blue", guide = "none") +
          labs(title=names(x$results)[p], subtitle=paste("QTL", q, "\n")) +
          facet_wrap(. ~ Parent, scales="free_x", ncol = 2, strip.position="bottom") +
          # facet_grid(Effects ~ Parent, scales="free_x", space="free_x") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.text.x.bottom = element_text(hjust = 1, vjust = 0.5))
        res = c(res, plot)
        ## if (is.null(names(res))){
        ##   names(res) = paste0(names(x$results)[p],"_",q)
        ## } else {
        ##   names(res) = c(names(res), paste0(names(x$results)[p],"_",q))          
        ## }
        print(plot)
      }
    }
  }
  invisible(res)
}

