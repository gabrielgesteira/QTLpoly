#' Fits multiple QTL models
#'
#' Fits alternative multiple QTL models by performing variance component estimation using REML.
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param model an object of class \code{qtlpoly.profile} or \code{qtlpoly.remim}.
#'
#' @param probs a character string indicating if either \code{"joint"} (genotypes) or \code{"marginal"} (parental gametes) conditional probabilities should be used.
#'
#' @param polygenes a character string indicating if either \code{"none"}, \code{"most"} or \code{"all"} QTL should be used as polygenes.
#'
#' @param keep if \code{TRUE} (default), stores all matrices and estimates from fitted model; if \code{FALSE}, nothing is stored.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced.
#' 
#' @param object an object of class \code{qtlpoly.fitted} to be summarized.
#'
#' @param pheno.col a numeric vector with the phenotype column numbers to be summarized; if \code{NULL} (default), all phenotypes from \code{'data'} will be included.
#'
#' @param ... currently ignored
#' 
#' @return An object of class \code{qtlpoly.fitted} which contains a list of \code{results} for each trait with the following components:
#'
#'     \item{pheno.col}{a phenotype column number.}
#'     \item{fitted}{a \pkg{sommer} object of class \code{mmer}.}
#'     \item{qtls}{a data frame with information from the mapped QTL.}
#'
#' @seealso \code{\link[qtlpoly]{read_data}}, \code{\link[qtlpoly]{remim}}
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
#'   fitted.mod = fit_model(data=data, model=remim.mod, probs="joint", polygenes="none")
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Covarrubias-Pazaran G (2016) Genome-assisted prediction of quantitative traits using the R package sommer. \emph{PLoS ONE} 11 (6): 1–15. \doi{10.1371/journal.pone.0156744}.
#' 
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \doi{10.1534/genetics.120.303080}.
#'
#' @useDynLib qtlpoly, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export fit_model2
fit_model2 <- function(data, model, probs="joint", polygenes="none", keep=TRUE, verbose=TRUE, pheno.col = NULL) {

  results <- vector("list", length(model$results))
  names(results) <- names(model$results)

  if(probs=="marginal") {
    
    ploidy <- data$ploidy
    Palleles <- letters[1:ploidy]
    Pgametes <- lapply(combn(Palleles, ploidy/2, simplify = FALSE), paste, collapse="")
    Qalleles <- letters[(ploidy+1):(2*ploidy)]
    Qgametes <- lapply(combn(Qalleles, ploidy/2, simplify = FALSE), paste, collapse="")
    sibs <- sapply( Pgametes, FUN=function(x) paste(x, Pgametes, sep="") )
    Pia <- matrix(data = NA, nrow = length(Pgametes), ncol = length(Pgametes))
    for(i in 1:ncol(sibs)) {
      for(j in 1:nrow(sibs)) {
        Pia[i,j] <- (ploidy-length(unique(strsplit(sibs[i,j], "")[[1]])))/(ploidy/2)
      }
    }
    Pib <- Pia # Pi1 = Pi2
    colnames(Pia) <- rownames(Pia) <- unlist(Pgametes)
    colnames(Pib) <- rownames(Pib) <- unlist(Qgametes)

    indnames <- dimnames(data$Z)[[3]]
    mrknames <- dimnames(data$Z)[[2]]
    nmrk = length(mrknames)
    nind = length(indnames)
    Za <- array(data = NA, dim = c(length(Pgametes), nmrk, nind), dimnames = list(c(Pgametes), c(mrknames), c(indnames))) # Za = Z1
    Zb <- array(data = NA, dim = c(length(Qgametes), nmrk, nind), dimnames = list(c(Qgametes), c(mrknames), c(indnames))) # Zb = Z2
    for(m in 1:nmrk) {
      for(i in 1:nind) {
        Za[,m,i] <- colSums(matrix(data$Z[,m,i], nrow=length(Pgametes)))
        Zb[,m,i] <- rowSums(matrix(data$Z[,m,i], nrow=length(Qgametes)))
      }
    }

  }
  ploidy <- data$ploidy
  for(p in 1:length(model$results)) {
    
    pheno.col <- model$results[[p]]$pheno.col
    qtl.mrk <- model$results[[p]]$qtl[,"Nmrk"]
    qtl.lgr <- model$results[[p]]$qtl[,"LG"]
    qtl.pos <- model$results[[p]]$qtl[,"Pos"]
    if(verbose) {
      if(length(qtl.mrk) == 1 && verbose) cat("There is ", length(qtl.mrk), " QTL in the model for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), ". Fitting model... ", sep="")
      if(length(qtl.mrk) >= 2 && verbose) cat("There are ", length(qtl.mrk), " QTL in the model for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), ". Fitting model... ", sep="")
    }
    
    if(!is.null(model$results[[p]]$qtls)) {
      nqtl <- dim(model$results[[p]]$qtls)[1]
      ind <- rownames(data$pheno)[which(!is.na(data$pheno[,pheno.col]))]
      Y <- data$pheno[ind,pheno.col]
      Zstar <- diag(length(ind))
      Zstar2 <- diag(ploidy*2)/4
      rownames(Zstar2) = colnames(Zstar2) = letters[1:(ploidy*2)]
      ETA <- NULL
      
      if(nqtl > 1) {
        
        if(probs=="joint" && polygenes=="none") {
          markers <- model$results[[p]]$qtls[,"Nmrk"]
          for(q in 1:nqtl) {
            ## eta <- list(list(Z=t(data$Z[,markers[q],ind]),K=data$Pi))
              ## eta <- list(list(Z=Zstar,K=data$G[ind,ind,markers[q]]))
              eta <- list(list(Z=t(data$Z[,markers[q],ind]),K=Zstar2))
            names(eta) <- paste("g", q, sep="")
            ETA <- c(ETA, eta)
          }
          fitted <- mmer_adapted(Y=Y, Z=ETA, silent=TRUE, date.warning=FALSE)
          qtl <- c()
          for(q in 1:nqtl) {
            H2 <- sum(c(fitted$sigma)[q])/sum(c(fitted$sigma))
            qtl <- c(qtl, c(model$results[[p]]$qtls[q,c(1:3)], NA, c(fitted$sigma)[q], NA, H2))
          }
          H2 <- sum(c(fitted$sigma)[1:nqtl])/sum(c(fitted$sigma))
          qtl <- c(qtl, NA, NA, NA, c(fitted$Beta), NA, c(fitted$sigma)[nqtl+1], H2)
          qtls <- as.data.frame(matrix(qtl, ncol=7, byrow=TRUE), stringsAsFactors=FALSE)
          colnames(qtls) <- c("LG", "Pos", "Nmrk", "Intercept", "Var(g)", "Var(e)", "h2")
        }
        
        if(probs=="joint" && polygenes=="most") {
          fitted <- vector("list", nqtl)
          qtl <- c()
          for(q in 1:nqtl) {
            markers <- model$results[[p]]$qtls[-q,"Nmrk"]
            Gstar <- apply(data$G[ind,ind,markers], MARGIN = c(1,2), sum)/length(markers); Gstar[1:5,1:5]
            m <- model$results[[p]]$qtls[q,"Nmrk"]
            ETA <- list(g=list(Z=t(data$Z[,m,ind]),K=data$Pi),gstar=list(Z=Zstar,K=Gstar))
            fitted[[q]] <- mmer_adapted(Y=Y, Z=ETA, silent=TRUE, date.warning=FALSE)
            H2 <- c(fitted[[q]]$sigma)[1]/sum(c(fitted[[q]]$sigma))
            qtl <- c(qtl, c(model$results[[p]]$qtls[q,c(1:3)], c(fitted[[q]]$Beta), c(fitted[[q]]$sigma), H2))
          }
          qtls <- as.data.frame(matrix(qtl, ncol=8, byrow=TRUE), stringsAsFactors=FALSE)
          colnames(qtls) <- c("LG", "Pos", "Nmrk", "Intercept", "Var(g)", "Var(g*)", "Var(e)", "h2")
        }
        
        if(probs=="joint" && polygenes=="all") {
          markers <- model$results[[p]]$qtls[,"Nmrk"]
          Gstar <- apply(data$G[ind,ind,markers], MARGIN = c(1,2), sum)/length(markers); Gstar[1:5,1:5]
          ETA <- list(gstar=list(Z=Zstar,K=Gstar))
          fitted <- mmer_adapted(Y=Y, Z=ETA, silent=TRUE, date.warning=FALSE)
          H2 <- c(fitted$sigma)[1]/sum(c(fitted$sigma))
          qtl <- c(rep(NA, 4*nqtl), c(fitted$Beta), c(fitted$sigma), H2)
          qtls <- as.data.frame(cbind(rbind(model$results[[p]]$qtls[,c(1:3)], rep(NA, 3)), matrix(qtl, ncol=4, byrow=TRUE), stringsAsFactors=FALSE))
          colnames(qtls) <- c("LG", "Pos", "Nmrk", "Intercept", "Var(g*)", "Var(e)", "h2")
        }
        
        if(probs=="marginal" && polygenes=="none") {
          markers <- model$results[[p]]$qtls[,"Nmrk"]
          for(q in 1:nqtl) {
            eta <- list(list(Z=t(Za[,markers[q],ind]),K=Pia),list(Z=t(Zb[,markers[q],ind]),K=Pib),list(Z=t(data$Z[,markers[q],ind]),K=data$Pi))
            names(eta) <- paste(c("a","b","c"), rep(q,3), sep="")
            ETA <- c(ETA, eta)
          }
          fitted <- mmer_adapted(Y=Y, Z=ETA, silent=TRUE, date.warning=FALSE)
          qtl <- c()
          for(q in 1:nqtl) {
            H2 <- sum(c(fitted$sigma)[(1:3)+((q-1)*3)])/sum(c(fitted$sigma))
            qtl <- c(qtl, c(model$results[[p]]$qtls[q,c(1:3)], NA, c(fitted$sigma)[(1:3)+((q-1)*3)], NA, H2))
          }
          H2 <- sum(c(fitted$sigma)[1:(nqtl*3)])/sum(c(fitted$sigma))
          qtl <- c(qtl, NA, NA, NA, c(fitted$Beta), rep(NA, 3), c(fitted$sigma)[(nqtl*3)+1], H2)
          qtls <- as.data.frame(matrix(qtl, ncol=9, byrow=TRUE), stringsAsFactors=FALSE)
          colnames(qtls) <- c("LG", "Pos", "Nmrk", "Intercept", "Var(a)", "Var(b)", "Var(d)", "Var(e)", "h2")
        }
        
        if(probs=="marginal" && polygenes=="most") {
          fitted <- vector("list", nqtl)
          qtl <- c()
          for(q in 1:nqtl) {
            markers <- model$results[[p]]$qtls[-q,"Nmrk"]
            Gstar <- apply(data$G[ind,ind,markers], MARGIN = c(1,2), sum)/length(markers); Gstar[1:5,1:5]
            m <- model$results[[p]]$qtls[q,"Nmrk"]
            ETA <- list(a=list(Z=t(Za[,m,ind]),K=Pia),b=list(Z=t(Zb[,m,ind]),K=Pib),d=list(Z=t(data$Z[,m,ind]),K=data$Pi),gstar=list(Z=Zstar,K=Gstar))
            fitted[[q]] <- mmer_adapted(Y=Y, Z=ETA, silent=TRUE, date.warning=FALSE)
            H2 <- sum(c(fitted[[q]]$sigma)[(1:3)])/sum(c(fitted[[q]]$sigma))
            qtl <- c(qtl, c(model$results[[p]]$qtls[q,c(1:3)], c(fitted[[q]]$Beta), c(fitted[[q]]$sigma), H2))
          }
          qtls <- as.data.frame(matrix(qtl, ncol=10, byrow=TRUE), stringsAsFactors=FALSE)
          colnames(qtls) <- c("LG", "Pos", "Nmrk", "Intercept", "Var(a)", "Var(b)", "Var(d)", "Var(g*)", "Var(e)", "h2")
        }
        
        if(probs=="marginal" && polygenes=="all") {
          markers <- model$results[[p]]$qtls[,"Nmrk"]
          A <- array(data = NA, dim = c(data$nind, data$nind, nqtl))
          B <- array(data = NA, dim = c(data$nind, data$nind, nqtl))
          for(q in 1:nqtl) {
            A[,,q] <- t(Za[,markers[q],])%*%Pia%*%Za[,markers[q],]
            B[,,q] <- t(Zb[,markers[q],])%*%Pib%*%Zb[,markers[q],]
          }
          Astar <- apply(A[ind,ind,], MARGIN = c(1,2), sum)/length(markers); Astar[1:5,1:5]
          Bstar <- apply(B[ind,ind,], MARGIN = c(1,2), sum)/length(markers); Bstar[1:5,1:5]
          Dstar <- apply(data$G[ind,ind,markers], MARGIN = c(1,2), sum)/length(markers); Dstar[1:5,1:5]
          ETA <- list(astar=list(Z=Zstar,K=Astar),bstar=list(Z=Zstar,K=Bstar),dstar=list(Z=Zstar,K=Dstar))
          fitted <- mmer_adapted(Y=Y, Z=ETA, silent=TRUE, date.warning=FALSE)
          H2 <- sum(c(fitted$sigma)[1:3])/sum(c(fitted$sigma))
          qtl <- c(rep(NA, 6*nqtl), c(fitted$Beta), c(fitted$sigma), H2)
          qtls <- as.data.frame(cbind(rbind(model$results[[p]]$qtls[,c(1:3)], rep(NA, 3)), matrix(qtl, ncol=6, byrow=TRUE), stringsAsFactors=FALSE))
          colnames(qtls) <- c("LG", "Pos", "Nmrk", "Intercept", "Var(a*)", "Var(b*)", "Var(d*)", "Var(e)", "h2")
        }
        
      }
      
      if(nqtl == 1) {
        
        if(probs=="joint") {
          markers <- model$results[[p]]$qtls[,"Nmrk"]
          ## ETA <- list(g=list(Z=t(data$Z[,markers,ind]),K=data$Pi))
          ## ETA <- list(g=list(Z=Zstar,K=data$G[ind,ind,markers]))
          ETA <- list(g=list(Z=t(data$Z[,markers,ind]),K=Zstar2))
          fitted <- mmer_adapted(Y=Y, Z=ETA, silent=TRUE, date.warning=FALSE)
          H2 <- sum(c(fitted$sigma)[1])/sum(c(fitted$sigma))
          qtl <- c(c(fitted$Beta), c(fitted$sigma), H2)
          qtls <- as.data.frame(c(model$results[[p]]$qtls[,c(1:3)], matrix(qtl, ncol=4, byrow=TRUE)), stringsAsFactors=FALSE)
          colnames(qtls) <- c("LG", "Pos", "Nmrk", "Intercept", "Var(g)", "Var(e)", "h2")
        }
        
        if(probs=="marginal") {
          markers <- model$results[[p]]$qtls[,"Nmrk"]
          ETA <- list(a=list(Z=t(Za[,markers,ind]),K=Pia),b=list(Z=t(Zb[,markers,ind]),K=Pib),d=list(Z=t(data$Z[,markers,ind]),K=data$Pi))
          fitted <- mmer_adapted(Y=Y, Z=ETA, silent=TRUE, date.warning=FALSE)
          H2 <- sum(c(fitted$sigma)[1:3])/sum(c(fitted$sigma))
          qtl <- c(c(fitted$Beta), c(fitted$sigma), H2)
          qtls <- as.data.frame(c(model$results[[p]]$qtls[,c(1:3)], matrix(qtl, ncol=6, byrow=TRUE)), stringsAsFactors=FALSE)
          colnames(qtls) <- c("LG", "Pos", "Nmrk", "Intercept", "Var(a)", "Var(b)", "Var(d)", "Var(e)", "h2")
        }
      }
      
      if(verbose) cat("Done! \n\n", sep="")
      
    }
    
    if(is.null(model$results[[p]]$qtls)) {
      if(verbose) cat("There are no QTL in the model for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), ". No model has been fitted! \n\n", sep="")
      fitted <- NULL
      qtls <- NULL
    }
    
    if(!keep) fitted <- NULL
    
    results[[p]] <- list(
      pheno.col=pheno.col,
      fitted=fitted,
      qtls=qtls)
    
  }

  structure(list(data=deparse(substitute(data)),
                 model=deparse(substitute(model)),
                 pheno.col=model$pheno.col,
                 probs=probs,
                 polygenes=polygenes,
                 results=results),
            class="qtlpoly.fitted")

}
