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
#' @param x an object of class \code{qtlpoly.fitted} to be summarized.
#'
#' @param pheno.col a numeric vector with the phenotype column numbers to be summarized; if \code{NULL}, all phenotypes from \code{'data'} will be included.
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
#'   \dontrun{
#'   # load raw data
#'   data(maps)
#'   data(pheno)
#'
#'   # estimate conditional probabilities using mappoly package
#'   library(mappoly)
#'   genoprob <- lapply(maps, calc_genoprob)
#'
#'   # prepare data
#'   data <- read_data(ploidy = 6, geno.prob = geno.prob, pheno = pheno, step = 1)
#'
#'   # perform remim
#'   remim.mod <- remim(data = data, w.size = 15, sig.fwd = 0.01, sig.bwd = 0.0001,
#'     d.sint = 1.5, n.clusters = 4, plot = "remim")
#'
#'   # fit model
#'   fitted.mod <- fit_remim(data=data, model=remim.mod, probs="joint", polygenes="none")
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Covarrubias-Pazaran G (2016) Genome-assisted prediction of quantitative traits using the R package sommer. \emph{PLoS ONE} 11 (6): 1–15. \url{http://doi.org/10.1371/journal.pone.0156744}.
#' 
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \url{http://doi.org/10.1534/genetics.120.303080}.
#'
#' @export fit_model
#' @importFrom sommer mmer

fit_model <- function(data, model, probs="joint", polygenes="none", keep=TRUE, verbose=TRUE) {

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

    indnames <- dimnames(data$Zd)[[3]]
    mrknames <- dimnames(data$Zd)[[2]]
    Za <- array(data = NA, dim = c(length(Pgametes), nmrk, nind), dimnames = list(c(Pgametes), c(mrknames), c(indnames))) # Za = Z1
    Zb <- array(data = NA, dim = c(length(Qgametes), nmrk, nind), dimnames = list(c(Qgametes), c(mrknames), c(indnames))) # Zb = Z2
    for(m in 1:nmrk) {
      for(i in 1:nind) {
        Za[,m,i] <- colSums(matrix(data$Z[,m,i], nrow=length(Pgametes)))
        Zb[,m,i] <- rowSums(matrix(data$Z[,m,i], nrow=length(Qgametes)))
      }
    }

  }
  for(p in 1:length(model$results)) {
    
    pheno.col <- model$results[[p]]$pheno.col
    qtl.mrk <- model$results[[p]]$qtl[,"Nmrk"]
    qtl.lgr <- model$results[[p]]$qtl[,"LG"]
    qtl.pos <- model$results[[p]]$qtl[,"Pos"]
    if(verbose) {
      if(length(qtl.mrk) == 1) cat("There is ", length(qtl.mrk), " QTL in the model for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), ". Fitting model... ", sep="")
      if(length(qtl.mrk) >= 2) cat("There are ", length(qtl.mrk), " QTL in the model for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), ". Fitting model... ", sep="")
    }
    
    if(!is.null(model$results[[p]]$qtls)) {
      nqtl <- dim(model$results[[p]]$qtls)[1]
      ind <- rownames(data$pheno)[which(!is.na(data$pheno[,pheno.col]))]
      Y <- data$pheno[ind,pheno.col]
      Zstar <- diag(length(ind))
      ETA <- NULL
      
      if(nqtl > 1) {
        
        if(probs=="joint" && polygenes=="none") {
          markers <- model$results[[p]]$qtls[,"Nmrk"]
          for(q in 1:nqtl) {
            eta <- list(list(Z=t(data$Z[,markers[q],ind]),K=data$Pi))
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
          ETA <- list(g=list(Z=t(data$Z[,markers,ind]),K=data$Pi))
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

#' @rdname fit_model
#' @export

summary.qtlpoly.fitted <- function(x, pheno.col = NULL) {
  if(any(class(x) == "qtlpoly.fitted")) cat("This is an object of class 'qtlpoly.fitted'\n")
  if(is.null(pheno.col)) {
    pheno.col <- 1:length(x$results)
  } else {
    pheno.col <- which(x$pheno.col %in% pheno.col)
  }
  for(p in pheno.col) {
    cat("\n* Trait", x$results[[p]]$pheno.col, sQuote(names(x$results)[[p]]), "\n")
    if(!is.null(x$results[[p]]$qtls)) print(x$results[[p]]$qtls)
    else cat("There are no QTL in the model \n")
  }
}

#' Adapt to sommer MNR core function
#'
#' Adapts genomic incidence and relationship (varcov) matrices to run using sommer's C++ core function (v. 4.0 or higher)
#' Function adapted from sommer v. 3.6 (Author: Giovanny Covarrubias-Pazaran)
#'
#' @param internal
#'
#' @author Gabriel de Siqueira Gesteira, \email{gdesiqu@@ncsu.edu}
#'
#' @references
#'     Covarrubias-Pazaran G (2016) Genome-assisted prediction of quantitative traits using the R package sommer. \emph{PLoS ONE} 11 (6): 1–15. \url{http://doi.org/10.1371/journal.pone.0156744}.
#' 
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \url{http://doi.org/10.1534/genetics.120.303080}.
mmer_adapted <- function(Y,X=NULL,Z=NULL,R=NULL,W=NULL,method="NR",init=NULL,iters=20,tolpar=1e-3,
                 tolparinv=1e-6,draw=FALSE,silent=FALSE, constraint=TRUE, 
                 EIGEND=FALSE, forced=NULL, IMP=FALSE, complete=TRUE, 
                 check.model=TRUE, restrained=NULL, REML=TRUE, init.equal=TRUE,
                 date.warning=TRUE, verbose=FALSE){
  #R=NULL
  if (inherits(Y, "formula")){
    stop("\nYou have to use mmer2 function for formula-based models.\n", call. = FALSE)
  }
  diso <- dim(as.data.frame(Y))[2]
  ## control for 2-level list structure
  if(!is.list(Z)){
    stop("Please provide the Z parameter as a 2 level list structure.\nFor example for 2 random effects 'A' and 'B' do:\n    ETA <- list( A=list( Z=myZ1, K=myK1 ) , B=list( Z=myZ2, K=myK2 ) )\n    mod <- mmer(Y=y, Z=ETA)\nwhere Z's and K's are the incidence and var-covar matrices respectively.\nIf any Z or K is not provided, an identity matrix will be assumed. ",call. = FALSE)
  }else{
    if(!is.list(Z[[1]])){
      stop("Please provide the Z parameter as a 2 level list structure.\nFor example for 2 random effects 'A' and 'B' do:\n    ETA <- list( A=list( Z=myZ1, K=myK1 ) , B=list( Z=myZ2, K=myK2 ) )\n    mod <- mmer(Y=y, Z=ETA)\nwhere Z's and K's are the incidence and var-covar matrices respectively.\nIf any Z or K is not provided, an identity matrix will be assumed. ",call. = FALSE)
    }
  }
  ## control for Z-K names
  zzzkkk <- unlist(lapply(Z,function(x){length(names(x))}))
  badRE <- which(zzzkkk==0) # BAD RE WITH NO NAMES
  if(length(badRE)>0){
    stop("Please when specifying a random effect use the names; \n'Z' for incidence and 'K' for variance-covariance matrices.\nFor example for 1 random effect (i.e. named 'A') model do:\n    ETA <- list( A=list( Z=M1, K=M2) )\n    mod <- mmer(Y=y, Z=ETA)\nSpecifying at least one; Z or K. You need to specify if is a 'Z' or 'K' \nsince this is the only way the program distinguishes between matrices.",call. = FALSE)
  }
  
  #########*****************************
  ## make sure user don't provide the same names for random effects
  his.names <- names(Z)
  if(!is.null(his.names)){
    badnames <- which(duplicated(his.names))
    if(length(badnames)>0){
      his.names[badnames] <- paste(his.names[badnames],1:(length(his.names[badnames])), sep=".")
      names(Z) <- his.names
    }
  }
  dZ <- unlist(lapply(Z,function(x){dim(x$Z)[1]}))
  
  if(is.null(dZ)){ #sometimes user don't specify the Z matrices
    dZ <- unlist(lapply(Z,function(x){dim(x$K)[1]}))
  }
  if(!is.null(X)){
    dZ <- c(dZ,dim(X)[1])
  }
  dall <- unlist(c(dZ,dim(as.matrix(Y))[1]))
  if(length(which(!duplicated(dall))) > 1){
    if(is.null(X)){
      stop("Matrices Y and Z's should have the same number of individuals. \nPlease check the dimensions of your matrices.", call. = FALSE)
    }else{
      stop("Matrices Y, X and Z's should have the same number of individuals. \nPlease check the dimensions of your matrices.", call. = FALSE)
    }
  }
  #########*****************************
  for(bb in 1:length(Z)){
    ss1 <- colnames(Z[[bb]]$Z) == colnames(Z[[bb]]$K)
    if(length(which(!ss1))>0){
      print(paste("Names of columns in matrices Z and K for the",bb,"th random effect do not match.")) 
      print("This can lead to incorrect estimation of variance components. Double check.")
    }
  }
  #########*****************************
  #   if(!is.null(W)){
  #     cat("Response is imputed for estimation of variance components in GWAS models.\n")
  #   }
  #########*****************************
  
  if(!is.null(X)){
    if(is.list(X)){
      stop("Multivariate models only accept one incidence matrix for fixed effects (X). Please modifiy your X argument.",call. = FALSE)
    }
  }
  #if(!silent){cat("Running multivariate model\n")}
  
  if(check.model){
    if(is.list(Z)){
      if(is.list(Z[[1]])){ ### -- if is a 2 level list -- ##
        provided <- lapply(Z, names)
        for(s in 1:length(provided)){ #for each random effect =============================
          provided2 <- names(Z[[s]])
          if(length(provided2) ==1){ #----the 's' random effect has one matrix only----
            if(provided2 == "K"){ #user only provided K
              #zz <- diag(length(y))#model.matrix(~rownames(Z[[s]][[1]]))
              zz <- diag(nrow(as.matrix(Y)))
              #colnames(zz) <- rownames(Z[[s]][[1]])
              Z[[s]] <- list(Z=zz, K=Z[[s]][[1]])
            }
            if(provided2 == "Z"){ # user only provided Z
              #kk <- diag(dim(Z[[s]][[1]])[2])
              kk <- diag(dim(Z[[s]][[1]])[2])
              attributes(kk)$diagon <- TRUE
              #rownames(kk) <- colnames(Z[[s]][[1]]); colnames(kk) <- rownames(kk)
              Z[[s]] <- list(Z=Z[[s]][[1]],K=kk) 
            }
          }else{ #----the 's' random effect has two matrices----
            dido<-lapply(Z[[s]], dim) # dimensions of Z and K
            condi<-(dido$Z[2] == dido$K[1] & dido$Z[2] == dido$K[2]) 
            # condition, column size on Z matches with a square matrix K
            if(!condi){
              cat(paste("ERROR! In the",s,"th random effect you have provided or created an incidence \nmatrix with dimensions:",dido$Z[1],"rows and",dido$Z[2],"columns. Therefore the \nvariance-covariance matrix(K) for this random effect expected was a \nsquare matrix with dimensions",dido$Z[2],"x",dido$Z[2]),", but you provided a",dido$K[1],"x",dido$K[2]," matrix \nas a variance-covariance matrix. Please double check your matrices. Might be that you didn't drop all levels while subsetting your data, if so try: yourdata <- droplevels(yourdata).")
              stop()
            }
          }#---------------------------------------------------------------------------
        } #for each random effect end =================================================
      }else{ # if is a one-level list !!!!!!!!!!!!!
        if(length(Z) == 1){ ## -- if the user only provided one matrix -- ##
          provided <- names(Z)
          if(provided == "K"){
            #zz <- diag(length(y))
            zz <- diag(nrow(as.matrix(Y)))
            Z <- list(Z=zz, K=Z[[1]])
          }
          if(provided == "Z"){
            #kk <- diag(dim(Z[[1]])[2])
            kk <- diag(dim(Z[[1]])[2])
            attributes(kk)$diagon <- TRUE
            #rownames(kk) <- colnames(Z[[1]]); colnames(kk) <- rownames(kk)
            Z <- list(Z=Z[[1]],K=kk) 
          }
        }else{ # there's 2 matrices in Z
          dido<-lapply(Z, dim) # dimensions of Z and K
          condi<-(dido$Z[2] == dido$K[1] & dido$Z[2] == dido$K[2]) 
          # condition, column size on Z matches with a square matrix K
          if(!condi){
            cat(paste("ERROR! In the",s,"th random effect you have provided or created an incidence \nmatrix with dimensions:",dido$Z[1],"rows and",dido$Z[2],"columns. Therefore the \nvariance-covariance matrix(K) for this random effect expected was a \nsquare matrix with dimensions",dido$Z[2],"x",dido$Z[2]),", but you provided a",dido$K[1],"x",dido$K[2]," matrix \nas a variance-covariance matrix. Please double check your matrices. Might be that you didn't drop all levels while subsetting your data, if so try: yourdata <- droplevels(yourdata).")
            stop()
          }else{Z=list(Z=Z)}
        }
      }
    }else{
      if(is.null(Z)){ # the user is not using the random part
        cat("Error. No random effects specified in the model. \nPlease use 'lm' or provide a diagonal matrix in Z\ni.e. Zu = list(A=list(Z=diag(length(y))))\n")
        stop()
      }else{
        #stop;
        cat("\nThe parameter 'Z' needs to be provided in a 2-level list structure. \n\nPlease see help typing ?mmer and look at the 'Arguments' section\n")
        cat("\nIf no random effects provided, the model will be fitted using the 'lm' function\n\n")
      }
    }
  }

  Y2 = matrix(Y)
  rownames(Y2) = names(Y)
  colnames(Y2) = "Yield"
  Y = Y2
  x = matrix(rep(1, nrow(Y))) # fixed effect (only intercept)
  colnames(x) = "(Intercept)"
  X = list(x)
  Gx = list(matrix(1)) # list of constraints for fixed effects
  ZETA = lapply(Z, function(x) return(x[[1]])) # list containing only incidence matrices of random effects (sparse form)
  ## names(ZETA) = NULL
  for (i in 1:length(ZETA)){
    attr(ZETA[[i]], "assign") = rep(1, ncol(ZETA[[i]]))
    attr(ZETA[[i]], "contrasts") = list("contr.treatment")
    names(attr(ZETA[[i]], "contrasts")) = paste0("QTL ", i)
  }
  ZETA <- lapply(ZETA,function(x){as(x, Class = "sparseMatrix")})
  K = lapply(Z, function(x) return(x[[2]])) # list of varcov matrices of random effects (sparse form)
  ## names(K) = NULL
  ## for (i in 1:length(K)){
  ##   dimnames(K[[i]]) = NULL
  ## }
  R = list(diag(1, nrow = nrow(Y), ncol = nrow(Y))) # list of residual matrices (sparse form)
  rownames(R[[1]]) = rownames(Y)
  colnames(R[[1]]) = paste0("units",sort(as.character(seq(1, nrow(Y), 1)), decreasing = F))
  attr(R[[1]], "assign") = rep(1, nrow(R[[1]]))
  attr(R[[1]], "constrasts") = list(units = "contr.treatment")
  R <- lapply(R,function(x){as(x, Class = "sparseMatrix")})
  GES = rep(list(matrix(0.15)), length(ZETA)) # list of initial values for variance components
  GES[[length(GES)+1]] = matrix(0.75)
  GESI = rep(list(matrix(1)), length(ZETA)+1) # list of constraints applied to random effects
  selected = FALSE
  getPEV = TRUE
  stepweight <- rep(0.9,iters); stepweight[1:2] <- c(0.5,0.7)
  emupdate <- rep(0, iters)
  ws = rep(1, nrow(Y))
    
  RES = .Call("_sommer_MNR",PACKAGE = "sommer",Y, X,Gx,ZETA,K,R,GES,GESI, ws, iters, tolpar, tolparinv, selected,getPEV,verbose, FALSE, stepweight, emupdate)
  RES$alleles = rownames(K[[1]])
  
  ## RES <- MNR(Y=Y,X=X,ZETA=Z,R=R,W=W,init=init,iters=iters,tolpar=tolpar,
  ##            tolparinv = tolparinv,draw=draw,silent=silent, 
  ##            constraint = constraint,EIGEND = EIGEND,
  ##            forced=forced,IMP=IMP,restrained=restrained, REML=REML, 
  ##            init.equal = init.equal)
  ## class(RES)<-c("MMERM")
  ## layout(matrix(1,1,1))
  return(RES)
}
