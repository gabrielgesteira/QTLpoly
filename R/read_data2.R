#' Read genotypic and phenotypic data
#'
#' Reads files in specific formats and creates a \code{qtlpoly.data} object to be used in subsequent analyses.
#'
#' @param ploidy a numeric value of ploidy level of the cross.
#'
#' @param geno.prob an object of class \code{mappoly.genoprob} from \pkg{mappoly}.
#'
#' @param geno.dose an object of class \code{mappoly.data} from \pkg{mappoly}.
#'
#' @param double.reduction if \code{TRUE}, double reduction genotypes are taken into account; if \code{FALSE}, no double reduction genotypes are considered.
#'
#' @param pheno a data frame of phenotypes (columns) with individual names (rows) identical to individual names in \code{geno.prob} and/or \code{geno.dose} object.
#'
#' @param weights a data frame of phenotype weights (columns) with individual names (rows) identical to individual names in \code{pheno} object.
#'
#' @param step a numeric value of step size (in centiMorgans) where tests will be performed, e.g. 1 (default); if \code{NULL}, tests will be performed at every marker.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @return An object of class \code{qtlpoly.data} which is a list containing the following components:
#'
#'     \item{ploidy}{a scalar with ploidy level.}
#'     \item{nlgs}{a scalar with the number of linkage groups.}
#'     \item{nind}{a scalar with the number of individuals.}
#'     \item{nmrk}{a scalar with the number of marker positions.}
#'     \item{nphe}{a scalar with the number of phenotypes.}
#'     \item{lgs.size}{a vector with linkage group sizes.}
#'     \item{cum.size}{a vector with cumulative linkage group sizes.}
#'     \item{lgs.nmrk}{a vector with number of marker positions per linkage group.}
#'     \item{cum.nmrk}{a vector with cumulative number of marker positions per linkage group.}
#'     \item{lgs}{a list with selected marker positions per linkage group.}
#'     \item{lgs.all}{a list with all marker positions per linkage group.}
#'     \item{step}{a scalar with the step size.}
#'     \item{pheno}{a data frame with phenotypes.}
#'     \item{G}{a list of relationship matrices for each marker position.}
#'     \item{Z}{a list of conditional probability matrices for each marker position for genotypes.}
#'     \item{X}{a list of conditional probability matrices for each marker position for alleles.}
#'     \item{Pi}{a matrix of identical-by-descent shared alleles among genotypes.}
#'
#' @seealso \code{\link[qtlpoly]{maps6x}}, \code{\link[qtlpoly]{pheno6x}}
#'
#' @examples
#'   \donttest{
#'   # Estimate conditional probabilities using mappoly package
#'   library(mappoly)
#'   library(qtlpoly)
#'   genoprob4x = lapply(maps4x[c(5)], calc_genoprob)
#'   data = read_data(ploidy = 4, geno.prob = genoprob4x, pheno = pheno4x, step = 1)
#'   }
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}, Gabriel de Siqueira Gesteira, \email{gdesiqu@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \doi{10.1534/genetics.120.303080}.
#'
#' @export read_data2
#' @importFrom abind abind
#' @importFrom gtools combinations
#' @importFrom mappoly calc_genoprob

read_data2 <- function(ploidy = 6, geno.prob, geno.dose = NULL, type=c("genome","mds","custom"), double.reduction = FALSE, pheno, weights = NULL, step = 1, verbose = TRUE) {

  if(inherits(geno.prob, "mappoly2.sequence")){
  ## if (class(geno.prob) == "mappoly2.sequence"){
        
  if(is.null(step)) step <- 1e-10

  homo.prob = geno.prob
   
  raw.individual.names = homo.prob$data$screened.data$ind.names
  
  type <- match.arg(type) # Ensures only "genome","mds", or  "custom" are allowed
  
  
  ## Converting object back to previous format
  if(type=="genome"){
    if(!exists(homo.prob$maps[[1]]$genome$p1p2$hmm.phase[[1]]$haploprob)){stop("Haploprobs ahve not been calculated. Please go back to mappoly2 to calculate halploprobs")}
    geno.prob = lapply(homo.prob$maps, function(x) {
      probs = x$genome$p1p2$hmm.phase[[1]]$haploprob
      a = split(1:nrow(probs), ceiling(seq_along(1:nrow(probs)) / (ploidy*2)))
      b = lapply(a, function(y) return(as.matrix(probs[y,-c(1:3)])))
      c = abind(b, along = 3)
      dimnames(c)[[1]] = letters[1:(ploidy*2)]
      dimnames(c)[[2]] = rownames(x$genome$p1p2$hmm.phase[[1]]$p1)
      dimnames(c)[[3]] = raw.individual.names
      mpgpt = calc_genoprob # to ensure mappoly's function is required in the package
      map = c(0, cumsum(imf_h(x$genome$p1p2$hmm.phase[[1]]$rf)))
      names(map) = rownames(x$genome$p1p2$hmm.phase[[1]]$p1)
      return(list(probs = c, map = map))
    })
  } else if(type =="mds"){
    if(!exists(homo.prob$maps[[1]]$mds$p1p2$hmm.phase[[1]]$haploprob)){stop("Haploprobs ahve not been calculated. Please go back to mappoly2 to calculate halploprobs")}
    geno.prob = lapply(homo.prob$maps, function(x) {
      probs = x$mds$p1p2$hmm.phase[[1]]$haploprob
      a = split(1:nrow(probs), ceiling(seq_along(1:nrow(probs)) / (ploidy*2)))
      b = lapply(a, function(y) return(as.matrix(probs[y,-c(1:3)])))
      c = abind(b, along = 3)
      dimnames(c)[[1]] = letters[1:(ploidy*2)]
      dimnames(c)[[2]] = rownames(x$mds$p1p2$hmm.phase[[1]]$p1)
      dimnames(c)[[3]] = raw.individual.names
      mpgpt = calc_genoprob # to ensure mappoly's function is required in the package
      map = c(0, cumsum(imf_h(x$mds$p1p2$hmm.phase[[1]]$rf)))
      names(map) = rownames(x$mds$p1p2$hmm.phase[[1]]$p1)
      return(list(probs = c, map = map))
    } else if(type =="custom"){
      if(!exists(homo.prob$maps[[1]]$custom$p1p2$hmm.phase[[1]]$haploprob)){stop("Haploprobs ahve not been calculated. Please go back to mappoly2 to calculate halploprobs")}
      geno.prob = lapply(homo.prob$maps, function(x) {
        probs = x$custom$p1p2$hmm.phase[[1]]$haploprob
        a = split(1:nrow(probs), ceiling(seq_along(1:nrow(probs)) / (ploidy*2)))
        b = lapply(a, function(y) return(as.matrix(probs[y,-c(1:3)])))
        c = abind(b, along = 3)
        dimnames(c)[[1]] = letters[1:(ploidy*2)]
        dimnames(c)[[2]] = rownames(x$custom$p1p2$hmm.phase[[1]]$p1)
        dimnames(c)[[3]] = raw.individual.names
        mpgpt = calc_genoprob # to ensure mappoly's function is required in the package
        map = c(0, cumsum(imf_h(x$custom$p1p2$hmm.phase[[1]]$rf)))
        names(map) = rownames(x$custom$p1p2$hmm.phase[[1]]$p1)
        return(list(probs = c, map = map))
      }
      
  
  ######### HAPLOTYPE DATA
  
  if(!is.null(geno.prob)) {
    nlgs <- length(geno.prob)
    probs <- vector("list", nlgs)
    lgs <- vector("list", nlgs)
    lgs.all <- vector("list", nlgs)
    lgs.size <- numeric(nlgs)
    lgs.nmrk <- numeric(nlgs)
    for(c in 1:nlgs) {
      sel.map <- which(!duplicated(geno.prob[[c]]$map))
      sel.probs <- geno.prob[[c]]$probs[,sel.map,]
      lgs.all[[c]] <- geno.prob[[c]]$map[sel.map]
      lgs[[c]] <- unique(c(lgs.all[[c]][!duplicated(floor(lgs.all[[c]]/step)*step)], last(lgs.all[[c]])))
      lgs.size[c] <- last(lgs[[c]])
      names(lgs[[c]]) <- names(lgs.all[[c]])[which(lgs.all[[c]] %in% lgs[[c]])]
      probs[[c]] <- sel.probs[,which(lgs.all[[c]] %in% lgs[[c]]),]
      lgs.nmrk[c] <- dim(probs[[c]])[2]
    }
    cum.size <- c(0, cumsum(lgs.size))
    cum.nmrk <- c(0, cumsum(lgs.nmrk))
    
    Z <- abind::abind(probs, along = 2); dim(Z)
    ind.names <- dimnames(Z)[[3]]
    mrk.names <- dimnames(Z)[[2]]
    nmrk <- dim(Z)[2]
    nind <- dim(Z)[3]
    
    if(double.reduction) {
      if(ploidy == 4) n.unique <- 1
      if(ploidy == 6) n.unique <- 2
      if(ploidy == 8) n.unique <- 3
      if(ploidy == 10) n.unique <- 4
      if(ploidy == 12) n.unique <- 5
    } else {
      n.unique <- ploidy/2
    }
    Palleles <- letters[1:ploidy]; length(Palleles)
    Pgametes <- gtools::combinations(length(Palleles), ploidy/2, Palleles, repeats.allowed = TRUE); dim(Pgametes)
    Punique <- apply(Pgametes, 1, unique); length(Punique)
    Pgametes <- apply(Pgametes[which(lapply(Punique, length) >= n.unique),], 1, paste, collapse=""); length(Pgametes)
    # Pgametes <- lapply(combn(Palleles, ploidy/2, simplify = FALSE), paste, collapse=""); length(Pgametes)
    Qalleles <- letters[(ploidy+1):(2*ploidy)]
    Qgametes <- gtools::combinations(length(Qalleles), ploidy/2, Qalleles, repeats.allowed = TRUE); dim(Qgametes)
    Qunique <- apply(Qgametes, 1, unique); length(Qunique)
    Qgametes <- apply(Qgametes[which(lapply(Qunique, length) >= n.unique),], 1, paste, collapse=""); length(Qgametes)
    # Qgametes <- lapply(combn(Qalleles, ploidy/2, simplify = FALSE), paste, collapse="")
    genotypes <- as.vector(t(outer(Pgametes, Qgametes, paste, sep="")))
    sibs <- sapply( genotypes, FUN=function(x) paste(x, genotypes, sep="") )
    Pi <- matrix(data = NA, nrow = length(Pgametes)^2, ncol = length(Pgametes)^2)
    for(i in 1:ncol(sibs)) {
      for(j in 1:nrow(sibs)) {
        Pi[i,j] <- ((2*ploidy)-length(unique(strsplit(sibs[i,j], "")[[1]])))/ploidy
      }
    }
    sib.names <- as.vector(t(outer(Pgametes, Qgametes, paste, sep=":")))
    colnames(Pi) <- rownames(Pi) <- sib.names
    
    # dimnames(Z)[[1]] <- sib.names #DOUBLE CHECK IF IT'S CORRECT
    
    G <- array(data = NA, dim = c(nind, nind, nmrk), dimnames = list(c(ind.names), c(ind.names), c(mrk.names)))
    for(m in 1:nmrk) {
      ## G[,,m] <- t(Z[sib.names,m,])%*%Pi%*%Z[sib.names,m,]
      G[,,m] <- (t(Z[,m,])%*%Z[,m,])/ploidy
    }
    
    nphe <- dim(pheno)[2]
    ## Check for matching individual names
    if (length(which(rownames(pheno) %in% dimnames(G)[[1]])) == 0) stop("Individual names between genotype and phenotype data do not match. Please check your datasets and try again.")
    pheno.new <- as.matrix(pheno[which(rownames(pheno) %in% dimnames(G)[[1]]),])
    rownames(pheno.new) <- rownames(pheno)[which(rownames(pheno) %in% dimnames(G)[[1]])]
    colnames(pheno.new) <-colnames(pheno)
    if(!is.null(weights)) {
      weights.new <- as.matrix(weights[which(rownames(weights) %in% dimnames(G)[[1]]),])
      rownames(weights.new) <- rownames(weights)[which(rownames(weights) %in% dimnames(G)[[1]])]
    } else {
      weights.new <- NULL
    }

    ## Copying Z to X since Z already contains the homolog probabilities
    X = Z
  } else G <- Pi <- Z <- X <- NULL
  
  ######### DOSAGE
  
  if(!is.null(geno.dose)) {
    
    nlgs <- sum(!is.na(unique(geno.dose$sequence)))
    lgs <- vector("list", nlgs)
    lgs.all <- vector("list", nlgs)
    sel.mrk <- vector("list", nlgs)
    lgs.size <- numeric(nlgs)
    lgs.nmrk <- numeric(nlgs)
    for(c in 1:nlgs) {
      sel.mrk[[c]] <- rownames(geno.dose$geno.dose)[which(geno.dose$sequence == c)]
      if (length(geno.dose$sequence.pos)>1) lgs.all[[c]] <- geno.dose$sequence.pos[sel.mrk[[c]]] else lgs.all[[c]] <- c(1:length(sel.mrk[[c]]))
      lgs[[c]] <- c(lgs.all[[c]][c(1:length(lgs.all[[c]]))]) #fixme: ta gerando tamanho menor no exemplo do russet
      lgs.size[c] <- last(lgs[[c]])
      names(lgs)[[c]] <- names(lgs.all)[[c]] <- unique(geno.dose$sequence)[c]
      lgs.nmrk[c] <- length(lgs[[c]])
    }
    cum.size <- c(0, cumsum(lgs.size))
    cum.nmrk <- c(0, cumsum(lgs.nmrk))
    nmrk <- length(unlist(sel.mrk))
    
    dsg.names <- rep(NA, ploidy+1)
    for(d in 0:ploidy) dsg.names[d+1] <- paste(paste(rep("A", ploidy-d), collapse = "", sep=""), paste(rep("B", d), collapse = "", sep=""), collapse = "", sep="")
    # ordered <- order(geno.dose$geno$ind)
    # geno.dose$geno <- geno.dose$geno[order(geno.dose$geno$ind),]
    Z.dose <- array(0, dim=c(nmrk, geno.dose$n.ind, ploidy+1), dimnames = list(unlist(sel.mrk), colnames(geno.dose$geno.dose), dsg.names))
    for(d in 0:ploidy) {
      Z.dose[,,d] <- as.matrix((geno.dose$geno.dose[unlist(sel.mrk),] == d) + 0)
    }
    dim(Z.dose)
    Z.dose <- aperm(Z.dose, c(3,1,2))
    dim(Z.dose)
    str(Z.dose)
    nmrk <- dim(Z.dose)[2]
    nind <- dim(Z.dose)[3]
    Pi.dose <- matrix(data = NA, nrow = length(dsg.names), ncol = length(dsg.names))
    colnames(Pi.dose) <- rownames(Pi.dose) <- dsg.names
    for(i in 1:ncol(Pi.dose)) {
      for(j in 1:nrow(Pi.dose)) {
        Pi.dose[i,j] <- sum(!is.na(pmatch(strsplit(colnames(Pi.dose)[i],"")[[1]], strsplit(rownames(Pi.dose)[j],"")[[1]])))/ploidy
      }
    }
    
    mrk.names <- dimnames(Z.dose)[[2]]
    ind.names <- dimnames(Z.dose)[[3]]
    
    G.dose <- array(data = NA, dim = c(nind, nind, nmrk), dimnames = list(c(dimnames(Z.dose)[[3]]), c(dimnames(Z.dose)[[3]]), c(dimnames(Z.dose)[[2]])))
    for(m in 1:dim(Z.dose)[[2]]) {
      G.dose[,,m] <- t(Z.dose[,m,])%*%Pi.dose%*%Z.dose[,m,]
    }
    # G.dose <- NULL
    
    nphe <- dim(pheno)[2]
    pheno.new <- as.matrix(pheno[which(rownames(pheno) %in% ind.names),])
    rownames(pheno.new) <- rownames(pheno)[which(rownames(pheno) %in% ind.names)]
    if(!is.null(weights)) {
      weights.new <- as.matrix(weights[which(rownames(weights) %in% ind.names),])
      rownames(weights.new) <- rownames(weights)[which(rownames(weights) %in% ind.names)]
    } else {
      weights.new <- NULL
    }
    
    X.dose <- as.matrix(geno.dose$geno.dose[unlist(sel.mrk),ind.names])
    X.dose[X.dose >= ploidy+1] <- NA
    
  } else G.dose <- Pi.dose <- Z.dose <- X.dose <- NULL
  
  ############
  if(verbose) {
    cat("Reading the following data: \n")
    cat("  Ploidy level:       ", ploidy, "\n", sep="")
    cat("  No. individuals:    ", nind, "\n", sep="")
    if (!is.null(G)) {
      cat("  No. linkage groups: ", nlgs, "\n", sep="")
      cat("  Step size:          ", step, " cM \n", sep="")
      cat("  Map size:           ", round(last(cum.size), 2), " cM (", last(cum.nmrk), " positions) \n", sep="")
    }
    if (!is.null(G.dose)) {
      cat("  No. chromosomes:    ", nlgs, "\n", sep="")
      cat("  Step size:          ", step, "\n", sep="")
      cat("  Genome size:        ", round(last(cum.size)/10e5, 2), " Mbp (", last(cum.nmrk), ") \n", sep="")
    }
    cat("  No. phenotypes:     ", nphe, "\n", sep="")
  }
  
  structure(list(ploidy = ploidy,
                 nlgs = nlgs,
                 nind = nind,
                 nmrk = nmrk,
                 nphe = nphe,
                 ind.names = ind.names,
                 mrk.names = mrk.names,
                 lgs.size = lgs.size,
                 cum.size = cum.size,
                 lgs.nmrk = lgs.nmrk,
                 cum.nmrk = cum.nmrk,
                 lgs = lgs,
                 lgs.all = lgs.all,
                 step = step,
                 pheno = matrix(as.numeric(pheno.new), nrow = nrow(pheno.new), ncol = ncol(pheno.new), dimnames = dimnames(pheno.new)),
                 weights = weights.new,
                 G = G,
                 Z = Z,
                 X = X,
                 Pi = Pi,
                 ## alleles.mat = alleles2,
                 G.dose = G.dose,
                 Z.dose = Z.dose,
                 X.dose = X.dose,
                 Pi.dose = Pi.dose
  ),
  class="qtlpoly.data")

  }
    
    else {

  # if (!is.null(pheno) && !is.null(geno.prob) || !is.null(geno.dose) ) {
  #   nphe <- dim(pheno)[2] 
  #   if (!is.null(geno.prob)) 
  #     if (length(which(rownames(pheno) %in% dimnames(geno.prob[[1]]$probs)[[3]])) == 0)
  #       stop("Names of individuals from both genotype and phenotype data do not match. Please, check data \n")
  # } 
  # else stop("Data is required \n")
  
  if(is.null(step)) step <- 1e-10
  
  ######### HAPLOTYPE DATA
  
  if(!is.null(geno.prob)) {
    nlgs <- length(geno.prob)
    probs <- vector("list", nlgs)
    lgs <- vector("list", nlgs)
    lgs.all <- vector("list", nlgs)
    lgs.size <- numeric(nlgs)
    lgs.nmrk <- numeric(nlgs)
    for(c in 1:nlgs) {
      sel.map <- which(!duplicated(geno.prob[[c]]$map))
      sel.probs <- geno.prob[[c]]$probs[,sel.map,]
      lgs.all[[c]] <- geno.prob[[c]]$map[sel.map]
      lgs[[c]] <- unique(c(lgs.all[[c]][!duplicated(floor(lgs.all[[c]]/step)*step)], last(lgs.all[[c]])))
      lgs.size[c] <- last(lgs[[c]])
      names(lgs[[c]]) <- names(lgs.all[[c]])[which(lgs.all[[c]] %in% lgs[[c]])]
      probs[[c]] <- sel.probs[,which(lgs.all[[c]] %in% lgs[[c]]),]
      lgs.nmrk[c] <- dim(probs[[c]])[2]
    }
    cum.size <- c(0, cumsum(lgs.size))
    cum.nmrk <- c(0, cumsum(lgs.nmrk))
    
    Z <- abind::abind(probs, along = 2); dim(Z)
    ind.names <- dimnames(Z)[[3]]
    mrk.names <- dimnames(Z)[[2]]
    nmrk <- dim(Z)[2]
    nind <- dim(Z)[3]
    
    if(double.reduction) {
      if(ploidy == 4) n.unique <- 1
      if(ploidy == 6) n.unique <- 2
      if(ploidy == 8) n.unique <- 3
      if(ploidy == 10) n.unique <- 4
      if(ploidy == 12) n.unique <- 5
    } else {
      n.unique <- ploidy/2
    }
    Palleles <- letters[1:ploidy]; length(Palleles)
    Pgametes <- gtools::combinations(length(Palleles), ploidy/2, Palleles, repeats.allowed = TRUE); dim(Pgametes)
    Punique <- apply(Pgametes, 1, unique); length(Punique)
    Pgametes <- apply(Pgametes[which(lapply(Punique, length) >= n.unique),], 1, paste, collapse=""); length(Pgametes)
    # Pgametes <- lapply(combn(Palleles, ploidy/2, simplify = FALSE), paste, collapse=""); length(Pgametes)
    Qalleles <- letters[(ploidy+1):(2*ploidy)]
    Qgametes <- gtools::combinations(length(Qalleles), ploidy/2, Qalleles, repeats.allowed = TRUE); dim(Qgametes)
    Qunique <- apply(Qgametes, 1, unique); length(Qunique)
    Qgametes <- apply(Qgametes[which(lapply(Qunique, length) >= n.unique),], 1, paste, collapse=""); length(Qgametes)
    # Qgametes <- lapply(combn(Qalleles, ploidy/2, simplify = FALSE), paste, collapse="")
    genotypes <- as.vector(t(outer(Pgametes, Qgametes, paste, sep="")))
    sibs <- sapply( genotypes, FUN=function(x) paste(x, genotypes, sep="") )
    Pi <- matrix(data = NA, nrow = length(Pgametes)^2, ncol = length(Pgametes)^2)
    for(i in 1:ncol(sibs)) {
      for(j in 1:nrow(sibs)) {
        Pi[i,j] <- ((2*ploidy)-length(unique(strsplit(sibs[i,j], "")[[1]])))/ploidy
      }
    }
    sib.names <- as.vector(t(outer(Pgametes, Qgametes, paste, sep=":")))
    colnames(Pi) <- rownames(Pi) <- sib.names
    
    # dimnames(Z)[[1]] <- sib.names #DOUBLE CHECK IF IT'S CORRECT
    
    G <- array(data = NA, dim = c(nind, nind, nmrk), dimnames = list(c(ind.names), c(ind.names), c(mrk.names)))
    for(m in 1:nmrk) {
      G[,,m] <- t(Z[sib.names,m,])%*%Pi%*%Z[sib.names,m,]
    }
    
    nphe <- dim(pheno)[2]
    ## Check for matching individual names
    if (length(which(rownames(pheno) %in% dimnames(G)[[1]])) == 0) stop("Individual names between genotype and phenotype data do not match. Please check your datasets and try again.")
    pheno.new <- as.matrix(pheno[which(rownames(pheno) %in% dimnames(G)[[1]]),])
    rownames(pheno.new) <- rownames(pheno)[which(rownames(pheno) %in% dimnames(G)[[1]])]
    
    if(!is.null(weights)) {
      weights.new <- as.matrix(weights[which(rownames(weights) %in% dimnames(G)[[1]]),])
      rownames(weights.new) <- rownames(weights)[which(rownames(weights) %in% dimnames(G)[[1]])]
    } else {
      weights.new <- NULL
    }
    
    alleles <- matrix(unlist(strsplit(dimnames(Z)[[1]], '')), ncol=(ploidy+1), byrow=TRUE)[,-c((ploidy/2)+1)]
    X <- array(data = NA, dim = c(nind, (ploidy*2), nmrk), dimnames = list(c(ind.names), letters[1:(ploidy*2)], c(mrk.names)))
    ## Replacing slow part with updated code
    ## Creating incidence matrix to transition from genotypes to homologs
    alleles2 = matrix(0, length(genotypes), ploidy*2)
    ## Creating dummy variables to associate genotypes with alleles
    for(i in 1:nrow(alleles)){
      for(j in 1:ncol(alleles)){
        alleles2[i,match(alleles[i,j], letters)]=1
      }
    }
    ## Getting homolog probabilities
    for(m in 1:nmrk) {
      X[,,m] <- t(Z[sib.names,m,])%*%alleles2
    }
    ## Replaced slow code starts here
    ## for(m in 1:nmrk) {
    ##   for(i in 1:nind) {
    ##     a <- vector("list", (ploidy*2))
    ##     names(a) <- letters[1:(ploidy*2)]
    ##     for(j in 1:(ploidy*2)) {
    ##       a[[j]] <- which(alleles == letters[j], arr.ind = TRUE)[,1]
    ##       a[[j]] <- sum(Z[,m,i][Reduce(intersect, list(a[[j]]))])
    ##     }
    ##     X[i,,m] <- unlist(a)
    ##   }
    ## }
    ## Finished replacing slow code
  } else G <- Pi <- Z <- X <- NULL
  
  ######### DOSAGE
  
  if(!is.null(geno.dose)) {
    
    nlgs <- sum(!is.na(unique(geno.dose$sequence)))
    lgs <- vector("list", nlgs)
    lgs.all <- vector("list", nlgs)
    sel.mrk <- vector("list", nlgs)
    lgs.size <- numeric(nlgs)
    lgs.nmrk <- numeric(nlgs)
    for(c in 1:nlgs) {
      sel.mrk[[c]] <- rownames(geno.dose$geno.dose)[which(geno.dose$sequence == c)]
      if (length(geno.dose$sequence.pos)>1) lgs.all[[c]] <- geno.dose$sequence.pos[sel.mrk[[c]]] else lgs.all[[c]] <- c(1:length(sel.mrk[[c]]))
      lgs[[c]] <- c(lgs.all[[c]][c(1:length(lgs.all[[c]]))]) #fixme: ta gerando tamanho menor no exemplo do russet
      lgs.size[c] <- last(lgs[[c]])
      names(lgs)[[c]] <- names(lgs.all)[[c]] <- unique(geno.dose$sequence)[c]
      lgs.nmrk[c] <- length(lgs[[c]])
    }
    cum.size <- c(0, cumsum(lgs.size))
    cum.nmrk <- c(0, cumsum(lgs.nmrk))
    nmrk <- length(unlist(sel.mrk))
    
    dsg.names <- rep(NA, ploidy+1)
    for(d in 0:ploidy) dsg.names[d+1] <- paste(paste(rep("A", ploidy-d), collapse = "", sep=""), paste(rep("B", d), collapse = "", sep=""), collapse = "", sep="")
    # ordered <- order(geno.dose$geno$ind)
    # geno.dose$geno <- geno.dose$geno[order(geno.dose$geno$ind),]
    Z.dose <- array(0, dim=c(nmrk, geno.dose$n.ind, ploidy+1), dimnames = list(unlist(sel.mrk), colnames(geno.dose$geno.dose), dsg.names))
    for(d in 0:ploidy) {
      Z.dose[,,d] <- as.matrix((geno.dose$geno.dose[unlist(sel.mrk),] == d) + 0)
    }
    dim(Z.dose)
    Z.dose <- aperm(Z.dose, c(3,1,2))
    dim(Z.dose)
    str(Z.dose)
    nmrk <- dim(Z.dose)[2]
    nind <- dim(Z.dose)[3]
    Pi.dose <- matrix(data = NA, nrow = length(dsg.names), ncol = length(dsg.names))
    colnames(Pi.dose) <- rownames(Pi.dose) <- dsg.names
    for(i in 1:ncol(Pi.dose)) {
      for(j in 1:nrow(Pi.dose)) {
        Pi.dose[i,j] <- sum(!is.na(pmatch(strsplit(colnames(Pi.dose)[i],"")[[1]], strsplit(rownames(Pi.dose)[j],"")[[1]])))/ploidy
      }
    }
    
    mrk.names <- dimnames(Z.dose)[[2]]
    ind.names <- dimnames(Z.dose)[[3]]
    
    G.dose <- array(data = NA, dim = c(nind, nind, nmrk), dimnames = list(c(dimnames(Z.dose)[[3]]), c(dimnames(Z.dose)[[3]]), c(dimnames(Z.dose)[[2]])))
    for(m in 1:dim(Z.dose)[[2]]) {
      G.dose[,,m] <- t(Z.dose[,m,])%*%Pi.dose%*%Z.dose[,m,]
    }
    # G.dose <- NULL
    
    nphe <- dim(pheno)[2]
    pheno.new <- as.matrix(pheno[which(rownames(pheno) %in% ind.names),])
    rownames(pheno.new) <- rownames(pheno)[which(rownames(pheno) %in% ind.names)]
    if(!is.null(weights)) {
      weights.new <- as.matrix(weights[which(rownames(weights) %in% ind.names),])
      rownames(weights.new) <- rownames(weights)[which(rownames(weights) %in% ind.names)]
    } else {
      weights.new <- NULL
    }
    
    X.dose <- as.matrix(geno.dose$geno.dose[unlist(sel.mrk),ind.names])
    X.dose[X.dose >= ploidy+1] <- NA
    
  } else G.dose <- Pi.dose <- Z.dose <- X.dose <- NULL
  
  ############
  if(verbose) {
    cat("Reading the following data: \n")
    cat("  Ploidy level:       ", ploidy, "\n", sep="")
    cat("  No. individuals:    ", nind, "\n", sep="")
    if (!is.null(G)) {
      cat("  No. linkage groups: ", nlgs, "\n", sep="")
      cat("  Step size:          ", step, " cM \n", sep="")
      cat("  Map size:           ", round(last(cum.size), 2), " cM (", last(cum.nmrk), " positions) \n", sep="")
    }
    if (!is.null(G.dose)) {
      cat("  No. chromosomes:    ", nlgs, "\n", sep="")
      cat("  Step size:          ", step, "\n", sep="")
      cat("  Genome size:        ", round(last(cum.size)/10e5, 2), " Mbp (", last(cum.nmrk), ") \n", sep="")
    }
    cat("  No. phenotypes:     ", nphe, "\n", sep="")
  }
  
  structure(list(ploidy = ploidy,
                 nlgs = nlgs,
                 nind = nind,
                 nmrk = nmrk,
                 nphe = nphe,
                 ind.names = ind.names,
                 mrk.names = mrk.names,
                 lgs.size = lgs.size,
                 cum.size = cum.size,
                 lgs.nmrk = lgs.nmrk,
                 cum.nmrk = cum.nmrk,
                 lgs = lgs,
                 lgs.all = lgs.all,
                 step = step,
                 pheno = matrix(as.numeric(pheno.new), nrow = nrow(pheno.new), ncol = ncol(pheno.new), dimnames = dimnames(pheno.new)),
                 weights = weights.new,
                 G = G,
                 Z = Z,
                 X = X,
                 Pi = Pi,
                 alleles.mat = alleles2,
                 G.dose = G.dose,
                 Z.dose = Z.dose,
                 X.dose = X.dose,
                 Pi.dose = Pi.dose
  ),
  class="qtlpoly.data")

 }
  
}

## Support function from mappoly
imf_h <- function(r) {
  r[r >= 0.5] <- 0.5 - 1e-14
  -50 * log(1 - 2 * r)
}
