# Internal functions from varComp package
# This set of functions was originally written by Long Qu and implemented in the 'varComp' package.
# Since the package was archived upon the author's request and is no longer available on CRAN (\url{https://cran.r-project.org/package=varComp}), we adapted part of its functions under the terms of the GPL-3 license.

#' @importFrom nlme fixef fixed.effects
anova.varComp <-
function (object, ..., test = c("KR", 'Satterthwaite'), L)
{
	cl=match.call(expand.dots=FALSE)
		diffLmat=function(qr0, X1)
		{
			z=(diag(1, nrow(X1))- tcrossprod(qr.Q(qr0)))%*%X1
			qrz=qr(zapsmall(z))    ### WARNING: this will depends on getOption("digits")
			 pivot <- qrz$pivot
			 oo <- order(pivot)
			 qr.R(qrz)[seq_len(qrz$rank), oo, drop=FALSE]
		}
	if('...'%in%names(cl)){
		if(!missing(L)) message("'L' will be ignored when '...' is given.")
		ddd.names=sapply(cl[['...']], deparse)
		ddd = list(...)
		nObj=length(ddd)
		nObj=nObj + 1L
		ddd[[nObj]] = object
		names(ddd) = c(ddd.names, deparse(substitute(object)))
		if(!all(sapply(ddd, inherits, what='varComp'))) stop("All `...` objects need to inherit the `varComp` class.")
		if(diff(range(sapply(ddd, nobs))) != 0 ) stop("The number of observations needs to be the same in all model fits.")
		n=nobs(ddd[[1L]])
		
		Xmats = lapply(ddd,  model.matrix, what='fixed')
		qrs = lapply(Xmats, qr)
		rks = sapply(qrs, '[[', 'rank')
		ork=order(rks)
		rks=rks[ork]
		
		ddd=ddd[ork]
		Xmats=Xmats[ork]
		qrs=qrs[ork]
		# (`raw.F`=drop(Fstat), `scale.F` = scale.overall, df1=rk, df2=F.ddf, `Pr(>F)`=drop(F.p))
		
		nas=rep(NA_real_, nObj)
		ans = data.frame(`F value`=nas, `Scale` = nas, numDF=nas, denDF=nas, `Pr(>F)` = nas, check.names=FALSE)
		rownames(ans) = names(ddd)
		
		if('model'%in%names(ddd[[1L]])){
			hasInt = attr(attr(ddd[[1L]]$model, 'terms'), 'intercept') == 1
		}else{
			ones = rep(1, n)
			hasInt = ( mean((tcrossprod(qr.Q(qrs[[1]]))%*%ones - ones)^2) < sqrt(.Machine$double.eps) )
		}
		if(hasInt && rks[1L] > 1L){
			Lmat = diffLmat(qr(rep(1,n)), Xmats[[1L]])
		}else{
			Lmat = diag(1, ncol(Xmats[[1L]]))
		}
		for(m in seq_len(nObj)){
			tmp = fixef(ddd[[m]], Lmat = Lmat, test=test)
			ans[m, ] = attr(attr(tmp, 'anova'), 'Overall')
			if(m < nObj)	Lmat = diffLmat(qrs[[m]], Xmats[[m+1L]])
		}
		return(ans)
	}
	
	### one object is given
	if(!missing(L)){
		tmp = fixef(object, Lmat=L, test=test)
	}else tmp = fixef(object, test=test)
	tmp.aov=attr(tmp, 'anova')
	ans = data.frame(`F value` = (tmp.aov[,'t value']*tmp.aov[,'Scale'])^2,
					 `Scale` = tmp.aov[,'Scale']^2, 
					 `numDF` = 1, 
					 `denDF` = tmp.aov[,'Df'], 
					 `Pr(>F)` = tmp.aov[, 'Pr(>|t|)'],
					 check.names=FALSE
					)
	if(missing(L)){
		message("Current implementation will test each fixed effect parameter separately when only one `varComp` object is provided.")
		X=model.matrix(object)
		n = nobs(object)
		if('model'%in%names(object)){
			hasInt = attr(attr(object$model, 'terms'), 'intercept') == 1
		}else{
			ones = rep(1, n)		
			hasInt = ( mean((tcrossprod(qr.Q(qr(X)))%*%ones - ones)^2) < sqrt(.Machine$double.eps) )
		}	
		if(hasInt){
			tmp = fixef(object, Lmat = diffLmat(qr(rep(1,n)), X), test=test)
			tmp.aov=attr(tmp, 'anova')
		}
	}
	ans = rbind(ans, attr(tmp.aov, 'Overall'))
	if(!missing(L) && !is.null(rownames(L))) rownames(ans) = c(rownames(L), 'Overall')
	return(ans)

}

#' @importFrom MASS ginv
coef.varComp <-
function(object, what=c('fixed','beta','random','varComp','var.ratio','tau'), ...)
{
# S3 method for getting coef's of class varComp object. 
# what:   'fixed' or 'beta': Return fixed effect parameters
#         'random' or 'varComp': Return REML estimates. 
#         'var.ratio' or 'tau': Return ratios of REML estimates to error variance. 
  what=match.arg(what)
  what=switch(what, fixed='beta', random='varComp', var.ratio='tau', what)

  if(what=='beta'){
	if(FALSE){
    call0=object$call
    call.args=names(call0)
    Y=if(is.null(object[['Y']])) eval(call0[['Y']], envir=object$frame) else object[['Y']]
	}
	if('fixef'%in%names(object)) return(object$fixef)
	if('Y'%in%names(object)) {
		Y = object$Y
	}else if ('model'%in%names(object)){
		Y = model.response(object$model)
	}else stop("response variable is not recored.")
    X=model.matrix(object, what='fixed')
    
	if(ncol(X)>0L){
		this.V=vcov(object, what='Y')
		this.Vbet=ginv(crossprod(X, solve(this.V,X)))
		this.bet=drop(this.Vbet%*%crossprod(X, solve(this.V, Y)))
		names(this.bet)=colnames(X)
	}else this.bet=numeric(0L)
    this.bet
  }else if (what=='varComp'){
    c(object$varComps, error=object$sigma2)
  }else if (what=='tau'){
    object$parms
  }
}

# fixef.varComp <-
# function(object, ...) coef(object, what='fixed', ...)
fixef.varComp <-
function(object, Lmat, alpha=.05, test=c('KR', 'Satterthwaite'), ...)
{
## S3 method for reporting fixed effect parameters from class varComp. 
# object: varComp object
# Lmat: The same as in satterth, default to all parameters (i.e., identity Lmat)
# alpha: Type I error level requested. 
	if(is.null(test[1L])) return(coef(object, 'fixed'))
	test=match.arg(test)

	if('Y'%in%names(object)){
		Y = object$Y
	}else if ('model'%in%names(object)){
		Y = model.response(object$model)
	}else stop("response variable is not recored.")

  
  X=model.matrix(object, what='fixed')
  K=model.matrix(object, what='K')[object$parms>0]
  if(missing(Lmat)) {Lmat=diag(1, ncol(X)); rownames(Lmat)=colnames(X)}
  if(!is.matrix(Lmat))	Lmat = matrix(Lmat, ncol=length(Lmat), dimnames = list("", names(Lmat)))
  
  if(!is.null(colnames(Lmat)) && !is.null(colnames(X))){
	if(ncol(Lmat) < ncol(X)){
		L=matrix(0, nrow(Lmat), ncol(X))
		colnames(L) = colnames(X)
		rownames(L) = rownames(Lmat)
		L[, colnames(Lmat)]=Lmat
		Lmat=L
	}else if (ncol(Lmat) == ncol(X)){
		Lmat = Lmat[, colnames(X), drop=FALSE]
	}else stop("`Lmat` has more columns than the number of fixed effect parameters.")	
  }
  if(ncol(Lmat) != ncol(X)) stop("`Lmat` has incorrect number of columns.")
  if(is.null(rownames(Lmat))) rownames(Lmat) = rep('', nrow(Lmat))
  
  this.V=vcov(object, what='Y')
  this.Vbet=if(ncol(X)>0L) solve(crossprod(X, solve(this.V,X))) else matrix(NA_real_, 0L, 0L)
  this.bet=this.Vbet%*%crossprod(X, solve(this.V, Y))
  
  est=drop(Lmat%*%this.bet)
  LVL=tcrossprod(Lmat%*%this.Vbet, Lmat)
  svd.LVL=if(length(LVL)>0L) svd(LVL) else list(d=numeric(0L), u=matrix(NA_real_, 0L,0L), v=matrix(NA_real_,0L,0L))
  rk=sum(svd.LVL$d>sqrt(.Machine$double.eps))
    if(test=='Satterthwaite'){
	  
	  ses=sqrt(diag(LVL))
	  t.dfs=numeric(nrow(Lmat))
	  for(i in seq_len(nrow(Lmat))) 
		t.dfs[i]=satterth.varComp(object=object, Lmat[i,,drop=FALSE], Vbet=this.Vbet, 
						  #svd.VLbet=list(d=svd.LVL$d[i], u=svd.LVL$u[i,,drop=FALSE], v=svd.LVL$v[i,,drop=FALSE]),
						  X=X, K=K, V=this.V )
	#  t.dfs=apply(Lmat, 1L, satterth, object=object, Vbet=this.Vbet, X, K=K, V=this.V)  ## X has naming conflict
		scaleF=rep(1, nrow(Lmat)); scale.overall = if(nrow(Lmat)>0L) 1 else numeric(0L)
	  #debug(satterth.varComp)
	  LVLI=tcrossprod(sweep(svd.LVL$v, 2, ifelse(svd.LVL$d>sqrt(.Machine$double.eps), 1/svd.LVL$d, 0), '*'), svd.LVL$u) 
	  Fstat=if(nrow(Lmat)>0L) crossprod(est, LVLI%*%est)/rk  else numeric(0L)
	  F.ddf=satterth.varComp(object, Lmat=Lmat, Vbet=this.Vbet, svd.VLbet=svd.LVL, K=K, V=this.V, X)
	  F.p=pf(Fstat, rk, F.ddf, lower.tail=FALSE)
	}else if (test=='KR') {
		scaleF = ses = t.dfs =numeric(nrow(Lmat))
		for(i in seq_len(nrow(Lmat))){
			tmp = KR.varComp(object=object, Lmat=Lmat[i,,drop=FALSE], Vbet=this.Vbet, X = X, K = K, V = this.V)
			t.dfs[i]=tmp[[1L]]
			ses[i] = sqrt(diag(Lmat[i,,drop=FALSE]%*%attr(tmp, 'vcov.beta')%*%t(Lmat[i,,drop=FALSE])))
			scaleF[i] = attr(tmp, 'Scale')
		}
		tmp=KR.varComp(object=object, Lmat=Lmat, Vbet=this.Vbet, svdVLbet = svd.LVL, X = X, K = K, V = this.V)
		scale.overall = attr(tmp, 'Scale')
		Fstat = attr(tmp, 'F value')
		F.ddf = if(length(tmp)>0L) tmp[[1L]] else numeric(0L)
		F.ndf = attr(tmp, 'numDF')
		if(length(F.ndf)>0L && F.ndf != rk) warning("numerical instability detected on the rank of `Lmat`")
		F.p = attr(tmp, 'Pr(>F)')
	}
	ll = est - ses/sqrt(scaleF) * qt(1-alpha/2, t.dfs) 
	ul = est + ses/sqrt(scaleF) * qt(1-alpha/2, t.dfs) 
	tstats = est / ses 
	t.pval = pt(abs(tstats)*sqrt(scaleF), t.dfs, lower.tail=FALSE)*2
	

  ans = est
  attrs = cbind(`Std. Error`=ses, Lower=ll, Upper=ul, `t value`=tstats*sqrt(scaleF), `Scale` = sqrt(scaleF), Df=t.dfs, `Pr(>|t|)`=t.pval)
  attr(attrs, 'Overall')=cbind(`F value`=drop(Fstat*scale.overall), `Scale` = scale.overall, numDF=if(nrow(Lmat)>0L)rk else integer(0L), denDF=F.ddf, `Pr(>F)`=drop(F.p))
  names(est) = rownames(attrs)=rownames(Lmat)
  rownames(attr(attrs, 'Overall')) = if(nrow(Lmat)>0L) 'Overall' else character(0L)
  attr(ans, 'anova') = attrs
  class(ans) = 'varCompFixEf'
  ans
}

print.varCompFixEf = function(x,...)
{
	cat("Individual fixef effect estimates:\n")
	tmp=cbind(Estimate=x, attr(x, 'anova'))
	print(tmp)
	cat("\nOverall fixed effect contrast:\n")
	print(attr(attr(x, 'anova'), 'Overall'))
	cat('\n')
	invisible(x)
}
colon = as.symbol(':')

callList2terms=function(cl)
{	lst=cl
	if(!is.list(lst)) return(deparse(lst))
	
	first=lst[[1L]]
	if(length(lst)==1L) return(as.character(lst[[1L]]))
	
	if(identical(first, colon) && length(lst)==3L ){
		return(c(Recall(as.list(lst[[2L]])), Recall(as.list(lst[[3L]]))))
	}
	tmp=sapply(lst[-1], as.character)
	tmpn = names(lst[-1])
	if(!is.null(tmpn)) for(i in which(nchar(tmpn)>0L)) tmp[i] = paste(tmpn[i], '=', tmp[i])
	paste0(as.character(first), "(", paste(tmp, collapse=', '), ')')
}

splitTerm=function(term)
{
	tcall=as.list(as.call(parse(text=term))[[1]])
	unique(callList2terms(tcall))
}

sortTerm=function(term, priority)
{
	singleTerms = splitTerm(term)
	if(missing(priority)) singleTerms=sort(singleTerms) else{
		idx=singleTerms%in%priority
		singleTerms=c(sort(singleTerms[idx]), sort(singleTerms[!idx]))
	}
	paste0(singleTerms, collapse=':')
}

normalizeTrace=function(x){x=as.matrix(x); x/mean(diag(x), na.rm=TRUE)}

Minkowski=function(x, p=1) 1-as.matrix(dist(x, method='minkowski', p=p)) * .5 / max(1, ncol(x))^(1/p) 

IBS=function(x)  1 - as.matrix(dist(x, method='manhattan') * .5 /max(1, ncol(x)) )  ## dist does scaling in the presence of missing values

Lin0=function(x) normalizeTrace(tcrossprod(x)/max(1,ncol(x)))
Quad1=function(x) normalizeTrace((base::tcrossprod(x)+1)^2/max(1,ncol(x)))
# Intxn2=function(K1, K2) normalizeTrace(K1*K2)
Polyk=function(x,c=0,d=1) normalizeTrace((base::tcrossprod(x)+c)^d)

cholRoot=function(x)
{
	R=suppressWarnings(chol(x, pivot=TRUE))
	oo <- order(attr(R, 'pivot'))
	rk = attr(R, 'rank')
	t(R[seq_len(rk), oo, drop=FALSE])
}

ibs=function(x)	cholRoot(IBS(x))
lin0=function(x)	cholRoot(Lin0(x))
quad1=function(x) cholRoot(Quad1(x))
polyk=function(x,c=0,d=1) cholRoot(Polyk(x,c,d))
minkowski=function(x,p=1) cholRoot(Minkowski(x,p))

#' @importFrom Matrix nearPD
KR.varComp=function(object, Lmat, Vbet, svd.VLbet, X, K, V, ...)
{
	bhat=coef(object, 'fixed'); 
	if(length(bhat)== 0L){
		return(structure(numeric(0L), names=character(0L), numDF = numeric(0L), Scale=numeric(0L), `F value`=numeric(0L), `Pr(>F)`=numeric(0L), vcov.beta = matrix(NA_real_, 0L, 0L)))
	}
	#if(missing(K)) 
		K= model.matrix(object, what='K')
	if(missing(X)) X= model.matrix(object, what='X')
  if(missing(Lmat)) {Lmat=diag(1, ncol(X)); rownames(Lmat)=colnames(X)}
  if(!is.matrix(Lmat))	Lmat = matrix(Lmat, 1L, dimnames = list("", names(Lmat)))
  
  if(!is.null(colnames(Lmat)) && !is.null(colnames(X))){
	if(ncol(Lmat) < ncol(X)){
		L=matrix(0, nrow(Lmat), ncol(X))
		colnames(L) = colnames(X)
		rownames(L) = rownames(Lmat)
		L[, colnames(Lmat)]=Lmat
		Lmat=L
	}else if (ncol(Lmat) == ncol(X)){
		Lmat = Lmat[, colnames(X), drop=FALSE]
	}else stop("`Lmat` has more columns than the number of fixed effect parameters.")	
  }
  if(ncol(Lmat) != ncol(X)) stop("`Lmat` has incorrect number of columns.")
  if(is.null(rownames(Lmat))) rownames(Lmat) = rep('', nrow(Lmat))
  if(missing(V)) V= vcov(object, 'Y'); VI=solve(V)
	if(missing(Vbet)) Vbet=vcov(object, what='beta', beta.correction=FALSE)
	Phi = Vbet
	p=length(bhat)
	ell=qr(Lmat)$rank
	r=length(coef(object, 'varComp'))
	J=r-1
	n=nrow(V)
	
	w=vcov.varComp(object, 'varComp', drop=FALSE)
	w=as.matrix(nearPD(w)$mat)
	tmp=which(object$parms==0)
	w[tmp,]=w[,tmp]=0

	
	VIX=solve(V, X)
	# Phi = solve(crossprod(X, VIX))

	Q=array(NA_real_, dim=c(p,p,r,r))
		for(i in seq_len(J)){
			for(j in seq_len(J))
				Q[,,i,j]=crossprod(VIX, K[[i]]%*%VI%*%K[[j]]%*%VIX)
			Q[,,i,r]=crossprod(VIX, K[[i]]%*%VI%*%VIX)
		}
		for(j in seq_len(J)) Q[,,r,j]=crossprod(VIX, VI%*%K[[j]]%*%VIX)
		Q[,,r,r]=crossprod(VIX, VI%*%VIX)
	P=array(NA_real_, dim=c(p,p,r))
		for(i in seq_len(J)) P[,,i]=-crossprod(VIX, K[[i]]%*%VIX)
		P[,,r]=-crossprod(VIX)

	LambdaTilde=matrix(0, p,p)
		for(i in seq_len(r))
			for(j in seq_len(r))
				LambdaTilde=LambdaTilde + w[i,j]*(Q[,,i,j]-P[,,i]%*%Phi%*%P[,,j])
		LambdaTilde=Phi%*%LambdaTilde%*%Phi

	Rstar=0
	Phi_A = Vbet + 2 * LambdaTilde - Rstar

	if(ell==0L) return(structure(NA_real_, names='denDF', numDF = ell, Scale=NA_real_, `F value`=NA_real_, `Pr(>F)`=NA_real_, vcov.beta = Phi_A))
	
	Lmat=t(Lmat)
	F=drop( crossprod(bhat, Lmat%*%solve(crossprod(Lmat, Phi_A%*%Lmat))%*%crossprod(Lmat, bhat))) / ell
	Theta = Lmat%*%solve(crossprod(Lmat, Phi%*%Lmat), t(Lmat))
	PTP=Phi%*%Theta%*%Phi
	
	A1=0;	for(i in seq_len(r)) for(j in seq_len(r)) A1=A1+w[i,j]*sum(PTP*P[,,i])*sum(PTP*P[,,j])
	A2=0;	for(i in seq_len(r)) for(j in seq_len(r)) A2=A2+w[i,j]*sum(diag(PTP%*%P[,,i]%*%PTP%*%P[,,j]))
	B=(A1+6*A2)/2/ell
	
	g=((ell+1)*A1-(ell+4)*A2)/(ell+2)/A2; c=3*ell+2*(1-g); 
	d1=g/c; d2=(ell-g)/c; d3=(ell-g+2)/c
	Estar=1/(1-A2/ell); Vstar=2/ell*(1+d1*B)/(1-d2*B)^2/(1-d3*B)
	rhoTilde=Vstar/2/Estar/Estar
	m=4+(2+ell)/(ell*rhoTilde-1); lambda=m/Estar/(m-2)
	
	structure(m, names='denDF', numDF = ell, Scale=lambda, `F value`=F*lambda, `Pr(>F)`=pf(lambda*F, ell, m, lower.tail=FALSE), vcov.beta = Phi_A)
}

logLik.varComp=function(object, ...)
{
	REML=object$control$REML
	if(!isTRUE(REML)) .NotYetImplemented()
	val=object$PREML
	vc=coef(object, "varComp")
	p=sum(vc>object$control$boundary.eps)
	
    N = attr(val, "nall") <- nobs(object)
    attr(val, "nobs") <- N - REML * length(coef(object, 'fixed'))
    attr(val, "df") <- p #+ length(coef(object[["modelStruct"]])) + 1
    class(val) <- "logLik"	
	val
}

#' @importFrom quadprog solve.QP
#' @importFrom Matrix nearPD
minque <-
function(y, varcov, start=rep(0, length(k)), lower.bound=-Inf, restricted=TRUE)
{ 
	k=varcov; varRatio=start
  if(!isTRUE(restricted)) .NotYetImplemented()
  stopifnot(is.numeric(y) && is.list(k) && is.numeric(varRatio))
  nK=length(k)
  varRatio=c(rep(varRatio, length=nK),1)
  k[[nK+1L]]=diag(1,nrow(k[[1L]]))  
  
  LI=t(backsolve(chol(Reduce(`+`,mapply(`*`,k,varRatio,SIMPLIFY=FALSE))), k[[nK+1]]))
  LIk=lapply(lapply(k, `%*%`, LI), crossprod, LI)
  LIy=LI%*%y

  S=outer(LIk,LIk,function(a,b)sapply(mapply(`*`,a,b,SIMPLIFY=FALSE),sum))
  u=sapply(lapply(LIk,'%*%',LIy),crossprod,LIy)

  if(lower.bound<0 && is.infinite(lower.bound)){
    ans0=solve(S, u) ## FIXME: add error handler
  }else{
    qp.rslt=tryCatch(solve.QP(Dmat=crossprod(S), dvec=crossprod(S,u), Amat=rbind(diag(1,nK),0), bvec=rep(lower.bound,nK), meq=0L, factorized=FALSE), 
      error=function(e)  {
        if(e$message=='matrix D in quadratic function is not positive definite!'){
          solve.QP(Dmat=as.matrix(nearPD(crossprod(S))$mat), dvec=crossprod(S,u), Amat=rbind(diag(1,nK),0), bvec=rep(lower.bound,nK), meq=0L, factorized=FALSE) 
        }else{
          e
        }
      }
    )
    ans0=qp.rslt$solution
  }
  ans=ans0[-nK-1L]/ans0[nK+1L]
  ans
}
model.matrix.varComp <-
function(object, what=c('fixed','random', 'varcov', 'X', 'K', 'Z'), ...)
{
# S3 method for getting design matrix of object of class varComp. 
# what: 'fixed' or 'X': Fixed effect design matrix. 
#       'random' or 'Z': random effect design matrix. NOTE: This is an equivalent design matrix, not necessarily the one supplied to varComp. 
#       'K': Returns the kernal matrix. 

  what=match.arg(what)
  what=switch(what, random='Z', fixed='X', varcov='K', what)
  if(what=='X'){
	if('X'%in%names(object)) object$X else stop("X matrix is not recored")
  }else if (what=='K' || what=='Z'){
	if(!('K'%in%names(object))) stop("K matrices are not recored")
	if(what=='Z') lapply(object$K, cholRoot) else object$K
  } 
}

pchibarsq=function(q, V, lower.tail=TRUE, log.p=FALSE)
{
  n=nrow(V)
  
  wts=wchibarsq(V) 
    
  ans=pchisq(q, 0, lower.tail=FALSE)*wts[1L]+
      pchisq(q, n, lower.tail=FALSE)*wts[n+1L]
  for(i in seq_len(n-1)){
    ans=ans+pchisq(q, i, lower.tail=FALSE)*wts[i+1L]
  }
  ans[q<=0]=1

  ans = if(isTRUE(lower.tail)) 1-ans else ans
  if(isTRUE(log.p)) log(ans) else ans
}

#' @importFrom mvtnorm pmvnorm
#' @importFrom utils combn
wchibarsq=function(V) ## weights
{## Ref: Kudo A. 1963. Biometrika, 50, pg 403  (page 414)
  seed=get.seed()
  stopifnot(is.matrix(V) && diff(dim(V))==0L && nrow(V)>0L)
  n=nrow(V)
  P=function(idx){## V, n, i
    pmvnorm(rep(0,n-i), rep(Inf, n-i), sigma=solve(V[-idx,-idx,drop=FALSE]))*
    pmvnorm(rep(0, i), rep(Inf, i), 
            sigma=V[idx,idx,drop=FALSE]-V[idx, -idx, drop=FALSE]%*%solve(V[-idx,-idx,drop=FALSE], V[-idx, idx, drop=FALSE])
           )
  }
  ans=numeric(n+1L)
  ans[1]=pmvnorm(rep(0,n),rep(Inf,n),sigma=solve(V))[[1]]
  ans[n+1L]=pmvnorm(rep(0,n), rep(Inf,n), sigma=V)[[1]]
  for(i in safeseq(1L, n-1L, by=1L)) ans[i+1]=sum(combn(n, i, P))
  attr(ans, 'seed')=seed
  ans
}

#' @importFrom quadprog solve.QP
#' @importFrom graphics plot abline
mchibarsq=function(V, order=1:2) # moments
{
  n=nrow(V)
  wts=wchibarsq(V)
  l2=log(2); lg=lgamma(1:n/2)
  sapply(order,function(ord) weighted.mean(c(0, exp(ord*l2+lgamma(ord+1:n/2)-lg)), wts))
}

print.varComp <-
function(x, ...)
{
	cat("Variance component model fit", '\n')
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
	if(inherits(x$fixef, 'varCompFixEf')) {
		print(x$fixef)
	}else{
		cat("\nFixed effect estimates:", '\n')
		print(coef(x, 'fixed'))
	}
	cat("\nVariance component estimates:", '\n')
	print(coef(x, 'varComp'))
	cat(sprintf("\nNumber of observations: %d\n", nobs(x)))
	invisible(x)
}

#' @importFrom nlme fixef fixed.effects
summary.varComp = function(object, ...)
{## FIXME: add more summary information
	object$fixef = fixef(object,...)
	class(object) = c('summary.varComp', class(object))
	object
}
print.varComp.test <-
function(x, ...)
{
#Printing method for the class `varComp.test`. 
  print(p.value(x))
}

p.value <-
function(x,...) UseMethod('p.value')

p.value.htest = function(x,...) x$p.value
p.value.default = p.value.htest

p.value.varComp.test <-
function(x,...)
{
#Extracting p-values from objects from varComp.test. 
#i.	x: an object of class varComp.test
  f=function(zz){
    if(inherits(zz, 'htest')) return(zz$p.value)
    if(is.list(zz)) return(unlist(sapply(zz, f)))
    NA
  }
  p=drop(unlist(sapply(x, f)))
  
  cn=function(zz){  
    if(inherits(zz, 'htest')) return(NULL)
    nn=names(zz)
    ans=character(0L)
    for(i in seq_along(nn)) {
      tmp=cn(zz[[i]])
      ans=c(ans, if(is.null(tmp)) nn[i] else paste(nn[i], tmp, sep='.'))
    }
    ans
  }
  
  names(p)=cn(x)
  p
}

p.value.varCompFixEf=function(x,...)
{
	Aov=attr(x, 'anova')
	ans=Aov[, 'Pr(>|t|)']
	attr(ans, 'Overall') = attr(Aov, 'Overall')[,'Pr(>F)']
	ans
}#satterth <-
#function(object, ...) UseMethod('satterth')

satterth.varComp <-
function(object, Lmat, Vbet, svd.VLbet, X, K, V, ...)
{
# S3 method for object of class varComp. 
# Lmat:   A matrix with each row being a linear combination of fixed effect parameters of interest. 
# Vbet:   vcov of fixed effect parameter estimates
# svd.VLbet: An optional svd object of vcov of Lmat%*%beta
# X:  Optional X design matrix
# K:  Optional kernal matrix
# V:  Optional vcov of response. 

	if(ncol(model.matrix(object, 'X')) == 0L) return(structure(numeric(0L), individual.df=numeric(0L)))
	
  VVC=vcov(object, what='varComp', drop=TRUE)
  if(missing(Vbet) ) Vbet=vcov(object, what='beta', beta.correction=FALSE)
  if(missing(svd.VLbet)){
    eig=eigen(tcrossprod(Lmat%*%Vbet, Lmat), TRUE)
  }else{
    eig=list(values=svd.VLbet$d, vectors=svd.VLbet$u)
  }
  if(missing(X)) X=model.matrix(object, what='fixed')
  if(missing(K)) K=model.matrix(object, what='K')[object$parms>0]
  if(!is.list(K)) K=list(K)
  if(missing(V)) V=vcov(object, what='Y')
  if(missing(Lmat)) {Lmat=diag(1, ncol(X)); rownames(Lmat)=colnames(X)}
  if(!is.matrix(Lmat))	Lmat = matrix(Lmat, 1L, dimnames = list("", names(Lmat)))
  
  if(!is.null(colnames(Lmat)) && !is.null(colnames(X))){
	if(ncol(Lmat) < ncol(X)){
		L=matrix(0, nrow(Lmat), ncol(X))
		colnames(L) = colnames(X)
		rownames(L) = rownames(Lmat)
		L[, colnames(Lmat)]=Lmat
		Lmat=L
	}else if (ncol(Lmat) == ncol(X)){
		Lmat = Lmat[, colnames(X), drop=FALSE]
	}else stop("`Lmat` has more columns than the number of fixed effect parameters.")	
  }
  if(ncol(Lmat) != ncol(X)) stop("`Lmat` has incorrect number of columns.")
  if(is.null(rownames(Lmat))) rownames(Lmat) = rep('', nrow(Lmat))
  
  q=sum(eig$values>=sqrt(.Machine$double.eps))
  l=crossprod(eig$vectors[,seq_len(q)], Lmat)
  
  nK=length(K)
  K[[nK+1L]]=diag(1, nrow(X))
  
  VIXUl=solve(V, tcrossprod(X%*%Vbet, l))
  
  ders=matrix(NA_real_, q, nK+1L)
  for(i in seq_len(nK+1L))
    ders[, i]=diag(crossprod(VIXUl, K[[i]]%*%VIXUl))
  
  vs=numeric(q)
  for(m in seq_len(q))
    vs[m]=2*eig$values[m]*eig$values[m]/crossprod(ders[m,], VVC%*%ders[m,])
  
  E=sum(vs/(vs-2)*(vs>2))
  ans=if(E>q) 2*E/(E-q) else 0
  attr(ans, 'individual.df')=vs  
  names(ans)='denDF'
  ans
}

# varCompSpecial=c('ibs','lin0','quad1','intxn2','am')
varComp.LinScore.LC12 <-
function(null.fit, X, Y, K, null, w, ...)
{ 
#Actual testing function for LC12. 
#i.	null.fit: object from varComp under the null. 
#ii.	X: the same as in varComp.test.
#iii.	Y: the same as in varComp.test. 
#iv.	K: the same as in varComp.test, after pre-/post-multiplying Z matrices, if applicable. 
#v.	null: the same as in varComp.test. 
#vi.	w: weights of scores (not used for the LC12 method). 
  if(!all(w==1)) warning('Weights are ignored in LC12 approximation method')  ## FIXME
  n=nrow(K[[1]])
  nNull=length(null)
  nK=length(K)
  if(nK!=1L+nNull){ # FIXME: implement
    return(structure(list(p.value=NA_real_),class='htest'))
  }
  sigma=diag(null.fit$sigma2,n)+Reduce('+', mapply('*', null.fit$varComps, K[null], SIMPLIFY=FALSE))
  sigma.ihalf=local({eg=eigen(sigma,TRUE); tcrossprod( sweep(eg$vec, 2L, sqrt(eg$val), '/'), eg$vec)})
  
  qrix=qr(sigma.ihalf%*%X)
  P01=sigma.ihalf%*%(diag(1,n)-tcrossprod(qr.Q(qrix)[,seq_len(qrix$rank),drop=FALSE]))%*%sigma.ihalf
  
  K3=K[-null][[1]]  ## FIXME: Check K3=Reduce('+', K[-null])
  SI=.5*crossprod(Y, P01%*%K3%*%P01%*%Y)

  Delta=matrix(NA_real_, nK, nK)
  Phi=numeric(nK)
    null1=c(null, nK+1)
    K[[nK+1]]=diag(1, n)
  for(i in 1:nK) {
    for(j in i:nK) Delta[i,j]=Delta[j,i]=sum(diag(P01%*%K[[null1[i]]]%*%P01%*%K[[null1[j]]]))
    Phi[i]=sum(diag(P01%*%K3%*%P01%*%K[[null1[i]]]))
  }
  nuI=.5*sum(diag(P01%*%K3%*%P01%*%K3))-.5*crossprod(Phi, solve(Delta, Phi))
  deltaI=.5*sum(diag(P01%*%K3))
  aI=nuI/2/deltaI
  gI=2*deltaI*deltaI/nuI
  pval=pchisq(SI/aI, gI, lower.tail=FALSE)
  ans=list(statistic=c(SI=SI), 
           p.value=pval, 
           alternative='greater', 
           parameter=c(scale=aI, df=gI), 
           null.value=structure(0, names='interaction variance component'), null.fit=null.fit, 
           method='Li and Cui (2012) Satterwaite-approximated Variance Component Test'
          )
  class(ans)='htest'
  ans
}

varComp.LinScore.LC12Boundary <-
function(null.fit, X, Y, K, null, ...)
{
#Boundary corrected LC12 method. Argument definitions are the same as in varComp.LinScore.LC12. 
  .Deprecated('varComp.test.aproximation.Satterthwaite')
  n=nrow(K[[1]])
  nNull=length(null)
  nK=length(K)
  stopifnot(nK==1L+nNull) ## FIXME
  sigma=diag(null.fit$sigma2,n)+Reduce('+', mapply('*', null.fit$varComps, K[null], SIMPLIFY=FALSE))
  sigma.ihalf=local({eg=eigen(sigma,TRUE); tcrossprod( sweep(eg$vec, 2L, sqrt(eg$val), '/'), eg$vec)})
  
  qrix=qr(sigma.ihalf%*%X)
  P01=sigma.ihalf%*%(diag(1,n)-tcrossprod(qr.Q(qrix)[,seq_len(qrix$rank),drop=FALSE]))%*%sigma.ihalf
  
  K3=K[-null][[1]]  ## FIXME: Check K3=Reduce('+', K[-null])
  SI=.5*crossprod(Y, P01%*%K3%*%P01%*%Y)

  Delta=matrix(NA_real_, nK, nK)
  Phi=mu12=S12=numeric(nK)
    null1=c(null, nK+1)
    K[[nK+1]]=diag(1, n)
  for(i in 1:nK) {
    for(j in i:nK) Delta[i,j]=Delta[j,i]=sum(diag(P01%*%K[[null1[i]]]%*%P01%*%K[[null1[j]]]))
    Phi[i]=sum(diag(P01%*%K3%*%P01%*%K[[null1[i]]]))
    mu12[i]=.5*sum(diag(P01%*%K[[null1[i]]]))
    S12[i]=.5*crossprod(Y, P01%*%K[[null1[i]]]%*%P01%*%Y)
  }
  nuI=.5*sum(diag(P01%*%K3%*%P01%*%K3))-.5*crossprod(Phi, solve(Delta, Phi))
  deltaI=.5*sum(diag(P01%*%K3))-crossprod(Phi, solve(Delta, mu12-S12))
  aI=nuI/2/deltaI
  gI=2*deltaI*deltaI/nuI
  pval=pchisq(SI/aI, gI, lower.tail=FALSE)
  ans=list(statistic=c(SI=SI), 
           p.value=pval, 
           alternative='greater', 
           parameter=c(scale=aI, df=gI), 
           null.value=structure(0, names='interaction variance component'), null.fit=null.fit, 
           method='Corrected Li and Cui (2012) Satterwaite-approximated Variance Component Test'
          )
  class(ans)='htest'
  ans
}


varComp.LinScore.approximation <-
function(approximation, ...)
{
#A Switcher that choose the approximation method for LinScore test when the null has additional variance components other than the error variance. 
#i.	approximation: A character scalar, specifying the mthod of approximation. 
#ii.	? Other arguments to be passed to the actual approximate testing function. 
  this.call=match.call()
  idx=which(names(this.call)=='approximation')
  stopifnot(length(idx)>0L && length(approximation)==1L)
  this.call=this.call[-idx]
  #this.call[[1]]=as.name(paste('varComp.LinScore', approximation, sep='.'))

  this.env=sys.frame(sys.parent(0L))
  parent.env(this.env)=parent.frame()

  fname=paste('varComp.LinScore', approximation, sep='.')
  toCall=get(fname, mode='function')
  eval(substitute({environment(zzz)=this.env}, list(zzz=as.name(fname))))

  this.call[[1]]=as.name(fname)
  eval(this.call)
  
#  idx=which(names(this.call)=='envir')
#  if(length(idx)>0L)this.call=this.call[-idx]
#  eval(this.call, envir=parent.frame()) # explicit envir=parent.frame() is needed!
#  this.env=sys.frame(sys.parent(0L))
#  bak.parent=parent.env(this.env)
#  parent.env(this.env)=parent.frame()
#  ans=eval(this.call)
#  parent.env(this.env)=bak.parent
#  ans
}

refugeEnvironment=new.env()
evalq(
	alt.fit <- control <- denom <- diag.1.n <- eigK <- k <- LI <- LIkLI <- LIy <- n <- negGrad <- negHess <- nK <- null <- null.fit <- nums <- nums2 <- numsPart <- tau.idx <- test <- tr1 <- tr2 <- y <- NULL
, refugeEnvironment)

evalq(
	infoFunc <- function()NULL
  , refugeEnvironment)
  
V=evalq(function(tau){ ## diag.1.n, k
diag.1.n + Reduce('+', mapply('*', tau, k, SIMPLIFY=FALSE))
}, refugeEnvironment)

updateLI=evalq(function(tau){
	if(isTRUE( attr(LI, 'tau')==tau )) return(LI)
	if(nK==1L){ 
		tmp= crossprod(eigK$tvec, 1/sqrt(eigK$val*tau+1) * eigK$tvec)
	}else
		tmp = t(backsolve(chol(V(tau)), diag.1.n))
	attr(tmp,'tau')=tau
	LI<<- tmp
}, refugeEnvironment)

updateLIkLI=evalq(function(){ #LI, k
	if(isTRUE( attr(LIkLI, 'tau')==attr(LI, 'tau'))) return(LIkLI)
	tmp=LIkLI
	for(j in 1:nK) tmp[[j]] = tcrossprod(LI%*%k[[j]], LI)
	attr(tmp, 'tau')=attr(LI, 'tau')
	LIkLI <<- tmp
}, refugeEnvironment)

updateLIy=evalq(function(){
	if(isTRUE( attr(LIy, 'tau')==attr(LI, 'tau'))) return(LIy)
	LIy<<- structure( LI%*%y, tau=attr(LI, 'tau'))
}, refugeEnvironment)

PREML=evalq(function(){ # LIy, LI, n
	drop(
		if(nK==1L){
			.5*(
				-n*log(crossprod(LIy)) - sum(log(eigK$val*attr(LI,'tau')+1)) -n -n*log(2*pi)-n*log(n)
			)
		}else
			.5*(
				-n*log(crossprod(LIy)) +sum(log(diag(LI)))*2 -n -n*log(2*pi)-n*log(n)
			)
	)
}, refugeEnvironment)

preprocPREML=evalq(function(tau){
	updateLI(tau)
	updateLIy(); 
}  , refugeEnvironment)

obj=evalq(function(tau){
	updateLI(tau)
	updateLIy()
	PREML()
}, refugeEnvironment)

obj2=evalq(function(ltau)obj(exp(ltau)), refugeEnvironment)

updateNumsPart=evalq(function(){ # LIkLI, LIy
	if(isTRUE( attr(numsPart, 'tau') == attr(LIkLI, 'tau') )) return( numsPart )
	numsPart<<- sapply(LIkLI, '%*%', LIy)
}, refugeEnvironment)

updateNums=evalq(function(){ # LIy, numsPart
	if(isTRUE( attr(nums, 'tau') == attr(numsPart, 'tau'))) return( nums )
	nums<<-drop(crossprod(LIy, numsPart))
}, refugeEnvironment)

updateDenom=evalq(function(){ # LIy
	if(isTRUE( attr(denom, 'tau') == attr(LIy, 'tau'))) return( denom )
	denom<<-drop(crossprod(LIy))
}, refugeEnvironment)

updateTr1=evalq(function(){ # LIkLI
	if(isTRUE( attr(tr1, 'tau') == attr(LIkLI, 'tau'))) return( tr1 )
	tr1<<-sapply(LIkLI, function(z)sum(diag(z)))
}, refugeEnvironment)

updateTr2=evalq(function(){ # LIkLI
	if(isTRUE( attr(tr2, 'tau') == attr(LIkLI, 'tau'))) return( tr2 )
	for(i in 1:nK) {
		for(j in i:nK){
			kij=LIkLI[[i]]%*%LIkLI[[j]]
			tr2[i,j]<<-tr2[j,i]<<-sum(diag(kij))
		}
	}
}, refugeEnvironment)

updateNums2=evalq(function(){ # numsPart
	if(isTRUE( attr(nums2, 'tau') == attr(numsPart, 'tau'))) return( nums2 )
	for(i in 1:nK) {
		for(j in i:nK){
			nums2[i,j]<<-nums2[j,i]<<-2*crossprod(numsPart[,i], numsPart[,j])
		}
	}
}, refugeEnvironment)

score=evalq(function(){ # LIy, LIkLI, n, denoms, nums
	.5*(n*nums/denom - tr1)
}, refugeEnvironment)

updateNegGrad=evalq(function(){
	negGrad <<- -score()
}, refugeEnvironment)

gradFunc=evalq(function(tau){ # this agrees with numerical gradient :)
	updateLI(tau)
	updateLIkLI()
	updateLIy(); 
	updateNumsPart();   
	updateNums(); updateDenom(); updateTr1()
	score()
}, refugeEnvironment)

OI=evalq(function(){ # nums2, nums, denom, n, tr2
	.5*(
		(nums2-outer(drop(nums), drop(nums))/denom ) /denom*n - tr2
	)
}, refugeEnvironment)

hess=evalq(function(tau){  # this agrees with numerical hessian :)
	updateLI(tau)
	updateLIkLI(); updateLIy(); updateNumsPart(); 
	updateNums(); updateDenom()
	updateNums2()
	updateTr2()
	-OI()
}, refugeEnvironment)

EI=evalq(function(){ # tr2, tr1, n
	.5/(n+2)*(n*tr2-outer(tr1,tr1))
}, refugeEnvironment)

AOI=evalq(function(){ # nums2, nums, denom
	.5*(
		(nums2/2-outer(nums, nums)/denom )/denom*n
	)
}, refugeEnvironment)

AEI=evalq(function(){ # nums2, tr1, denom
	.5/(n+2)*(
		n*n*nums2/denom/2 - outer(tr1, tr1)
	)
}, refugeEnvironment)

WAI=evalq(function(){ # nums2, denom, nums, tr1, n
	.5*.5/(n+1)*(
		n*n*nums2/denom - n*n*outer(nums, nums)/denom/denom - outer(tr1, tr1)
	)
}, refugeEnvironment)

updateNegHess=evalq(function(){
	negHess <<- infoFunc()
}, refugeEnvironment)

varComp.LinScore.Normal <-
function(	# null.fit, 
	all.scores, lin.form, infoMat, null, w, tr1, n, LIkLI, tau.idx,  ...)
{
#Actual normal approximation for LinScore. Argument definitions are the same as in varComp.LinScore.SSAS155, except the non.pd is not used. 
  nK=length(all.scores)
  nonNull=seq_len(nK)[-null]
  Phi=infoMat[null, nonNull, drop=FALSE]  
  Delta=infoMat[null, null, drop=FALSE]   

  var.nonNull.score=infoMat[nonNull, nonNull, drop=FALSE] - crossprod(Phi, solve(Delta, Phi))
  mean.nonNull.score= crossprod(Phi, solve(Delta, all.scores[null]))
  mean.obs.score=sum(mean.nonNull.score*w)
  var.obs.score=drop(crossprod(w, var.nonNull.score%*%w))
  
  mean.theo.score=0 
  var.theo.score=drop(crossprod(w, infoMat[nonNull, nonNull, drop=FALSE]%*%w))
  
  adjusted.obs.score=sqrt(1/var.obs.score)*(lin.form-mean.obs.score)
  
  pval=pnorm(adjusted.obs.score, lower.tail=FALSE)
  pval0=pnorm(lin.form/sqrt(var.obs.score), lower.tail=FALSE)
  structure(list(statistic0=c(z=lin.form/sqrt(var.obs.score)), 
                 p.value0=pval0, 
                 statistic=c(z=adjusted.obs.score), 
                 p.value=pval,
                 parameter=c(lin.form=lin.form, sd=sqrt(var.obs.score)), 
                 method='Profiled Variance Component Test Based on Asymptotic Normality', 
                 alternative='greater',
                 null.value=structure(numeric(length(nonNull)), names=sprintf('variance component %d', tau.idx)) 	#, null.fit=null.fit
                ),
            class='htest'
           )
}
varComp.LinScore.Satterthwaite <-
function( # null.fit,
	all.scores, infoMat, null, w, tr1, n, ...)
{
#Actual Satterthwaite approximation for LinScore. Argument definitions are the same as in varComp.LinScore.SSAS155, except the LIkLI, tau.idx, and non.pd are not used. 
  nK=length(all.scores)
  rqforms=(all.scores+tr1/2)  # n/2 * ratio of qforms
#  
  nonNull=seq_len(nK)[-null]
#  
  Phi=infoMat[null, nonNull, drop=FALSE]  
  Delta=infoMat[null, null, drop=FALSE]   
#  
  var.score=infoMat[nonNull, nonNull, drop=FALSE] - crossprod(Phi, solve(Delta, Phi))
  mean.score=tr1[nonNull]*.5 - crossprod(Phi, solve(Delta, tr1[null]*.5-rqforms[null]))
#  
  sum.score=sum(w*mean.score)
  var.sum=drop(crossprod(w, var.score%*%w))
#  
  scale=var.sum/2/sum.score
  df=2*sum.score^2/var.sum
#
  obs.stat=sum(rqforms[-null]*w)
  pval=pchisq(obs.stat/scale, df, lower.tail=FALSE)
#  
  ans=list(statistic=c(`sum of ratio of quadratic forms`=obs.stat), 
           p.value=pval, 
           alternative='greater', 
           parameter=c(scale=scale, df=df), 
           null.value=structure(numeric(length(nonNull)), names=sprintf('variance component %d', nonNull)), #null.fit=null.fit, 
           method='Satterthwaite Approximation of Profiled Variance Component Test'
          )
  class(ans)='htest'
  ans  
}

#' @importFrom CompQuadForm davies
varComp.LinScore.SSAS155 <-
function(	#null.fit, 
	all.scores, lin.form, infoMat, null, w, tr1, n, LIkLI, tau.idx, non.pd, control, ...)
{
#Actual shifted and scaled AS155
#i.	null.fit: object from varComp under the null. 
#ii.	all.scores: All scores.
#iii.	lin.form: test statistic. 
#iv.	infoMat: information matrix. 
#v.	null: the same as in varComp.test.
#vi.	w: vector of weights. 
#vii.	tr1: vector of traces.
#viii.	n: residual sample size
#ix.	LIkLI: list of matrices, each being sqrt(V^{-1})' K_i sqrt(V ^{-1})
#x.	tau.idx: complement of null.
#xi.	non.pd: A logical indicating that the infoMat is not positive definite. This only applies when the information argument of varComp.test is not `EI`. 
#control: Optional control object, primarily for getting extremely small p-values
#xii.	...: place holder for unused arguments.
  nK=length(all.scores)
  nonNull=seq_len(nK)[-null]
  Phi=infoMat[null, nonNull, drop=FALSE]  
  Delta=infoMat[null, null, drop=FALSE]   

  var.nonNull.score=infoMat[nonNull, nonNull, drop=FALSE] - crossprod(Phi, solve(Delta, Phi))
  mean.nonNull.score= crossprod(Phi, solve(Delta, all.scores[null]))
  mean.obs.score=sum(mean.nonNull.score*w)
  var.obs.score=drop(crossprod(w, var.nonNull.score%*%w))
  
  mean.theo.score=0 
  var.theo.score=drop(crossprod(w, infoMat[nonNull, nonNull, drop=FALSE]%*%w))
  
  adjusted.obs.score=sqrt(var.theo.score/var.obs.score)*(lin.form-mean.obs.score) + mean.theo.score
  
  pval=davies(0, eigen(n*Reduce('+', mapply('*',w,LIkLI[tau.idx],SIMPLIFY=FALSE)),TRUE,TRUE)$val - sum(tr1[tau.idx]*w)-adjusted.obs.score*2, acc=control$acc, lim=control$lim)
  pval0=davies(0, eigen(n*Reduce('+', mapply('*',w,LIkLI[tau.idx],SIMPLIFY=FALSE)),TRUE,TRUE)$val - sum(tr1[tau.idx]*w)-lin.form*2,acc=control$acc, lim=control$lim)
  structure(list(statistic0=c(`linear form`=lin.form), 
                 p.value0=pval0$Qq, 
                 statistic=c(`adjusted linear form`=adjusted.obs.score), 
                 p.value=pval$Qq,
                 parameter=c(w=w, 
                                 acc=control$acc, 
                                 lim=control$lim, 
                                 ifault=pval$ifault, 
                                 trace=pval$trace, 
                                 non.pd=non.pd
                                ), 
                 method='Profiled Variance Component Test through Shifted & Scaled AS155', 
                 alternative='greater',
                 null.value=structure(numeric(length(nonNull)), names=sprintf('variance component %d', tau.idx)) #, null.fit=null.fit
                ),
            class='htest'
           )
}

#' @importFrom quadprog solve.QP
#' @importFrom Matrix nearPD
#' @importFrom CompQuadForm davies
varComp.LinScore.test <-
function(control, n, tau.idx, LIkLI, tr1, infoMat, all.scores, non.pd)  #, ...)
{
#### external ref: null.fit, X, Y, K, non.pd

#Compute weights for LinScore; do test when the null has no other variance components; pass to varComp.LinScore.approximation when the null has other variance components. 
#i.	control: the control object from varCompTest.control.
#ii.	null: the same as in varComp.test. 
#iii.	n: residual sample size.
#v.	tau.idx: the complement of null. 
#vi.	LIkLI: the same as in varComp.VM03.test.
#vii.	tr1: the same as varComp.VM03.test
#viii.	obs.score: observed score vector for the alternative\null
#ix.	infoMat: the same as varComp.VM03.test
#x.	all.scores: all scores
#xi.	...: place holder for unused arguments. 
  null = seq_along(all.scores)[-tau.idx]
  nNull=length(null)
  wts=control$wt
  obs.score = all.scores[tau.idx]
  if(length(null)==0L) {
    methods=control$method[[1]] 
    if(!any(methods%in%c('AS155', 'exact', 'Davies', 'SSAS155')))
      message('LinScore test method was changed to the exact method "AS155" when null=integer(0L). ')
    methods='AS155'
  }else{
    methods=control$method[[2]]
    tmp=intersect(methods, c('AS155', 'exact', 'Davies'))
    if(length(tmp)>0L){
      warning(sprintf('LinScore test method(s) [%s] was changed to "SSAS155" when null!=integer(0L).', paste(tmp,collapse=', ')))
      methods=unique(c('SSAS155', setdiff(methods, tmp)))
    }    
  }
  n.tau=length(tau.idx); nK=length(LIkLI)
  ans=vector('list', length(wts)); names(ans)=wts
  
  if(nNull>0L){
    partNegHess=infoMat[tau.idx,tau.idx,drop=FALSE]-
      infoMat[tau.idx,-tau.idx,drop=FALSE]%*%solve(infoMat[-tau.idx,-tau.idx,drop=FALSE],infoMat[-tau.idx,tau.idx,drop=FALSE])
    invNegHess=tryCatch(solve(partNegHess), error=function(e){partNegHess<<-as.matrix(nearPD(partNegHess)$mat);solve(partNegHess)})
  }else{
    partNegHess=infoMat; invNegHess=solve(infoMat)
  }
  if(any(grepl('^Wd', wts))){ ## Wald variant
    inv.info=solve(infoMat)
    mle=inv.info%*%all.scores
  }
    
  for(w.i in seq_along(wts)){
    tmpwts=gsub('^Wd','', wts[w.i])
    if(tmpwts!=wts[w.i]){  ## Wald
      bak.partNegHess=partNegHess
      partNegHess=inv.info[tau.idx, tau.idx, drop=FALSE]
    }
    w=switch(tmpwts,
        EqWt=rep(1/n.tau, n.tau),
        InvSTD=local({tmp=1/sqrt(diag(as.matrix(partNegHess))); tmp/sum(tmp)}),
        InvSqrtV=local({tmp=pmax(0, base::colSums(t(backsolve(chol(partNegHess), diag(1,nK))))); tmp/sum(tmp)}),
        MinVar=solve.QP(Dmat=partNegHess, dvec=numeric(n.tau), Amat=cbind(1, diag(1,n.tau)), bvec=rep(1:0,c(1,n.tau)), meq=1, factorized=FALSE)$solution 
      )
    if(tmpwts!=wts[w.i]){  ## Wald
      w=drop(w%*%inv.info[tau.idx, , drop=FALSE])
      partNegHess=bak.partNegHess
      lin.form=sum(w*all.scores)
    }else
      lin.form=sum(w*obs.score)

    this.ans=vector('list', length(methods)); names(this.ans)=methods
    for(m.i in seq_along(methods)){
      if(nNull==0L){
        stopifnot(methods[m.i]=='AS155')
        #pval=davies(0, eigen(n*Reduce('+', mapply('*',w,LIkLI[tau.idx],SIMPLIFY=FALSE)),TRUE,TRUE)$val - sum(tr1[tau.idx]*w)-lin.form*2, acc=control$acc, lim=control$lim) ## when length(null)==0L, tau.idx subsetting is not needed. 
        pval=davies(0, eigen(n*Reduce('+', mapply('*',w,LIkLI,SIMPLIFY=FALSE)),TRUE,TRUE)$val - sum(tr1*w)-lin.form*2, acc=control$acc, lim=control$lim)
        this.ans[[m.i]]=structure(list(statistic=c(`linear form`=lin.form), 
                                       p.value=pval$Qq, 
                                       parameter=c(w=w,
                                                   acc=control$acc, lim=control$lim,
                                                   trace=pval$trace, 
                                                   ifault=pval$ifault
                                                      ), 
                                       method='Exact variance component test', 
                                       alternative='greater',
                                       null.value=structure(numeric(n.tau), names=sprintf('variance component %d', tau.idx))
                                      ),
                                  class='htest'
                                 )
      }else{
        this.ans[[m.i]]=do.call(
			what = paste0('varComp.LinScore.', methods[m.i]),
			args=list(null=null, 
				w=w, lin.form=lin.form, LIkLI=LIkLI, tr1=tr1
				,infoMat=infoMat, all.scores=all.scores, n=n, tau.idx=tau.idx, non.pd=non.pd, control=control
				)
			)
      }
      ans[[w.i]]=this.ans
    }
  }
  class(ans)='varComp.LinScore.test'
  ans
}

varComp.LinScore.beta <-
function(...) .NotYetImplemented()

varComp.LinScore.cumulant3 <-
function(...) .NotYetImplemented()

varComp.RWD88.test <-
function(control,...)
{
  .NotYetImplemented()
}

varComp.Wald.test <-
function(control,...)
{ 
  .NotYetImplemented()
}

varComp.LinScore.pboot <-
function(...) .NotYetImplemented()

varComp.LinScore.saddlepoint <-
function(...) .NotYetImplemented()

nlminb.control=function(
	eval.max=200L, 
	iter.max=150L,
	trace = 0L, 
	abs.tol = 0, 
	rel.tol = 1e-10, 
	x.tol = 1.5e-8, 
	xf.tol = 2.2e-14, 
	step.min = 1, 
	step.max = 1, 
	sing.tol =rel.tol,
	...)
{structure(list(
	eval.max=eval.max, 
	iter.max=iter.max,
	trace = trace, 
	abs.tol = abs.tol, 
	rel.tol = rel.tol, 
	x.tol = x.tol, 
	xf.tol = xf.tol, 
	step.min = step.min, 
	step.max = step.max, 
	sing.tol =sing.tol, 
	...) 
	, class = 'nlminb.control')
}

#' @importFrom graphics plot abline
varComp.control=function(verbose = FALSE, start=NULL, REML = TRUE, 
information=c('AOI', 'WAI', 'AEI', 'OI', 'EI'), boundary.eps=5e-4, 
  nlminb=nlminb.control(iter.max=200L, eval.max=500L), 
  plot.it=FALSE, keepXYK=TRUE)
{	optMethods=c('nlminb', 'optim', 'NRGD')
	optMethod=optMethods[match('nlminb', optMethods)]
	structure(
	  list(optMethod=optMethod,
		 verbose=verbose, 
		 starts=start, 
		 REML = REML, 
		 information = match.arg(information, c('AOI', 'WAI', 'AEI', 'OI', 'EI')), 
		 boundary.eps= boundary.eps, 
		 nlminb = nlminb, 
		 plot.it=plot.it, 
		 keepXYK=keepXYK) #nStepHalving
	  , class = 'varComp.control'
	)
}

varComp=function(fixed, data, random, varcov, weights, subset, family = stats::gaussian('identity'), na.action, offset, control = varComp.control(...), 
     doFit = TRUE, normalizeTrace = TRUE, 
     contrasts = NULL, model = TRUE, X = TRUE, Y = TRUE, K = TRUE, ...)
{
	if(missing(fixed) || !is.formula(fixed) || length(as.list(fixed))!=3L)	stop('fixed needs to be a two-sided formula')
	if(!missing(random)) {
		if(!is.formula(random) || length(as.list(random))!=2L)  stop('random needs to be a right-sided formula, when non-missing')
		if(!identical(environment(fixed) , environment(random)) )warning("environment(random) does not match environment(fixed).")
		fixedRandom = update.formula(fixed, as.formula(paste0('.~.+',random[2],collapse=''))) # random[2] is a call to the RHS
		random = update.formula(random, as.formula('~.')) ## this will expand the formula (removing * shorthand)
	}else fixedRandom = fixed 
	
	if(!missing(varcov) && !is.matrix(varcov) && !is.list(varcov)) stop("varcov needs to be a matrix or a list, when non-missing")
	if(!missing(varcov) && is.matrix(varcov)) varcov=list(varcov)
    if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family))  family <- family()
	if(!all.equal(stats::gaussian('identity'), family)) { 
		warning("Currently only Gaussian family with identity link is supported")
		family=stats::gaussian(link='identity')
	}
	if (missing(data)) data <- environment(formula)
	ret.x = X; ret.y=Y; ret.k=K; ret.mod=model
	
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)

    m <- match(c("fixed", "data", "random", "weights" ,"subset", "na.action", "offset"
			), names(mf), 0L)
    mf <- mf[c(1L, m)]  ## mf[1] is varComp()
	
	if(m[5L]>0 && is.logical(subset)) {
		subset=which(subset)
		mf$subset=subset
	}
	if(m[4L]>0 && any(weights<=0)) {	## handling non-positive weights to subset argument
		if(m[5L]>0) subset = unique(setdiff(subset, which(weights<=0))) else subset = which(weights>0)
		mf$subset=subset
	}
	if(m[3L]>0L) mf['random']=NULL
	mf['fixed']=NULL
	mf$formula=fixedRandom		
	
	mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mfAll <- eval(mf, parent.frame())	## removing obs with missing data or not in subset
	
	if(!is.null( na.vec <- attr(mfAll, 'na.action'))){	## missing data removed. 
		if(m[5L]> 0 ){ ## there exists subset argument
			subset = sort(setdiff(subset, na.vec))
		}else   subset = sort(setdiff(seq_len(nrow(mfAll)+length(na.vec)), na.vec))
		mf$subset = subset
		mfAll = eval(mf, parent.frame())
	}
	mtAll = attr(mfAll, 'terms')
	
	# mf$data = as.name('mfAll')
	mf$formula=fixed
	eframe=sys.frame(sys.nframe())
	mf.fixed=eval(mf, parent.frame())
	
    mt.fixed <- attr(mf.fixed, "terms")
	fixedTerms = attr(mt.fixed, 'term.labels')
	fixedTerms = sapply(fixedTerms, sortTerm)
	fixedVars  = sapply(attr(mt.fixed, 'variables'), as.character)[-1L]
    Y <- model.response(mf.fixed, "numeric")
    w <- as.vector(model.weights(mf.fixed))
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf.fixed))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
		Y = Y - offset  ## CHECKME when family()$link is not identity
    }
    if (is.empty.model(mt.fixed)) {
        X <- matrix(0, length(Y), 0L)
    } else {
        X <- model.matrix(mt.fixed, mf.fixed, contrasts)
    }
	
	if(missing(random)){	
		if(missing(varcov)) varcov=list()
		K=varcov
		nK=length(K)
	}else{
	## preparing a vector of all random terms
		mf$formula=random
		# if(m[2L]>0L) mf$data = cl$data
		mf.random=eval(mf, parent.frame())
		mt.random=attr(mf.random, 'terms')
		rterms = attr(mt.random, 'term.labels')
		# rterms = unlist(strsplit(as.character(random[2]), ' *\\+ *'))
		rterms = setdiff(rterms, -1:1) ## removing all terms invoving intercept(s)
		rterms = sapply(rterms, sortTerm, priority=fixedVars)
		rterms = setdiff(rterms, fixedTerms)
		rterms = local({  ## expand all fixed-by-random interactions to ensure that later when fixed-by-random interactions are encountered, the marginal random effect does not need to be added. 
			ans=rterms
			for(iz in seq_along(rterms)){
				this.form=update.formula(random, as.formula(paste('~', rterms[iz], '+0', collapse='')))
				mf$formula = this.form
				mf.this <- eval(mf, parent.frame())
				# isFact <- vapply(mf.this, function(x) is.factor(x) || is.logical(x), NA)
				isFixed = names(mf.this)%in%fixedVars
				# idx = (isFact & isFixed)
				idx = isFixed
				if(any(idx) && !all(idx)) ans[iz] = gsub(":", "*", ans[iz], fixed=TRUE)
			}
			ans
		})
		mf$formula = update.formula(random, as.formula(paste0("~", paste0(rterms, collapse='+'))))
		mf.random=eval(mf, parent.frame())
		mt.random=attr(mf.random, 'terms')
		rterms = attr(mt.random, 'term.labels')
		rterms = setdiff(rterms, -1:1) ## removing all terms involving intercept(s)
		rterms = sapply(rterms, sortTerm, priority=fixedVars)
		rterms = setdiff(rterms, fixedTerms)

		#### processing each random term
		nRterms = length(rterms)
		Z=vector('list', nRterms)
		j=1L
		for(iz in seq_len(nRterms)){
			this.form=update.formula(random, as.formula(paste('~', rterms[iz], '+0', collapse='')))
			mf$formula = this.form
			mf.this <- eval(mf, parent.frame())
			mt.this <- attr(mf.this, "terms")
			tmp.model.matrix = suppressWarnings(model.matrix(mt.this, mf.this, contrasts))
			isFact <- vapply(mf.this, function(x) is.factor(x) || is.logical(x), NA)
			isFixed = names(mf.this)%in%fixedVars
			idx = which(isFact & isFixed)
			if(length(idx)>0L){ ## random interactions involving fixed variables: previously, it has been ensured that random marginal effects have already been included
				fterm=paste0(sort(names(mf.this)[idx]), collapse=':')
				if(fterm%in%fixedTerms){ ## truly fixed-by-random interaction
					this.form=update.formula(this.form, paste0("~",fterm))
					mf$formula=this.form
					mf.tmp=eval(mf, parent.frame())
					mt.tmp=attr(mf.tmp, 'terms')
					oldContr=getOption('contrasts')
					options(contrasts=c('contr.sum','contr.sum'))
					fixedX = suppressWarnings(model.matrix(mt.tmp, mf.tmp, contrasts)[,-1L,drop=FALSE])
					options(contrasts=oldContr)
					
					tmpRterm = substr(rterms[iz], nchar(fterm)+2L, nchar(rterms[iz]))
					if(tmpRterm=='') {warning("DEBUG ME: tmpRterm should not be empty"); browser()}
					this.form=update.formula(this.form, paste0("~", tmpRterm, '+0'))
					mf$formula=this.form
					mf.tmp=eval(mf, parent.frame())
					mt.tmp=attr(mf.tmp, 'terms')
					rdmZ=suppressWarnings(model.matrix(mt.tmp, mf.tmp, contrasts))
					
					for(tmpi in seq_len(ncol(fixedX))){
						Z[[j]] = fixedX[, tmpi] * rdmZ
						names(Z)[j] = paste(colnames(fixedX)[tmpi], tmpRterm, sep=':')
						j=j+1L
					}
				}else{	## the fixed interaction is actually random
					Z[[j]] = tmp.model.matrix
					names(Z)[j] = rterms[iz]
					j=j+1L
				}				
			}else{	## no fixed terms involved
				Z[[j]] = tmp.model.matrix
				names(Z)[j] = rterms[iz]
				j = j+1L
			}
		}
		nRterms = length(Z)
		if(missing(varcov)){
		  nK=nRterms
		  K=vector('list', nK)
		  for(j in seq_len(nK))  K[[j]]=tcrossprod(Z[[j]])
		}else{
			K=varcov
			nK=length(K)
		  if(nRterms!=nK) stop('The number of matrices in "varcov" needs to equal the number of random effect terms in "random", when both are provided.')
		  ## match the order of elements in Z with that in K the same way as matching arguments in function calls
			tmpf=function(){mget(ls())}
			tmpz=seq_along(Z); namesZ=names(tmpz)=names(Z)
			formals(tmpf)=as.pairlist(tmpz)
			tmpk=as.list(seq_along(K)); names(tmpk)=names(K)
			zk.idx=do.call('tmpf', tmpk)
			K.bak=K
		  for(j in 1:nK)  K[[j]]=tcrossprod(Z[[j]]%*%K.bak[[ zk.idx[[ namesZ[j] ]] ]], Z[[j]])
		}
 	    names(K)=names(Z)
 	}
	
	if(isTRUE(normalizeTrace)) K=lapply(K, normalizeTrace)
	
	if(!missing(weights)){
		keep.weights= which(w>0)
		w.5 = sqrt(w[keep.weights])
		rw.5 = w.5[rep(seq_along(keep.weights), each=length(keep.weights))]
		Y=Y[keep.weights]*w.5
		X=w.5*X[keep.weights, , drop=FALSE]
		for(ik in seq_len(nK)) K[[ik]] = w.5*K[[ik]][keep.weights, keep.weights, drop=FALSE]*rw.5
	}
	
	if(any(ret.x, ret.y, ret.k)) control$keepXYK=TRUE
	ansCall = call('varComp.fit', Y=Y, X=X, K=K, control=control)
	
	
	ans = list(doFit = ansCall)
	ans$fixef = matrix(NA_real_, ncol(X), 0L, dimnames=list(colnames(X),NULL))
	
    class(ans) <- "varComp"
    ans['na.action'] <- list(attr(mfAll, "na.action"))  ## could be null
    ans['offset'] <- list(offset) ## could be null
    ans['contrasts'] <- list(attr(X, "contrasts"))  ## if X is empty, this could be null
    ans$xzlevels <- .getXlevels(mtAll, mfAll)
    ans$call <- cl
    ans$terms <- mtAll
	ans$nobs = length(Y)
	ans$control=control
    ## CHEKME: 
	if (ret.mod)  ans$model <- mfAll
    if (!ret.x)     ans$X <- NULL else ans$X=X
    if (!ret.y)     ans$Y = NULL else ans$Y = Y
	if (!ret.k) 	 ans$K = NULL else ans$K = K
	ans$random.labels = if(is.null(names(K))) {if(nK>0L) paste("varComp", seq_along(K), sep='.') else character(0L)} else names(K)
	if(missing(weights)) ans$weights = NULL else ans$weights = w
	
	if(isTRUE(doFit)) doFit.varComp(ans) else ans
}

#' @importFrom utils combn
#' @importFrom graphics plot abline
varComp.fit = function(Y, X=matrix(0,length(Y),0L), K, control=varComp.control())
{
	{
	#                  a.	varComp: The main fitting function using REML criterion
	#                      i.	Y: response
	#                      ii.	X: fixed effect design matrix. Treated as intercept if missing. A zero matrix can be supplied for REML. 
	#                      iii.	Z: random effect design matrix, or a list of such matrices. Treated as identity if missing. 
	#                      iv.	K: A matrix or a list of matrices, each determining the correlations among observations. This has to be at least p.s.d. in this implementation. 
	#                      v.	information: A character specifying the information matrix to be used during optimization. For doing score tests, this has to be EI; for fitting models, others are OK and the default AOI is usually fast enough.
	#                      vi.	optMethod: A character specifying the optimization method. The default should work well in most cases. 'NRGD' is only a naive implementation and has not been tested extensively. 
	#                      vii.	boundary.eps: A small numeric, below which estimates of tau will be checked for hitting boundary (0). This is a simple attempt at differentiating small but non-zero variance components vs truly zero variance components. Another way of doing this is to divide the K matrices by some large number. 
	#                      viii.	conv: A small numeric, the convergence criterion, only used when optMethod='NRGD'.
	#                      ix.	nStepHalving: A positive integer, giving the max number of step halving during NRGD. 
	#                      x.	nIter: A positive integer, giving the max number of iterations allowed for NRGD.
	#                      xi.	starts: A vector of nonnegative doubles of the same length as the number of variance components. This is the starting value for tau, i.e., the ratio of each variance components to the error variance. 
	#                      xii.	plot.it: A logical scalar, indicating if PREML surface will be plotted (for one variance components only). Please send Long Qu an Email with your data set if you find multiple local maxima, for possible improvements of this function. 
	#                      xiii.	verbose: A logical scalar, indicating if extra information is printed. 
	#                             keepXYK: Logical, indication if original X, Y and K matrices are stored in the results. 
	}
  #if(!is.matrix(X)) X=as.matrix(X)  ## does it matter?
  #if(diff(range(X))>0) X=model.matrix(~X)  ##  dealing with constant X including zero X
  #Y=as.numeric(as.vector(Y)) ## plan to move to varComp formula interface for efficiency
	{
		 optMethod=control$optMethod
		 verbose= control$verbose 
		 starts= control$starts 
		 REML = control$REML 
			if(!isTRUE(REML)) stop("Currently only REML method is implemented")
		 information = control$information 
		 boundary.eps= control$boundary.eps 
		 nlminb.control = control$nlminb 
		 plot.it= control$plot.it 
		 keepXYK= control$keepXYK
		 
		 environment(V)=
		 environment(updateLI)=
		 environment(updateLIkLI)=
		 environment(updateLIy)=
		 environment(PREML)=
		 environment(preprocPREML)=
		 environment(obj)=
		 environment(obj2)=
		 environment(updateNumsPart)=
		 environment(updateNums)=
		 environment(updateDenom)=
		 environment(updateTr1)=
		 environment(updateTr2)=
		 environment(updateNums2)=
		 environment(score)=
		 environment(updateNegGrad)=
		 environment(gradFunc)=
		 environment(OI)=
		 environment(hess)=
		 environment(EI)=
		 environment(AOI)=
		 environment(AEI)=
		 environment(WAI)=
		 environment(updateNegHess)=sys.frame(sys.nframe())
	}

  
  qrx=qr(X)
  # * this block checks if an intercept is in the model or not; 
  # * don't needed this anymore
	 # Q1=qr.Q(qrx)
	 # if(qrx$rank>0 && max(abs(Q1%*%crossprod(Q1, rep(1,length(Y))) - rep(1,length(Y)))) > sqrt(.Machine$double.eps)){  ## no intercept
	   # X=cbind(`(Intercept)`=1, X)
	   # qrx=qr(X)
	 # }
  
  Q2=if(qrx$rank>0) qr.Q(qrx, complete=TRUE)[, -seq_len(qrx$rank),drop=FALSE] else diag(1, length(Y)) #  qr.Q(qrx, complete=TRUE) == IdMat
  y=crossprod(Q2, Y)
  n=length(y)
  
  if(missing(K) || length(K)==0) {  ## linear model fit
		null.sig2=drop(crossprod(y))  / n   
		null.preml=.5*(
		  -n*log(crossprod(y)) -n -n*log(2*pi)-n*log(n)
		)

	  ans=list(
	    ## varComp.fit specific block
		parms=numeric(0L),
		gradients=numeric(0L), 
		hessian=matrix(numeric(0L), 0L, 0L),
		sigma2=drop(null.sig2), 
		varComps=numeric(0),
		n.iter=0L, PREML=drop(null.preml),
		X.Q2=Q2, 
		residual.contrast=y, 
		working.cor=vector('list', 0L),
		
		## varComp common block
		# na.action=NULL,
		# offset = NULL,
		# contrasts=NULL,
		# xzlevels = NULL,
		# terms = NULL, 
		call=match.call(), 
		nobs = length(Y), 
		control=control, 
		random.labels = character(0L), 
		doFit= TRUE,
		
		# frame=if(isTRUE(keepXYK)) NULL else parent.frame(), 
		X=if(isTRUE(keepXYK)) X else NULL,
		# qrx = if(isTRUE(keepXYK)) qrx else NULL, 
		Y=if(isTRUE(keepXYK)) Y else NULL, 
		K=if(isTRUE(keepXYK)) vector('list', 0L) else NULL
	  )
	  class(ans)='varComp'	
	  return(ans)
  }
  
  nK=length(K)
  if(nK>1) plot.it=FALSE
  
  k=vector('list',nK)
  for(j in 1:nK) k[[j]]=crossprod(Q2, K[[j]]%*%Q2)
  
  LIkLI=k; LIy=numeric(n); numsPart=matrix(NA_real_, n, nK); 
  tr2=negHess=nums2=matrix(NA_real_, nK, nK)
  tr1=nums=negGrad=numeric(nK); denom=numeric(1L)
  LI=matrix(NA_real_, n, n)
  diag.1.n=diag(1,n)
  if(nK==1) {
	eigK=eigen(k[[1]],TRUE)
	eigK$tvector=t(eigK$vector)
  }
  
  infoFunc=get(information)
  preprocScore=expression({
	updateLIkLI()
	updateNumsPart();   
	updateNums(); updateDenom(); updateTr1()
  })
  preprocInfo=switch(information,
	OI=expression({ updateNums2(); updateTr2() }),
	EI=expression({ updateTr2();  }),
	AOI=expression({updateNums2();    }),
	AEI=expression({updateNums2();   }),
	WAI=expression({updateNums2();   })
  )
  
  
  last.tau=tau=if(is.null(starts)) {
	minque(y, k, rep(0,nK), lower.bound=.Machine$double.eps^.5, restricted=TRUE)
  }else rep(starts, length=nK)
  ltau=log(max(.Machine$double.eps, tau))
  
  if(optMethod=='nlminb'){
	objNeg=function(tau)-obj(tau)
	gradNeg=function(tau)-gradFunc(tau)
	hessNeg=function(tau)-hess(tau)
	nlminb.fit=nlminb(tau, objNeg, gradNeg, hessNeg, lower=rep(0,nK), control=nlminb.control)
	tau=nlminb.fit$par
	n.nr=nlminb.fit$iterations
  }else if(optMethod=='optim'){
	optim.fit=optim(tau, obj, gradFunc, method='L-BFGS-B', lower=rep(0, nK), control=list(fnscale=-1))
	tau=optim.fit$par
	n.nr=optim.fit$counts
  }else if(optMethod=='NRGD'){
	nStepHalving = 20L
	preprocPREML(tau)
	last.func = PREML()
	n.nr=0L
	repeat{
	  if(verbose) cat('PREML=',last.func, 'tau=',tau,'\n')
	  eval(preprocScore)
	  eval(preprocInfo)
	  updateNegGrad()
	  updateNegHess()
	  
	  negGrad=tau*negGrad
	  negHess=negHess*outer(tau,tau)
	  diag(negHess)=diag(negHess)+negGrad
	  
	  adj.ltau=solve(negHess, negGrad)
	  
	  doNewton=TRUE
	  repeat{
		step.size=1
		n.step=0L
		ltau.new=ltau
		repeat{
		  ltau.new=ltau - step.size * adj.ltau
		  preprocPREML(exp(ltau.new))
		  this.func=PREML()
		  if( ( (this.func>last.func) & doNewton) | ((!doNewton) & (this.func>=last.func)) ) {
			doNewton=TRUE
			break
		  }
		  step.size=step.size/2
		  n.step=n.step+1L
		  if((n.step > nStepHalving) & doNewton ) {  ## try gradient descent
			#warning('Step halving failed')
			adj.ltau=negGrad
			doNewton=!doNewton
			break
		  }
		}
		if(isTRUE(doNewton)) break
	  }
	  
	  last.tau=tau
	  ltau=ltau.new
	  tau=exp(ltau)
	  
	  if(max(abs(last.tau-tau)) / max(abs(last.tau)+abs(tau)) < control$nlminb$x.tol){
		break
	  }
	  n.nr=n.nr+1L
	  last.func=this.func
	  if(n.nr >= control$nlminb$iter.max){
		warning(paste(control$nlminb$iter.max, 'iterations reached'))
		break
	  }
	}    
  }else stop('Method not implemented')
  

  bd.idx=which(tau<boundary.eps)
  if(FALSE){
	  if(length(bd.idx)>0L){#browser()
		preprocPREML(tau)
		eval(preprocScore)
		all.score=score()
		bd.idx=bd.idx[ all.score[bd.idx] < -boundary.eps ]
		
		if(length(bd.idx)==length(tau)) {
		  tau=rep(0,length(tau))
		}else if(length(bd.idx)>0L) {
		  bd.idx0= bd.idx [ all.score[bd.idx]==0 ]
		  if(length(bd.idx0) != length(bd.idx)){
			new.fit=Recall(y, rep(0,length(y)), , k[-bd.idx], information, method, boundary.eps, conv, nStepHalving, control$nlminb$iter.max, starts=tau[-bd.idx], plot.it, verbose)
			if(new.fit$PREML >= PREML()){
			  tau[-bd.idx]=new.fit$parms
			  tau[bd.idx]=0
			}
		  }
		}
	  }
  }

  if((nearZero = length(bd.idx))>0L){  # check boundary
	bd.idx.all= if(nearZero == 1) list(matrix(bd.idx)) else lapply(seq_len(nearZero), combn, x=bd.idx)
	cur.obj=obj(tau)
	tau.bak = tau
	for(i.0 in bd.idx.all){
		if(nrow(i.0)<nK){
			for(j.0 in seq_len(ncol(i.0))){
				this.0 = i.0[,j.0]
				setZero=function(tau){tau0=tau.bak; tau0[this.0]=0; tau0[-this.0]=tau; tau0}  # this.0
				objNeg0=function(tau) objNeg(setZero(tau))
				gradNeg0=function(tau) gradNeg(setZero(tau))[-this.0]
				hessNeg0=function(tau) hessNeg(setZero(tau))[-this.0, -this.0, drop=FALSE]
				nlminb.fit=nlminb(tau[-this.0], objNeg0, gradNeg0, hessNeg0, lower=rep(0,nK-length(this.0)), control=nlminb.control)
				if( -nlminb.fit$objective > cur.obj) {tau = setZero(nlminb.fit$par); cur.obj=-nlminb.fit$objective}
			}
		}else{
			if( (tmp = obj(rep(0, nK))) > cur.obj){ tau = rep(0, nK); cur.obj = tmp}
		}
	}
  }

  preprocPREML(tau)  
  eval(preprocScore)
  eval(preprocInfo)
  updateNegGrad()
  updateNegHess()
  sigma2=crossprod(LIy)/n
	
  if(isTRUE(plot.it)){
	taus=exp(seq(log(1e-9), log(2*tau), length=500))
	taus=sort(unique(c(taus, seq(0, 2*tau, length=500))))
	objs=sapply(taus, obj)
	# x11()
	plot(taus, objs, xlab='tau', ylab='objective', type='o',main='Profiled residual log likelihood')
	abline(v=tau)
  }
  
  nm=names(K)
  if(is.null(nm)) nm = if(nK>0L) paste('varComp', seq_along(K), sep='.') else character(0L)
  ans=list(
	## varComp.fit specific block
	parms=structure(tau, names=nm),
	gradients=structure(-negGrad, names=nm), 
	hessian=structure(hess(tau), dimnames=list(nm,nm)),
	sigma2=drop(sigma2), 
	varComps=structure(drop(tau*sigma2), names=nm),
	n.iter=n.nr, PREML=PREML(),
	X.Q2=Q2, 
	residual.contrast=y, 
	working.cor=structure(k, names=nm),
	
	## varComp common block
	# na.action=NULL,
	# offset = NULL,
	# contrasts=NULL,
	# xzlevels = NULL,
	# terms = NULL, 
	call=match.call(), 
	nobs = length(Y), 
	control=control, 
	random.labels = nm, 
	doFit= TRUE,
	# frame=if(isTRUE(keepXYK)) NULL else parent.frame(), 
	
	X=if(isTRUE(keepXYK)) X else NULL,
	# qrx = if(isTRUE(keepXYK)) qrx else NULL, 
	Y=if(isTRUE(keepXYK)) Y else NULL, 
	K=if(isTRUE(keepXYK)) K else NULL
  )
  class(ans)='varComp'
  ans
  
}

#' @importFrom MASS ginv
doFit.varComp=function(object)
{
	if(isTRUE(object$doFit)) return(object)
	ans=eval(object$doFit)

	if(ncol(ans$X)>0L){
		this.V=vcov(ans, what='Y')
		this.Vbet=ginv(crossprod(ans$X, solve(this.V,ans$X)))
		this.bet=drop(this.Vbet%*%crossprod(ans$X, solve(this.V, ans$Y)))
		names(this.bet)=colnames(ans$X)
		ans$fixef = this.bet
	}else ans$fixef = numeric(0L)
	ans$na.action = object$na.action
	ans$offset = ans$offset
	ans$contrasts = object$contrasts
	ans$xzlevels = object$xzlevels
	ans$call = object$call
	ans$terms = object$terms
	ans$model = object$model
	ans$X = object$X
	ans$Y = object$Y
	ans$K = object$K
	ans$random.labels = object$random.labels
	ans$weights = object$weights
	ans	
}

#' @importFrom RLRsim RLRTSim
varComp.RLRT.test <-
function(control,negHess, alt.fit, null.fit, tau.idx)  #, k,null,...)
{
#Actual testing function for RLRT. 
#i.	control. the control object from varCompTest.control. 
#ii.	negHess: Negative expected Hessian matrix. 
#iii.	alt.fit: Object from varComp under alternative hypothesis. 
#iv.	null.fit: Object form varComp under null hypothesis. 
#v.	tau.idx: complement of null. 
#vi.	k: list of corrected correlation matrices. 
#vii.	null: the same as the null in varComp.test. 
#viii.	...: place holder for unused arguments. 
	
  nNull=length(null.fit$random.labels); 
  null = seq_along(alt.fit$random.labels)[-tau.idx]
  n.tau=length(tau.idx)
  k = alt.fit$working.cor
  if(n.tau==1L) methods=control$method[[1]] else methods=control$method[[2]]
#  if(nNull>0L) { ## FIXME: implement pboot
#    null.htest=structure(list(p.value=NA_real_),class='htest')
#    ans=vector('list', length(methods)); names(ans)=methods
#    for(m in methods) ans[[m]]=null.htest
#    class(ans)='varComp.RLRT.test'; return(ans) 
#  }
  LRstat0=max(0, 2*(alt.fit$PREML-null.fit$PREML))
  
  nsim=control$nsim
  ans=vector('list', length(methods)); names(ans)=methods
  for(m.i in seq_along(methods)){
    if(methods[m.i]=='pboot' && nNull==0L && n.tau==1L)  methods[m.i]='exact' ## C & R (2004) situation
    if(methods[m.i]=='pboot'){
      .NotYetImplemented()
    }else if( methods[m.i]=='exact'){
      if(n.tau!=1L){
        #warning('RLRT is not implemented for multiple variance components')
        ans[[m.i]]=structure(list(statistic=c(RLR=LRstat0), p.value=NA_real_, parameters=c(nsim=nsim), null.stats=NA_real_), class='htest')
      }else{
        K.half=suppressWarnings(t(chol(k[[tau.idx]], pivot=TRUE)))
        if(nNull>0L)  warning('RLRT when the null model has other variance component than the error variance is only approximate!')
        
        oo=order(attr(K.half, 'pivot'))                                       ## fixed by LQ on 050412
        K.des=K.half[oo, seq_len(attr(K.half,'rank'))]
        LRstats=RLRTSim(matrix(NA_real_,nrow(K.des),0L), Z=K.des, sqrt.Sigma=diag(1,ncol(K.des)), nsim=nsim-1L)
        pval=(1+sum(LRstats>=LRstat0))/(nsim)
        ans[[m.i]]=structure(list(statistic=c(RLR=LRstat0), p.value=pval, parameters=c(nsim=nsim),null.stats=LRstats
                                  ,null.fit=null.fit, alt.fit=alt.fit
                                  ), class='htest')
      }
    }else stop('Unrecognized RLRT method')
  }
  class(ans)='varComp.RLRT.test'
  ans
}

#' @importFrom quadprog solve.QP
#' @importFrom Matrix nearPD
varComp.SS95.test <-
function(control,infoMat, tau.idx, LIkLI, LIy, tr1, n, all.scores)# , ...)
{
#Actual testing function for SS95 test. Arguments are the same as in varComp.VM03.test. 
	null = seq_along(all.scores)[-tau.idx]
  nNull=length(null); n.tau=length(tau.idx)
  if(length(null)==0L) methods=control$method[[1]] else methods=control$method[[2]]
  if(nNull>0L && length(setdiff(methods, 'pboot'))==0L) { ## FIXME: implement pboot
	warning("parametric bootstrap has not been implemented yet for the situation when the null model is not a linear model")
    null.htest=structure(list(p.value=NA_real_),class='htest')
    ans=vector('list', length(methods)); names(ans)=methods
    for(m in methods) ans[[m]]=null.htest
    class(ans)='varComp.SS95.test'; return(ans) 
  }
  #sumLIkLI=mapply('*',w,LIkLI[tau.idx])
  if(nNull>0L){
    Phi=infoMat[-tau.idx, tau.idx, drop=FALSE]
    Delta=infoMat[-tau.idx, -tau.idx, drop=FALSE]
    partNegHess=infoMat[tau.idx,tau.idx,drop=FALSE]-
      infoMat[tau.idx,-tau.idx,drop=FALSE]%*%solve(infoMat[-tau.idx,-tau.idx,drop=FALSE],infoMat[-tau.idx,tau.idx,drop=FALSE])
    invNegHess=tryCatch(solve(partNegHess), error=function(e){partNegHess<<-as.matrix(nearPD(partNegHess)$mat);solve(partNegHess)})
    mean.nonNull.score=drop(crossprod(Phi, solve(Delta, all.scores[-tau.idx])))
    adj.nonNull.score=all.scores[tau.idx] - mean.nonNull.score
    SS95stat0=max(0, -2*solve.QP(partNegHess, drop(adj.nonNull.score), diag(1,n.tau), rep(0,n.tau))$value)
  }else{
    partNegHess=infoMat; invNegHess=solve(infoMat)
    Amat=diag(1, n.tau); bvec=numeric(n.tau)
    SS95stat=function(z){
      S=.5* (sapply(LIkLI[tau.idx], function(lik) drop(crossprod(z, lik%*%z)/crossprod(z)))*n-tr1[tau.idx] )
      tvec=solve.QP(partNegHess, drop(S), Amat, bvec)$solution
      max(0, crossprod(tvec, partNegHess%*%tvec))
    }
    SS95stat0=drop(SS95stat(LIy))
  }
  nsim=control$nsim
  ans=vector('list', length(methods)); names(ans)=methods
  for(m.i in seq_along(methods)){
    if(methods[m.i]=='pboot'){
      if(nNull==0L){
        seed=get.seed()
        SS95stats=replicate(nsim-1L, SS95stat(rnorm(n)))
        pval=(1+sum(SS95stats>=SS95stat0))/(nsim)
        ans[[m.i]]=structure(list(statistic=c(SS95=SS95stat0), p.value=pval, null.stats=SS95stats, seed=seed, parameters=c( nsim=nsim)), class='htest')
      }else{
        .NotYetImplemented()
      }
    }else if(methods[m.i]=='ChiBarSq'){
      seed=get.seed()
      pval=pchibarsq(SS95stat0, invNegHess, lower.tail=FALSE)
      ans[[m.i]]=structure(list(statistic=c(SS95=SS95stat0), p.value=pval, seed=seed, parameters=invNegHess), class='htest')
    }else stop('Undefined methods for SS95')
  }
  class(ans)='varComp.SS95.test'
  ans
}

varCompTest.control <-
function(
  test="LinScore", 
  LinScore.wt="InvSTD", 
  LinScore.acc=1e-8, LinScore.lim=1e6L,
  LinScore.method=c('AS155', 'SSAS155'), 
  VM03.method=c('SSChiBarSq', 'pboot'),
  VM03.nsim=1e4L,
  SS95.method=c('pboot', 'pboot'),
  SS95.nsim=1e4L,
  RLRT.method=c('exact', 'pboot'), 
  RLRT.nsim=1e4L, 
  information = 'EI'
  # , Wald.method=list('pboot', 'pboot'),
  # RWD88.method=list('pboot', 'pboot')
)
{
#a.	varCompTest.control: Returns a list that determines the testing method for varComp.test
#i.	test: A character vector specifying the test to be performed. 
#1.	LinScore: Linear score tests
#2.	VM03: Projected, quadratic score test of Verbeke & Molenberghs (2003, Biometrics, 59, 254)
#3.	HP01: Projected, quadratic score test of Hall & Praestgaard (2001, Biometrika, 88, 739), which is the same as SS95. 
#4.	SS95: Projected, quadratic score test of Silvapulle & Silvapulle (1995, JASA, 90, 342)
#5.	RLRT: The restricted likelihood ratio test of Crainiceanu & Ruppert (2003, JRSSB, 66, 165) or its pseudo-likelihood heuristic extension by Greven et al. (2008, JCGS, 17, 870)
#6.	Wald: The Wald test (not implemented)
#7.	RWD88: The global score test suggested by Robertson, Wright, and Dykstra (1988, Order-Restricted Statistical Inferenc, p. 321) (not implemented).
#ii.	LinScore.wt: Character vector giving the method of finding weights of scores:
#1.	EqWt: equal weights
#2.	InvSTD: 1/standard deviation of scores
#3.	InvSqrtV: colSums of inverse square root of variance matrix of scores. 
#4.	MinVar: Minimizing the variance of convex combination of scores. 
#5.	WdEqWt: not implemented
#6.	WdInvSTD: not implemented
#7.	WdInvSqrtV: not implemented
#8.	WdMinVar: not implemented
#iii.	LinScore.acc: The same as the acc in CompQuadForm:::davies.
#iv.	LinScore.lim: The same as the lim in CompQuadForm:::davies.
#v.	LinScore.method: A list of two components for the LinScore test, each component being a character string. The first component specifies the method for obtaining null distributions when the null hypothesis contains no variance components other than the error variance; the second component specifies the method for obtaining null distributions when the null hypothesis contains additional variance parameters other than the error variance. 
#1.	AS155: Applied Statistics algorithm 155, i.e., in CompQuadForm:::davies.
#2.	exact: An alias of AS155.
#3.	Davies: An alias of AS155.
#4.	SSAS155: Shifted and scaled AS155 (for the 2nd component only)
#5.	Satterthwaite: scaled chi-square approximation (for the 2nd component only)
#6.	Normal: Normal approximation (for the 2nd component only)
#7.	LC12: Incorrect scaled chi-square approximation of Li and Cui (2012, AoAS) (for the 2nd component only)
#8.	LC12Boundary: Boundary corrected scaled chi-square approximation (for the 2nd component only)
#9.	beta: not implemented
#10.	cumulant3: not implemented
#11.	pboot: not implemented
#12.	saddlepoint: not implemented
#vi.	VM03.method: A list of two components for the VM03 test, each component being a character string. The first component specifies the method for obtaining null distributions when the null hypothesis contains no variance components other than the error variance; the second component specifies the method for obtaining null distributions when the null hypothesis contains additional variance parameters other than the error variance. Currently, the 2nd component is discarded. 
#1.	SSChiBarSq: shifted and scaled Chi-bar-square approximation (not used)
#2.	ChiBarSq: Chi-bar-square asymptotic null. 
#3.	pboot: Monte Carlo null.
#vii.	SS95.method: A list of two components for the SS95 test, each component being a character string. The first component specifies the method for obtaining null distributions when the null hypothesis contains no variance components other than the error variance; the second component specifies the method for obtaining null distributions when the null hypothesis contains additional variance parameters other than the error variance.
#1.	ChiBarSq: Chi-bar-square asymptotic null. 
#2.	pboot: Monte Carlo null (not implemented for the 2nd component)
#viii.	RLRT.method: A list of two components for the RLRT test, each component being a character string. The first component specifies the method for obtaining null distributions when the null hypothesis contains no variance components other than the error variance; the second component specifies the method for obtaining null distributions when the null hypothesis contains additional variance parameters other than the error variance.
#1.	exact: exact or pseudo-likelihood heuristic
#2.	pboot: For the 1st component, this is an alias of exact; for the 2nd component, this is not implemented. 
#ix.	Wald.method: Not implemented. 
  #x.	RWD88.method: Not implemented.
  varCompScoreTests = c(score='LinScore',score='VM03',score='SS95')
  varCompTests=c(varCompScoreTests,score='HP01',RLRT='RLRT') #,'Wald', 'RWD88')
  LinScoreWeightingMethods=c('EqWt', 'InvSTD', 'InvSqrtV', 'MinVar')
  test=match.arg(test, varCompTests, several.ok=TRUE)
  test[test=='HP01']='SS95'; test=unique(test)
  LinScore.wt=match.arg(LinScore.wt, LinScoreWeightingMethods, several.ok=TRUE)
  information = match.arg(information, c('AOI', 'WAI', 'AEI', 'OI', 'EI'))

  stopifnot(length(LinScore.method)==2L)                       
  LinScore.method=lapply(LinScore.method, function(zzz, tab) tab[pmatch(zzz, tab,duplicates.ok=TRUE)],  
                         tab=c('AS155', 'exact', 'Davies', 'SSAS155', 'Satterthwaite', 'Normal', 'LC12', 'LC12Boundary', 'beta', 'cumulant3', 'pboot', 'saddlepoint')
                        )

  stopifnot(length(VM03.method)==2L)                       
  VM03.method=lapply(VM03.method, function(zzz, tab) tab[pmatch(zzz, tab,duplicates.ok=TRUE)],  
                       tab=c('SSChiBarSq', 'ChiBarSq', 'pboot'))

  stopifnot(length(SS95.method)==2L)                       
  SS95.method=lapply(SS95.method, function(zzz, tab) tab[pmatch(zzz, tab,duplicates.ok=TRUE)],  
                       tab=c('ChiBarSq', 'pboot'))
                       
  stopifnot(length(RLRT.method)==2L)                       
  RLRT.method=lapply(RLRT.method, function(zzz, tab) tab[pmatch(zzz, tab,duplicates.ok=TRUE)],  
                       tab=c('exact', 'CR04', 'pboot'))

  # stopifnot(length(Wald.method)==2L)                       
  # Wald.method=lapply(LinScore.method, function(zzz, tab) tab[pmatch(zzz, tab,duplicates.ok=TRUE)],  
                       # tab=c('pboot'))

  # stopifnot(length(RWD88.method)==2L)                       
  # RWD88.method=lapply(RWD88.method, function(zzz, tab) tab[pmatch(zzz, tab,duplicates.ok=TRUE)],  
                       # tab=c('pboot'))

  ntest=length(test)
  ans=vector('list', ntest); names(ans)=test

  for(m in test){
    if(m=='LinScore') {
      ans[[m]]=list(wt=LinScore.wt, acc=LinScore.acc, lim=LinScore.lim, method=LinScore.method)
    }else if(m=='VM03') {
      ans[[m]]=list(nsim=VM03.nsim, method=VM03.method)
    }else if(m=='SS95') {
      ans[[m]]=list(nsim=SS95.nsim, method=SS95.method)
    }else if(m=='RLRT') {
      ans[[m]]=list(nsim=RLRT.nsim, method=RLRT.method)
    }else{
      ans[[m]]=list(method=get(paste(m,'method',sep='.')))
    }
  }
  varCompScoreTests = c(score='LinScore',score='VM03',score='SS95')
  if(any(test%in%varCompScoreTests)) ans$information=information
  class(ans) = 'varCompTest.control'
  ans
}

varComp.test = function(object, ...) UseMethod("varComp.test", object)

varComp.test.formula = function(object, data, random1, varcov1, random2, varcov2, fit.control, test.control, ...)
{
	this.call=match.call()
	this.call[[1L]] = as.name('varComp')
	m = match(c('random1', 'varcov1', 'random2', 'varcov2','fit.control','object'), names(this.call), 0L)
	if(m[1L]>0L) names(this.call)[m[1L]] = 'random'
	if(m[2L]>0L) names(this.call)[m[2L]] = 'varcov'
	this.call$test.control=this.call$test=NULL
	if(m[5L]>0L) names(this.call)[m[5L]] = 'control'
	if(m[6L]>0L) names(this.call)[m[6L]] = 'fixed'
	this.call$random2 = this.call$varcov2 = NULL
	this.call$doFit = FALSE
	nul=eval.parent(this.call)
	
	this.call=match.call()
	this.call[[1L]] = as.name('varComp')
	m = match( c('random1', 'varcov1', 'random2', 'varcov2','fit.control','object'), names(this.call), 0L)
	if(m[3L]>0L) names(this.call)[m[3L]] = 'random'
	if(m[4L]>0L) names(this.call)[m[4L]] = 'varcov'
	this.call$test.control=this.call$test=NULL
	if(m[5L]>0L) names(this.call)[m[5L]] = 'control'
	if(m[6L]>0L) names(this.call)[m[6L]] = 'fixed'
	this.call$random1 = this.call$varcov1 = NULL
	this.call$doFit = FALSE
	alt=eval.parent(this.call)
	
	if(missing(test.control)) varComp.test(nul, alt, ...) else varComp.test(nul, alt, control=test.control, ...)
}  

varComp.test.varComp = function(object, object2, 
	additional.varcov, null, test='LinScore', 
  control=varCompTest.control(test),  ...)
{
  # information=match.arg(information, informationTypes)
  varCompScoreTests = c(score='LinScore',score='VM03',score='SS95')
  varCompTests=c(varCompScoreTests,score='HP01',RLRT='RLRT') #,'Wald', 'RWD88')
  test=match.arg(test, varCompTests, several.ok=TRUE)
  test[test=='HP01']='SS95'; test=unique(test)
  ntest=length(test)
  if(ntest == 0L) stop(sprintf("Test needs to be given among [%s].", paste0(varCompTests, collapse=',')))

	if(!missing(object2)){  ## two object comparison
		if(!missing(additional.varcov) || !missing(null)) warning("'additional.varcov' and/or 'null' will be ignored when 'object2' is provided.")
		
		if(length(object2$random.labels) > length(object$random.labels)){
			nul=object
			alt=object2
		}else{
			nul=object2
			alt=object
		}
		if(!all(nul$random.labels %in% alt$random.labels)) stop("Two model ares not nested.")
		varComp.test.2modelDoTest(null.fit=nul, alt.fit=alt, test=test, control=control, ...)
	}else if(!missing(additional.varcov)){  ## null fit is available
		if(!missing(null)) warning("'null' is ignored when 'additional.varcov' is provided.")
		varComp.test.nulDoTest(null.fit=object, additional.varcov = additional.varcov, test=test, control=control, ...)
	}else if(!missing(null)){ ##  alt fit is available
		varComp.test.altDoTest(alt.fit=object, null=null, test=test,  control=control, ...)
	}else {  ## default testing all components of object
		varComp.test.altDoTest(alt.fit=object, null=integer(0L), test=test,  control=control, ...)
	}
}
varComp.test.2modelDoTest = function(null.fit, alt.fit, test='LinScore', control=varCompTest.control(test), ...)
{
  # information=match.arg(information, informationTypes)
  varCompScoreTests = c(score='LinScore',score='VM03',score='SS95')
  varCompTests=c(varCompScoreTests,score='HP01',RLRT='RLRT') #,'Wald', 'RWD88')
  test=match.arg(test, varCompTests, several.ok=TRUE)
  test[test=='HP01']='SS95'; test=unique(test)
  ntest=length(test)
  if(ntest == 0L) stop(sprintf("Test needs to be given among [%s].", paste0(varCompTests, collapse=',')))
  
	if(!isTRUE(null.fit$doFit)) {
		if(any(test%in%c(varCompScoreTests, 'RLRT'))){ ## needs a full fitting of null model
			null.fit=doFit.varComp(null.fit)
		}else { ## needs a partial fitting only
			### FIXME: consider Wald and RWD88 separately? 
			bak.control = tmp.control = null.fit$control
			tmp.control$nlminb$iter.max=0L
			null.fit$control = tmp.control
			null.fit=doFit.varComp(null.fit)
			null.fit$control=bak.control
		}
	}
	if((!isTRUE(alt.fit$doFit)) && any(test %in%c('RLRT','Wald'))){	## needs a full alternative fit
		alt.fit = doFit.varComp(alt.fit) 
	}else{ ## needs a partial alternative fit
		tmp.control = bak.control = alt.fit$control
		tmp.control$nlminb$iter.max=0L
		alt.fit$control = tmp.control
		alt.fit = doFit.varComp(alt.fit)
		alt.fit$control = bak.control
	}
  if(!all(null.fit$random.labels %in% alt.fit$random.labels)) stop("Two models are not nested.")
  null=which(alt.fit$random.labels %in% null.fit$random.labels)
  nNull = length(null)
  nK = length(alt.fit$working.cor) 
  
  tau.idx=seq_len(nK)
  if(nNull>0L) tau.idx=tau.idx[-null]
  n.tau=length(tau.idx)	
  
	environment(varComp.test.Common)  = sys.frame(sys.nframe())
	varComp.test.Common()
}
varComp.test.nulDoTest = function(null.fit, additional.varcov, test='LinScore', control=varCompTest.control(test), alt.fit=NULL, ...)
{
  # information=match.arg(information, informationTypes)
  varCompScoreTests = c(score='LinScore',score='VM03',score='SS95')
  varCompTests=c(varCompScoreTests,score='HP01',RLRT='RLRT') #,'Wald', 'RWD88')
  test=match.arg(test, varCompTests, several.ok=TRUE)
  test[test=='HP01']='SS95'; test=unique(test)
  ntest=length(test)
  if(ntest == 0L) stop(sprintf("Test needs to be given among [%s].", paste0(varCompTests, collapse=',')))

	if(!isTRUE(null.fit$doFit)) {
		if(any(test%in%c(varCompScoreTests, 'RLRT'))){ ## needs a full fitting of null model
			null.fit=doFit.varComp(null.fit)
		}else { ## needs a partial fitting only
			### FIXME: consider Wald and RWD88 separately? 
			bak.control = tmp.control = null.fit$control
			tmp.control$nlminb$iter.max=0L
			null.fit$control = tmp.control
			null.fit=doFit.varComp(null.fit)
			null.fit$control=bak.control
		}
	}
	
  null=seq_along(null.fit$working.cor)
  nNull = length(null)
  if(is.matrix(additional.varcov)) additional.varcov = list(additional.varcov)
  nK = length(additional.varcov) + nNull

  tau.idx=seq_len(nK)
  if(nNull>0L) tau.idx=tau.idx[-null]
  n.tau=length(tau.idx)
  
	if(!isTRUE(alt.fit$doFit)){
		if(is.null(alt.fit)){
			k = vector('list', nK)
			k[seq_len(nNull)] = null.fit$working.cor
			Q2 = null.fit$X.Q2
			for(i.k in (nNull+1L):nK) k[[i.k]] = crossprod(Q2, additional.varcov[[i.k-nNull]]%*%Q2)
		}
		if(any(test %in%c('RLRT','Wald'))){	## needs a full alternative fit
			if(!is.null(alt.fit)) {
				alt.fit = doFit.varComp(alt.fit) 
			}else {
				alt.fit = varComp.fit(Y=null.fit$residual.contrast, K = k, control = null.fit$control)
			}
		}else if(any(test%in%varCompScoreTests)){ ## needs a partial fit
			if(!is.null(alt.fit)){
				tmp.control=bak.control=alt.fit$control
				tmp.control$nlminb$iter.max=0L
				alt.fit$control=tmp.control
				alt.fit=doFit.varComp(alt.fit)
				alt.fit$control=bak.control
			}else{
				tmp.control=null.fit$control
				tmp.control$nlminb$iter.max=0L
				alt.fit=varComp.fit(Y=null.fit$residual.contrast, K=k, control=tmp.control)
			}
		}else stop("should not reach here")
	}
	environment(varComp.test.Common)  = sys.frame(sys.nframe())
	varComp.test.Common()
}

varComp.test.altDoTest = function(alt.fit, null=integer(0L), test='LinScore', control=varCompTest.control(test), null.fit=NULL, ...)
{
  # information=match.arg(information, informationTypes)
  varCompScoreTests = c(score='LinScore',score='VM03',score='SS95')
  varCompTests=c(varCompScoreTests,score='HP01',RLRT='RLRT') #,'Wald', 'RWD88')
  test=match.arg(test, varCompTests, several.ok=TRUE)
  test[test=='HP01']='SS95'; test=unique(test)
  ntest=length(test)
  if(ntest == 0L) stop(sprintf("Test needs to be given among [%s].", paste0(varCompTests, collapse=',')))
  null=as.integer(round(null))
  nNull=length(null)
	
    if( !isTRUE(alt.fit$doFit) ) {
		if( any(test %in% c('RLRT','RWD88','Wald')) ) { ## needs a full fit of alternative
			alt.fit = doFit.varComp(alt.fit)
		}else if( any(test %in%varCompScoreTests)){ ## needs a quick fit of alternative
			tmp.control=bak.control=alt.fit$control
			tmp.control$nlminb$iter.max=0L
			tmp.call=alt.fit$doFit
			tmp.call[['control']] = tmp.control
			alt.fit=eval(tmp.call)
			alt.fit$control = bak.control
		}
	}

	nK = length(alt.fit$working.cor)
	if(any(null<1L | null>nK) | nNull>=nK) {
		stop('The number of null components should be less than the number of non-error variance components')		
	}
	tau.idx=seq_len(nK)
	if(nNull>0L) tau.idx=tau.idx[-null]
	n.tau=length(tau.idx)
  
	if( !isTRUE(null.fit$doFit)) {
		if(any(test%in%c(varCompScoreTests,'RLRT') ) ) { ## full null fit
			if(is.null(null.fit)) {
				null.fit = varComp.fit(Y=alt.fit$residual.contrast, K=alt.fit$working.cor[null], control=alt.fit$control)
			}else null.fit=doFit.varComp(null.fit)
		}else{
			null.fit=NULL
		}
	}

	environment(varComp.test.Common) = sys.frame(sys.nframe())
	varComp.test.Common()
}

#' @importFrom Matrix nearPD
varComp.test.Common = evalq(function()  ## not to be called directly! 
{	
	# This function requires external objects: 
	#	"test", "varCompScoreTests", "control", "alt.fit",  "null.fit", 
	#	"null", "updateNegGrad",  "updateNegHess", "updateTr2", "tau.idx"

	call.lists=list(LinScore = c('infoMat', 'all.scores','tr1','n','LIkLI','non.pd'),
				 SS95 = c('infoMat', 'all.scores', 'tr1','n','LIkLI','LIy'),
				 RLRT = c('alt.fit', 'null.fit'),
				 VM03 = c('infoMat', 'all.scores', 'tr1','n','LIkLI','LIy'))
	#environment(char2list) = sys.frame(sys.nframe())
	varCompScoreTests = c(score='LinScore',score='VM03',score='SS95')
  
	if(any(test%in%varCompScoreTests)){
		{
			 environment(V)=
			 environment(updateLI)=
			 # environment(updateLIk)=
			 environment(updateLIkLI)=
			 environment(updateLIy)=
			 environment(PREML)=
			 environment(preprocPREML)=
			 environment(obj)=
			 environment(obj2)=
			 environment(updateNumsPart)=
			 environment(updateNums)=
			 environment(updateDenom)=
			 environment(updateTr1)=
			 environment(updateTr2)=
			 environment(updateNums2)=
			 environment(score)=
			 environment(updateNegGrad)=
			 environment(gradFunc)=
			 environment(OI)=
			 environment(hess)=
			 environment(EI)=
			 environment(AOI)=
			 environment(AEI)=
			 environment(WAI)=
			 environment(updateNegHess)=sys.frame(sys.nframe())
			infoFunc=get(control$information)
			preprocScore=expression({
				updateLIkLI()
				updateNumsPart();   
				updateNums(); updateDenom(); updateTr1()
			})
			preprocInfo=switch(control$information,
				OI=expression({ updateNums2(); updateTr2() }),
				EI=expression({ updateTr2();  }),
				AOI=expression({updateNums2();    }),
				AEI=expression({updateNums2();   }),
				WAI=expression({updateNums2();   })
			)			 
		}

		y=alt.fit$residual.contrast
		n=length(y); diag.1.n = diag(1, n)
		k=alt.fit$working.cor
		nK=length(k)
		
		LIkLI=k; LIy=numeric(n); numsPart=matrix(NA_real_, n, nK); 
		tr2=negHess=nums2=matrix(NA_real_, nK, nK)
		tr1=nums=negGrad=numeric(nK); denom=numeric(1L)
		LI=matrix(NA_real_, n, n)
		if(nK==1) {
			eigK=eigen(k[[1]],TRUE)
			eigK$tvector=t(eigK$vector)
		}
		
		tau=numeric(nK); tau[null]=null.fit$parms
		preprocPREML(tau)
		eval(preprocScore)
		eval(preprocInfo)
		updateNegGrad()
		updateNegHess()
		non.pd=FALSE
		if(control$information!='EI' && (min(eigen(negHess,TRUE,TRUE)$value)< -1e-6 || any(diag(as.matrix(negHess))<0) ) ){
			non.pd=TRUE		
			updateTr2();
			negHess=EI()
		}
		infoMat = as.matrix(nearPD(negHess)$mat)    ## added line on 062212
		sigma2=null.fit$sigma2   #sigma2=crossprod(LIy)/n
		all.scores= - negGrad    #score()
	}else if("RLRT"%in%test){
	}else stop("should not reach this line")
  
  test.list=vector('list', length(test)); names(test.list)=test
  for(i.test in seq_along(test)){
    # this.control=control[[test[i.test]]] ## do NOT move this into the call below, because the name 'control' has conflict.
	# this.fun=paste("varComp", test[i.test], "test", sep='.')
	# thisToCall=call(this.fun, control=this.control
          # ,null=null, null.fit=null.fit, alt.fit=alt.fit, X=X, Y=Y, K=K, k=k, lin.form=lin.form, LIkLI=LIkLI, LIy=LIy, tr1=tr1
          # ,infoMat=negHess, negHess=negHess, all.scores=all.scores, n=n, tau.idx=tau.idx, non.pd=non.pd
          # ,obs.score=obs.score)
		  
    # test.list[[i.test]]=eval(thisToCall)
	test.list[[i.test]] = do.call(
		what = paste('varComp', test[i.test], 'test', sep='.'),
		args = c(mget(call.lists[[test[i.test]]], inherits=TRUE), list(tau.idx = tau.idx, control=control[[test[i.test]]]))
    )
  }
  class(test.list)='varComp.test'
  test.list
}, refugeEnvironment)

is.formula=function(x)inherits(x, 'formula')

get.seed <-
function(){
#Retrieve the current random number seed the RNGkind() for reproducibility. 
  if(!exists('.Random.seed', envir=globalenv(), mode='numeric', inherits=FALSE)) runif(1L)
  seed=get('.Random.seed', envir=globalenv(), mode='numeric', inherits=FALSE)
  attr(seed, 'RNGkind')=RNGkind()
  seed
}

safeseq=function(from=1L, to=1L, by=1L,...)
{
  disc=by*(from-to)
  if(disc>0){
    vector(class(disc), 0L)
  }else seq(from=from, to=to, by=by, ...)
}

if(FALSE){# mget replaces this
char2list=function(Names)
{
	pfm=parent.frame(); pfm
	ans=lapply(Names, get, envir=pfm)
	names(ans) = Names
	ans
}

}

#' @importFrom quadprog solve.QP
#' @importFrom Matrix nearPD
varComp.VM03.test <-
function(control, infoMat, tau.idx, LIkLI, tr1, n, LIy, all.scores)# , ...)
{#	varComp.VM03.test: Actual testing function for VM03 test
#i.	control: the control object from varCompTest.control. 
#ii.	infoMat: Negative expected Hessian matrix.
#iii.	tau.idx: Integer vector that indexes the additional variance components of the alternative hypothesis compared to the null. That is, this is the complement of the null argument of varComp.test.
#iv.	LIkLI: A list of matrices being  sqrt(V^{-1})' K_i sqrt(V ^{-1})
#v.	tr1: Vector of trace terms. 
#vi.	n: residual sample size. 
#vii.	LIy: sqrt(V^{-1})'y
#viii.	null: the same as in varComp.test. 
#ix.	all.scores: vector of all score statistics. 
#x.	...: place holder for unused arguments. 
  null = seq_along(all.scores)[-tau.idx]
  nNull=length(null)
  n.tau=length(tau.idx)
  if(length(null)==0L) methods=control$method[[1]] else methods=control$method[[2]]
  if(nNull>0L && length(setdiff(methods, 'pboot'))==0L) { ## FIXME: implement pboot
 	warning("parametric bootstrap has not been implemented yet for the situation when the null model is not a linear model")
   null.htest=structure(list(p.value=NA_real_),class='htest')
    ans=vector('list', length(methods)); names(ans)=methods
    for(m in methods) ans[[m]]=null.htest
    class(ans)='varComp.VM03.test'; return(ans) 
  }
  #sumLIkLI=mapply('*',w,LIkLI[tau.idx])
#  invNegHess=solve(infoMat)[tau.idx, tau.idx,drop=FALSE]  ## this should be the same as below
  if(nNull>0L){ #warning('VM03 is not seriously tested when null!=integer(0) and should be used with extreme caution!')
    Phi=infoMat[-tau.idx, tau.idx, drop=FALSE]
    Delta=infoMat[-tau.idx, -tau.idx, drop=FALSE]
    partNegHess=infoMat[tau.idx,tau.idx,drop=FALSE]- crossprod(Phi, solve(Delta, Phi))
    invNegHess=tryCatch(solve(partNegHess), error=function(e){partNegHess<<-as.matrix(nearPD(partNegHess)$mat);solve(partNegHess)})
    
    mean.nonNull.score=drop(crossprod(Phi, solve(Delta, all.scores[-tau.idx])))

    adj.nonNull.score=all.scores[tau.idx] - mean.nonNull.score
    Amat=diag(1,n.tau); bvec=rep(0,n.tau)
    fact=1; i=0L
    repeat{
      i=i+1
      tmp=try( solve.QP(Dmat=invNegHess/fact, dvec=drop(invNegHess%*%all.scores[tau.idx])/fact, Amat=Amat, bvec=bvec, meq=0L, factorized=FALSE)$value )
      if(class(tmp)!='try-error'){
        qp.term=tmp*fact
        break
      }else if(i==20L){
        null.htest=structure(list(p.value=NA_real_),class='htest')
        ans=vector('list', length(methods)); names(ans)=methods
        for(m in methods) ans[[m]]=null.htest
        class(ans)='varComp.VM03.test'; 
        warning('Quadratic solver error')
        return(ans) 
      }else fact=fact*2
    }
    SSstat0=max(0, -2*qp.term      )
  }else{
    partNegHess=infoMat; invNegHess=solve(infoMat)
    Amat=diag(1,n.tau); bvec=rep(0,n.tau)
    SSstat=function(z){
      S=.5* (sapply(LIkLI, function(lik) drop(crossprod(z, lik%*%z)/crossprod(z)))*n-tr1 )
      SInvNegHess=crossprod(invNegHess,S)
      #max(0, -2*solve.QP(Dmat=invNegHess, dvec=drop(SInvNegHess), Amat=Amat, bvec=bvec, meq=0L, factorized=FALSE)$value )
      sol=drop(solve.QP(Dmat=invNegHess, dvec=drop(SInvNegHess), Amat=Amat, bvec=bvec, meq=0L, factorized=FALSE)$solution)
      max(0, crossprod(sol, invNegHess%*%sol))
    }
    SSstat0=drop(SSstat(LIy))
    nsim=control$nsim
  }
  ans=vector('list', length(methods)); names(ans)=methods
  for(m.i in seq_along(methods)){
    if(methods[m.i]=='pboot'){
      if(nNull==0L){
        seed=get.seed()
        SSstats=replicate(nsim-1L, SSstat(rnorm(n)))
        pval=(1+sum(SSstats>=SSstat0))/(nsim)
        ans[[m.i]]=structure(list(statistic=c(VM03=SSstat0), p.value=pval, null.stats=SSstats, seed=seed, parameters=c(nsim=nsim)), class='htest')      
      }else{
        .NotYetImplemented()
      }
    }else if(methods[m.i]=='ChiBarSq'){
      if(nNull==0L){
        seed=get.seed()
        pval=pchibarsq(SSstat0, partNegHess, lower.tail=FALSE)
        ans[[m.i]]=structure(list(statistic=c(VM03=SSstat0), p.value=pval, seed=seed, parameters=partNegHess), class='htest')
      }else{
        if(length(tau.idx)==1L){  ## this is a misnomer, as the mixture is not a chi-bar-square
          pval=if(SSstat0==0) 1 else pnorm(sqrt(SSstat0), mean.nonNull.score*sqrt(invNegHess), lower.tail=FALSE)
          ans[[m.i]]=structure(list(statistic=c(VM03=SSstat0), p.value=pval, parameters=c(mean=mean.nonNull.score, var=partNegHess)), class='htest')
        }else{
          .NotYetImplemented()
        }
      }
    }else if(methods[m.i]=='Bound'){
      if(nNull==0L){ ## the same as ChiBarSq
        seed=get.seed()
        pval=pchibarsq(SSstat0, partNegHess, lower.tail=FALSE)
        ans[[m.i]]=structure(list(statistic=SSstat0, p.value=pval, seed=seed, parameters=partNegHess), class='htest')
      }else{
        .NotYetImplemented()
      }
    }else if(methods[m.i]=='SSChiBarSq'){
      if(nNull==0L){
        seed=get.seed()
        SSstats=replicate(nsim, SSstat(rnorm(n)))
        m=mean(SSstats); s=sd(SSstats)
        theo.mom=mchibarsq(partNegHess, order=1:2)
        theo.sd=sqrt(theo.mom[2]-theo.mom[1]^2)
        SSstat0.adj=(SSstat0-m)/s*theo.sd+theo.mom[1]
        pval=pchibarsq(SSstat0.adj, partNegHess, lower.tail=FALSE)
        ans[[m.i]]=structure(list(statistic=c(VM03=SSstat0, VM03.adj=SSstat0.adj), p.value=pval, 
                                  null.stats=SSstats, seed=seed, 
                                  parameters=c(nsim=nsim, emp.mean=m, emp.sd=s, 
                                                theo.mean=theo.mom[1], theo.sd=theo.sd)), 
                                  
                             class='htest')    
      }else{
        .NotYetImplemented()
        seed=get.seed()
        pval=pchibarsq(SSstat0, partNegHess, lower.tail=FALSE)
        ans[[m.i]]=structure(list(statistic=SSstat0, p.value=pval, seed=seed, parameters=partNegHess), class='htest')
      }      
    }else stop('Undefined methods for VM03')
  }
  class(ans)='varComp.VM03.test'
  ans
}

#' @importFrom MASS ginv
vcov.varComp <-
function(object, what=c('fixed','beta','random','varComp','var.ratio','tau','response','Y'), drop=TRUE, beta.correction=TRUE, ...)
{
# S3 method for returning variance-covariance matrix of interesting quantities. 
# what: 'fixed' or 'beta': vcov for fixed effect parameter estimates
#       'random' or 'varComp': vcov for REML estimates
#       'var.ratio' or 'tau': vcov for ratios of variances
#       'response' or 'Y': vcov of response variable. 
  what=match.arg(what)
  what=switch(what, fixed='beta', random='varComp', var.ratio='tau', response='Y', what)
  
  if(what=='Y'){
    K=model.matrix(object, what='K')
    V=Reduce('+', mapply('*', object$varComps, K, SIMPLIFY=FALSE))
	if(is.null(V)){  ## fixed effect only model
		if('weights'%in%names(object)) V = diag(object$sigma2 / object$weights, nobs(object)) else V=diag(object$sigma2, nobs(object))
	}else{
		if('weights'%in%names(object)) V = V + diag(object$sigma2 / object$weights, nobs(object)) else V = V + diag(object$sigma2, nobs(object))
	}
    return(V)    
  }else if(what=='varComp'){  
    if(isTRUE(drop)) idx=which(object$parms>0) else idx=seq_along(object$parms)
  
    this.Q2Y=object$residual.contrast
    k=object$working.cor[idx]
    nk=length(k)
    k[[nk+1L]]=diag(1,length(this.Q2Y))
    V=Reduce('+', mapply('*', c(object$parms[idx],1)*object$sigma2, k, SIMPLIFY=FALSE))
    # LI=solve(t(chol(V)))
	LI = t(backsolve(chol(V), diag(1, nrow(V))))
    lik=lapply(k, function(kk) tcrossprod(LI%*%kk, LI))
    liy=LI%*%this.Q2Y
    
    ans=matrix(NA_real_, nk+1L, nk+1L)
		rownames(ans)=colnames(ans) = c(object$random.labels[idx], 'error')
    for(i in seq_len(nk+1)){
      for(j in i:(nk+1)){
        tmp=lik[[i]]%*%lik[[j]]
        ans[i,j]=ans[j,i]=-crossprod(liy, tmp%*%liy)+.5*sum(diag(tmp))
      }
    }
    ans=-solve(ans)
    #attr(ans, 'V')=V
    ans
  }else if(what=='tau'){
    if(isTRUE(drop)) idx=which(object$parms>0) else idx=seq_along(object$parms)
    if(length(idx)>0L) ans=-solve(object$hessian[idx,idx, drop=FALSE]) else ans = matrix(NA_real_, 0L, 0L)
	rownames(ans)=colnames(ans)=object$random.labels[idx]
    ans 
  }else if(what=='beta'){
    X=model.matrix(object, what='X')
	if(ncol(X)==0L) return(matrix(0,0L,0L))
    ans = ginv(crossprod(X, solve(vcov(object, what='Y'),X)))
	if(isTRUE(beta.correction)){
		ans = attr( KR.varComp(object, Lmat=matrix(0,0L,ncol(X)), Vbet=ans), 'vcov.beta' )
	}
	if(!is.null(colnames(X))) rownames(ans) = colnames(ans) = colnames(X)
	ans
  }
}
