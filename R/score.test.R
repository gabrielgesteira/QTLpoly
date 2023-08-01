score.test = function(y, X, K, tau){
  nK = nonNull = length(K) # n of components and component r index
  nNull = nK-1 # n of null components
  null = seq(nK)[-nonNull] # null components index
  tau = c(tau,0)
  qr = qr(X)
  A = t(qr.Q(qr, complete = T)[,-seq(qr$rank)])
  y = A%*%y
  n = nrow(y)
  diag.1.n = diag(1, n)
  K = lapply(K, function(x) A%*%x%*%t(A))
  V = diag.1.n + Reduce("+", mapply("*", tau, K, SIMPLIFY = FALSE))
  
  # varComp:::varComp.test.varComp()
  # varComp:::varComp.test.altDoTest()
  # varComp:::varComp.test.Common()
  # varComp:::varComp.LinScore.test()
  # varComp:::updateLI()
  if (nK == 1L) {
    eigK = eigen(K[[1]], TRUE)
    eigK$tvector = t(eigK$vector)
    LI = crossprod(eigK$tvec, 1/sqrt(eigK$val * tau + 1) * 
                      eigK$tvec)
  } else LI = t(backsolve(chol(V), diag.1.n))
  
  LIy = LI %*% y
  # varComp:::updateLIkLI()
  LIKLI = lapply(K, function(x) LI %*% x %*% t(LI))
  
  # varComp:::updateNumsPart()
  numsPart = sapply(LIKLI, "%*%", LIy)
  
  # varComp:::updateNums()
  nums <- drop(crossprod(LIy, numsPart))
  
  # varComp:::updateDenom()
  denom <- drop(crossprod(LIy))
  
  # varComp:::updateTr1()
  tr1 <- sapply(LIKLI, function(z) sum(diag(z)))
  
  # varComp:::updateTr2()
  tr2 = matrix(NA_real_, nK, nK)
  
  for (i in 1:nK) {
    for (j in i:nK) {
      Kij = LIKLI[[i]] %*% LIKLI[[j]]
      tr2[i, j] <- tr2[j, i] <- sum(diag(Kij))
    }
  }
  
  # varComp:::updateNegHess
  negHess = 0.5/(n + 2) * (n * tr2 - outer(tr1, tr1))
  
  infoMat = as.matrix(Matrix::nearPD(negHess)$mat)
  all.scores = 0.5 * (n * nums/denom - tr1)
  
  # varComp:::varComp.LinScore.test() ####
  
  obs.score = all.scores[nonNull]
  
  if (nNull > 0L) {
    partNegHess = infoMat[nonNull, nonNull, drop = FALSE] -
      infoMat[nonNull, -nonNull, drop = FALSE] %*% solve(infoMat[-nonNull,
                                                                 -nonNull, drop = FALSE], infoMat[-nonNull, nonNull,
                                                                                                  drop = FALSE])
    invNegHess = tryCatch(solve(partNegHess), error = function(e) {
      partNegHess <- as.matrix(nearPD(partNegHess)$mat)
      solve(partNegHess)
    })
  } else {
    partNegHess = infoMat
    invNegHess = solve(infoMat)
  }
  w = (1/sqrt(diag(as.matrix(partNegHess))))/sum(1/sqrt(diag(as.matrix(partNegHess))))
  
  lin.form = sum(w * obs.score)
  
  if (nNull == 0L) {
    pval = davies(0, eigen(n * Reduce("+", mapply("*", 
                                                  w, LIKLI, SIMPLIFY = FALSE)), TRUE, TRUE)$val - 
                    sum(tr1 * w) - lin.form * 2, acc = 1e-08, 
                  lim = 1e6)
    ans = c(st = lin.form, pv = pval$Qq)
    
  } else {
    # varComp.LinScore.SSAS155()
    
    Phi = infoMat[null, nonNull, drop = FALSE]
    Delta = infoMat[null, null, drop = FALSE]
    var.nonNull.score = infoMat[nonNull, nonNull, drop = FALSE] - 
      crossprod(Phi, solve(Delta, Phi))
    mean.nonNull.score = crossprod(Phi, solve(Delta, all.scores[null]))
    mean.obs.score = sum(mean.nonNull.score * w)
    var.obs.score = drop(crossprod(w, var.nonNull.score %*% 
                                     w))
    mean.theo.score = 0
    var.theo.score = drop(crossprod(w, infoMat[nonNull, nonNull, 
                                               drop = FALSE] %*% w))
    adjusted.obs.score = sqrt(var.theo.score/var.obs.score) * 
      (lin.form - mean.obs.score) + mean.theo.score
    pval = davies(0, eigen(n * Reduce("+", mapply("*", w, LIKLI[nonNull], 
                                                  SIMPLIFY = FALSE)), TRUE, TRUE)$val - sum(tr1[nonNull] * 
                                                                                              w) - adjusted.obs.score * 2, acc = 1e-08, lim = 1e6)
    pval0 = davies(0, eigen(n * Reduce("+", mapply("*", w, LIKLI[nonNull], 
                                                   SIMPLIFY = FALSE)), TRUE, TRUE)$val - sum(tr1[nonNull] * 
                                                                                               w) - lin.form * 2, acc = 1e-08, lim = 1e6)
    
    ans = c(st = adjusted.obs.score, pv = pval$Qq)
  }
  as.matrix(ans)
}