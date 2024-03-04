#' Function to carry out irregular MFPCA
#'
#' @param components restrict number of components to use
#' @param split whether to split the output
#' @param ... list of fpca results
#'
#' @return eigenfunctions and vectors from mfpca
#' @export
#'
irregMFPCA = function(components = NULL, split = F, ...){
  J = length(list(...))

  eigs = lapply(list(...), function(x) x$phi)
  xis = lapply(list(...), function(x) x$xiEst)
  eigs_all = do.call(cbind, eigs)
  XiEst = do.call(cbind, xis)

  ms = sapply(eigs, function(x) ncol(x))
  Mplus = sum(ms)

  N = nrow(XiEst)

  Zhat = (N-1)^(-1)*t(XiEst) %*% (XiEst)
  eigen_analysis = eigen(Zhat)

  cm = eigen_analysis$vectors

  cmat = matrix(0, Mplus, J*Mplus)

  for(j in 1:J){
    if(j == 1){
      cmat[1:ms[j], 1:Mplus] = cm[1:ms[j], ]
    } else{
      cmat[(ms[(j-1)]+1):(ms[(j-1)] + ms[j]), ((j-1)*Mplus+1):(j*Mplus)] = cm[(ms[(j-1)]+1):(ms[(j-1)] + ms[j]), ]
    }
  }

  psi = eigs_all %*% cmat

  spl = split(1:ncol(psi), ceiling(1:ncol(psi)/Mplus))
  list_psi = lapply(spl, function(c) psi[,c])
  stack = do.call(rbind, list_psi)

  rho = XiEst %*% cm


  if(is.null(components)){
    Xhat = rho %*% t(stack)
    if(split == T){
      spl = split(1:ncol(Xhat), ceiling(1:ncol(Xhat)/51))
      list_Xhat = lapply(spl, function(c) Xhat[,c])
      return(list("unstackpsi" = psi,
                  "stackpsi" = stack,
                  "rho" = rho,
                  "Dhat" = t(rho) %*% rho / N,
                  "Xhat" = list_Xhat))
    }else{
      return(list("unstackpsi" = psi,
                  "stackpsi" = stack,
                  "rho" = rho,
                  "Dhat" = t(rho) %*% rho / N,
                  "Xhat" = Xhat))
    }
  } else{
    cols = 1:components
    for(j in 2:J){cols = append(cols, 1:components + Mplus)}
    psi = psi[, cols]
    stack = stack[, 1:components]
    rho = rho[, 1:components]
    Xhat = rho %*% t(stack)
    if(split == T){
      spl = split(1:ncol(Xhat), ceiling(1:ncol(Xhat)/51))
      list_Xhat = lapply(spl, function(c) Xhat[,c])
      return(list("unstackpsi" = psi,
                  "stackpsi" = stack,
                  "rho" = rho,
                  "Dhat" = t(rho) %*% rho / N,
                  "Xhat" = list_Xhat))
    }else{
      return(list("unstackpsi" = psi,
                  "stackpsi" = stack,
                  "rho" = rho,
                  "Dsquarehat" = t(rho) %*% rho / N,
                  "Xhat" = Xhat))
    }
    # TODO: Standard deviations are already built into the scores
  }
}
