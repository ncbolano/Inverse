# Constructs variance estimator


#############################################
# Construction of J matrix
# Input x vector; ensure length n is a multiple of 2

# JJ0 = mvfft(x)/sqrt(nrow(x))
# n = nrow(x)
# backward = c(((n/2)+1):n)
# JJ = rbind(JJ0[backward,],JJ0[c(1:(n/2)),])

## JJ rearranges the frequencies. In particular,
#  the zero frequency is shifted to n/2.

#############################################
# extract J-nu vector

#' Uniform kernel
#'
#' @param x Input index
#'
#' @return Returns 1 for values ( in absolute terms) less or equal 1, else 0.
#' @export
#'
#' @examples
Kernel_Uniform=function(x)
{
  return(as.numeric(abs(x)<=1))
}

#' Triangular kernel
#'
#' @param x Input index
#'
#' @return Returns 1-|x| for values ( in absolute terms) less or equal 1, else 0.
#' @export
#'
#' @examples
Kernel_Triangular=function(x)
{
  return(as.numeric(abs(x)<=1)*(1-abs(x)))
}


' Triangular kernel
#'
#' @param x Input index
#'
#' @return Returns 1-|x| for values ( in absolute terms) less or equal 1, else 0.
#' @export
#'
#' @examples
Kernel_Triangular=function(x)
{
  return(as.numeric(abs(x)<=1)*(1-abs(x)))
}

Kernel_Quadtratic=function(x)
{
  return(as.numeric(abs(x)<=1)*(1-abs(x)^2))
}

J.to.J.k.nu <- function(J, k, nu) {
  n <- nrow(J)
  J1 <- J[(k + (-nu:nu)), ]
  J2 <- c(t(J1))
  return(J2)
}

J.to.J.k.nu.a <- function(J, k, nu, a) {
  n <- nrow(J)
  p <- ncol(J)
  Jtemp <- c(t(J[(k + (-nu:nu)), ]))
  Jtemp2 <- Jtemp[-(p * nu + a)] # is this right, seems to take all entries except one , rather than one
  return(Jtemp2)
}

# Make R-hat vector

Hat.R <- function(J, k, nu, W, M) {
  n <- nrow(J)
  p <- ncol(J)

  R_nu_M <- matrix(0, p * (2 * nu + 1), p * (2 * nu + 1))

  for (s in (-M:M))
  {
    J_k_S <- J.to.J.k.nu(J, k + s, nu)
    R_nu_M <- R_nu_M + W(s/M)*J_k_S %*% t(Conj(J_k_S))
  }
  norm=sum(sapply((-M:M),function(s) W(s/M)))
  return(R_nu_M / norm)
}

Hat.R.reduced <- function(J, k, nu, W,M1, a) {
  n <- nrow(J)
  p <- ncol(J)
  R_nu_M <- matrix(0, (p * (2 * nu + 1) - 1), (p * (2 * nu + 1) - 1))

  for (s in (-M1:M1))
  {
    J_k_S <- J.to.J.k.nu.a(J, k + s, nu, a)
    R_nu_M <- R_nu_M + W(s/M1)*J_k_S %*% t(Conj(J_k_S))
  }
  norm=sum(sapply((-M1:M1),function(s) W(s/M1)))
  return(R_nu_M / norm)
}

Hat.r <- function(J, k, nu,W, M, a) {
  n <- nrow(J)
  p <- ncol(J)

  r_nu_M <- rep(0, (p * (2 * nu + 1) - 1))

  for (s in (-M:M))
  {
    J_k_S <- J.to.J.k.nu.a(J, k + s, nu, a)
    r_nu_M <- r_nu_M + W(s/M)* J_k_S * Conj(J[k + s, a])
  }
  norm=sum(sapply((-M:M),function(s) W(s/M)))
  return(r_nu_M / norm)
}

beta <- function(J, k, nu,W, M, a, delta) {
  HatR <- Hat.R.reduced(J, k, nu, W,M, a)
  if(delta>0)
  {
    dim1 <- ncol(HatR)
    HatR <-HatR+ diag(rep(delta, dim1))
    
  }
  hatr <- Hat.r(J, k, nu, W,M, a)
  betaC <- solve(a = HatR,b = hatr)
  return(beta_beta_r(betaC,nu,ncol(J)))
  # return(betaC)
}

#(beta_r+beta_{-r})/2
beta_beta_r=function(beta,nu,p)
{
  tmp1=as.vector(unlist(sapply((-nu):nu,function(i) 
  {
    if(i==0)
      return(
        (nu)*p+(1:(p-1))
      )
    if(i<0)
      return(
        (abs(i)+nu)*p+1:p-1
      )
    if(i>0)
      return(
        (nu-abs(i))*p+1:p
      )
  }
  )))
  tmp2=seq_along(beta)
  beta_erg=beta
  beta_erg[tmp2!=tmp1]=((beta+Conj(beta[tmp1]))[tmp2!=tmp1])/2
  return(beta_erg)
}


# beta_cv=function(J, k, nu, M_grid, a, delta)
# {
#   n_points=floor(max(M_grid)/(2*nu+1)/3)
#   tmp=sapply(M_grid,function(M)
#   {
#     J_tmp=J
#     #Computational less demanding; one beta for all left out frequencies
#     #more demanding, each left ouf beta gets its own frequency
#     J_tmp[k+(-n_points):n_points*(2*nu+1),a]=0
#     beta_tmp=beta(J_tmp,k,nu,M,a,delta)
#     beta_tmp2=beta(J,k,nu,M,a,delta)
#     Eval=sapply(k+(-n_points):n_points*(2*nu+1),function(k_i)
#     {
#       J[k_i, a]-J.to.J.k.nu.a(J, k_i, nu, a)%*%Conj(beta_tmp)
#     })
#     # Eval2=sapply(k+(-n_points):n_points*(2*nu+1),function(k_i)
#     # {
#     #   J[k_i, a]-J.to.J.k.nu.a(J, k_i, nu, a)%*%Conj(beta_tmp2)
#     # })
#     # return(mean(abs(Eval)^2)/mean(abs(Eval2)^2))
#     return(mean(abs(Eval)^2))
#   })
#   names(tmp)=M_grid
#   return(tmp)
# }


beta_cv2=function(J, k, nu,W, M_grid, a, delta)
{
  n_points=floor(median(M_grid)/(2*nu+1))
  tmp=sapply(M_grid,function(M)
  {
    J_tmp=J
    #Computational less demanding; one beta for all left out frequencies
    #more demanding, each left ouf beta gets its own frequency
    Eval=sapply(k+(-n_points):n_points*(2*nu+1),function(k_i)
    {
      J_tmp[k_i,a]=0
      HatR <- Hat.R.reduced(J, k_i, nu,W, M, a)
      hatr <- Hat.r(J_tmp, k_i, nu,W, M, a)
      beta_tmp <- solve(a = HatR,b = hatr)
      # =beta(J_tmp,k_i,nu,M,a,delta)
      J[k_i, a]-J.to.J.k.nu.a(J, k_i, nu, a)%*%Conj(beta_tmp)
    })
    return(mean(abs(Eval)^2))
  })
  names(tmp)=M_grid
  return(tmp)
}


# beta_cv3=function(J, k, nu, M_grid, a, delta)
# {
#   n_points=floor(median(M_grid)/(2*nu+1))
#   tmp=sapply(M_grid,function(M)
#   {
#     # beta_tmp2=beta(J,k,nu,M,a,delta)
#     Eval2=sapply(k+(-n_points):n_points*(2*nu+1),function(k_i)
#     {
#       beta_tmp2=beta(J,k_i,nu,M,a,delta)
#       J[k_i, a]-J.to.J.k.nu.a(J, k_i, nu, a)%*%Conj(beta_tmp2)
#     })
#     return(mean(abs(Eval2)^2))
#   })
#   names(tmp)=M_grid
#   return(tmp)
# }



#########################################################################
############################################################################
# Real coefficient ###
#########################################################################

variance.estimator.Re.v2 <- function(J, k, nu, M, a, L, delta) {
  betaC <- beta(J, k, nu, M, a, delta)
  dim1 <- length(betaC)
  HatRred <- Hat.R.reduced(J, k, nu, M, a)
  ridge <- diag(rep(delta, dim1))
  invHatRreduced <- solve(HatRred + ridge)
  # invHatRreduced = solve(Hat.R.reduced(J,k,nu,M1=M,a)+diag(rep(delta,dim1)))
  #          dim1 = ncol(invHatRreduced)
  Y <- matrix(rep(0, (2 * M + 1) * dim1), nrow = dim1)
  p <- ncol(J)

  for (s in (-M:M))
  {
    J_k_S <- J.to.J.k.nu.a(J, k + s, nu, a)
    temp1 <- (Conj(J[k + s, a]) - sum(betaC * Conj(J_k_S)))
    temp2 <- invHatRreduced %*% J_k_S
    Y[, s + M + 1] <- Re(temp1 * temp2 / (2 * M + 1))
  }

  sigma1 <- matrix(rep(0, dim1**2), ncol = dim1)

  for (s1 in (1:(2 * M + 1))) {
    for (s2 in (1:(2 * M + 1))) {
      if (abs(s1 - s2) < L + 1) sigma1 <- sigma1 + Y[, s1] %*% t(Y[, s2])
    }
  }

  sigma2 <- matrix(rep(0, dim1**2), ncol = dim1)
  k1 <- (k - n / 2 - 1)
  for (s1 in (-M:M)) {
    for (s2 in (-M:M)) {
      if (abs(s1 + s2 + 2 * k1) < L + 1) sigma2 <- sigma2 + Y[, s1 + M + 1] %*% t(Y[, s2 + M + 1])
    }
  }

  sigmaRe <- diag(sigma1 + sigma2)

  return(sigmaRe)
}

variance.estimator.Im.v2 <- function(J, k, nu, M, a, L, delta) {
  betaC <- beta(J, k, nu, M, a, delta)
  dim1 <- length(betaC)
  HatRred <- Hat.R.reduced(J, k, nu, M, a)
  ridge <- diag(rep(delta, dim1))
  invHatRreduced <- solve(HatRred + ridge)
  # invHatRreduced = solve(Hat.R.reduced(J,k,nu,M1=M,a)+diag(rep(delta,dim1)))
  #          dim1 = ncol(invHatRreduced)
  Y <- matrix(rep(0, (2 * M + 1) * dim1), nrow = dim1)
  p <- ncol(J)

  for (s in (-M:M))
  {
    J_k_S <- J.to.J.k.nu.a(J, k + s, nu, a)
    temp1 <- (Conj(J[k + s, a]) - sum(betaC * Conj(J_k_S)))
    temp2 <- invHatRreduced %*% J_k_S
    Y[, s + M + 1] <- Im(temp1 * temp2 / (2 * M + 1))
  }

  sigma1 <- matrix(rep(0, dim1**2), ncol = dim1)

  for (s1 in (1:(2 * M + 1))) {
    for (s2 in (1:(2 * M + 1))) {
      if (abs(s1 - s2) < L+1) sigma1 <- sigma1 + Y[, s1] %*% t(Y[, s2])
    }
  }

  sigma2 <- matrix(rep(0, dim1**2), ncol = dim1)
  k1 <- (k - n / 2 - 1)
  for (s1 in (-M:M)) {
    for (s2 in (-M:M)) {
      if (abs(s1 + s2 + 2 * k1) < L + 1) sigma2 <- sigma2 + Y[, s1 + M + 1] %*% t(Y[, s2 + M + 1])
    }
  }

  sigmaIm <- diag(sigma1 + sigma2)

  return(sigmaIm)
}


variance.estimator.v2 <- function(J, k, nu,W, M, a, L, delta) {
  betaC <- beta(J, k, nu, W,M, a, delta)
  dim1 <- length(betaC)
  HatRred <- Hat.R.reduced(J, k, nu,W, M, a)
  ridge <- diag(rep(delta, dim1))
  invHatRreduced <- solve(HatRred + ridge)
  norm=sum(sapply((-M:M),function(s) W(s/M)))
  # invHatRreduced = solve(Hat.R.reduced(J,k,nu,M1=M,a)+diag(rep(delta,dim1)))
  #          dim1 = ncol(invHatRreduced)
  Y <- matrix(rep(0, 2*(2 * M + 1) * dim1), nrow = 2*dim1)
  p <- ncol(J)
  
  for (s in (-M:M))
  {
    J_k_S <- J.to.J.k.nu.a(J, k + s, nu, a)
    temp1 <- (Conj(J[k + s, a]) - sum(betaC * Conj(J_k_S)))
    temp2 <- invHatRreduced %*% J_k_S
    Y[, s + M + 1] <- W(s/M)*c(Re(temp1 * temp2 / norm),Im(temp1 * temp2 / norm))
  }
  
  sigma1 <- matrix(rep(0, (2*dim1)**2), ncol = 2*dim1)
  
  for (s1 in (1:(2 * M + 1))) {
    for (s2 in (1:(2 * M + 1))) {
      if (abs(s1 - s2) < L + 1) sigma1 <- sigma1 + Y[, s1] %*% t(Y[, s2])
    }
  }
  
  sigma2 <- matrix(rep(0, (2*dim1)**2), ncol = 2*dim1)
  k1 <- (k - n / 2 - 1)
  for (s1 in (-M:M)) {
    for (s2 in (-M:M)) {
      if (abs(s1 + s2 + 2 * k1) < L + 1) sigma2 <- sigma2 + Y[, s1 + M + 1] %*% t(Y[, s2 + M + 1])
    }
  }
  
  return(sigma1+sigma2)
  # sigmaRe <- diag(sigma1 + sigma2)
  
  # return(sigmaRe)
}



###############################################################################
