
#' @useDynLib SIVCMTest
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import Rcpp
#' @importFrom MASS mvrnorm 
#' @importFrom stats rnorm rt runif sd
#' @importFrom nleqslv nleqslv

# ======================================================
# Generate data
# ======================================================

# Example 1
####   g and eta functions 
eta1 <- function(s){ sqrt(2)*sin(2*pi*s) }
eta2 <- function(s){ sqrt(2)*cos(2*pi*s) }
g <- function(xb){ sin(2*xb) + 2*cos(2+xb)}
g2 <- function(xb){xb}


#' @export
GenData.Sa <- function(n, m, a){
  
  p <- 3
  rho <- 0.6   # X, cov
  lambda <- c(1, 0.5)  # eta, std
  # sigma <- 0.3  # error standard variance
  
  # mean and covariance matrix for X
  muX <- matrix(0, 1, p)
  SigmaX <- matrix(0, p, p)
  SigmaX <- SigmaXC(rho, SigmaX)
  
  # true coefficients
  beta <- matrix(0, p, m)
  tm <- sort(runif(m, min = 0, max = 1))
  for (s in 1:m) {
    beta[1, s] <- 1 + tm[s]^2
    beta[2, s] <- -(1 - tm[s])^2
    beta[3, s] <- -1 + 4*(tm[s]-0.5)^2
    beta[ , s] <- beta[,s]/sqrt(sum(beta[,s]^2))
  }
  
  
  beta2 <- matrix(0, p, m)
  for (s in 1:m) {
    beta2[1, s] <- -(tm[s]-0.5)^2
    beta2[2, s] <- sqrt(tm[s])
    beta2[3, s] <- cos(3*pi*tm[s])
    beta2[ , s] <- beta2[,s]/sqrt(sum(beta2[,s]^2))
    
  }
  
  x <- mvrnorm(n, muX, SigmaX)
  eta <- matrix(sqrt(lambda[1])*rnorm(n), ncol=1) %*%  matrix(eta1(tm), nrow=1) +
    matrix(sqrt(lambda[2])*rnorm(n), ncol=1) %*%  matrix(eta2(tm), nrow=1)
  
  simu.g <- matrix(0, nrow=n, ncol=m)
  ally <- matrix(0, nrow=n, ncol=m)
  #for (s in 1:m) {
  #  xb <- x %*% beta[,s]
  #  xb2 <- cos(x) %*% beta2[,s]
  #  simu.g[,s] <- g(xb)
  #  ally[,s] <- simu.g[,s] + eta[,s] + a*exp(xb2)
  #}
  
  for (s in 1:m) {
    xb <- x %*% beta[,s]
    xb2 <- exp(x) %*% beta2[,s]
    simu.g[,s] <- g(xb)
    ally[,s] <- simu.g[,s] + eta[,s] + a*xb2
  }
  
  data = list(x=x, ally=ally, tm=tm, beta=beta)
  return(data)
  
}


GenData.Sb <- function(n, m, a){
  
  p <- 3
  rho <- 0.6   # X, cov
  lambda <- c(1, 0.5)  # eta, std
  
  # mean and covariance matrix for X
  muX <- matrix(0, 1, p)
  SigmaX <- matrix(0, p, p)
  SigmaX <- SigmaXC(rho, SigmaX)
  
  # true coefficients
  beta <- matrix(0, p, m)
  tm <- sort(runif(m, min = 0, max = 1))
  for (s in 1:m) {
    beta[1, s] <- 1 + tm[s]^2
    beta[2, s] <- -(1 - tm[s])^2
    beta[3, s] <- -1 + 4*(tm[s]-0.5)^2
    beta[ , s] <- beta[,s]/sqrt(sum(beta[,s]^2))
  }
  
  beta2 <- matrix(0, p, m)
  for (s in 1:m) {
    beta2[1, s] <- (tm[s]-0.5)^2-0.1
    beta2[2, s] <- sqrt(tm[s])
    beta2[3, s] <- sin(2*pi*tm[s])
    beta2[ , s] <- beta2[,s]/sqrt(sum(beta2[,s]^2))
  }  
  
  x <- mvrnorm(n, muX, SigmaX)
  eta <- matrix(sqrt(lambda[1])*rnorm(n), ncol=1) %*%  matrix(eta1(tm), nrow=1) +
    matrix(sqrt(lambda[2])*rnorm(n), ncol=1) %*%  matrix(eta2(tm), nrow=1)
  
  simu.g <- matrix(0, nrow=n, ncol=m)
  ally <- matrix(0, nrow=n, ncol=m)
  for (s in 1:m) {
    xb <- x %*% beta[,s]
    xb2 <- cos(0.4*pi*x) %*% beta2[,s]
    simu.g[,s] <- g(xb)
    ally[,s] <- simu.g[,s] + eta[,s] + a*xb2
  }
  
  data = list(x=x, ally=ally, simu.g=simu.g, tm=tm, beta=beta, eta=eta)
  return(data)
  
}



GenData.Sc <- function(n, m, a){
  
  p <- 3
  rho <- 0.6   # X, cov
  lambda <- c(1, 0.5)  # eta, std
  
  # mean and covariance matrix for X
  muX <- matrix(0, 1, p)
  SigmaX <- matrix(0, p, p)
  SigmaX <- SigmaXC(rho, SigmaX)
  
  # true coefficients
  beta <- matrix(0, p, m)
  tm <- sort(runif(m, min = 0, max = 1))
  for (s in 1:m) {
    beta[1, s] <- 1 + tm[s]^2
    beta[2, s] <- -(1 - tm[s])^2
    beta[3, s] <- -1 + 4*(tm[s]-0.5)^2
    beta[ , s] <- beta[,s]/sqrt(sum(beta[,s]^2))
  }
  
  beta2 <- matrix(0, p, m)
  for (s in 1:m) {
    beta2[1, s] <- -(tm[s]-0.5)^2
    beta2[2, s] <- sqrt(tm[s])
    beta2[3, s] <- cos(3*pi*tm[s])
    beta2[ , s] <- beta2[,s]/sqrt(sum(beta2[,s]^2))
  }
  x <- mvrnorm(n, muX, SigmaX)
  eta <- matrix(sqrt(lambda[1])*rnorm(n), ncol=1) %*%  matrix(eta1(tm), nrow=1) +
    matrix(sqrt(lambda[2])*rnorm(n), ncol=1) %*%  matrix(eta2(tm), nrow=1)
  
  simu.g <- matrix(0, nrow=n, ncol=m)
  ally <- matrix(0, nrow=n, ncol=m)
  for (s in 1:m) {
    xb <- x %*% beta[,s]
    xb2 <- exp(x) %*% beta2[,s]
    simu.g[,s] <- g2(xb)
    ally[,s] <- simu.g[,s] + eta[,s] + a*xb2
  }
  
  data = list(x=x, ally=ally, simu.g=simu.g, tm=tm, beta=beta, eta=eta)
  return(data)
  
}



GenData.Sd <- function(n, m, a){
  
  p <- 3
  rho <- 0.6   # X, cov
  lambda <- c(1, 0.5)  # eta, std
  
  # mean and covariance matrix for X
  muX <- matrix(0, 1, p)
  SigmaX <- matrix(0, p, p)
  SigmaX <- SigmaXC(rho, SigmaX)
  
  # true coefficients
  beta <- matrix(0, p, m)
  tm <- sort(runif(m, min = 0, max = 1))
  for (s in 1:m) {
    beta[1, s] <- 1 + tm[s]^2
    beta[2, s] <- -(1 - tm[s])^2
    beta[3, s] <- -1 + 4*(tm[s]-0.5)^2
    beta[ , s] <- beta[,s]/sqrt(sum(beta[,s]^2))
  }
  
  beta2 <- matrix(0, p, m)
  for (s in 1:m) {
    beta2[1, s] <- (tm[s]-0.5)^2-0.1
    beta2[2, s] <- sqrt(tm[s])
    beta2[3, s] <- sin(2*pi*tm[s])
    beta2[ , s] <- beta2[,s]/sqrt(sum(beta2[,s]^2))
  }  
  
  x <- mvrnorm(n, muX, SigmaX)
  eta <- matrix(sqrt(lambda[1])*rnorm(n), ncol=1) %*%  matrix(eta1(tm), nrow=1) +
    matrix(sqrt(lambda[2])*rnorm(n), ncol=1) %*%  matrix(eta2(tm), nrow=1)
  
  simu.g <- matrix(0, nrow=n, ncol=m)
  ally <- matrix(0, nrow=n, ncol=m)
  for (s in 1:m) {
    xb <- x %*% beta[,s]
    xb2 <- cos(0.4*pi*x) %*% beta2[,s]
    simu.g[,s] <- g2(xb)
    ally[,s] <- simu.g[,s] + eta[,s] + a*xb2
  }
  
  data = list(x=x, ally=ally, simu.g=simu.g, tm=tm, beta=beta, eta=eta)
  return(data)
  
}


################################ t dist ###################################

GenData.Sa.t <- function(n, m, a){
  
  p <- 3
  rho <- 0.6   # X, cov
  # lambda <- c(1, 0.5)  # eta, std
  # sigma <- 0.3  # error standard variance
  
  # mean and covariance matrix for X
  muX <- matrix(0, 1, p)
  SigmaX <- matrix(0, p, p)
  SigmaX <- SigmaXC(rho, SigmaX)
  
  # true coefficients
  beta <- matrix(0, p, m)
  tm <- sort(runif(m, min = 0, max = 1))
  for (s in 1:m) {
    beta[1, s] <- 1 + tm[s]^2
    beta[2, s] <- -(1 - tm[s])^2
    beta[3, s] <- -1 + 4*(tm[s]-0.5)^2
    beta[ , s] <- beta[,s]/sqrt(sum(beta[,s]^2))
  }
  
  
  beta2 <- matrix(0, p, m)
  for (s in 1:m) {
    beta2[1, s] <- -(tm[s]-0.5)^2
    beta2[2, s] <- sqrt(tm[s])
    beta2[3, s] <- cos(3*pi*tm[s])
    beta2[ , s] <- beta2[,s]/sqrt(sum(beta2[,s]^2))
    
  }
  
  x <- mvrnorm(n, muX, SigmaX)
  eta <- matrix(rt(n, df=3), ncol=1) %*%  matrix(eta1(tm), nrow=1) +
    matrix(rt(n, df=5), ncol=1) %*%  matrix(eta2(tm), nrow=1)
  
  simu.g <- matrix(0, nrow=n, ncol=m)
  ally <- matrix(0, nrow=n, ncol=m)
  
  
  for (s in 1:m) {
    xb <- x %*% beta[,s]
    xb2 <- exp(x) %*% beta2[,s]
    simu.g[,s] <- g(xb)
    ally[,s] <- simu.g[,s] + eta[,s] + a*xb2
  }
  
  data = list(x=x, ally=ally, simu.g=simu.g, tm=tm, beta=beta, eta=eta)
  return(data)
  
}


GenData.Sb.t <- function(n, m, a){
  
  p <- 3
  rho <- 0.6   # X, cov
  # lambda <- c(1, 0.5)  # eta, std
  
  # mean and covariance matrix for X
  muX <- matrix(0, 1, p)
  SigmaX <- matrix(0, p, p)
  SigmaX <- SigmaXC(rho, SigmaX)
  
  # true coefficients
  beta <- matrix(0, p, m)
  tm <- sort(runif(m, min = 0, max = 1))
  for (s in 1:m) {
    beta[1, s] <- 1 + tm[s]^2
    beta[2, s] <- -(1 - tm[s])^2
    beta[3, s] <- -1 + 4*(tm[s]-0.5)^2
    beta[ , s] <- beta[,s]/sqrt(sum(beta[,s]^2))
  }
  
  beta2 <- matrix(0, p, m)
  for (s in 1:m) {
    beta2[1, s] <- (tm[s]-0.5)^2-0.1
    beta2[2, s] <- sqrt(tm[s])
    beta2[3, s] <- sin(2*pi*tm[s])
    beta2[ , s] <- beta2[,s]/sqrt(sum(beta2[,s]^2))
  }  
  
  x <- mvrnorm(n, muX, SigmaX)
  eta <- matrix(rt(n, df=3), ncol=1) %*%  matrix(eta1(tm), nrow=1) +
    matrix(rt(n, df=5), ncol=1) %*%  matrix(eta2(tm), nrow=1)
  
  simu.g <- matrix(0, nrow=n, ncol=m)
  ally <- matrix(0, nrow=n, ncol=m)
  for (s in 1:m) {
    xb <- x %*% beta[,s]
    xb2 <- cos(0.4*pi*x) %*% beta2[,s]
    simu.g[,s] <- g(xb)
    ally[,s] <- simu.g[,s] + eta[,s] + a*xb2
  }
  
  data = list(x=x, ally=ally, simu.g=simu.g, tm=tm, beta=beta, eta=eta)
  return(data)
  
}



GenData.Sc.t <- function(n, m, a){
  
  p <- 3
  rho <- 0.6   # X, cov
  # lambda <- c(1, 0.5)  # eta, std
  
  # mean and covariance matrix for X
  muX <- matrix(0, 1, p)
  SigmaX <- matrix(0, p, p)
  SigmaX <- SigmaXC(rho, SigmaX)
  
  # true coefficients
  beta <- matrix(0, p, m)
  tm <- sort(runif(m, min = 0, max = 1))
  for (s in 1:m) {
    beta[1, s] <- 1 + tm[s]^2
    beta[2, s] <- -(1 - tm[s])^2
    beta[3, s] <- -1 + 4*(tm[s]-0.5)^2
    beta[ , s] <- beta[,s]/sqrt(sum(beta[,s]^2))
  }
  
  beta2 <- matrix(0, p, m)
  for (s in 1:m) {
    beta2[1, s] <- -(tm[s]-0.5)^2
    beta2[2, s] <- sqrt(tm[s])
    beta2[3, s] <- cos(3*pi*tm[s])
    beta2[ , s] <- beta2[,s]/sqrt(sum(beta2[,s]^2))
  }
  x <- mvrnorm(n, muX, SigmaX)
  eta <- matrix(rt(n, df=3), ncol=1) %*%  matrix(eta1(tm), nrow=1) +
    matrix(rt(n, df=5), ncol=1) %*%  matrix(eta2(tm), nrow=1)
  
  simu.g <- matrix(0, nrow=n, ncol=m)
  ally <- matrix(0, nrow=n, ncol=m)
  for (s in 1:m) {
    xb <- x %*% beta[,s]
    xb2 <- exp(x) %*% beta2[,s]
    simu.g[,s] <- g2(xb)
    ally[,s] <- simu.g[,s] + eta[,s] + a*xb2
  }
  
  data = list(x=x, ally=ally, simu.g=simu.g, tm=tm, beta=beta, eta=eta)
  return(data)
  
}



GenData.Sd.t <- function(n, m, a){
  
  p <- 3
  rho <- 0.6   # X, cov
  # lambda <- c(1, 0.5)  # eta, std
  
  # mean and covariance matrix for X
  muX <- matrix(0, 1, p)
  SigmaX <- matrix(0, p, p)
  SigmaX <- SigmaXC(rho, SigmaX)
  
  # true coefficients
  beta <- matrix(0, p, m)
  tm <- sort(runif(m, min = 0, max = 1))
  for (s in 1:m) {
    beta[1, s] <- 1 + tm[s]^2
    beta[2, s] <- -(1 - tm[s])^2
    beta[3, s] <- -1 + 4*(tm[s]-0.5)^2
    beta[ , s] <- beta[,s]/sqrt(sum(beta[,s]^2))
  }
  
  beta2 <- matrix(0, p, m)
  for (s in 1:m) {
    beta2[1, s] <- (tm[s]-0.5)^2-0.1
    beta2[2, s] <- sqrt(tm[s])
    beta2[3, s] <- sin(2*pi*tm[s])
    beta2[ , s] <- beta2[,s]/sqrt(sum(beta2[,s]^2))
  }  
  
  x <- mvrnorm(n, muX, SigmaX)
  eta <- matrix(rt(n, df=3), ncol=1) %*%  matrix(eta1(tm), nrow=1) +
    matrix(rt(n, df=5), ncol=1) %*%  matrix(eta2(tm), nrow=1)
  
  simu.g <- matrix(0, nrow=n, ncol=m)
  ally <- matrix(0, nrow=n, ncol=m)
  for (s in 1:m) {
    xb <- x %*% beta[,s]
    xb2 <- cos(0.4*pi*x) %*% beta2[,s]
    simu.g[,s] <- g2(xb)
    ally[,s] <- simu.g[,s] + eta[,s] + a*xb2
  }
  
  data = list(x=x, ally=ally, simu.g=simu.g, tm=tm, beta=beta, eta=eta)
  return(data)
  
}


# ======================================================
# Estimate beta
# ======================================================
# Input:
#   x: n*p matrix
#   ally: n*m matrix
#   s, m, hx, hy, h: scalar
#   tm: 1*m vector
#   beta0: p*1 matrix
# Output:
#   betaest: p*1 matrix
#

getBeta <- function(x, ally, m, tm, s, hx, hy, h, beta0){

  model <- function(beta){
    sumseffC(x, ally, m, tm, s, beta, hx, hy, h)
  }

  betaest <- nleqslv(beta0, model, method = 'Broyden', global = 'dbldog',
                     control=list(xtol = 0.01, ftol = 0.1))$x
  
  betaest <- matrix(betaest, ncol=1)
  betaest <- betaest / sqrt(sum(betaest^2))    # standarization
  betaest <- betaest * sign(betaest[1,1])

  return(betaest)

}

# ======================================================
# Compute test statistic
# ======================================================
anglefun <- function(n, p, X){
  
  angle0 <- matrix(0, nrow = n, ncol = n)
  if (p == 1){ 
    for (k in 1:n) {
      tmp.X <- matrix( rep( X[k,], n), nrow = n, ncol = p, byrow = T) 
      angle0.k <- (X - tmp.X <= 0) %*% t(X - tmp.X <= 0)
      angle0 <- angle0 + angle0.k
    }
    angle <- angle0 # B: n*n
  } else if (p > 1){
    for (k in 1:n) {
      tmp.X <- matrix( rep( X[k,], n), nrow = n, ncol = p, byrow = T)   # X : n*p
      nume <- (X - tmp.X) %*% t(X - tmp.X)
      nume[k,k] <- 1    # special cases
      deno <- sqrt( matrix(diag(nume), nrow = n) %*% matrix(diag(nume), ncol = n) )
      angle0.k <- abs(pi - acos(nume / deno))
      angle0.k[k,] <- pi      # special cases
      angle0.k[,k] <- pi
      angle0.k[k,k] <- 2 * pi
      angle0 <- angle0 + angle0.k
    }
    angle <- angle0 * pi^(p / 2 - 1) / gamma(p / 2 + 1) # B: n*n
  } 
  
  return(angle)
}


# ======================================================
# CV for h1 
# ======================================================
# Input:
#   x: n*p matrix
#   betaest: p*m matrix
#   allxb, ally: n*m matrix
# Output:
#   h1: scalar
#

cvh1 <- function(x, betaest, allxb, ally){
  
  n <- nrow(allxb)
  m <- ncol(allxb)
  srange <- max(allxb) - min(allxb)
  nh <- 20
  hmin <- srange/m
  hmax <- srange/20  
  vh <- seq(hmin, hmax, length.out = nh)
  
  vec_allxb <- matrix(allxb, ncol=1)
  vec_y <- matrix(ally, ncol=1)
  vec_cv <- matrix(0, nrow=n*m-1, ncol=1)
  
  CV <- matrix(0, nh, 1)
  gest <- matrix(0, n, m)
  CV <- cvC(vh, x, ally, betaest, vec_allxb, vec_y, vec_cv, gest, CV)
  
  flag <- which.min(CV)
  h1 <- vh[flag]
  
  return(h1)
  
}


# ======================================================
# Main function
# ======================================================
#' @export
SIVCMTest <- function(n,p,m,tm,x,ally,B){
  
  #### beta0 <- beta  # initial beta
  beta0 <- solve(t(x) %*% x) %*% (t(x) %*% ally)
  
  # bandwidth
  c1 <- mean(apply((x%*%beta0), 2, sd))
  hx <- n^(-1/3)*c1    # Eqn.(8)
  hy <- n^(-1/5)*c1    # Eqn.(9)
  h1 <- n^(-1/4)*m^(-1/5)*c1   # overfit
  h <- m^(-1/5)* sd(tm)  # Eqn.(11), w(s_m,s)
  
  # Estimate beta
  betaest <- matrix(0, nrow=p, ncol=m)
  for (s in 1:m) {
    betaest[, s] <- getBeta(x, ally, m, tm, s, hx, hy, h, beta0[,s,drop=F])
    betaest[, s] <- betaest[,s]*sign(betaest[1,s])
  }
  
  # Estimate g
  allxb <- x %*% betaest
  #### h1 <- cvh1(x, betaest, allxb, ally)   # cv for h1
  gest <- matrix(0, nrow=n, ncol=m)
  gest <- gestC(h1, x, betaest, matrix(allxb, ncol=1), matrix(ally, ncol=1), gest)
  
  # residual
  res <- ally - gest
  
  #====================================
  #   Step 3: Compute test statistic 
  #====================================
  
  angle <- anglefun(n, p, x)
  Tn <- TnC(res, angle, n, m)
  
  #=====================================
  #   Step 4: Resampling
  #=====================================
  
  NewTn <- matrix(NA, nrow=1, ncol=B)    # Bootstrap Tn
  
  for(b in 1:B){
    
    ei <- matrix(rep(rnorm(n),m), nrow = n, byrow = F)
    allyboot <- gest + ei * res  # Bootsrap response y
    
    # Estimate beta
    betaestboot <- matrix(0, nrow=p, ncol=m)
    beta0 <- solve(t(x) %*% x) %*% (t(x) %*% allyboot)
    for (s in 1:m) {
      betaestboot[, s] <- getBeta(x, allyboot, m, tm, s, hx, hy, h, beta0[,s,drop=F])
      betaestboot[, s] <- betaestboot[,s]*sign(betaestboot[1,s])
    }
    
    # Estimate g
    ### h1boot <- cvh1(x, betaestboot, x%*%betaestboot, allyboot)  # cv for h1, overfit, thus we do not use cv method
    gestboot <- matrix(0, nrow=n, ncol=m)
    gestboot <- gestC(h1, x, betaestboot, matrix(x%*%betaestboot, ncol=1), matrix(allyboot, ncol=1), gestboot)
    resboot <- allyboot - gestboot # bootstrap residual
    
    NewTn[b] <- TnC(resboot, angle, n, m) # bootstrap Tn
    
  } # for(b in 1:B)
  
  result <- list(TestStat = Tn, Pvalue = mean(NewTn>Tn), Cri95 = NewTn[floor(B*0.95)])
  
  return(result)
}







