##############################################################################
## A fast algorithm for the accelerated failure time model 
## with high-dimensional time-to-event data
## By Taehwa Choi and Sangbum Choi
## Update: August 28, 2020
##############################################################################
library(dplyr)
library(mvtnorm)

glike = function(x,h=0.1) h*log(1+exp(-x/h))
gkern = function(x,h=0.1) 1/(1+exp(x/h)) 

# soft-thresholding
S = function(z, lambda) {
  case_when(z > lambda ~ z - lambda, z < -lambda ~ z + lambda, TRUE ~ 0)
}

# Ln: objective function
fn = function(y, z, delta, beta, h=0.1){
  n = nrow(z)
  e = log(y) - c(z %*% beta)
  di = rep(delta, each = n)
  ei = rep(e, each = n)
  ej = rep(e, times = n)
  res = di * glike((ei - ej),h=h)
  
  return(sum(res))
} 

# Un: estimating equation
G = function(y, z, delta, beta, h = 0.1){
  z = as.matrix(z)
  n = nrow(z); p = ncol(z)
  e = log(y) - c(z %*% beta)
  
  di = rep(delta, each = n)
  ei = rep(e, each = n)
  ej = rep(e, times = n)
  
  idx1 = rep(1:n, each = n)
  idx2 = rep(1:n, times = n)
  
  zi = z[idx1,]
  zj = z[idx2,]
  
  pi = gkern(ei-ej, h=h)
  
  sp = split(data.frame(di * (zi - zj) * pi), idx2)
  g = t(sapply(sp,colSums))
  
  return(colSums(g))
}

G.vec = function(y, z, delta, beta, h = 0.1){
  z = as.matrix(z)
  n = nrow(z); p = ncol(z)
  e = log(y) - c(z %*% beta)
  
  di = rep(delta, each = n)
  ei = rep(e, each = n)
  ej = rep(e, times = n)
  
  idx1 = rep(1:n, each = n)
  idx2 = rep(1:n, times = n)
  
  zi = z[idx1,]
  zj = z[idx2,]
  
  pi = gkern(ei-ej, h=h)
  
  sp = split(data.frame(di * (zi - zj) * pi), idx2)
  
  return(t(sapply(sp,colSums)))
}

# Hn: Hessian
H = function(y, z, delta, beta, h = 0.1){
  z = as.matrix(z)
  n = nrow(z); p = ncol(z)
  e = log(y) - c(z %*% beta)
  
  di = rep(delta, each = n)
  ei = rep(e, each = n)
  ej = rep(e, times = n)
  
  idx1 = rep(1:n, each = n)
  idx2 = rep(1:n, times = n)
  
  zi = z[idx1,]
  zj = z[idx2,]
  
  pi = gkern(ei-ej, h=h)
  wi = di * pi * (1 - pi)/h
  
  return(crossprod(wi * (zi - zj), (zi - zj)))
}

# Simulation data 
simudata=function(n, beta0, c){
  require(mvtnorm)
  p = length(beta0)
  rho = 0.5; sigma = matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i,j] = rho^(abs(i-j))
    }
  }  
  z = rmvnorm(n, mean = rep(0, p), sigma = sigma)
  T = exp(z %*% beta0 + rchisq(n, 1))
  C = runif(n, 0, c)
  y = c(pmin(T, C))
  delta = c(ifelse(T <= C, 1, 0))
  
  return(data.frame(y, delta, z))
}

# Unpenalized smoothed AFT with Newton-Raphson
aft_nr = function(y, z, delta, init, bandwidth = NULL, eps = 1e-8, iter.max = 100)
{
  iter = 0
  if(is.null(bandwidth)) bandwidth = sd(resid(lm(log(y)[delta == 1] ~ z[delta == 1,])))*length(y)^(-0.26)
  if(sum(init == rep(0, ncol(z))) == ncol(z)) warning("Improper Initial value for NR !!") & break
  old.beta = init
  for (i in 1:iter.max)
  {
    old.f = fn(y, z, delta, old.beta, h = bandwidth)
    old.G = G(y, z, delta, old.beta, h = bandwidth)
    old.H = H(y, z, delta, old.beta, h = bandwidth)
    
    iter = iter + 1
    if (iter == iter.max) break & print("Maximum iteration reached !!")
    new.beta = old.beta - (solve(old.H) %*% old.G)
    if (max(abs(new.beta - old.beta)) < eps) break
    old.beta = new.beta
  }
  vec = G.vec(y, z, delta, new.beta, h = bandwidth)
  A = solve(old.H); B = t(vec) %*% vec
  ase = c(sqrt(diag(A %*% B %*% A)))
  
  return(list(beta = c(new.beta),
              hess = old.H, grad = c(old.G),
              like = c(old.f), ase = ase))
}

# Penalized smoothed AFT
aft_pen =  function(y, z, delta, beta0, lambda, bandwidth = NULL, r = 3.7, type, max.iter = 100, eps = 1e-8)
{
  fit = aft_nr(y, z, delta, beta0, bandwidth = bandwidth)
  init = fit$beta
  # Weights for (adaptive) lasso 
  if(type == "ALASSO") tilde.beta = init
  if(type == "LASSO") tilde.beta = rep(1, length(beta0))
  
  h = fit$hess
  g = fit$grad
  # Step 1
  x = chol(h)
  y =  forwardsolve(t(x), h%*%init - g)
  
  n =  length(y)
  p =  ncol(x)
  
  # standardize 
  u =  y - mean(y)
  z =  t(t(x) - apply(x, 2, mean))
  norm.z =  apply(z^2, 2, mean)
  z =  t(t(z)/sqrt(norm.z)) #scale(x)
  init = coef(lm(u ~ -1 + z))
  init = ifelse(is.na(init), eps, init)
  
  # initialize beta
  beta =  init
  
  # residual update
  resid =  (u - z %*% beta)
  
  # start update
  for (t in 1:max.iter)
  {
    new.beta =  beta
    for (j in 1:p)                                  # Step 4
    {
      zj = crossprod(z[,j], resid)/n + beta[j]     # Step 2
      
      new.beta[j] = if (type == "SCAD") {          # Step 3
        S(zj, lambda)*(abs(zj) <= 2*lambda) + 
          S(zj, r*lambda/(r-1))/(1-1/(r-1))*(2*lambda < abs(zj) & abs(zj) <= r*lambda) +
          zj*(abs(zj) > r*lambda)
      } else {
        S(zj, lambda/abs(tilde.beta[j])) 
      }
      
      resid =  resid - z[,j] * (new.beta[j] - beta[j]) 
    }     
    # Step 5
    if (max(abs(beta - new.beta)) < eps) break
    beta =  new.beta
  }
  
  # transform back
  return(beta / sqrt(norm.z))
}

n = 300; c = 10
beta0 = c(0.7, 0.7, 0, 0, 0, 0.7, 0, 0, 0)
set.seed(1)
dt = simudata(n, beta0, c)
y = dt$y; delta = dt$delta; z = as.matrix(dt[,-c(1,2)])

# Plug-in bandwdith (Heller, 2007)
bandwidth = sd(resid(lm(log(y)[delta == 1] ~ z[delta == 1,])))*length(y)^(-0.26)

# Initialize with uncensored data
init = coef(lm(log(y)[delta==1] ~ z[delta==1,]))

# LASSO 
lambda = 5
lasso_coef = aft_pen(y, z, delta, beta0, lambda = lambda, type = "LASSO")
lasso_id = which(lasso_coef != 0)
A = lambda * diag(1 / (abs(lasso_coef[lasso_id])))
hess = H(y, z[,lasso_id], delta, lasso_coef[lasso_id], h = bandwidth)
vec = G.vec(y, z[,lasso_id], delta, lasso_coef[lasso_id], h = bandwidth)
A = solve(hess + n*diag(A))
B = t(vec) %*% vec
lasso_ase = sqrt(diag(A %*% B %*% A))

#SCAD
lambda = 10
scad_coef = aft_pen(y, z, delta, beta0, lambda = lambda, type = "SCAD")
scad_id = which(scad_coef != 0)
A = lambda * diag(1 / (abs(scad_coef[scad_id])))
hess = H(y, z[,scad_id], delta, scad_coef[scad_id], h = bandwidth)
vec = G.vec(y, z[,scad_id], delta, scad_coef[scad_id], h = bandwidth)
A = solve(hess + n*diag(A))
B = t(vec) %*% vec
scad_ase = sqrt(diag(A %*% B %*% A))

# Coefficients
cbind(lasso_coef, scad_coef)

# Standard error estimation
lasso_ase
scad_ase
