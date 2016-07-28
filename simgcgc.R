library(Rcpp)
library(RcppArmadillo)
library(geoR)
library(VGAM)
library(MASS)
library(snowfall)

#### Simulating geoStatistical Count data (isotropic) in Gaussian copula models
#### With choice of marginals: NB2, ZIP (parameterized in terms of mean and over-dispersion)
#### Poisson and Binomial
#set.seed(1234)
#simtmp = sim.gcgc(locs = cbind(xloc,yloc), mu = 4, od = 10, sim.n = 3, 
#                  family = "NB2", range = 0.3, nugget = 0.1)
#set.seed(1234)
#simtmp.y = sim.gcgc(locs = cbind(xloc, yloc), mu = exp(0.5+0.5*xloc+0.5*yloc),
#                    od = 0.5, sim.n =1, family = "NB2", range = 0.3, nugget = 0)

sim.gcgc <- function(locs, mu, od, sim.n = 1, binomial.N = NULL, p = NULL, family, 
                     matern.kappa = 0.5, range, nugget = NULL){
  if(!is.matrix(locs) & !is.data.frame(locs))
    stop(" Input 'locs' must be a matrix or data frame!")
  D <- as.matrix(dist(locs, method = "euclidean", diag = T, upper = T))
  if(is.null(nugget)){
    R <- geoR::matern(D, range, matern.kappa)
  }else{
    R <- (1-nugget)*geoR::matern(D, range, matern.kappa) + nugget*diag(ncol(D))
  }
  if (requireNamespace("MASS", quietly = TRUE)) {
    sim.grf <- MASS::mvrnorm(sim.n+1, mu = rep(0, nrow(locs)), Sigma = R)
  }else{
    stop("Please install and load {MASS} first!")
  }
  if(family == "NB2"){
 #   binomial.N <- NULL; p <- NULL
    if(length(mu) == 1) mu <- rep(mu, nrow(locs))
    if(length(mu)!= nrow(locs)) stop("length of mu must be equal to number of locations")
    ans <- t(apply(sim.grf, 1, function (x)
      qnbinom(pnorm(x), size= rep(1/od, nrow(locs)), mu = mu)))[-(sim.n+1),]
  }
  if(family == "ZIP"){
  #  binomial.N <- NULL; p <- NULL
    if(length(mu) == 1) mu <- rep(mu, nrow(locs))
    if(length(mu)!= nrow(locs)) stop("length of mu must be equal to number of locations")
    ans <- t(apply(sim.grf,1,function(x)
      qzip2(pnorm(x), size = rep(1/od, nrow(locs)), mu = mu)))[-(sim.n+1),]
  }
  if(family == "Poisson"){
    od <- NULL
    if(length(mu) == 1) mu = rep(mu, nrow(locs))
    if(length(mu)!= nrow(locs)) stop("length of mu must be equal to number of locations")
    ans <- t(apply(sim.grf,1,function(x)
      qpois(pnorm(x), lambda = mu)))[-(sim.n+1),]
  }
  if(family == "Binomial"){
    mu <- NULL; od <- NULL
    if(length(p) == 1) p = rep(p, nrow(locs))
    if(length(p)!= nrow(locs)) stop("length of p must be equal to number of locations")
    ans <- t(apply(sim.grf,1,function(x)
      qbinom(pnorm(x), size = binomial.N, prob = p)))[-(sim.n+1),]
  }
  if(!family %in% c("NB2", "Poisson", "ZIP", "Binomial"))  
     stop("'family' not recognized. Must be one of: 'NB2', 'Poisson', 'ZIP' or 'Binomial'.")
  return(ans)
}

#######################################################################################################
#### OUTPUT FUNCTION: GENERATE random samples (Truncated Correlated MVN)
#### Using GHK simulator. 
#######################################################################################################

rtnormGHK <- function(mu, lower, upper, sigma, n.sample = 1L){
  if(n.sample <= 0)
    stop("Input n.sample must be > 0!" )
  if(!is.matrix(sigma)) 
    stop("Input 'sigma' must be of form matrix!")
  if(!all(lower <= upper)) 
    stop("Elements in 'lower' must be <=  the corresponding elements in 'upper'!")
  if(!isSymmetric(sigma)) 
    stop("Input covariance matrix 'sigma' must be symmetric!")
  L <- try(t(chol(sigma)),silent=TRUE)
  if( inherits(L,"try-error") ) 
    stop("Cholesky Decomposition failed. Input matrix sigma is not a valid covariance matrix!")
  if(!all.equal(length(mu), nrow(sigma), length(lower), length(upper)))
    stop("Input 'mu', lower' and 'upper' must have same length as dimension of the sigma!")
  ans <- ghksimM(mu = mu, R = sigma, lower = lower, upper = upper, M = n.sample)
  if(any(apply(cbind(ans$simulation, lower), 2, function(x) (x-lower)*(x-upper)) > 0))
     stop("Algorithm Failed. The input vector 'mu' is too far from bounds. 
           Sample is not valid due to numeric issue!")
  return(ans)
}


