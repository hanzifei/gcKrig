library(EQL)

#### Calculate Frechet Hoeffding Upper Bound 
#### Must source FHUB.cpp first
###########################################################################################
####  Using Nelsen's summation rule: E(xy) = sum(S(xy)), S(xy) is joint survival function
###########################################################################################

FHUB.discrete.nelsen <- function(marg1, marg2, mu1, mu2, od1 = 0, od2 = 0,
                                 binomial.n1 = 1, binomial.n2 = 1)
{
  if (is.function(marg1) | is.function(marg2)){
    stop("Use method 'integral' or 'mc' !")
  }else{  
    if (!marg1 %in% c("Poisson", "ZIP", "NB", "binomial"))
      stop("'marg1' must be one of the following: 'Poisson', 'ZIP', 'NB' or 'binomial'. ")
    if (!marg2 %in% c("Poisson", "ZIP", "NB", "binomial"))
      stop("'marg2' must be one of the following: 'Poisson', 'ZIP', 'NB' or 'binomial'. ")
    
    if(marg1 == 'Poisson') od1 = 0; if(marg2 == 'Poisson') od2 = 0
    
    if( (marg1 == 'ZIP' & marg2 == 'ZIP') | (marg1 == 'Poisson' & marg2 == 'Poisson')
        | (marg1 == 'Poisson' & marg2 == 'ZIP') | (marg1 == 'ZIP' & marg2 == 'Poisson') )
      corr <- try(FHUBZIP(m1 = mu1, m2 = mu2, od1 = od1, od2 = od2), silent = T)
    
    if(marg1 == 'NB' & marg2 == 'NB')
      corr <- try(FHUBNB2(m1 = mu1, m2 = mu2, od1 = od1, od2 = od2), silent = T)
    
    if(marg1 == 'binomial' & marg2 == 'binomial')
      corr <- try(FHUBbinom(m1 = mu1, m2 = mu2, n1 = binomial.n1, n2 = binomial.n2), silent = T)
    
    
    if( (marg1 == 'NB' & marg2 == 'ZIP') | (marg1 == 'NB' & marg2 == 'Poisson') )
      corr <- try(FHUBZIPNB2(zipmu = mu2, nbmu = mu1, zipod = od2, nbod = od1), silent = T)
    if( (marg2 == 'NB' & marg1 == 'ZIP') | (marg2 == 'NB' & marg1 == 'Poisson') )
      corr <- try(FHUBZIPNB2(zipmu = mu1, nbmu = mu2, zipod = od1, nbod = od2), silent = T)
    
    
    if( (marg1 == 'binomial' & marg2 == 'ZIP') | (marg1 == 'binomial' & marg2 == 'Poisson') )
      corr <- try(FHUBZIPbinomial(zipmu = mu2, bmu = mu1, zipod = od2, bn = binomial.n1), silent = T)
    if( (marg2 == 'binomial' & marg1 == 'ZIP') | (marg2 == 'binomial' & marg1 == 'Poisson') )
      corr <- try(FHUBZIPbinomial(zipmu = mu1, bmu = mu2, zipod = od1, bn = binomial.n2), silent = T)
    
    
    if( (marg1 == 'NB' & marg2 == 'binomial'))
      corr <- try(FHUBNB2binomial(nbmu = mu1, bmu = mu2, nbod = od1, bn = binomial.n2), silent = T)
    if( (marg2 == 'NB' & marg1 == 'binomial'))
      corr <- try(FHUBNB2binomial(nbmu = mu2, bmu = mu1, nbod = od2, bn = binomial.n1), silent = T)
    
    if(is(corr,"try-error"))
      stop("Marginal mean too large or n too large in binomial: 
           computation of the summation is time consuming! Try use other method.")
    
    return(corr)
  }
}


###########################################################################################################
#### Specify a class of marginal functions with the potential of user-defined marginals
###########################################################################################################

#### Continuous Distributions...

gaussian.fhub <- function(mean = 0, sd = 1){
  q <- function(p) qnorm(p = p, mean = mean, sd = sd)
  ans <- list()
  ans$margvar <- sd^2
  
  ans$int.marg <- function (order) {
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
      order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- rnorm(nrep, mean = mean, sd = sd)
  class(ans) <- c("marginal.fhub")
  return(ans)
}


gamma.fhub <- function(shape = 1, rate = 1){
  
  q <- function(p) qgamma(p = p, shape = shape, rate = rate)
  ans <- list()
  ans$margvar <- shape/(rate^2)
  ans$int.marg <- function (order) {
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
              order = order, -8, 8, 
      subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- rgamma(nrep, shape = shape, rate = rate)
  class(ans) <- c("marginal.fhub")
  return(ans)
}



beta.fhub <- function(shape1 = 1, shape2 = 1){
  q <- function(p) qbeta(p = p, shape1 = shape1, shape2 = shape2)
  ans <- list()
  ans$margvar <- shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1))
  
  ans$int.marg <- function (order) {
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
      order = order, -8, 8, 
      subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- rbeta(nrep, shape1 = shape1, shape2 = shape2)
  class(ans) <- c("marginal.fhub")
  return(ans)
}

weibull.fhub <- function(shape = 1, scale = 1){
  q <- function(p) qweibull(p = p, shape = shape, scale = scale)
  ans <- list()
  ans$margvar <- scale^2*(base::gamma(1+2/shape) - base::gamma(1+1/shape)^2)
  ans$int.marg <- function(order){ 
    integrate(function(x, order) 
             ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
            q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
    order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- rweibull(nrep, shape = shape, scale = scale)
  class(ans) <- c("marginal.fhub")
  return(ans)
}

  
#### Some Discrete Distributions...

poisson.fhub <- function(lambda = 1){
  q <- function(p) qpois(p = p, lambda = lambda)
  ans <- list()
  ans$margvar <- lambda
  ans$int.marg <- function(order){ 
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
      order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- rpois(nrep, lambda = lambda)
  class(ans) <- c("marginal.fhub")
  return(ans)
}

  

negbin.fhub <- function(mu = 1, od = 1){
  q <- function(p) qnbinom(p = p, size = 1/od, mu = mu)
  ans <- list()
  ans$margvar <- mu*(1+mu*od)
  ans$int.marg <- function(order){ 
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
      order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- rnbinom(nrep, size = 1/od, mu = mu)
  class(ans) <- c("marginal.fhub")
  return(ans)
}
  
  
 
binomial.fhub <- function(size = 1, prob = 0.5){
  q <- function(p) qbinom(p = p, size = size, prob = prob)
  ans <- list()
  ans$margvar <- size*prob*(1-prob)
  ans$int.marg <- function(order){ 
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
      order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- rbinom(nrep, size = size, prob = prob)
  class(ans) <- c("marginal.fhub")
  return(ans)
} 
  
  
zip.fhub <- function(mu = 1, od = 1){
  q <- function(p) qpois(pmax( 0, (p-od/(1+od))/(1-od/(1+od)) ), (1+od)*mu)
  ans <- list()
  ans$margvar <- mu*(1+mu*od)
  ans$int.marg <- function(order){ 
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
      order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  
  ans$rsample <- ifelse(rbinom(nrep, 1, od/(1+od)), 0, rpois(nrep, (1+od)*mu))
  class(ans) <- c("marginal.fhub")
  return(ans)
}


###########################################################################################
####  The function that return the FHUB, with: 
####  Three methods
####  Four discrete marginals: Poisson, nb, binomial and zero-inflated poisson
####  Flexible function that can be implemented for arbitrary marginals.
###########################################################################################

FHUB <- function(marg1, marg2, nrep = 1000, method = "integral"){
 
 if(!inherits(marg1, "marginal.fhub"))  stop("'marg1' must be a function of the class marg.fhub")
 if(!inherits(marg2, "marginal.fhub")) stop("'marg2' must be a function of the class marg.fhub")
  
  if(method == "integral"){
    assign("nrep", 1, envir=.GlobalEnv); atmp <- c()
    for(u in 1:50){
      atmp[u] = marg1$int.marg(order = u)[[1]]*marg2$int.marg(order = u)[[1]]/factorial(u)
      if (abs(atmp[u]) < .Machine$double.eps^0.25*sqrt(marg1$margvar*marg2$margvar) )
        break
    }
    corr <- sum(atmp)/sqrt(marg1$margvar*marg2$margvar)
  }
  if(method == "mc"){
    assign("nrep", nrep, envir=.GlobalEnv)
    corr <- cor(sort(marg1$rsample), sort(marg2$rsample))
  }
  corr <- ifelse(corr>1, 1, corr)
  return(corr) 
}

#### Some Examples 
FHUB(marg1 = gamma.fhub(shape = 10, rate = 1), 
     marg2 = beta.fhub(shape1 = 20, shape2 = 3), method = "integral")

FHUB(marg1 = gamma.fhub(shape = 10, rate = 1), 
     marg2 = beta.fhub(shape1 = 20, shape2 = 3), nrep = 100000,method = "mc")




FHUB(marg1 = binomial.fhub(size = 10, prob = 0.5), marg2 = zip.fhub(mu = 4, od = 0.3), 
     method = "integral")

FHUB(marg1 = binomial.fhub(size = 10, prob = 0.5), marg2 = zip.fhub(mu = 4, od = 0.3), 
     nrep = 100000, method = "mc")

FHUB.discrete.nelsen(marg1 = 'binomial', marg2 = 'ZIP', 
                     mu1 = 5, mu2 = 4, od1 = 0, od2 = 0.3, binomial.n1 = 10)
