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



###########################################################################################
####  Integrate Directly
###########################################################################################

library(EQL)

intnb <- function(mu, od, order){
  integrate(function(x, mu, od, order) 
              ifelse((qnbinom(pnorm(x), size=1/od, mu=mu)==Inf | dnorm(x) < .Machine$double.eps),0,
             qnbinom(pnorm(x), size=1/od, mu=mu)*dnorm(x)*EQL::hermite(x, order, prob = T))
             , mu = mu, od = od, order = order, -8, 8, 
            subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
}


intzip <- function (mu, od, order) {
  integrate(function(x, mu, od, order) 
    ifelse( (qpois(pmax(0, (pnorm(x)-od/(1+od))/(1-od/(1+od))), (1+od)*mu) ==Inf | 
               dnorm(x) < .Machine$double.eps), 0, 
            qpois(pmax(0, (pnorm(x)-od/(1+od))/(1-od/(1+od))), 
                  (1+od)*mu)*dnorm(x)*EQL::hermite(x, order, prob = T)),
    mu = mu, od = od, order = order, -8, 8, 
    subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  
}


intbinomial <- function (mu, n.binomial, order) {
  integrate(function(x, mu, n.binomial, order) 
    ifelse((qbinom(pnorm(x), size = n.binomial, prob = mu/n.binomial)==Inf | 
              dnorm(x) < .Machine$double.eps), 0, qbinom(pnorm(x), size= n.binomial, 
              prob = mu/n.binomial)*dnorm(x)*EQL::hermite(x, order, prob = T)), 
            mu = mu, n.binomial = n.binomial, order = order, -8, 8, 
            subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
}



self.marg.fhub1 <- function(){
  q <- function(p) qgamma(p = p, shape = 4, rate = 1)
  ans <- list()
  ans$margvar <- 4
  ans$int.marg <- function (order) {
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), order = order, -8, 8, 
      subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
# if(method == 'mc'){
 ans$rsample <- function(rep) rgamma(rep, shape = 4, rate = 1)
# }
  class(ans) <- c("marg.fhub")
  return(ans)
}
  

self.marg.fhub2 <- function(){
  q <- function(p) qgamma(p = p, shape = 3, rate = 1)
  ans <- list()
  ans$margvar <- 3
  ans$int.marg <- function (order) {
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), order = order, -8, 8, 
      subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
#  if(method == 'mc'){
    ans$rsample <- function(rep) rgamma(rep, shape = 3, rate = 1)
 #  
  class(ans) <- c("marg.fhub")
  return(ans)
}


FHUB.discrete.integral <- function(marg1, marg2, mu1, mu2, od1 = 0, od2 = 0, 
                                   binomial.n1 = 1, binomial.n2 = 1)
{
  atmp <- c()
  
  if (is.function(marg1) & is.function(marg2)){
    funtmp1 <- marg1()
    funtmp2 <- marg2()
    
  if(!inherits(funtmp1, "marg.fhub"))
    stop("'marg1' must be one of the following: 'Poisson', 'ZIP', 'NB', 'binomial'
         or self-defined function of the class marg.fhub")
    
  if(!inherits(funtmp2, "marg.fhub"))
    stop("'marg2' must be one of the following: 'Poisson', 'ZIP', 'NB', 'binomial'
         or self-defined function of the class marg.fhub")
  
  for(u in 1:30){
    atmp[u] = funtmp1$int.marg(order = u)[[1]]*funtmp2$int.marg(order = u)[[1]]/factorial(u)
    if (abs(atmp[u]) < .Machine$double.eps^0.25*sqrt(funtmp1$margvar*funtmp2$margvar) )
      break
  }
    corr <- sum(atmp)/sqrt(funtmp1$margvar*funtmp2$margvar)
  }else{
    if (!marg1 %in% c("Poisson", "ZIP", "NB", "binomial"))
      stop("'marg1' must be one of the following: 'Poisson', 'ZIP', 'NB', 'binomial'
         or self-defined function of the class marg.fhub")
    
    if (!marg2 %in% c("Poisson", "ZIP", "NB", "binomial"))
      stop("'marg2' must be one of the following: 'Poisson', 'ZIP', 'NB', 'binomial'
           or self-defined function of the class marg.fhub")
    
    if(marg1 == 'Poisson') od1 = 0; if(marg2 == 'Poisson') od2 = 0
    
  if( (marg1 == 'ZIP' & marg2 == 'ZIP') | (marg1 == 'Poisson' & marg2 == 'Poisson')
      | (marg1 == 'Poisson' & marg2 == 'ZIP') | (marg1 == 'ZIP' & marg2 == 'Poisson') )
  { 
    for(u in 1:30){
    atmp[u] <- intzip(mu = mu1, od = od1, order = u)[[1]]*intzip(mu = mu2, od = od2, order = u)[[1]]/factorial(u)
    if (abs(atmp[u]) < .Machine$double.eps^0.25*sqrt(mu1*mu2*(1+od1*mu1)*(1+od2*mu2)) )
      break}
    corr <- sum(atmp)/sqrt(mu1*mu2*(1+od1*mu1)*(1+od2*mu2))
  }
  
  if(marg1 == 'NB' & marg2 == 'NB')
  {
    for(u in 1:30){
      atmp[u] <- intnb(mu = mu1, od = od1, order = u)[[1]]*intnb(mu = mu2, od = od2, order = u)[[1]]/factorial(u)
      if (abs(atmp[u]) < .Machine$double.eps^0.25*sqrt(mu1*mu2*(1+od1*mu1)*(1+od2*mu2)))
        break}
    corr <- sum(atmp)/sqrt(mu1*mu2*(1+od1*mu1)*(1+od2*mu2))
  }
  
  if(marg1 == 'binomial' & marg2 == 'binomial')
  {
    for(u in 1:30){
      atmp[u] <- intbinomial(mu = mu1, n.binomial = binomial.n1, order = u)[[1]]*
                 intbinomial(mu = mu2, n.binomial = binomial.n2, order = u)[[1]]/factorial(u)
      if (abs(atmp[u]) < .Machine$double.eps^0.25*sqrt(mu1*mu2*(1-mu1/binomial.n1)*(1-mu2/binomial.n2)))
        break}
    corr <- sum(atmp)/sqrt(mu1*mu2*(1-mu1/binomial.n1)*(1-mu2/binomial.n2))
  }
  
  if( (marg1 == 'NB' & marg2 == 'ZIP') | (marg1 == 'NB' & marg2 == 'Poisson') )
  {
    for(u in 1:30){
      atmp[u] <- intnb(mu = mu1, od = od1, order = u)[[1]]*
                 intzip(mu = mu2, od = od2, order = u)[[1]]/factorial(u)
      if (abs(atmp[u]) < .Machine$double.eps^0.25*sqrt(mu1*mu2*(1+od1*mu1)*(1+od2*mu2)))
        break}
    corr <- sum(atmp)/sqrt(mu1*mu2*(1+od1*mu1)*(1+od2*mu2))
  }  
    
  if( (marg2 == 'NB' & marg1 == 'ZIP') | (marg2 == 'NB' & marg1 == 'Poisson') )
  {
    for(u in 1:30){
      atmp[u] <- intzip(mu = mu1, od = od1, order = u)[[1]]*
                 intnb(mu = mu2, od = od2, order = u)[[1]]/factorial(u)
      if (abs(atmp[u]) < .Machine$double.eps^0.25*sqrt(mu1*mu2*(1+od1*mu1)*(1+od2*mu2)))
        break}
    corr <- sum(atmp)/sqrt(mu1*mu2*(1+od1*mu1)*(1+od2*mu2))
  }
  
  if( (marg1 == 'binomial' & marg2 == 'ZIP') | (marg1 == 'binomial' & marg2 == 'Poisson') )
  {
    for(u in 1:30){
      atmp[u] <- intbinomial(mu = mu1, n.binomial = binomial.n1, order = u)[[1]]*
                 intzip(mu = mu2, od = od2, order = u)[[1]]/factorial(u)
      if (abs(atmp[u]) < .Machine$double.eps^0.25*sqrt(mu1*mu2*(1-mu1/binomial.n1)*(1+od2*mu2)))
        break}
    corr <- sum(atmp)/sqrt(mu1*mu2*(1-mu1/binomial.n1)*(1+od2*mu2))
  }
  
  if( (marg2 == 'binomial' & marg1 == 'ZIP') | (marg2 == 'binomial' & marg1 == 'Poisson') )
  {
    for(u in 1:30){
      atmp[u] <- intbinomial(mu = mu2, n.binomial = binomial.n2, order = u)[[1]]*
                 intzip(mu = mu1, od = od1, order = u)[[1]]/factorial(u)
      if (abs(atmp[u]) < .Machine$double.eps^0.25*sqrt(mu1*mu2*(1-mu2/binomial.n2)*(1+od1*mu1)))
        break}
    corr <- sum(atmp)/sqrt(mu1*mu2*(1-mu2/binomial.n2)*(1+od1*mu1))
  }
  
  if( (marg1 == 'NB' & marg2 == 'binomial'))
  {
    for(u in 1:30){
      atmp[u] <- intbinomial(mu = mu2, n.binomial = binomial.n2, order = u)[[1]]*
                 intnb(mu = mu1, od = od1, order = u)[[1]]/factorial(u)
      if (abs(atmp[u]) < .Machine$double.eps^0.25*sqrt(mu1*mu2*(1-mu2/binomial.n2)*(1+od1*mu1)))
        break}
    corr <- sum(atmp)/sqrt(mu1*mu2*(1-mu2/binomial.n2)*(1+od1*mu1))
  }
    
  if( (marg2 == 'NB' & marg1 == 'binomial'))
  {
    for(u in 1:30){
      atmp[u] <- intbinomial(mu = mu1, n.binomial = binomial.n1, order = u)[[1]]*
                 intnb(mu = mu2, od = od2, order = u)[[1]]/factorial(u)
      if (abs(atmp[u]) < .Machine$double.eps^0.25*sqrt(mu1*mu2*(1-mu1/binomial.n1)*(1+od2*mu2)))
        break}
    corr <- sum(atmp)/sqrt(mu1*mu2*(1-mu1/binomial.n1)*(1+od2*mu2))
  }
  
}
   corr <- ifelse(corr>1, 1, corr)
    return(corr)
}



###########################################################################################
####  Using Monte Carlo
###########################################################################################

FHUB.discrete.mc<- function(marg1, marg2, mu1, mu2, od1 = 0, od2 = 0, 
                            binomial.n1 = 1, binomial.n2 = 1, nrep = 100000)
{
  if (is.function(marg1) & is.function(marg2)){
    funtmp1 <- marg1()
    funtmp2 <- marg2()
    
    if(!inherits(funtmp1, "marg.fhub"))
      stop("'marg1' must be one of the following: 'Poisson', 'ZIP', 'NB', 'binomial'
           or self-defined function of the class marg.fhub")
    
    if(!inherits(funtmp2, "marg.fhub"))
      stop("'marg2' must be one of the following: 'Poisson', 'ZIP', 'NB', 'binomial'
           or self-defined function of the class marg.fhub")
    
    sim1 <- funtmp1$rsample(rep = nrep); sim2 <- funtmp2$rsample(rep = nrep)
    
  }else{
    
  if (!marg1 %in% c("Poisson", "ZIP", "NB", "binomial"))
    stop("'marg1' must be one of the following: 'Poisson', 'ZIP', 'NB' or 'binomial'. ")
  if (!marg2 %in% c("Poisson", "ZIP", "NB", "binomial"))
    stop("'marg2' must be one of the following: 'Poisson', 'ZIP', 'NB' or 'binomial'. ")
    
  if(marg1 == 'Poisson') od1 = 0; if(marg2 == 'Poisson') od2 = 0
  
  if( marg1 == 'ZIP' | marg1 == 'Poisson')
    sim1 <- ifelse(rbinom(nrep, 1, od1/(1+od1)), 0, rpois(nrep, (1+od1)*mu1))
  if( marg1 == 'NB')
    sim1 <- rnbinom(n = nrep, size = 1/od1, mu = mu1)
  if( marg1 == 'binomial')
    sim1 <- rbinom(n = nrep, size = binomial.n1, prob = mu1/binomial.n1)
  
  if( marg2 == 'ZIP' | marg2 == 'Poisson')
    sim2 <- ifelse(rbinom(nrep, 1, od2/(1+od2)), 0, rpois(nrep, (1+od2)*mu2))
  if( marg2 == 'NB')
    sim2 <- rnbinom(n = nrep, size = 1/od2, mu = mu2)
  if( marg2 == 'binomial')
    sim2 <- rbinom(n = nrep, size = binomial.n2, prob = mu2/binomial.n2)
  }
    corr <- cor(sort(sim1), sort(sim2))  
  return(corr)
}


###########################################################################################
####  The function that return the FHUB, with: 
####  Three methods
####  Four discrete marginals: Poisson, nb, binomial and zero-inflated poisson
####  Flexible function that can be implemented for arbitrary marginals.
###########################################################################################

FHUB <- function(marg1, marg2, mu1, mu2, od1 = 0, od2 = 0, binomial.n1 = 1, binomial.n2 = 1, 
                 nrep = 10000, method = "nelsen")
{
  if(method == "nelsen"){
    ans <- FHUB.discrete.nelsen(marg1 = marg1, marg2 = marg2, mu1 = mu1, mu2 = mu2, od1 = od1, 
                                od2 = od2, binomial.n1 = binomial.n1, binomial.n2 = binomial.n2)
  }
  if(method == "integral"){
    ans <- FHUB.discrete.integral(marg1 = marg1, marg2 = marg2, mu1 = mu1, mu2 = mu2, od1 = od1, 
                                  od2 = od2, binomial.n1 = binomial.n1, binomial.n2 = binomial.n2)
  }
  if(method == "mc"){
    ans <- FHUB.discrete.mc(marg1 = marg1, marg2 = marg2, mu1 = mu1, mu2 = mu2, od1 = od1, od2 = od2,
                            binomial.n1 = binomial.n1, binomial.n2 = binomial.n2, nrep = nrep)
  }
  return(ans) 
}

unix.time(
  FHUB(marg1 = 'NB', marg2 = 'NB', 
       mu1 = 20, mu2 = 1, od1 = 1, od2 = 1, binomial.n1 = 1, binomial.n2 = 10, 
       nrep = 10000, method = "nelsen")
)

unix.time(
FHUB(marg1 = self.marg.fhub1, marg2 = self.marg.fhub2, 
    mu1, mu2, od1 = 0, od2 = 0, binomial.n1 = 1, binomial.n2 = 1, nrep = 100000, method = "mc")
)

FHUB(marg1 = 'NB', marg2 = 'binomial',mu1 = 20, mu2 = 1, od1 = 1, od2 = 1, 
     binomial.n1 = 1, binomial.n2 = 10, nrep = 10000, method = "nelsen")
