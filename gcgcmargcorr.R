
###########################################################################################################
#### Specify a class of marginal functions with the potential of user-defined marginals
#### gcgc1: parameter values are given by users, for simulation/computing FHUB purposes
###########################################################################################################
#### Continuous Distributions...

gaussian.gcgc1 <- function(mean = 0, sd = 1){
  q <- function(p) qnorm(p = p, mean = mean, sd = sd)
  ans <- list()
  ans$margvar <- sd^2
  ans$int.marg <- function (order) {
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
      order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- function(nrep) rnorm(n = nrep, mean = mean, sd = sd)
  ans$q <- q
  class(ans) <- c("marginal.gcgc1")
  return(ans)
}


gamma.gcgc1 <- function(shape = 1, rate = 1){
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
  ans$rsample <- function(nrep) rgamma(n = nrep, shape = shape, rate = rate)
  ans$q <- q
  class(ans) <- c("marginal.gcgc1")
  return(ans)
}

beta.gcgc1 <- function(shape1 = 1, shape2 = 1){
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
  ans$rsample <- function(nrep) rbeta(n = nrep, shape1 = shape1, shape2 = shape2)
  ans$q <- q
  class(ans) <- c("marginal.gcgc1")
  return(ans)
}

weibull.gcgc1 <- function(shape = 1, scale = 1){
  q <- function(p) qweibull(p = p, shape = shape, scale = scale)
  ans <- list()
  ans$margvar <- scale^2*(base::gamma(1+2/shape) - base::gamma(1+1/shape)^2)
  ans$int.marg <- function(order){ 
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
      order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- function(nrep) rweibull(n = nrep, shape = shape, scale = scale)
  ans$q <- q
  class(ans) <- c("marginal.gcgc1")
  return(ans)
}

#### Some Discrete Distributions...

poisson.gcgc1 <- function(lambda = 1){
  q <- function(p) qpois(p = p, lambda = lambda)
  ans <- list()
  ans$margvar <- lambda
  ans$int.marg <- function(order){ 
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
      order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- function(nrep) rpois(n = nrep, lambda = lambda)
  ans$q <- q
  class(ans) <- c("marginal.gcgc1")
  return(ans)
}



negbin.gcgc1 <- function(mu = 1, od = 1){
  q <- function(p) qnbinom(p = p, size = 1/od, mu = mu)
  ans <- list()
  ans$margvar <- mu*(1+mu*od)
  ans$int.marg <- function(order){ 
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
      order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- function(nrep) rnbinom(n = nrep, size = 1/od, mu = mu)
  ans$q <- q
  class(ans) <- c("marginal.gcgc1")
  return(ans)
}



binomial.gcgc1 <- function(size = 1, prob = 0.5){
  q <- function(p) qbinom(p = p, size = size, prob = prob)
  ans <- list()
  ans$margvar <- size*prob*(1-prob)
  ans$int.marg <- function(order){ 
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
      order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- function(nrep) rbinom(n = nrep, size = size, prob = prob)
  ans$q <- q
  class(ans) <- c("marginal.gcgc1")
  return(ans)
} 


zip.gcgc1 <- function(mu = 1, od = 1){
  q <- function(p) qpois(pmax( 0, (p-od/(1+od))/(1-od/(1+od)) ), (1+od)*mu)
  ans <- list()
  ans$margvar <- mu*(1+mu*od)
  ans$int.marg <- function(order){ 
    integrate(function(x, order) 
      ifelse( (q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0, 
              q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = T)), 
      order = order, -8, 8, subdivisions = 1500, rel.tol = 0.01, stop.on.error = FALSE)
  }
  ans$rsample <- function(nrep) ifelse(rbinom(n=nrep, 1, od/(1+od)), 0, rpois(n=nrep, (1+od)*mu))
  ans$q <- q
  class(ans) <- c("marginal.gcgc1") 
  return(ans)
}


####### Correlaiton Functions in the class of gcgc1 for simulation purpose


matern.gcgc1 <- function(range = 0, kappa = 0.5, nugget = 0){
  ans <- list()
  ans$S <- function(D) (1-nugget)*geoR::matern(D, range, kappa) + nugget*diag(NROW(D))
  class(ans) <- c("corr.gcgc1")
  return(ans)
}


powerexp.gcgc1 <- function(range = 0, order = 1, nugget = 0){
  ans <- list()
  ans$S <- function(D) (1-nugget)*exp(-abs( (D/range)^(order) )) + nugget*diag(NROW(D))
  class(ans) <- c("corr.gcgc1")
  return(ans)
}


spherical.gcgc1 <- function(range = 0, nugget = 0){
  ans <- list()
  ans$S <- function(D) (1-nugget)*(1 - 1.5*D/range + 0.5*(D/range)^3)*(D<=range)+nugget*diag(NROW(D))
  class(ans) <- c("corr.gcgc1")
  return(ans)
}



###########################################################################################################
#### Specify a class of marginal functions with the potential of user-defined marginals
#### gcgc2: parameter values are estimated given data, for computing MLE/Prediction
###########################################################################################################
#### Output Objects: 
# ans$npar.cor :  Number of correlation parameters
# ans$start: Starting points
# ans$chol: Cholesky decomposition of correlation matrix
# ans$matD: Distance Matrix

matern.gcgc2 <- function(kappa = 0.5, nugget = T){
    ans <- list()
  if(nugget == TRUE){
    ans$nug <- 1
    ans$npar.cor <- 2
    ans$start <- function(matD) c(median(matD), 0.2)
    ans$corr <- function(corrpar, matD, not.na ){
      S <- (1 - corrpar[2]) * geoR::matern(matD, corrpar[1], kappa) + corrpar[2]*diag(NROW(matD));S
      }
    }else{
    ans$nug <- 0
    ans$npar.cor <- 1
    ans$start <- function(matD) median(matD)
    ans$corr <- function(corrpar, matD, not.na ){
      S <- geoR::matern(matD, corrpar, kappa);S 
    }
  }
  return(ans)
}


powerexp.gcgc2 <- function(order = 1, nugget = T){
  ans <- list()
  if(nugget == TRUE){
    ans$nug <- 1
    ans$npar.cor <- 2
    ans$start <- function(matD) c(median(matD), 0.2)
    ans$corr <- function(corrpar, matD, not.na ){
      S <- (1 - corrpar[2])*exp(-abs( (matD/corrpar[1])^(order) )) + corrpar[2]*diag(NROW(matD));S
    }
  }else{
    ans$nug <- 0
    ans$npar.cor <- 1
    ans$start <- function(matD) median(matD)
    ans$corr <- function(corrpar, matD, not.na ){
      S <- exp(-abs((matD/corrpar[1])^(order))); S
    }
  }
  return(ans)
}



spherical.gcgc2 <- function(nugget = T){
  ans <- list()
  if(nugget == TRUE){
    ans$nug <- 1
    ans$npar.cor <- 2
    ans$start <- function(matD) c(median(matD), 0.2)
    ans$corr <- function(corrpar, matD, not.na ){
      S <- (1 - corrpar[2])*(1 - 1.5*matD/corrpar[1] + 0.5*(matD/corrpar[1])^3)*(matD<=corrpar[1]) + 
            corrpar[2]*diag(NROW(matD)); S
    }
  }else{
    ans$nug <- 0
    ans$npar.cor <- 1
    ans$start <- function(matD) median(matD)
    ans$corr <- function(corrpar, matD, not.na ){
      S <- (1 - 1.5*matD/corrpar[1] + 0.5*(matD/corrpar[1])^3)*(matD<=corrpar[1]); S
    }
  }
  return(ans)
}


#### Now Try to Create a class of Marginals starting from negative binomial

poisson.gcgc2 <- function(link = "log"){
  fm <- poisson( substitute( link ) ); ans <- list()
  ans$start <- function(y, x, effort) {
     mfit <- suppressWarnings(glm.fit(x, y/effort, family = fm))
     mu0 <- fitted(mfit)
     reg0 <- coef(mfit)
    return(list(start = c(reg0), res = qnorm(ppois(q = y, lambda = mu0)) ))
  }
  ans$nod <- 0
  ans$bounds <- function(y, x, pars, effort) {
      M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
      a <- qnorm(ppois( y-1, lambda = M));  b <- qnorm(ppois( y, lambda = M))
    return(list(lower = a, upper = b))
  }
  ans$pdf <- function(y, x, pars, effort){
      M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
      pdf <- dpois(y, lambda = M, log = T)
    return(pdf)
  }
  ans$cdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    cdf <- ppois( y, lambda = M)
    return(cdf)
  }
  ans
}


negbin.gcgc2 <- function(link = "log"){
  fm <- poisson( substitute( link ) ); ans <- list()
  
  ans$start <- function(y, x, effort) {
    mfit <- suppressWarnings(glm.fit(x, y/effort, family = fm))
    reg0 <- coef(mfit); mu <- fitted(mfit)
    od0 <- max(10*.Machine$double.eps, mean(((y-mu)^2-mu)/mu^2))
     return(list(start = c(reg0, od0), res = qnorm(pnbinom(y, size = 1/od0, mu = mu))))
  }
  ans$nod <- 1
  ans$bounds <- function(y, x, pars, effort) {
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    S <- 1/pars[ncol(x)+1]
    a <- qnorm(pnbinom( y-1, size = S, mu = M));  b <- qnorm(pnbinom( y, size = S, mu = M))
    return(list(lower = a, upper = b))
  }
  ans$pdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    S <- 1/pars[ncol(x)+1]
    pdf <- dnbinom(y, size = S, mu = M, log = T)
    return(pdf)
  }
  ans$cdf <- function(y, x, pars, effort){
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort
    S <- 1/pars[ncol(x)+1]
    cdf <- pnbinom(y, size = S, mu = M)
    return(cdf)
  }
 ans
}

  
binomial.gcgc2 <- function(link = "logit"){  
  
   fm <- binomial( substitute( link ) )
   ans <- list()
   ans$start <- function(y, x, effort) {
     mfit <- suppressWarnings(glm.fit(x, y/(effort), family = fm))
     reg0 <- coef(mfit)
     return(list(start = reg0, res = qnorm(pbinom(y, size = effort, prob = fitted(mfit)))))
  }  
   ans$nod <- 0
   ans$bounds <- function(y, x, pars, effort) {
     p <- fm$linkinv(pars[1:ncol(x)]%*%t(x))
     a <- qnorm(pbinom( y-1, size = effort, prob = p))
     b <- qnorm(pbinom( y, size = effort, prob = p))
     return(list(lower = a, upper = b))
  }
   ans$pdf <- function(y, x, pars, effort){
     p <- fm$linkinv(pars[1:ncol(x)]%*%t(x))
     pdf <- dbinom(y, size = effort, prob = p, log = T)
     return(pdf)
   }
   ans$cdf <- function(y, x, pars, effort){
     p <- fm$linkinv(pars[1:ncol(x)]%*%t(x))
     cdf <- pbinom(y, size = effort, prob = p)
     return(cdf)
   }
  ans
}


zip.gcgc2 <- function(link = "log"){
  fm <- poisson( substitute( link ) ); ans <- list()
  ans$start <- function(y, x, effort) {
    mfit <- suppressWarnings(glm.fit(x, y/effort, family = fm))
    reg0 <- coef(mfit); mu <- fitted(mfit)
    od0 <- max(10*.Machine$double.eps, mean(((y-mu)^2-mu)/mu^2))
    return(list(start = c(reg0, od0), 
                res = qnorm((y >= 0)*od0/(1+od0)+ppois(y, (1+od0)*mu)/(1+od0))))
  }
  ans$nod <- 1
  ans$bounds <- function(y, x, pars, effort) {
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort;  od <- pars[ncol(x)+1]
    a <- qnorm( (y>=1)*od/(1+od) + ppois(y-1, lambda = (1+od)*M)/(1+od) )
    b <- qnorm( (y>=0)*od/(1+od) + ppois(y, lambda = (1+od)*M)/(1+od) )
    return(list(lower = a, upper = b))
  }
  ans$pdf <- function(y, x, pars, effort){ 
   M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort;  od <- pars[ncol(x)+1]
   pdf <- log( (y==0)*od/(1+od) + dpois(y, (1+od)*M)/(1+od) )
   return(pdf)
  }
  ans$cdf <- function(y, x, pars, effort){ 
    M <- fm$linkinv(pars[1:ncol(x)]%*%t(x))*effort;  od <- pars[ncol(x)+1]
    cdf <- (y >= 0)*od/(1+od)+ppois(y, (1+od)*M)/(1+od)
    return(cdf)
  }
  ans 
}






# A special Example of specifying your own link function

binomial_t.gcgc2 <- function(df.t = 6){ 
   ans <- list()
   ans$start <- function(y, x, effort) {
     mfit <- suppressWarnings(glm.fit(x, y/effort, family = binomial(link = "logit")))
     reg0 <- coef(mfit)
     return(list(start = reg0, res = qnorm(pbinom(y, size = effort, prob = fitted(mfit)))))
  }
  ans$nod <- 0
  ans$bounds <- function(y, x, pars, effort) {
    p <- pt(pars[1:ncol(x)]%*%t(x), df = df.t)
    a <- qnorm(pbinom( y-1, size = effort, prob = p))
    b <- qnorm(pbinom( y, size = effort, prob = p))
    return(list(lower = a, upper = b))
  }
  ans$pdf <- function(y, x, pars, effort){
    p <- pt(pars[1:ncol(x)]%*%t(x), df = df.t)
    pdf <- dbinom(y, size = effort, prob = p, log = T)
    return(pdf)
  } 
  ans$cdf <- function(y, x, pars, effort){
    p <- pt(pars[1:ncol(x)]%*%t(x), df = df.t)
    cdf <- pbinom(y, size = effort, prob = p, log = T)
    return(cdf)
  }
  ans
}