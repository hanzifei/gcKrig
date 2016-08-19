library(Rcpp)
library(RcppArmadillo)
library(geoR)
library(VGAM)
library(MASS)
library(snowfall)

#### Create function likGQT which returns the negative (approximate) log-likelihood using GQT
#### C++ is implemented with editor GQTLik.cpp
#######################################################################################################################
#######################################################################################################################
#### Inner function only used to be called by another functions inside package
#######################################################################################################################
#######################################################################################################################
likGQT <- function(pars, y, x = NULL, locs, family, corr, effort){
  matD <- as.matrix(dist(locs, method = "euclidean", diag = T, upper = T))
  npars <- length(pars)
  nparcorr <- corr$npar.cor
  R <- corr$corr(pars[(npars-nparcorr+1):npars], matD = matD)
  logindpdf <- sum(family$pdf(y = y, x = x, pars = pars[1:(npars-nparcorr)], effort = effort))
  u1 <- family$cdf(y = y-1, x = x, pars = pars, effort = effort)
  u2 <- family$cdf(y = y, x = x, pars = pars, effort = effort)
  u <- (u1+u2)*0.5
  loglik <- ldgc(u = u, sigma = R)+logindpdf
  if(is.nan(loglik)) loglik = -1e6 
  if (loglik == Inf) loglik = 1e6
  return(-loglik)
}

likGQT(pars = c(0.7,1.1,0.2,0.5,0.2), y = simtmp.y, x = cbind(1,xloc,yloc), locs = cbind(xloc,yloc),
       family = poisson.gcgc2(link = 'log'), 
       corr = matern.gcgc2(matern.kappa = 0.5, nugget = T), effort = 1)

#mleGQT( y = sim.y, obs.locs = cbind(xloc,yloc), family = 'NB')

mleGQT <- function(y, x = NULL, obs.locs, family, corr, effort = 1, max.range = NULL)
{
  #### Create initial values in optimization 
  matD <- as.matrix(dist(obs.locs, method = "euclidean", diag = T, upper = T))
  if(is.null(max.range)) max.range <- 2*max(matD)
  marg0 <- family$start(y = y, x = x, effort = effort)$start
  n.nugget0 <- corr$nug; MIN <- .Machine$double.eps^0.25
  z0 <- family$start(y = y, x = x, effort = effort)$res
  
  gauslik <- function(corrpar) -mvtnorm::dmvnorm(z0, mean = rep(0, nrow(matD)), 
                                                 sigma = corr$corr(corrpar, matD), log = TRUE)
  
  corr0 <- try(optim(par =  corr$start(matD) , fn = gauslik, method = "L-BFGS-B", 
                     lower = c(MIN, rep(0,n.nugget0)), upper = c(max.range, rep(1-MIN, n.nugget0))), silent = T)
  
  if(inherits(corr0, "try-error")){
    corpar0 <- corr$start(matD)}
  else{
    corpar0 <- corr0$par;
  }
  
  est <- c(family$start(y = y, x = x, effort = effort)$start, corpar0)
  n.reg0 <- ncol(x); n.od0 <- family$nod 
  
  fit <- optim(par = est, fn = likGQT, y = y, x = x, locs = obs.locs, family = family, corr = corr, 
                 effort = effort, method = "L-BFGS-B", 
                 lower = c(rep(-Inf,n.reg0), rep(MIN, n.od0), MIN, rep(0, n.nugget0)),
                 upper = c(rep(Inf,n.reg0), rep(Inf, n.od0), max.range, rep(1-MIN, n.nugget0)))
  est <- fit$par
  
  if (fit$par[n.reg0+n.od0+1] == max.range)
    warning("Maximum Range Value Reached, MLE may not exist. Try to input a larger One")
  
  if (is.null(fit$convergence) || fit$convergence != 0)   
    stop("Maximum likelihood estimation failed. Algorithm does not converge") 
  
  result.list <- list()
  result.list$MLE <- est; result.list$log.lik <- -fit$value
  k <- length(est); N <- length(y)
  result.list$AIC <- 2*k+2*fit$value
  result.list$AICc <- 2*k+2*fit$value + 2*k*(k+1)/(N-k-1)
  result.list$BIC <- 2*fit$value + k*log(N)
 # class(result.list) <- c( "pred.ghk")
  return(result.list)
}

mleGQT(y = simtmp.y, x = cbind(1, xloc,yloc), obs.locs = cbind(xloc,yloc), 
       family = negbin.gcgc2(link = 'log'), corr = matern.gcgc2( matern.kappa = 0.5, nugget = F), 
       effort = 1, max.range = NULL)

mleGHK(y = simtmp.y, x = cbind(1, xloc,yloc), obs.locs = cbind(xloc,yloc), 
       family = negbin.gcgc2(link = 'log'), corr = matern.gcgc2( matern.kappa = 0.5, nugget = T), 
       effort = 1, max.range = 100, nrep = c(100, 1000), seed = 123)

####################################################################################################
#######################################################################################################################
#### OUTPUT FUNCTION: Prediction using GQT: serial version
#######################################################################################################################
#######################################################################################################################

#### Output: n(number of prediction locations)*2 matrix
#### First Column: predicting value
#### Second Column: Estimated MSPE

predGQT <- function(y, x = NULL, pred.x = NULL, obs.locs, pred.locs, family, corr,
                    sample.effort = 1, pred.effort = 1, max.range = NULL, max.count = NULL)
{
  x <- cbind(rep(1,length(y)), x)
  
  if(!is.matrix(pred.locs) & !is.data.frame(pred.locs)) 
    stop("Input 'pred.locs' must be a data frame or matrix!")
  
  if(!is.matrix(obs.locs) & !is.data.frame(obs.locs)) 
    stop("Input 'obs.locs' must be a data frame or matrix!")
  
  if(length(sample.effort) == 1) sample.effort <- rep(sample.effort, nrow(obs.locs))
  if(!length(sample.effort) == nrow(obs.locs)) 
    stop("Sampling Effort must be equal to the number of sampling locations!")
  
  if(length(pred.effort) == 1) pred.effort <- rep(pred.effort, nrow(pred.locs))
  
  if(!length(pred.effort) == nrow(pred.locs)) 
    stop("Prediction Effort must be equal to the number of prediction locations!")
  
  #### First calculate the MLE and output the log-likelihood as denominator
  
  MLE.est <- mleGQT(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                    effort = sample.effort, max.range = max.range)
  
  loglik <- MLE.est$log.lik; estpar <- MLE.est$MLE
  
  if(!is.matrix(pred.x) & !is.data.frame(pred.x) & !is.null(pred.x))
    stop("Input pred.x must be a data frame or matrix")
  
  if(is.null(pred.x)) pred.x <- matrix(1, nrow = nrow(pred.locs), ncol = 1)
  else {pred.x <- cbind(rep(1, nrow(pred.x)) , pred.x)}
  
  if(!is.matrix(pred.x) & !is.data.frame(pred.x))
    stop("Input pred.x must be a data frame or matrix")
  
  if(nrow(pred.x)!= nrow(pred.locs))
    stop("Number of prediction locations did not match the number of covariates") 
  
  if(nrow(pred.locs) == 1){
    indexloc <- which.min(FNN::get.knnx(pred.locs, obs.locs, 1)$nn.dist)
    m0 <- n0 <- round(unique(pred.effort)*y[indexloc]/sample.effort[indexloc])+1
  }else{
    m0 <- n0 <- round(pred.effort*apply(pred.locs, 1, function(x) y[which.min(FNN::get.knnx(t(as.matrix(x)), 
                obs.locs, 1)$nn.dist)]/sample.effort[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]))+1
  }
  
  NPL <- length(m0)
  #### m0 and n0: initial prediction values for probability search. Scalor or Vector
  
  if(is.null(max.count)) max.count <- ceiling(50*max(y))
  
  #### A for loop which cannot be avoided
  
  ans <- matrix(NA, nrow = NPL, ncol = 2)
  
  for(j in 1:NPL){  
    
    tmpfun2 <- function(xtmp){
      exp(-likGQT(pars = estpar, y = c(y, xtmp), x = rbind(as.matrix(x), pred.x[j,]),
                  locs = rbind(obs.locs, pred.locs[j,]), family = family, corr = corr,
                  effort = c(sample.effort, pred.effort[j])) - loglik)
    }
    
    p.m0 <- p.n0 <- tmpfun2(m0[j]); mu.m0 <- mu.n0 <- p.m0*m0[j];  mu2.m0 <- mu2.n0 <- p.m0*m0[j]^2
    MM1 <- matrix(0, nrow = 2, ncol = 3); MM2 <- matrix(0, nrow = 2, ncol = 3)
    MM1[1,] <- c(p.m0, mu.m0, mu2.m0);  MM2[1,] <- c(p.n0, mu.n0, mu2.n0)
    
    p.m0 <- tmpfun2(m0[j]-1); mu.m0 <- p.m0*(m0[j]-1); 
    mu2.m0 <- p.m0*(m0[j]-1)^2;  MM1[2,] <- c(p.m0, mu.m0, mu2.m0)
    
    while( (p.m0 > sqrt(.Machine$double.eps) | MM1[nrow(MM1), 2] > MM1[nrow(MM1)-1, 2]) & m0[j] > 1)
    {
      p.m0 <- tmpfun2(m0[j]-2); mu.m0 <- p.m0*(m0[j]-2); mu2.m0 <- p.m0*(m0[j]-2)^2
      MM1 <- rbind(MM1, c(p.m0, mu.m0, mu2.m0));   m0[j] <- m0[j]-1
    }
    #### Search from n0 to the right
    p.n0 <- tmpfun2(n0[j]+1); mu.n0 <- p.n0*(n0[j]+1);  mu2.n0 <- p.n0*(n0[j]+1)^2
    MM2[2, ] <- c(p.n0, mu.n0, mu2.n0)
    
    while( (p.n0 > sqrt(.Machine$double.eps) | MM2[nrow(MM2), 2] > MM2[nrow(MM2)-1, 2]) & n0[j] < max.count)
    {
      p.n0 <- tmpfun2(n0[j]+2); mu.n0 <- p.n0*(n0[j]+2);   mu2.n0 <- p.n0*(n0[j]+2)^2
      MM2 <- rbind(MM2, c(p.n0, mu.n0, mu2.n0));  n0[j] <- n0[j]+1
    }
    
    MM2 <- MM2[-1, ];  MM.all <- rbind(MM1, MM2)
    
    #### Due to some Monte Carlo Error, put sum to 1 constraint
    weight <- 1/sum(MM.all[,1])
    ans[j, ] <- c(sum(MM.all[,2])*weight, sum(MM.all[,3]*weight)-(sum(MM.all[,2])*weight)^2)
  } 
  return(ans)
}


## Example 
predGQT(y = simtmp.y, x = cbind(xloc,yloc), obs.locs = cbind(xloc, yloc), 
        pred.locs = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),
        family = binomial_t.gcgc2(df.t = 5), 
        corr = matern.gcgc2(matern.kappa = 0.5, nugget = T), 
        sample.effort = rep(c(11,19),50), pred.effort = c(10,12,12,19), 
        pred.x = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T), max.range = NULL, max.count = NULL)


#######################################################################################################################
#######################################################################################################################
#### OUTPUT FUNCTION: Prediction using GQT: parallel version via snowfall
#######################################################################################################################
#######################################################################################################################
predGQT.sf <- function(y, x = NULL, pred.x = NULL, obs.locs, pred.locs, family, corr, 
                       sample.effort = 1, pred.effort = 1, max.range = NULL, max.count = NULL,
                       n.cores = 2, cluster.type="SOCK")
{
  x <- cbind(rep(1,length(y)), x)
  
  if(!is.matrix(pred.locs) & !is.data.frame(pred.locs)) 
    stop("Input 'pred.locs' must be a data frame or matrix!")
  
  if(!is.matrix(obs.locs) & !is.data.frame(obs.locs)) 
    stop("Input 'obs.locs' must be a data frame or matrix!")
  
  if(length(sample.effort) == 1) sample.effort <- rep(sample.effort, nrow(obs.locs))
  
  if(!length(sample.effort) == nrow(obs.locs)) 
    stop("Sampling Effort must be equal to the number of sampling locations!")
  
  if(length(pred.effort) == 1) pred.effort <- rep(pred.effort, nrow(pred.locs))
  
  if(!length(pred.effort) == nrow(pred.locs)) 
    stop("Prediction Effort must be equal to the number of prediction locations!")
  
  #### First calculate the MLE and output the log-likelihood as denominator
  
  MLE.est <- mleGQT(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                    effort = sample.effort, max.range = max.range)
  
  loglik <- MLE.est$log.lik; estpar <- MLE.est$MLE
  
  if(!is.matrix(pred.x) & !is.data.frame(pred.x) & !is.null(pred.x))
    stop("Input pred.x must be a data frame or matrix")
  
  if(is.null(pred.x)) pred.x <- matrix(1, nrow = nrow(pred.locs), ncol = 1)
  else {pred.x <- cbind(rep(1, nrow(pred.x)) , pred.x)}
  
  if(!is.matrix(pred.x) & !is.data.frame(pred.x))
    stop("Input pred.x must be a data frame or matrix")
  
  if(nrow(pred.x)!= nrow(pred.locs))
    stop("Number of prediction locations did not match the number of covariates") 
  
  #### Begin to parallel

  if(nrow(pred.locs) == 1){
    indexloc <- which.min(FNN::get.knnx(pred.locs, obs.locs, 1)$nn.dist)
    m0 <- n0 <- round(unique(pred.effort)*y[indexloc]/sample.effort[indexloc])+1
  }else{
    m0 <- n0 <- round(pred.effort*apply(pred.locs, 1, function(x) y[which.min(FNN::get.knnx(t(as.matrix(x)),
         obs.locs, 1)$nn.dist)]/sample.effort[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]))+1
  }
  NPL <- length(m0)
  if(is.null(max.count)) max.count <- ceiling(50*max(y))
  #  ans <- matrix(NA, nrow = NPL, ncol = 2); 

  #### Begin to parallel
  
    if (requireNamespace("snowfall", quietly = TRUE)) { 
    snowfall::sfInit(parallel=TRUE, cpus= n.cores, type= cluster.type)
    suppressMessages(snowfall::sfExportAll( except=NULL, debug=FALSE ))
    suppressMessages(snowfall::sfLibrary(Rcpp))
    suppressMessages(snowfall::sfLibrary(RcppArmadillo))
    suppressMessages(snowfall::sfLibrary(geoR))
    suppressMessages(snowfall::sfLibrary(VGAM)) 
 
    #  snowfall::sfClusterSetupRNG()
    #  sfClusterEvalQ( ls() )

    par.pred.inner <- function(j){
      #### Q1: How to avoid function reloading and namespace problem. 
      #### Loading the own package directly may avoid this problem.
      ########################################################################################
      ########################################################################################
      ###################### FUNCTION RELOADING ############################################
      ########################################################################################
      cppFunction(depends = 'RcppArmadillo',   code = '
          double ldgc(arma::rowvec u, arma::mat sigma) { 
                  int dimn = sigma.n_rows; 
                  arma::mat RI = arma::inv(sigma)-arma::speye(dimn, dimn);
                  arma::rowvec phi(dimn);
                  for(int i=0; i<dimn; i++){ 
                  phi[i] = Rcpp::stats::qnorm_0(u[i], 1, 0);
                  }
                  arma::mat outtmp = phi*RI*phi.t();
                  double logdet = sum(arma::log(arma::eig_sym(sigma)));
                  double out = -0.5*(logdet+outtmp(0));
                  if(std::isinf(out) == true || std::isnan(out) == true){
                  out = -9.532493e+14;
                  }
                  return(out);
    }')

      likGQT <- function(pars, y, x = NULL, locs, family, corr, effort){
        matD <- as.matrix(dist(locs, method = "euclidean", diag = T, upper = T))
        npars <- length(pars)
        nparcorr <- corr$npar.cor
        R <- corr$corr(pars[(npars-nparcorr+1):npars], matD = matD)
        logindpdf <- sum(family$pdf(y = y, x = x, pars = pars[1:(npars-nparcorr)], effort = effort))
        u1 <- family$cdf(y = y-1, x = x, pars = pars, effort = effort)
        u2 <- family$cdf(y = y, x = x, pars = pars, effort = effort)
        u <- (u1+u2)*0.5
        loglik <- ldgc(u = u, sigma = R)+logindpdf
        if(is.nan(loglik)) loglik = -1e6 
        if (loglik == Inf) loglik = 1e6
        return(-loglik)
      }
 #########################################################################################################
 #########################################################################################################
 ###################### FUNCTION RELOADING END ###########################################################
 #########################################################################################################
      
      tmpfun2 <- function(xtmp){
        exp(-likGQT(pars = estpar, y = c(y, xtmp), x = rbind(as.matrix(x), pred.x[j,]),
                    locs = rbind(obs.locs, pred.locs[j,]), family = family, corr = corr,
                    effort = c(sample.effort, pred.effort[j])) - loglik)
      }
      
      p.m0 <- p.n0 <- tmpfun2(m0[j]); mu.m0 <- mu.n0 <- p.m0*m0[j];  mu2.m0 <- mu2.n0 <- p.m0*m0[j]^2
      MM1 <- matrix(0, nrow = 2, ncol = 3); MM2 <- matrix(0, nrow = 2, ncol = 3)
      MM1[1,] <- c(p.m0, mu.m0, mu2.m0);  MM2[1,] <- c(p.n0, mu.n0, mu2.n0)
      
      # To avoid 0 length, repeated here
      
      p.m0 <- tmpfun2(m0[j]-1); mu.m0 <- p.m0*(m0[j]-1)
      mu2.m0 <- p.m0*(m0[j]-1)^2;  MM1[2,] <- c(p.m0, mu.m0, mu2.m0)
      
      while( (p.m0 > sqrt(.Machine$double.eps) |  MM1[nrow(MM1), 2] > MM1[nrow(MM1)-1, 2]) & m0[j] > 1)
      {
        p.m0 <- tmpfun2(m0[j]-2); mu.m0 <- p.m0*(m0[j]-2); mu2.m0 <- p.m0*(m0[j]-2)^2
        MM1 <- rbind(MM1, c(p.m0, mu.m0, mu2.m0));  m0[j] <- m0[j]-1
      } 
      
      #### Search from n0 to the right
      
      p.n0 <- tmpfun2(n0[j]+1); mu.n0 <- p.n0*(n0[j]+1);  mu2.n0 <- p.n0*(n0[j]+1)^2
      MM2[2, ] <- c(p.n0, mu.n0, mu2.n0)
      
      while( (p.n0 > sqrt(.Machine$double.eps) | MM2[nrow(MM2), 2] > MM2[nrow(MM2)-1, 2]) & n0[j] < max.count)
      {
        p.n0 <- tmpfun2(n0[j]+2); mu.n0 <- p.n0*(n0[j]+2); mu2.n0 <- p.n0*(n0[j]+2)^2
        MM2 <- rbind(MM2, c(p.n0, mu.n0, mu2.n0));  n0[j] <- n0[j]+1
      }
      MM2 <- MM2[-1, ];  MM.all <- rbind(MM1, MM2)
      
      #### Due to some Monte Carlo Error, put sum to 1 constraint
      weight <- 1/sum(MM.all[,1])
      return( c(sum(MM.all[,2])*weight, sum(MM.all[,3]*weight)-(sum(MM.all[,2])*weight)^2) )
    }
    out = sfLapply(1:NPL, par.pred.inner)
    sfStop()
 }else{
   stop("Please install {snowfall} first before using this function!")
 }
   return(out)
}


#### Example 
predGQT.sf(y = simtmp.y, x = cbind(xloc,yloc), obs.locs = cbind(xloc, yloc), 
        pred.locs = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),
        family = binomial_t.gcgc2(df.t = 5), 
        corr = matern.gcgc2(matern.kappa = 0.5, nugget = T), 
        sample.effort = rep(c(11,19),50), pred.effort = c(10,12,12,19), 
        pred.x = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T), 
        max.range = 5, max.count = 15, n.cores = 4)
