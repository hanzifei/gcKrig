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
likGQT <- function(param, y, x = NULL, obs.locs, family, nugget = FALSE, matern.kappa = 0.5){
  nparam <- length(param)
  D <- as.matrix(dist(obs.locs, method = "euclidean", diag = T, upper = T))
  if(nugget == FALSE | missing(nugget)){
    R <- geoR::matern(D, param[nparam], matern.kappa)
  }else{
    R <- (1-param[nparam])*geoR::matern(D, param[nparam-1], matern.kappa) + param[nparam]*diag(ncol(D))
  }
  if( is.null(x) ){
    M <- rep(exp(param[1]), length(y));  S <- param[2]
  }else{
    M <- exp(param[1] + param[2:(ncol(x)+1)]%*%t(x));  S <- param[ncol(x)+2]
  }
  if (family == "Poisson")  loglik <- gqtlikPois(count = y, mu = M, sigma = R)
  if (family == "ZIP")   loglik <- gqtlikZIP(count = y, mu = M, od = S, sigma = R)
  if (family == "NB")   loglik <- gqtlikNB(count = y, mu = M, od = S, sigma = R)
  if(is.nan(loglik)) loglik = -1e6
  if (loglik == Inf) loglik = 1e6
  return(-loglik)
}

#mleGQT( y = sim.y, obs.locs = cbind(xloc,yloc), family = 'NB')

mleGQT <- function(y, x = NULL, obs.locs, family, nugget = FALSE, 
                    matern.kappa = 0.5, max.range = NULL)
{
  #### Create initial values in optimization 
  D <- as.matrix(dist(obs.locs, method = "euclidean", diag = T, upper = T))
  if(is.null(x))  mfit <- glm.fit(rep(1,length(y)), y, family = poisson(link='log'))
    if(!is.matrix(x) & !is.data.frame(x) &!is.null(x)){
    stop("Input 'x' must be a matrix or data frame!")} else{
      mfit <- glm.fit(cbind(rep(1,length(y)), x), y, family = poisson(link='log'))
    }
  reg0 <- coef(mfit); mu <- fitted(mfit); od0 <- NULL
  if (family %in% c("NB", "ZIP")) od0 <- max(10*.Machine$double.eps, mean(((y-mu)^2-mu)/mu^2))
  if(is.null(max.range)) max.range <- 5*max(D)
  range0 <- median(D)
  if(nugget == FALSE | missing(nugget)){
    nugget0 <- NULL
  }else{
    nugget0 <- 0.2
  }
  par0 <- c(reg0, od0, range0, nugget0)
  n.reg0 <- length(reg0); n.od0 <- length(od0); n.range0 <- length(range0)
  n.nugget0 <- length(nugget0); MIN <- .Machine$double.eps^0.25
  fit <- optim(par = par0, fn = likGQT, y = y, x = x, obs.locs = obs.locs, family = family, 
                  nugget = nugget, matern.kappa = matern.kappa, method = "L-BFGS-B",
                  lower = c(rep(-Inf,n.reg0), rep(MIN, n.od0), rep(MIN, n.range0), rep(0, n.nugget0)),
                  upper = c(rep(Inf,n.reg0), rep(Inf, n.od0), rep(max.range, n.range0), rep(1-MIN, n.nugget0)))
  if (fit$par[n.reg0+n.od0+1] == max.range)
    warning("Maximum Range Value Reached, MLE may not exist. Try to input a larger One")
  if (is.null(fit$convergence) || fit$convergence != 0)   
    stop("Maximum likelihood estimation failed. Algorithm does not converge") 
  result.list <- list()
  result.list$MLE <- fit$par
  result.list$log.lik <- -fit$value
 # class(result.list) <- c( "pred.ghk")
  return(result.list)
}


#######################################################################################################################
#######################################################################################################################
#### OUTPUT FUNCTION: Prediction using GQT: serial version
#######################################################################################################################
#######################################################################################################################

#### Output: n(number of prediction locations)*2 matrix
#### First Column: predicting value
#### Second Column: Estimated MSPE

predGQT <- function(y, x = NULL, pred.x = NULL, obs.locs, pred.locs, family, nugget = FALSE, 
                      matern.kappa = 0.5, max.range = NULL, max.count = NULL)
{
  #  pred.locs = matrix(c(0.4,0.4,0.5,0.5),2,2,byrow = T)
  #  y = sim.nbnew1[1,1:200]; obs.locs = cbind(xloc,yloc)
  #  family = 'NB'; nugget = FALSE; matern.kappa = 0.5; nrep = 1000;
  #  seed = 123; max.range = 4; max.count = 20
  
  #### Nearest Neighbour Search and output the index.
  #### if pred.locs only contains one location info
  
  #  if((!is.matrix(pred.locs) & !is.data.frame(pred.locs)) | 
  #     (is.matrix(pred.locs) & nrow(pred.locs) == 1)) 
  #     { pred.locs = matrix(pred.locs, nrow = 1)
  #      m0 = n0 = y[which.min(FNN::get.knnx(pred.locs, obs.locs, 1)$nn.dist)] + 1
  
  if(!is.matrix(pred.locs) & !is.data.frame(pred.locs)) 
    stop("Input 'pred.locs' must be a data frame or matrix!")
  if(!is.matrix(obs.locs) & !is.data.frame(obs.locs)) 
    stop("Input 'obs.locs' must be a data frame or matrix!")
  if(nrow(pred.locs) == 1){
    m0 <- n0 <- y[which.min(FNN::get.knnx(pred.locs, obs.locs, 1)$nn.dist)] + 1
  }else{
    m0 <- n0 <- apply(pred.locs, 1, function(x) y[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]+1 )
  }
  NPL <- length(m0)  
  #### m0 and n0: initial prediction values for probability search. Scalor or Vector
  
  if(is.null(max.count)) max.count <- ceiling(50*max(y))
  MLE.est <- mleGQT(y = y, x = x, obs.locs = obs.locs, family = family, nugget = nugget, 
                     matern.kappa = matern.kappa, max.range = max.range)
  loglik <- MLE.est$log.lik; estpar <- MLE.est$MLE
  
  if(is.null(x))  x <- rep(1, length(y))
  if(is.null(pred.x)) pred.x <- matrix(1, nrow = NPL, ncol = 1)
  if(!is.matrix(pred.x) & !is.data.frame(pred.x))
    stop("Input pred.x must be a data frame or matrix")
  if(nrow(pred.x)!= nrow(pred.locs))
    stop("Number of prediction locations did not match number of covariates") 
  
  ans <- matrix(NA, nrow = NPL, ncol = 2)
  
  for(j in 1:NPL){  
    p.m0 <- p.n0 <- exp(-likGQT(param = estpar, y = c(y, m0[j]), x = rbind(as.matrix(x), pred.x[j,]),
                                 obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                                 matern.kappa = matern.kappa) - loglik)
    
    mu.m0 <- mu.n0 <- p.m0*m0[j];  mu2.m0 <- mu2.n0 <- p.m0*m0[j]^2
    MM1 <- matrix(0, nrow = 2, ncol = 3); MM2 <- matrix(0, nrow = 2, ncol = 3)
    MM1[1,] <- c(p.m0, mu.m0, mu2.m0);  MM2[1,] <- c(p.n0, mu.n0, mu2.n0)
    
    p.m0 <- exp(-likGQT(param = estpar, y = c(y, m0[j] -1), x = rbind(as.matrix(x), pred.x[j,]),
                         obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                         matern.kappa = matern.kappa) - loglik)
    mu.m0 <- p.m0*(m0[j]-1); mu2.m0 <- p.m0*(m0[j]-1)^2;  MM1[2,] <- c(p.m0, mu.m0, mu2.m0)
    
    while( (p.m0 > sqrt(.Machine$double.eps) | MM1[nrow(MM1), 2] > MM1[nrow(MM1)-1, 2]) & m0[j] > 1)
    {
      p.m0 <- exp(-likGQT(param = estpar, y = c(y, m0[j]-2), x = rbind(as.matrix(x), pred.x[j,]),
                           obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                           matern.kappa = matern.kappa) - loglik)
      mu.m0 <- p.m0*(m0[j]-2); mu2.m0 <- p.m0*(m0[j]-2)^2
      MM1 <- rbind(MM1, c(p.m0, mu.m0, mu2.m0));   m0[j] <- m0[j]-1
    }
    #### Search from n0 to the right
    p.n0 <- exp(-likGQT(param = estpar, y = c(y, n0[j]+1), x = rbind(as.matrix(x), pred.x[j,]),
                         obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                         matern.kappa = matern.kappa) - loglik)
    mu.n0 <- p.n0*(n0[j]+1);  mu2.n0 <- p.n0*(n0[j]+1)^2
    MM2[2, ] <- c(p.n0, mu.n0, mu2.n0)
    
    while( (p.n0 > sqrt(.Machine$double.eps) | MM2[nrow(MM2), 2] > MM2[nrow(MM2)-1, 2]) & n0[j] < max.count)
    {
      p.n0 <- exp(-likGQT(param = estpar, y = c(y, n0[j]+2), x = rbind(as.matrix(x), pred.x[j,]),
                           obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                           matern.kappa = matern.kappa) - loglik)
      mu.n0 <- p.n0*(n0[j]+2);   mu2.n0 <- p.n0*(n0[j]+2)^2
      MM2 <- rbind(MM2, c(p.n0, mu.n0, mu2.n0));  n0[j] <- n0[j]+1
    }
    
    MM2 <- MM2[-1, ];  MM.all <- rbind(MM1, MM2)
    
    #### Due to some Monte Carlo Error, put sum to 1 constraint
    weight <- 1/sum(MM.all[,1])
    ans[j, ] <- c(sum(MM.all[,2])*weight, sum(MM.all[,3]*weight)-(sum(MM.all[,2])*weight)^2)
  } 
  return(ans)
}



#######################################################################################################################
#######################################################################################################################
#### OUTPUT FUNCTION: Prediction using GQT: parallel version via snowfall
#######################################################################################################################
#######################################################################################################################
predGQT.sf <- function(y, x = NULL, pred.x = NULL, obs.locs, pred.locs, family, nugget = FALSE, 
                          matern.kappa = 0.5, max.range = NULL, max.count = NULL,
                          n.cores = 2, cluster.type="SOCK")
{
  #### Nearest Neighbour Search and output the index.
  #### if pred.locs only contains one location info
  #  if((!is.matrix(pred.locs) & !is.data.frame(pred.locs)) | 
  #     (is.matrix(pred.locs) & nrow(pred.locs) == 1)) 
  #     { pred.locs = matrix(pred.locs, nrow = 1)
  #      m0 = n0 = y[which.min(FNN::get.knnx(pred.locs, obs.locs, 1)$nn.dist)] + 1
  
  if(!is.data.frame(pred.locs) & !is.matrix(pred.locs))
    stop(" 'pred.locs' must be of a data frame or matrix")
  if(!is.data.frame(obs.locs) & !is.matrix(obs.locs))
    stop("'obs.locs' must be a data frame or matrix")
  
  m0 <- n0 <- apply(pred.locs, 1, function(x) y[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]+1 )
  NPL <- length(m0)  
  
  #### m0 and n0: initial prediction values for probability search. Scalor or Vector
  
  if(is.null(max.count)) max.count <- ceiling(50*max(y))
  
  #### First calculate the MLE and output the log-likelihood as denominator
  
  MLE.est <- mleGQT(y = y, x = x, obs.locs = obs.locs, family = family, nugget = nugget, 
                     matern.kappa = matern.kappa, max.range = max.range)
  
  loglik <- MLE.est$log.lik;  estpar <- MLE.est$MLE
  
  if(is.null(x))  x <- rep(1, length(y))
  if(is.null(pred.x)) pred.x <- matrix(1, nrow = NPL, ncol = 1)
  
  if(!is.matrix(pred.x) & !is.data.frame(pred.x))
    stop("Input pred.x must be a data frame or matrix")
  if(nrow(pred.x)!= nrow(pred.locs))
    stop("Number of prediction locations did not match number of covariates") 
  
  #### Begin to parallel
  
  # ans <- matrix(NA, nrow = NPL, ncol = 2)
  if (requireNamespace("snowfall", quietly = TRUE)) {    
    snowfall::sfInit(parallel=TRUE, cpus= n.cores, type= cluster.type)
    snowfall::sfExportAll( except=NULL, debug=FALSE )
    snowfall::sfLibrary(Rcpp);snowfall::sfLibrary(RcppArmadillo)
    snowfall::sfLibrary(geoR);snowfall::sfLibrary(VGAM)
    #snowfall::sfLibrary("geoCount", character.only= TRUE)
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
          double gqtlikNB(arma::vec count, arma::vec mu, double od, arma::mat sigma) { 
                  int dimn = sigma.n_rows; 
                  double SZ = 1/od; arma::vec PB = SZ/(SZ+mu);
                  arma::mat RI = arma::inv(sigma)-arma::speye(dimn, dimn);
                  arma::rowvec phi(dimn); arma::rowvec u(dimn);
                  
                  double loglik = 0;  double logdet = sum(arma::log(arma::eig_sym(sigma)));
                  
                  for(int i=0; i<dimn; i++){ 
                  u[i] = (R::pnbinom(count[i], SZ, PB[i], 1, 0) + R::pnbinom(count[i]-1, SZ, PB[i], 1, 0))/2;
                  phi[i] = Rcpp::stats::qnorm_0(u[i], 1, 0);
                  loglik += R::dnbinom(count[i], SZ, PB[i], 1);
                  }
                  arma::mat outtmp = phi*RI*phi.t();
                  double out = -0.5*(logdet+outtmp(0))+loglik;
                  if(std::isinf(out) == true || std::isnan(out) == true){
                  out = -9.532493e+14;
                  }
                  return(out);
    }')

      cppFunction(depends = 'RcppArmadillo',   code = '
         double gqtlikPois(arma::vec count, arma::vec mu, arma::mat sigma) { 
  
                  int dimn = sigma.n_rows; 
                  arma::mat RI = arma::inv(sigma)-arma::speye(dimn, dimn);
                  arma::rowvec phi(dimn); arma::rowvec u(dimn);
                  
                  double loglik = 0;  double logdet = sum(arma::log(arma::eig_sym(sigma)));
                  
                  for(int i=0; i<dimn; i++){ 
                  u[i] = (R::ppois(count[i], mu[i], 1, 0) + R::ppois(count[i]-1, mu[i], 1, 0))/2;
                  phi[i] = Rcpp::stats::qnorm_0(u[i], 1, 0);
                  loglik += R::dpois(count[i], mu[i], 1);
                  }
                  arma::mat outtmp = phi*RI*phi.t();
                  double out = -0.5*(logdet+outtmp(0))+loglik;
                  if(std::isinf(out) == true || std::isnan(out) == true){
                  out = -9.532493e+14;
                  }
                  return(out);
    }')
      
      cppFunction(depends = 'RcppArmadillo',   code = '
       double gqtlikZIP(arma::vec count, arma::vec mu, double od, arma::mat sigma) { 
  
                  int dimn = sigma.n_rows; 
                  arma::mat RI = arma::inv(sigma)-arma::speye(dimn, dimn);
                  arma::rowvec phi(dimn); arma::rowvec u(dimn);
                  
                  double loglik = 0;  double logdet = sum(arma::log(arma::eig_sym(sigma)));
                  double tmp = 0;
                  
                  for(int i=0; i<dimn; i++){ 
                  
                  if(count[i] == 0){
                  u[i] = (od/(1+od) + exp(-(1+od)*mu[i])/(1+od))/2;
                  tmp = log(od/(1+od) + exp(-(1+od)*mu[i])/(1+od));
                  }
                  if(count[i] == 1){
                  u[i] = (exp(-(1+od)*mu[i])/(1+od)+R::ppois(1, (1+od)*mu[i], 1, 0)/(1+od))/2 + od/(1+od);
                  tmp = R::dpois(1, (1+od)*mu[i], 1) - log(1+od);
                  }
                  if(count[i] > 1){
                  u[i] = (R::ppois(count[i], (1+od)*mu[i], 1, 0)/(1+od)+
                  R::ppois(count[i]-1, (1+od)*mu[i], 1, 0)/(1+od))/2+od/(1+od);
                  tmp = R::dpois(count[i], (1+od)*mu[i], 1) - log(1+od);
                  }
                  phi[i] = Rcpp::stats::qnorm_0(u[i], 1, 0);
                  
                  loglik += tmp;
                  }
                  arma::mat outtmp = phi*RI*phi.t();
                  double out = -0.5*(logdet+outtmp(0))+loglik;
                  if(std::isinf(out) == true || std::isnan(out) == true){
                  out = -9.532493e+14;
                  }
                  return(out);
    }')
      
      cppFunction(depends = 'RcppArmadillo',   code = '
      double gqtlikBin(arma::vec count, arma::vec p, double N, arma::mat sigma) { 

                  int dimn = sigma.n_rows; 
                  //  double SZ = 1/od; arma::vec PB = SZ/(SZ+mu);
                  arma::mat RI = arma::inv(sigma)-arma::speye(dimn, dimn);
                  arma::rowvec phi(dimn); arma::rowvec u(dimn);
                  
                  double loglik = 0;  double logdet = sum(arma::log(arma::eig_sym(sigma)));
                  
                  for(int i=0; i<dimn; i++){ 
                  u[i] = (R::pbinom(count[i], N, p[i], 1, 0) + R::pbinom(count[i]-1, N, p[i], 1, 0))/2;
                  phi[i] = Rcpp::stats::qnorm_0(u[i], 1, 0);
                  loglik += R::dbinom(count[i], N, p[i], 1);
                  }
                  arma::mat outtmp = phi*RI*phi.t();
                  double out = -0.5*(logdet+outtmp(0))+loglik;
                  if(std::isinf(out) == true || std::isnan(out) == true){
                  out = -9.532493e+14;
                  }
                  return(out);
    }')
      
      likGQT <- function(param, y, x = NULL, obs.locs, family, nugget = FALSE, 
                          matern.kappa = 0.5){
        nparam <- length(param)
        D <- as.matrix(dist(obs.locs, method = "euclidean", diag = T, upper = T))
        if(nugget == FALSE | missing(nugget)){
          R <- geoR::matern(D, param[nparam], matern.kappa)
        }else{
          R <- (1-param[nparam])*geoR::matern(D, param[nparam-1], matern.kappa) + param[nparam]*diag(ncol(D))
        }
        if(is.null(x)){
          M <- exp(param[1]) ;S <- 1/param[2]
        }else{
          M <- exp(param[1] + param[2:(ncol(x)+1)]%*%t(x));  S <- 1/param[ncol(x)+2]
        }
        if (family == "Poisson")  loglik <- gqtlikPois(count = y, mu = M, sigma = R)
        if (family == "ZIP")  loglik <- gqtlikZIP(count = y, mu = M, od = S, sigma = R)
        if (family == "NB")  loglik <- gqtlikNB(count = y, mu = M, od = S, sigma = R)
        if (loglik == Inf) loglik = 1e6
        return(-loglik)
      }
      
      ########################################################################################
      ########################################################################################
      ###################### FUNCTION RELOADING END ############################################
      ########################################################################################
      
      p.m0 <- p.n0 <- exp(-likGQT(param = estpar, y = c(y, m0[j]), x = rbind(as.matrix(x), pred.x[j,]), 
                                   obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                                   matern.kappa = matern.kappa) - loglik)
      mu.m0 <- mu.n0 <- p.m0*m0[j];  mu2.m0 <- mu2.n0 <- p.m0*m0[j]^2
      
      MM1 <- matrix(0, nrow = 2, ncol = 3); MM2 <- matrix(0, nrow = 2, ncol = 3)
      MM1[1,] <- c(p.m0, mu.m0, mu2.m0);  MM2[1,] <- c(p.n0, mu.n0, mu2.n0)
      # To avoid 0 length, repeated here
      
      p.m0 <- exp(-likGQT(param = estpar, y = c(y, m0[j] -1), x = rbind(as.matrix(x), pred.x[j,]), 
                           obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                           matern.kappa = matern.kappa) - loglik)
      mu.m0 <- p.m0*(m0[j]-1);  mu2.m0 <- p.m0*(m0[j]-1)^2;  MM1[2,] <- c(p.m0, mu.m0, mu2.m0)
      
      while( (p.m0 > sqrt(.Machine$double.eps) |  MM1[nrow(MM1), 2] > MM1[nrow(MM1)-1, 2]) & m0[j] > 1)
      {
        p.m0 <- exp(-likGQT(param = estpar, y = c(y, m0[j]-2), x = rbind(as.matrix(x), pred.x[j,]), 
                             obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                             matern.kappa = matern.kappa) - loglik)
        
        mu.m0 <- p.m0*(m0[j]-2);   mu2.m0 <- p.m0*(m0[j]-2)^2
        MM1 <- rbind(MM1, c(p.m0, mu.m0, mu2.m0));  m0[j] <- m0[j]-1
      } 
      
      #### Search from n0 to the right
      
      p.n0 <- exp(-likGQT(param = estpar, y = c(y, n0[j]+1),x = rbind(as.matrix(x), pred.x[j,]), 
                           obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                           matern.kappa = matern.kappa) - loglik)
      mu.n0 <- p.n0*(n0[j]+1);  mu2.n0 <- p.n0*(n0[j]+1)^2
      MM2[2, ] <- c(p.n0, mu.n0, mu2.n0)
      
      while( (p.n0 > sqrt(.Machine$double.eps) | MM2[nrow(MM2), 2] > MM2[nrow(MM2)-1, 2]) & n0[j] < max.count)
      {
        p.n0 <- exp(-likGQT(param = estpar, y = c(y, n0[j]+2), x = rbind(as.matrix(x), pred.x[j,]), 
                             obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                             matern.kappa = matern.kappa) - loglik)
        mu.n0 <- p.n0*(n0[j]+2);    mu2.n0 <- p.n0*(n0[j]+2)^2
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







#predGQT(y = sim.y, obs.locs = cbind(xloc,yloc), pred.locs = matrix(c(0.4, 0.4),ncol = 2),
#         family = 'NB', nugget = FALSE, matern.kappa = 0.5, max.range = 4, max.count = 20)

#predGQT.sf(y = sim.y, obs.locs = cbind(xloc, yloc), 
#              pred.locs = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),family = 'ZIP', 
#              nugget = F, max.range = 2, max.count = 15, n.cores = 4, cluster.type = "SOCK")

#predGQT.sf(y = simtmp.y, x = cbind(xloc,yloc), obs.locs = cbind(xloc, yloc), 
#              pred.x = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),
#              pred.locs = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),family = 'NB', 
#              nugget = T, max.range = 2, max.count = 15, n.cores = 4, cluster.type = "SOCK")

predGHK.sf(y = simtmp.y, x = cbind(xloc,yloc), obs.locs = cbind(xloc, yloc), 
           pred.x = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),
           pred.locs = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),family = 'NB', 
           nugget = T, max.range = 2, max.count = 15, n.cores = 4, cluster.type = "SOCK")
