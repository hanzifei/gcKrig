library(Rcpp)
library(RcppArmadillo)
library(geoR)
library(VGAM)
library(MASS)
library(snowfall)

#### Likelihood Expression 
#### L: lower cholesky decomposition of correlation matrix
#### a: lower bound of the integration
#### b: upper bound of the intergration

cppFunction(depends = 'RcppArmadillo',   code = '
                  double ghkLik(arma::mat L, arma::vec lower, arma::vec upper, int nrep) {
            
            int dim = L.n_cols;double z_[dim];
            double res_ = 0, prod, mu_, eta_, gamma_, ans_;
            const double EPS = std::numeric_limits<double>::min(); 
            
            double* eta = & eta_; double* res = & res_; double* ans = &ans_;
            double* gamma = & gamma_; double* mu = & mu_; double* z = &z_[0];
            
            for(int i=0;i<nrep;i++) {
            prod = 1.0;
            for(int j=0;j<dim;j++) {     
            *mu=0; 
            for(int k=0;k<j;k++)  *mu += L(j,k)*(*(z+k));
            *eta = Rcpp::stats::pnorm_1((lower(j)-*mu)/L(j,j),0,1,0);
            *gamma = Rcpp::stats::pnorm_1((upper(j)-*mu)/L(j,j),0,1,0);
            prod *= *gamma-*eta; 
            double u = Rf_runif(0, 1);
            *ans=u*(*gamma)+(1-u)*(*eta);
            if(*ans > 1-EPS) *ans= 1-EPS;
            if(*ans < EPS) *ans= EPS ;
            *(z+j) = Rcpp::stats::qnorm_1(*ans,0,1,0);
            }
            *res += prod;
            }
            *res /= nrep;
            if(*res == 0){
            return -1e6;
            }
            else{
            return log(*res);
            }
            }')


#######################################################################################################################
#######################################################################################################################
#### OUTPUT FUNCTION: Likelihood Evaluation (multivariate normal rectange probability) ####
#######################################################################################################################
#######################################################################################################################


cppFunction(depends = 'RcppArmadillo',   code = '
            double ghkLiknew(arma::vec mu, arma::mat L, arma::vec lower, arma::vec upper, int nrep) {
            
            int dim = L.n_cols;double z_[dim];
            double res_ = 0, prod, MU_, eta_, gamma_, ans_;
            const double EPS = std::numeric_limits<double>::min(); 
            
            double* eta = & eta_; double* res = & res_; double* ans = &ans_;
            double* gamma = & gamma_; double* MU = & MU_; double* z = &z_[0];
            
            for(int i=0;i<nrep;i++) {
            prod = 1.0;
            for(int j=0;j<dim;j++) {     
            *MU=0; 
            for(int k=0;k<j;k++)  *MU += L(j,k)*(*(z+k));
            
            *eta = Rcpp::stats::pnorm_1((lower(j)-*MU-mu(j))/L(j,j),0,1,0);
            *gamma = Rcpp::stats::pnorm_1((upper(j)-*MU-mu(j))/L(j,j),0,1,0);
            
            prod *= *gamma-*eta; 
            double u = Rf_runif(0, 1);
            *ans=u*(*gamma)+(1-u)*(*eta);
            if(*ans > 1-EPS) *ans= 1-EPS;
            if(*ans < EPS) *ans= EPS ;
            *(z+j) = Rcpp::stats::qnorm_1(*ans,0,1,0);
            }
            *res += prod;
            }
            *res /= nrep;
            if(*res == 0){
            return -1e6;
            }
            else{
            return log(*res);
            }
            }')

#### mu: mean of the truncated normal distribution
#### sigma: Covariace Matrix
#### lower: lower bound in truncation
#### upper: upper bound in truncation


mvnintGHK <- function(mu, sigma, lower, upper, nrep = 1000, log = T){
  
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
  if(log == T){
    ans <- ghkLiknew(mu = mu, L = L, lower = lower, upper = upper,nrep = nrep)
  }else{
    ans <- exp(ghkLiknew(mu = mu, L = L, lower = lower, upper = upper,nrep = nrep))
  }
  return(ans)
}


#set.seed(1234)
mvnintGHK(mu = rep(0.1, 200), sigma = exp(-D/0.2), lower = rep(0.1, 200), upper = rep(0.5,200),
          nrep =10000)

#set.seed(1234)
#log(pmvnorm(lower = rep(-1, 200), upper = rep(0.5,200), mean = rep(0.2,200), 
#        sigma = exp(-D/0.2))[1])



#######################################################################################################################
#######################################################################################################################
#### Inner function only used to be called by another functions inside package
#######################################################################################################################
#######################################################################################################################

#### param: Last one (or two) elements are: Range and Nugget (if applicable)
#### The first (1 + number of covariates) parameters are: regression parameters
#### Before Range parameter: overdispersion parameter (if applicable)
#### y: counts.
#### x: covariates
#### matern.kappa: smoothness parameter in matern family of correlations

likGHK <- function(param, y, x = NULL, obs.locs, family, nugget = FALSE, matern.kappa = 0.5, 
                   effort = sample.effort, nrep = 1000, seed = 123)
  {
  nparam <- length(param)
  D <- as.matrix(dist(obs.locs, method = "euclidean", diag = T, upper = T))

  if(nugget == FALSE){
    R <- geoR::matern(D, param[nparam], matern.kappa)
  }else{
    R <- (1-param[nparam])*geoR::matern(D, param[nparam-1], matern.kappa) + param[nparam]*diag(ncol(D))
  }
  
  L <- try(t(chol(R)),silent=TRUE)
   if( inherits(L,"try-error") ) stop("Cholesky Decomposition failed in likelihood evaluation. 
                                      Possible ill-conditioned matrix involved.")

  if(is.null(x)){
    M <- exp(param[1])*effort;  S <- 1/param[2]
  }else{
    M <- exp(param[1] + param[2:(ncol(x)+1)]%*%t(x))*effort;  S <- 1/param[ncol(x)+2]
  }
  
  if (family == "Poisson")
#### Even though the "size" parameter S is created, in this case, it is not used. 
  {
    a <- qnorm(ppois( y-1, lambda = M));  b <- qnorm(ppois( y, lambda = M))
  }
  if (family == "ZIP")
  {
    a <- qnorm(VGAM::pzipois(q = y-1, lambda = M+M/S, pstr0 = 1/(S+1)))
    b <- qnorm(VGAM::pzipois(q = y, lambda = M+M/S, pstr0 = 1/(S+1)))
  }
  if (family == "NB")
  {
    a <- qnorm(pnbinom( y-1, size = S, mu = M));  b <- qnorm(pnbinom( y, size = S, mu = M))
  }
  set.seed(seed)
  loglik <- ghkLik(L = L, lower = a, upper = b, nrep = nrep)
  if(is.nan(loglik)) loglik = -1e6
  if (loglik == Inf) loglik = 1e6
  return(-loglik)
}


likGHK(param = c(0,0,1,0.3,0.2), y = simtmp.y, x = cbind(xloc,yloc), 
       obs.locs = cbind(xloc,yloc), family = 'NB', nugget = FALSE, 
       matern.kappa = 0.5, nrep = 1000, seed = 123)
  



#######################################################################################################################
#######################################################################################################################
#### Inner function only used to be called by another functions inside package
#######################################################################################################################
#######################################################################################################################
#### Now Consider Add optimization procedure
#### Output: MLE and Value of the simulated log likelihood evaluated at the MLEs

mleGHK <- function(y, x = NULL, obs.locs, family, nugget = FALSE, matern.kappa = 0.5,
                     binomial.n = NULL, sample.effort = 1, nrep = 1000, seed = 123, max.range = NULL)
{
#  y = simtmp.y; x = cbind(xloc,yloc); obs.locs = cbind(xloc,yloc);
#  family = "NB"; nugget = FALSE; matern.kappa = 0.5; 
#  nrep = 1000; seed = 123; max.range = 10
#### Create initial values in optimization 
  D <- as.matrix(dist(obs.locs, method = "euclidean", diag = T, upper = T))
  if(is.null(x)) mfit <- glm.fit(rep(1,length(y)), y, family = poisson(link='log'))

  if(!is.matrix(x) & !is.data.frame(x) & !is.null(x)){
    stop("Input 'x' must be a matrix or data frame!")} else{
    mfit <- glm.fit(cbind(rep(1,length(y)), x), y, family = poisson(link='log'))
  }
  reg0 <- coef(mfit); mu <- fitted(mfit); od0 <- NULL
  
  if (family %in% c("NB", "ZIP")) od0 <- max(10*.Machine$double.eps, mean(((y-mu)^2-mu)/mu^2))
  if(is.null(max.range)) max.range <- 2*max(D)
  range0 <- median(D)
  
  if(nugget == FALSE){
    nugget0 <- NULL
  }else{
    nugget0 <- 0.2
  }
  
  par0 <- c(reg0, od0, range0, nugget0)
  n.reg0 <- length(reg0); n.od0 <- length(od0); n.range0 <- length(range0)
  n.nugget0 <- length(nugget0); MIN <- .Machine$double.eps^0.25
  
  fittmp <- optim(par = par0, fn = likGHK, y = y, x = x, obs.locs = obs.locs, family = family, nugget = nugget, 
                  matern.kappa = matern.kappa, effort = sample.effort, nrep = floor(nrep/10), seed = seed, method = "L-BFGS-B", 
              lower = c(rep(-Inf,n.reg0), rep(MIN, n.od0), rep(MIN, n.range0), rep(0, n.nugget0)),
              upper = c(rep(Inf,n.reg0), rep(Inf, n.od0), rep(max.range, n.range0), rep(1-MIN, n.nugget0)))
  
  fit <- optim(par = fittmp$par, fn = likGHK, y = y, x = x, obs.locs = obs.locs, family = family, nugget = nugget, 
               matern.kappa = matern.kappa, effort = sample.effort, nrep = nrep, seed = seed, method = "L-BFGS-B",
                 lower = c(rep(-Inf,n.reg0), rep(MIN, n.od0), rep(MIN, n.range0), rep(0, n.nugget0)),
                 upper = c(rep(Inf,n.reg0), rep(Inf, n.od0),rep(max.range, n.range0), rep(1-MIN, n.nugget0)))
                
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

#### Some Test Runs
#mleGHK( y = sim.y, obs.locs = cbind(xloc,yloc), family = 'NB', nrep = 1000, seed = 1234)
predGHK(y = sim.y, obs.locs = cbind(xloc,yloc), pred.locs = matrix(c(0.9, 0.8, 0.5,0.5),ncol = 2),
         family = 'NB', nugget = FALSE, matern.kappa = 0.5, nrep = 1000, seed = 123, max.range = 4, max.count = 20)
####

predGHK(y = simtmp.y, x = cbind(xloc,yloc), obs.locs = cbind(xloc, yloc), 
        pred.x = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),
      pred.locs = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),family = 'NB', 
     nugget = T, max.range = 2, max.count = 15)

mleGHK(y = simtmp.y, x = cbind(xloc,yloc), obs.locs = cbind(xloc, yloc),
       family = 'NB', nugget = TRUE, matern.kappa = 0.5, nrep = 1000, seed = 1234, max.range = NULL)

#######################################################################################################################
#######################################################################################################################
#### OUTPUT FUNCTION: Prediction using GHK simulator: serial version
#######################################################################################################################
#######################################################################################################################

#### Output: n(number of prediction locations)*2 matrix
#### First Column: predicting value
#### Second Column: Estimated MSPE

predGHK <- function(y, x = NULL, pred.x = NULL, obs.locs, pred.locs, family, nugget = FALSE, 
                    matern.kappa = 0.5, sample.effort = 1,  pred.effort = 1, binomial.n = 1, 
                    binomial.pred.n = 1, nrep = 1000, seed = 123, max.range = NULL, max.count = NULL)
{
  #  pred.locs = matrix(c(0.4,0.4,0.5,0.5),2,2,byrow = T)
  #  y = sim.nbnew1[1,1:200]; obs.locs = cbind(xloc,yloc)
  #  family = 'NB'; nugget = FALSE; matern.kappa = 0.5; nrep = 1000;
  #  seed = 123; max.range = 4; max.count = 20
  #  y = simtmp.y; x = cbind(xloc,yloc); pred.x = matrix(c(0.2, 0.3), nrow = 1)
  #  obs.locs = cbind(xloc, yloc); pred.locs = matrix(c(0.2,0.3),nrow = 1);family = 'ZIP' 
  #  nugget = F; max.range = 2; max.count = 15
  
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

  if(length(sample.effort) == 1) sample.effort <- rep(sample.effort, nrow(obs.locs))
  
  if(!length(sample.effort) == nrow(obs.locs)) 
    stop("Sampling Effort must be equal to the number of sampling locations!")
  
  if(length(pred.effort) == 1) pred.effort <- rep(pred.effort, nrow(pred.locs))
  
  if(!length(pred.effort) == nrow(pred.locs)) 
    stop("Prediction Effort must be equal to the number of prediction locations!")
  
  indexloc <- which.min(FNN::get.knnx(pred.locs, obs.locs, 1)$nn.dist)
  
  if(nrow(pred.locs) == 1){
    m0 <- n0 <- round(pred.effort*y[indexloc]/sample.effort[indexloc])+ 1
  }else{
    m0 <- n0 <- round(apply(pred.locs, 1, function(x) y[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]*pred.effort/
                        sample.effort[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]))+1
  }
  
  NPL <- length(m0)  
  
  #### m0 and n0: initial prediction values for probability search. Scalor or Vector
  
  if(is.null(max.count)) max.count <- ceiling(50*max(y))
  
  #### First calculate the MLE and output the log-likelihood as denominator
  
  MLE.est <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, nugget = nugget, 
                    matern.kappa = matern.kappa, binomial.n = binomial.n, 
                    sample.effort = sample.effort, nrep = nrep, seed = seed, max.range = max.range)
  
  loglik <- MLE.est$log.lik; estpar <- MLE.est$MLE
  
  if(is.null(x))  x <- rep(1, length(y))
  if(is.null(pred.x)) pred.x <- matrix(1, nrow = NPL, ncol = 1)
  if(!is.matrix(pred.x) & !is.data.frame(pred.x))
    stop("Input pred.x must be a data frame or matrix")
  
  if(nrow(pred.x)!= nrow(pred.locs))
    stop("Number of prediction locations did not match the number of covariates") 
  
  #### A for loop which cannot be avoided
  
  ans <- matrix(NA, nrow = NPL, ncol = 2)
  
  for(j in 1:NPL){  
    p.m0 <- p.n0 <- exp(-likGHK(param = estpar, y = c(y, m0[j]), x = rbind(as.matrix(x), pred.x[j,]), 
                 obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget, matern.kappa = matern.kappa,
                 effort = c(sample.effort, pred.effort[j]),nrep = nrep, seed = seed) - loglik)
    
    mu.m0 <- mu.n0 <- p.m0*m0[j];  mu2.m0 <- mu2.n0 <- p.m0*m0[j]^2
    MM1 <- matrix(0, nrow = 2, ncol = 3); MM2 <- matrix(0, nrow = 2, ncol = 3)
    MM1[1,] <- c(p.m0, mu.m0, mu2.m0);  MM2[1,] <- c(p.n0, mu.n0, mu2.n0)
    
  # To avoid 0 length, repeated here
    p.m0 <- exp(-likGHK(param = estpar, y = c(y, m0[j]-1), x = rbind(as.matrix(x), pred.x[j,]), 
                 obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget, matern.kappa = matern.kappa,
                 effort = c(sample.effort, pred.effort[j]), nrep = nrep, seed = seed) - loglik)
    
    mu.m0 <- p.m0*(m0[j]-1); mu2.m0 <- p.m0*(m0[j]-1)^2; MM1[2,] <- c(p.m0, mu.m0, mu2.m0)
    
    while( (p.m0 > sqrt(.Machine$double.eps) | MM1[nrow(MM1), 2] > MM1[nrow(MM1)-1, 2]) & m0[j] > 1)
    {
      p.m0 <- exp(-likGHK(param = estpar, y = c(y, m0[j]-2), x = rbind(as.matrix(x), pred.x[j,]), 
                          matern.kappa = matern.kappa, obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, 
                          nugget = nugget, effort = c(sample.effort, pred.effort[j]), nrep = nrep, seed = seed) - loglik)
      mu.m0 <- p.m0*(m0[j]-2); mu2.m0 <- p.m0*(m0[j]-2)^2
      MM1 <- rbind(MM1, c(p.m0, mu.m0, mu2.m0));  m0[j] <- m0[j]-1
    } 
    
    #### Search from n0 to the right
    
    p.n0 <- exp(-likGHK(param = estpar, y = c(y, n0[j]+1), x = rbind(as.matrix(x), pred.x[j,]), 
                         obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                         matern.kappa = matern.kappa, effort = c(sample.effort, pred.effort[j]), 
                        nrep = nrep, seed = seed) - loglik)
    
    mu.n0 <- p.n0*(n0[j]+1);  mu2.n0 <- p.n0*(n0[j]+1)^2
    MM2[2, ] <- c(p.n0, mu.n0, mu2.n0)
    
    while( (p.n0 > sqrt(.Machine$double.eps) | MM2[nrow(MM2), 2] > MM2[nrow(MM2)-1, 2]) & n0[j] < max.count)
    {
      p.n0 <- exp(-likGHK(param = estpar, y = c(y, n0[j]+2), x = rbind(as.matrix(x), pred.x[j,]), 
                           obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                           matern.kappa = matern.kappa, effort = c(sample.effort, pred.effort[j]), 
                          nrep = nrep, seed = seed) - loglik)
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
#### OUTPUT FUNCTION: Prediction using GHK simulator: parallel version via snowfall.
#######################################################################################################################
#######################################################################################################################

#### Output: n(number of prediction locations)*2 matrix
#### First Column: predicting value
#### Second Column: Estimated MSPE

predGHK.sf <- function(y, x = NULL, pred.x = NULL, obs.locs, pred.locs, family, nugget = FALSE, 
                       matern.kappa = 0.5, sample.effort = 1,  pred.effort = 1, binomial.n = 1, 
                       binomial.pred.n = 1, nrep = 1000, seed = 123, max.range = NULL, 
                       max.count = NULL, n.cores = 2, cluster.type="SOCK")
{
  
  #    pred.locs = matrix(c(0,0.4,0.5,0.5),2,2,byrow = T);
  #    x = NULL; pred.x = NULL
  #    y = sim.y; obs.locs = cbind(xloc,yloc)
  #    family = 'ZIP'; nugget = FALSE; matern.kappa = 0.5; nrep = 100;
  #    seed = 123; max.range = 4; max.count = 20
  
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
  
  if(length(sample.effort) == 1) sample.effort <- rep(sample.effort, nrow(obs.locs))
  if(!length(sample.effort) == nrow(obs.locs)) 
    stop("Sampling Effort must be equal to the number of sampling locations!")
  
  if(length(pred.effort) == 1) pred.effort <- rep(pred.effort, nrow(pred.locs))
  if(!length(pred.effort) == nrow(pred.locs)) 
    stop("Prediction Effort must be equal to the number of prediction locations!")
  
  m0 <- n0 <- round(apply(pred.locs, 1, function(x) y[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]*pred.effort/
                            sample.effort[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]))+1

   NPL <- length(m0)  
  
  #### m0 and n0: initial prediction values for probability search. Scalor or Vector
  
  if(is.null(max.count)) max.count <- ceiling(50*max(y))
  
  MLE.est <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, nugget = nugget, 
                    matern.kappa = matern.kappa, binomial.n = binomial.n, 
                    sample.effort = sample.effort, nrep = nrep, seed = seed, max.range = max.range)
  
  loglik <- MLE.est$log.lik;  estpar <- MLE.est$MLE
  if(is.null(x))  x <- rep(1, length(y))
  if(is.null(pred.x)) pred.x <- matrix(1, nrow = NPL, ncol = 1)
  if(!is.matrix(pred.x) & !is.data.frame(pred.x))
    stop("Input pred.x must be a data frame or matrix")
  if(nrow(pred.x)!= nrow(pred.locs))
    stop("Number of prediction locations did not match number of covariates") 
  
  #### Begin to parallel
  
    if (requireNamespace("snowfall", quietly = TRUE)) {    
    snowfall::sfInit(parallel=TRUE, cpus= n.cores, type= cluster.type)
    snowfall::sfExportAll( except=NULL, debug=FALSE )
    snowfall::sfLibrary(Rcpp);snowfall::sfLibrary(RcppArmadillo)
    snowfall::sfLibrary(geoR);snowfall::sfLibrary(VGAM)
    
    #snowfall::sfLibrary("geoCount", character.only= TRUE)
    #  snowfall::sfClusterSetupRNG()
    #  sfClusterEvalQ( ls() )
    #  sfInit(parallel=T,cpus = 4,type="SOCK")
    #  sfLibrary(MASS);sfLibrary(Rcpp);sfLibrary(RcppArmadillo)
    #  sfExport(list=list("sim.ber1","GQTestber1","xloc","yloc","gqtllmle.ber1" ))
    
    par.pred.inner <- function(j){
      
      #### Q1: How to avoid function reloading and namespace problem. 
      #### Loading the own package directly may avoid this problem.
      
      ########################################################################################
      ########################################################################################
      ###################### FUNCTION RELOADING ############################################
      ########################################################################################
      
      cppFunction(depends = 'RcppArmadillo',   code = '
                  double ghkLik(arma::mat L, arma::vec lower, arma::vec upper, int nrep) {
                  
                  int dim = L.n_cols;double z_[dim];
                  double res_ = 0, prod, mu_, eta_, gamma_, ans_;
                  const double EPS = std::numeric_limits<double>::min();  
                  
                  double* eta = & eta_; double* res = & res_; double* ans = &ans_;
                  double* gamma = & gamma_; double* mu = & mu_; double* z = &z_[0];
                  
                  for(int i=0;i<nrep;i++) {
                  prod = 1.0;
                  for(int j=0;j<dim;j++) {     
                  *mu=0; 
                  for(int k=0;k<j;k++)  *mu += L(j,k)*(*(z+k));
                  *eta = Rcpp::stats::pnorm_1((lower(j)-*mu)/L(j,j),0,1,0);
                  *gamma = Rcpp::stats::pnorm_1((upper(j)-*mu)/L(j,j),0,1,0);
                  prod *= *gamma-*eta; 
                  double u = Rf_runif(0, 1);
                  *ans=u*(*gamma)+(1-u)*(*eta);
                  if(*ans > 1-EPS) *ans= 1-EPS;
                  if(*ans < EPS) *ans= EPS ;
                  *(z+j) = Rcpp::stats::qnorm_1(*ans,0,1,0);
                  }
                  *res += prod;
                  }
                  *res /= nrep;
                  if(*res == 0){
                  return -1000;
                  }
                  else{
                  return log(*res);
                  }
                  }')

      
      
      likGHK <- function(param, y, x = NULL, obs.locs, family, nugget = FALSE, matern.kappa = 0.5, 
                         effort = sample.effort, nrep = 1000, seed = 123)
      {
        nparam <- length(param)
        D <- as.matrix(dist(obs.locs, method = "euclidean", diag = T, upper = T))
        
        if(nugget == FALSE){
          R <- geoR::matern(D, param[nparam], matern.kappa)
        }else{
          R <- (1-param[nparam])*geoR::matern(D, param[nparam-1], matern.kappa) + param[nparam]*diag(ncol(D))
        }
        
        L <- try(t(chol(R)),silent=TRUE)
        if( inherits(L,"try-error") ) stop("Cholesky Decomposition failed in likelihood evaluation. 
                                           Possible ill-conditioned matrix involved.")
        
        if(is.null(x)){
          M <- exp(param[1])*effort;  S <- 1/param[2]
        }else{
          M <- exp(param[1] + param[2:(ncol(x)+1)]%*%t(x))*effort;  S <- 1/param[ncol(x)+2]
        }
        
        if (family == "Poisson")
          #### Even though the "size" parameter S is created, in this case, it is not used. 
        {
          a <- qnorm(ppois( y-1, lambda = M));  b <- qnorm(ppois( y, lambda = M))
        }
        if (family == "ZIP")
        {
          a <- qnorm(VGAM::pzipois(q = y-1, lambda = M+M/S, pstr0 = 1/(S+1)))
          b <- qnorm(VGAM::pzipois(q = y, lambda = M+M/S, pstr0 = 1/(S+1)))
        }
        if (family == "NB")
        {
          a <- qnorm(pnbinom( y-1, size = S, mu = M));  b <- qnorm(pnbinom( y, size = S, mu = M))
        }
        set.seed(seed)
        loglik <- ghkLik(L = L, lower = a, upper = b, nrep = nrep)
        if(is.nan(loglik)) loglik = -1e6
        if (loglik == Inf) loglik = 1e6
        return(-loglik)
      }
      
      
      
      ########################################################################################
      ########################################################################################
      ###################### FUNCTION RELOADING END ############################################
      ########################################################################################
      
      p.m0 <- p.n0 <- exp(-likGHK(param = estpar, y = c(y, m0[j]), x = rbind(as.matrix(x), pred.x[j,]), 
                                   obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                                   matern.kappa = matern.kappa, effort = c(sample.effort, pred.effort[j]), 
                                  nrep = nrep, seed = seed) - loglik)
      
      mu.m0 <- mu.n0 <- p.m0*m0[j];  mu2.m0 <- mu2.n0 <- p.m0*m0[j]^2
      
      MM1 <- matrix(0, nrow = 2, ncol = 3); MM2 <- matrix(0, nrow = 2, ncol = 3)
      MM1[1,] <- c(p.m0, mu.m0, mu2.m0);  MM2[1,] <- c(p.n0, mu.n0, mu2.n0)
      
      p.m0 <- exp(-likGHK(param = estpar, y = c(y, m0[j] -1), x = rbind(as.matrix(x), pred.x[j,]), 
                           obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                           matern.kappa = matern.kappa, effort = c(sample.effort, pred.effort[j]), 
                          nrep = nrep, seed = seed) - loglik)
      
      mu.m0 <- p.m0*(m0[j]-1); mu2.m0 <- p.m0*(m0[j]-1)^2; MM1[2,] <- c(p.m0, mu.m0, mu2.m0)
      
      while( (p.m0 > sqrt(.Machine$double.eps) |  MM1[nrow(MM1), 2] > MM1[nrow(MM1)-1, 2]) & m0[j] > 1)
      {
        p.m0 <- exp(-likGHK(param = estpar, y = c(y, m0[j]-2), x = rbind(as.matrix(x), pred.x[j,]), 
                             obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                             matern.kappa = matern.kappa, effort = c(sample.effort, pred.effort[j]), 
                            nrep = nrep, seed = seed) - loglik)
        
        mu.m0 <- p.m0*(m0[j]-2);   mu2.m0 <- p.m0*(m0[j]-2)^2
        MM1 <- rbind(MM1, c(p.m0, mu.m0, mu2.m0)); m0[j] <- m0[j]-1
      } 
      
      p.n0 <- exp(-likGHK(param = estpar, y = c(y, n0[j]+1),x = rbind(as.matrix(x), pred.x[j,]), 
                           obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                           matern.kappa = matern.kappa, effort = c(sample.effort, pred.effort[j]), 
                          nrep = nrep, seed = seed) - loglik)
      
      mu.n0 <- p.n0*(n0[j]+1);  mu2.n0 <- p.n0*(n0[j]+1)^2
      MM2[2, ] <- c(p.n0, mu.n0, mu2.n0)
      
      while( (p.n0 > sqrt(.Machine$double.eps) | MM2[nrow(MM2), 2] > MM2[nrow(MM2)-1, 2]) & n0[j] < max.count)
      {
        p.n0 <- exp(-likGHK(param = estpar, y = c(y, n0[j]+2), x = rbind(as.matrix(x), pred.x[j,]), 
                             obs.locs = rbind(obs.locs, pred.locs[j,]), family = family, nugget = nugget,
                             matern.kappa = matern.kappa, effort = c(sample.effort, pred.effort[j]), 
                            nrep = nrep, seed = seed) - loglik)
        
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


predGHK.sf(y = simtmp.y, x = cbind(xloc,yloc), obs.locs = cbind(xloc, yloc), 
           pred.x = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),
           pred.locs = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),family = 'NB', 
           nugget = T, max.range = 2, max.count = 15, n.cores = 4, cluster.type = "SOCK")
