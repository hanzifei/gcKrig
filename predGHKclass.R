library(Rcpp)
library(RcppArmadillo)
library(geoR)
library(VGAM)
library(MASS)
library(snowfall)
library(FNN)
library(mvtnorm)

#### Likelihood Expression 
#### mu: vector of mean (in the truncated gaussian)
#### R: correlation matrix
#### lower: lower bound of the integration
#### upper: upper bound of the intergration
#### nrep: Monte Carlo size
#### reorder: if TRUE, using reorder algorithm; otherwise not. 
#######################################################################################################################
#######################################################################################################################
#### OUTPUT FUNCTION: Likelihood Evaluation (multivariate normal rectange probability) ####
#######################################################################################################################
#######################################################################################################################

cppFunction(depends = 'RcppArmadillo',   code = '
           double ghkLiko(arma::vec mu, arma::mat R, arma::vec lower, arma::vec upper, int nrep){
            const double eps = std::numeric_limits<double>::epsilon();  
            const double EPS = std::numeric_limits<double>::min(); 
            int n = R.n_cols;
            arma::mat L = R; arma::vec ap = lower; arma::vec bp = upper; 
            arma::vec y;  y.zeros(n);
            double ctmp_, vtmp_, dem_, am_, bm_, de_, ai_, bi_, s_; 
            double res_ = 0, prod, MU_, eta_, gamma_, ans_;
            double z_[n];
            double* ctmp = &ctmp_;  double* vtmp = &vtmp_; double* dem = &dem_;
            double* am = &am_; double* bm = &bm_; double* de = &de_; 
            double* ai = &ai_; double* bi = &bi_; double* s = &s_;
            double* eta = & eta_; double* res = & res_; double* ans = &ans_;
            double* gamma = & gamma_; double* MU = & MU_; double* z = &z_[0];
            int im; arma::mat tmp;  arma::mat mattmp;
            arma::vec d = arma::sqrt(arma::diagvec(L));
            for(int i = 0; i<n; i++){
            if(d(i) > 0){
            L(arma::span::all, i) = L(arma::span::all, i)/d(i); 
            L(i, arma::span::all) = L(i, arma::span::all)/d(i);
            ap(i) = ap(i)/d(i); bp(i) = bp(i)/d(i);
            }
            }
            for(int k=0; k<n; k++){
            im = k;  *vtmp = 0;  *dem = 1.0; *s = 0;
            for(int i=k; i<n; i++){
            if(L(i,i) > eps){
            *ctmp = sqrt(R::fmax2(L(i,i), 0 ));
            if(k == 0){*s = 0;}
            if(i > 0 && k > 0){ 
            mattmp = L(i, arma::span(0,k-1))*y(arma::span(0,k-1)); *s = mattmp(0,0);
            }
            *ai = (ap(i)-*s)/(*ctmp);  *bi = (bp(i)-*s)/(*ctmp); 
            *de = R::pnorm(*bi,0,1,1,0) - R::pnorm(*ai,0,1,1,0);
            if(*de<=*dem){ *vtmp = *ctmp;  *dem = *de;  *am = *ai;  *bm = *bi; im = i;}
            }
            }
            if(im > k){
            ap.swap_rows(im,k); bp.swap_rows(im,k);
            L(im, im) = L(k,k);
            if(k > 0){
            L(arma::span::all,arma::span(0,k-1)).swap_rows(im,k);
            }
            if( n-im >= 2){
            L(arma::span(im+1,n-1),arma::span::all).swap_cols(im,k);
            }
            if(im - k >= 2){ 
            tmp = L(arma::span(k+1,im-1),k); L(arma::span(k+1,im-1),k) = L(im,arma::span(k+1,im-1)).t();
            L(im,arma::span(k+1,im-1)) = tmp.t();
            }
            }
            if(n - k >= 2)  L(k,arma::span(k+1,n-1)).zeros();
            if(*vtmp > eps*k){
            L(k,k) = *vtmp; 
            for(int i=k+1; i<n; i++){
            L(i,k) = L(i,k)/(*vtmp);
            L(i, arma::span(k+1,i)) = L(i, arma::span(k+1,i))-L(i,k)*L(arma::span(k+1,i),k).t();
            }
            if(std::abs(*dem) > eps){
            y(k) = (R::dnorm(*am,0,1,0) - R::dnorm(*bm,0,1,0) )/ (*dem);
            }else{
            if(*am < -10) {y(k) = *bm;}
            else if(*bm > 10) {y(k) = *am;}
            else{y(k) = (*am+*bm)/2;}
            }
            }
            else{L(arma::span(k,n-1), k).zeros(); y(k) = 0; }
            }
            for(int ii=0; ii<nrep; ii++) {
            prod = 1.0;
            for(int jj=0; jj<n; jj++) {    
            *MU=0; 
            for(int kk=0; kk<jj; kk++)  *MU += L(jj,kk)*(*(z+kk));
            *eta = Rcpp::stats::pnorm_1((ap(jj)-*MU-mu(jj))/L(jj,jj),0,1,0);
            *gamma = Rcpp::stats::pnorm_1((bp(jj)-*MU-mu(jj))/L(jj,jj),0,1,0);
            prod *= *gamma-*eta; 
            double u = Rf_runif(0, 1);
            *ans=u*(*gamma)+(1-u)*(*eta);
            if(*ans > 1-EPS) *ans = 1-EPS;
            if(*ans < EPS) *ans = EPS ;
            *(z+jj) = Rcpp::stats::qnorm_1(*ans,0,1,0);
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


cppFunction(depends = 'RcppArmadillo',   code = '
            double ghkLik(arma::vec mu, arma::mat R, 
            arma::vec lower, arma::vec upper, int nrep) {
            int dim = R.n_cols; double z_[dim];
            arma::mat L = arma::chol(R).t();
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
#### ReOrdered Version of the GHK simulator

mvnintGHK <- function(mu, sigma, lower, upper, nrep = 1000, log = T, reorder = T){
  
  if(!is.matrix(sigma)) 
    stop("Input 'sigma' must be of form matrix!")
  if(!all(lower <= upper)) 
    stop("Elements in 'lower' must be <=  the corresponding elements in 'upper'!")
  if(!isSymmetric(sigma)) 
    stop("Input covariance matrix 'sigma' must be symmetric!")
  if( inherits(try(chol(sigma),silent=TRUE),"try-error") ) 
    stop("Cholesky Decomposition failed. Input matrix sigma is not a valid covariance matrix!")
  if(!all.equal(length(mu), nrow(sigma), length(lower), length(upper)))
    stop("Input 'mu', lower' and 'upper' must have same length as dimension of the sigma!")
  
  if(reorder == T){
    ans <- ghkLiko(mu = mu, R = sigma, lower = lower, upper = upper, nrep = nrep)
  }else{
    ans <- ghkLik(mu = mu, R = sigma, lower = lower, upper = upper, nrep = nrep)
  }
  if(log == F) ans <- exp(ans)
  return(ans)
}


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

likGHK <- function(pars, y, x = NULL, locs, family, corr, effort, nrep = 1000, seed = 123)
{
  matD <- as.matrix(dist(locs, method = "euclidean", diag = T, upper = T))
  npars <- length(pars)
  nparcorr <- corr$npar.cor
  R <- corr$corr(pars[(npars-nparcorr+1):npars], matD)
  bounds <- family$bounds(y = y, x = x, pars = pars, effort = effort)
  set.seed(seed)
  loglik <- ghkLiko(mu = rep(0, nrow(R)), R = R, lower = bounds$lower, upper = bounds$upper, nrep = nrep)
  if(is.nan(loglik)) loglik = -1e6;  if (loglik == Inf) loglik = 1e6
  return(-loglik)
}


unix.time(
likGHK(pars = c(0.7,1.1,0.2,0.3,0.5,0.2), y = simtmp.y, x = cbind(1,xloc,yloc), locs = cbind(xloc,yloc),
       family = negbin.gcgc2(link = 'log'), 
       corr = matern.gcgc2(matern.kappa = 0.5, nugget = T), effort = 1,
       nrep = 1000, seed = 123)
)


mleGHK(y = simtmp.y, x = cbind(1, xloc,yloc), obs.locs = cbind(xloc,yloc), 
       family = poisson.gcgc2(link = 'log'),
       corr = matern.gcgc2( matern.kappa = 0.5, nugget = T), 
       effort = 1,  nrep = c(100,1000), seed = 123, max.range = 100)

mleGQT(y = simtmp.y, x = cbind(1, xloc,yloc), obs.locs = cbind(xloc,yloc), 
       family = negbin.gcgc2(link = 'log'),
       corr = spherical.gcgc2(nugget = F), 
       effort = 1, max.range = 100)

mleGHK(y = simtmp.y, x = cbind(1, xloc,yloc), obs.locs = cbind(xloc,yloc), 
       family = poisson.gcgc2(link = 'log'),
       corr = spherical.gcgc2( nugget = F), 
       effort = 1,  nrep = c(100,1000), seed = 1234, max.range = 100)


#######################################################################################################################
#######################################################################################################################
#### Inner function only used to be called by another functions inside package
#######################################################################################################################
#######################################################################################################################
#### Now Consider Add optimization procedure
#### Output: MLE and Value of the simulated log likelihood evaluated at the MLEs

mleGHK <- function(y, x = NULL, obs.locs, family, corr, effort = 1, nrep = c(100, 1000), 
                   seed = 123, max.range = NULL)
{ 
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
    corpar0 <- corr$start(matD)
  }else{
    corpar0 <- corr0$par
  }
  
  est <- c(family$start(y = y, x = x, effort = effort)$start, corpar0)
  n.reg0 <- ncol(x); n.od0 <- family$nod 
  
  for(i in 1:length(nrep)){
    fit <- optim(par = est, fn = likGHK, y = y, x = x, locs = obs.locs, family = family, corr = corr, 
                    effort = effort, nrep = nrep[i], seed = seed, method = "L-BFGS-B", 
                    lower = c(rep(-Inf,n.reg0), rep(MIN, n.od0), MIN, rep(0, n.nugget0)),
                    upper = c(rep(Inf,n.reg0), rep(Inf, n.od0), max.range, rep(1-MIN, n.nugget0)))
    est <- fit$par
  }
  if (est[n.reg0+n.od0+1] == max.range)
    warning("Maximum Range Value Reached, MLE may not exist. Try to input a larger One")
  
  if (is.null(fit$convergence) || fit$convergence != 0)   
    warnings("Algorithm might converge bad in MLE") 
  result.list <- list()
  result.list$MLE <- est; result.list$log.lik <- -fit$value
  k <- length(est); N <- length(y)
  result.list$AIC <- 2*k+2*fit$value
  result.list$AICc <- 2*k+2*fit$value + 2*k*(k+1)/(N-k-1)
  result.list$BIC <- 2*fit$value + k*log(N)
  # class(result.list) <- c( "pred.ghk")
  return(result.list)
}




#### Some Test Runs
#mleGHK( y = sim.y, obs.locs = cbind(xloc,yloc), family = 'NB', nrep = 1000, seed = 1234)


predGHK(y = simtmp.y, x = cbind(xloc,yloc), obs.locs = cbind(xloc, yloc), 
        pred.locs = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),
        family = negbin.gcgc2(link = 'log'), 
        corr = matern.gcgc2( matern.kappa = 0.5, nugget = F), 
        pred.x = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),
         max.range = 5, max.count = 15)


#######################################################################################################################
#######################################################################################################################
#### OUTPUT FUNCTION: Prediction using GHK simulator: serial version
#######################################################################################################################
#######################################################################################################################

#### Output: n(number of prediction locations)*2 matrix
#### First Column: predicting value
#### Second Column: Estimated MSPE


predGHK(y = simtmp.y, x = cbind(xloc,yloc), obs.locs = cbind(xloc, yloc), 
        pred.locs = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),
        family = binomial_t.gcgc2(df.t = 5), 
        corr = matern.gcgc2(matern.kappa = 0.5, nugget = T), 
        sample.effort = rep(c(11,19),50), pred.effort = c(10,12,12,19), 
        pred.x = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),
        nrep = c(100,1000),seed = 124, max.range = 5, max.count = 15)


predGHK <- function(y, x = NULL, pred.x = NULL, obs.locs, pred.locs, family, corr, sample.effort = 1, 
                    pred.effort = 1, nrep = c(100, 1000), seed = 123, max.range = NULL, max.count = NULL)
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
  
  MLE.est <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                    effort = sample.effort, nrep = nrep, seed = seed, max.range = max.range)
  
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
  #    m0 <- n0 <- round(apply(pred.locs, 1, function(x) 
  #      y[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]*pred.effort/
  #        sample.effort[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]))+1
  
  NPL <- length(m0)
  #### m0 and n0: initial prediction values for probability search. Scalor or Vector
  if(is.null(max.count)) max.count <- ceiling(50*max(y))
  
  #### A for loop which cannot be avoided
  
  ans <- matrix(NA, nrow = NPL, ncol = 2); nnrep <- length(nrep)
  for(j in 1:NPL){  
    tmpfun <- function(xtmp) {
      exp(-likGHK(pars = estpar, y = c(y, xtmp), x = rbind(as.matrix(x), pred.x[j,]), 
                  locs = rbind(obs.locs, pred.locs[j,]), family = family, corr = corr,
                  effort = c(sample.effort, pred.effort[j]),nrep = nrep[nnrep], seed = seed) - loglik)
    }
    
    p.m0 <- p.n0 <- tmpfun(m0[j]); mu.m0 <- mu.n0 <- p.m0*m0[j]; mu2.m0 <- mu2.n0 <- p.m0*m0[j]^2
    MM1 <- matrix(0, nrow = 2, ncol = 3); MM2 <- matrix(0, nrow = 2, ncol = 3)
    MM1[1,] <- c(p.m0, mu.m0, mu2.m0);  MM2[1,] <- c(p.n0, mu.n0, mu2.n0)
    
    # To avoid 0 length, repeat here
    p.m0 <- tmpfun(m0[j]-1); mu.m0 <- p.m0*(m0[j]-1)
    mu2.m0 <- p.m0*(m0[j]-1)^2; MM1[2,] <- c(p.m0, mu.m0, mu2.m0)
    
    while( (p.m0 > sqrt(.Machine$double.eps) | MM1[nrow(MM1), 2] > MM1[nrow(MM1)-1, 2]) & m0[j] > 1)
    {
      p.m0 <- tmpfun(m0[j]-2); mu.m0 <- p.m0*(m0[j]-2); mu2.m0 <- p.m0*(m0[j]-2)^2
      MM1 <- rbind(MM1, c(p.m0, mu.m0, mu2.m0));  m0[j] <- m0[j]-1
    } 
    #### Search from n0 to the right
    
    p.n0 <- tmpfun(n0[j]+1); mu.n0 <- p.n0*(n0[j]+1);  mu2.n0 <- p.n0*(n0[j]+1)^2
    MM2[2, ] <- c(p.n0, mu.n0, mu2.n0)
    
    while( (p.n0 > sqrt(.Machine$double.eps) | MM2[nrow(MM2), 2] > MM2[nrow(MM2)-1, 2]) & n0[j] < max.count)
    {
      p.n0 <- tmpfun(n0[j]+2); mu.n0 <- p.n0*(n0[j]+2); mu2.n0 <- p.n0*(n0[j]+2)^2
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

predGHK.sf <- function(y, x = NULL, pred.x = NULL, obs.locs, pred.locs, family, corr, 
                       sample.effort = 1, pred.effort = 1, nrep = 1000, seed = 123, max.range = NULL, 
                       max.count = NULL, n.cores = 2, cluster.type="SOCK")
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
  
  MLE.est <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                    effort = sample.effort, nrep = nrep, seed = seed, max.range = max.range)
  
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
  #    m0 <- n0 <- round(apply(pred.locs, 1, function(x) 
  #      y[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]*pred.effort/
  #        sample.effort[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]))+1
  
  NPL <- length(m0)
  if(is.null(max.count)) max.count <- ceiling(50*max(y));nnrep <- length(nrep)
#  ans <- matrix(NA, nrow = NPL, ncol = 2); 
  
  if (requireNamespace("snowfall", quietly = TRUE)) {    
    snowfall::sfInit(parallel =TRUE, cpus = n.cores, type = cluster.type)
    suppressMessages(snowfall::sfExportAll(except = NULL, debug = FALSE))
    suppressMessages(snowfall::sfLibrary(Rcpp))
    suppressMessages(snowfall::sfLibrary(RcppArmadillo))
    suppressMessages(snowfall::sfLibrary(geoR))
    suppressMessages(snowfall::sfLibrary(VGAM))
    
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
                double ghkLiko(arma::vec mu, arma::mat R, arma::vec lower, arma::vec upper, int nrep){
                const double eps = std::numeric_limits<double>::epsilon();  
                const double EPS = std::numeric_limits<double>::min(); 
                int n = R.n_cols;
                arma::mat L = R; arma::vec ap = lower; arma::vec bp = upper; 
                arma::vec y;  y.zeros(n);
                double ctmp_, vtmp_, dem_, am_, bm_, de_, ai_, bi_, s_; 
                double res_ = 0, prod, MU_, eta_, gamma_, ans_;
                double z_[n];
                double* ctmp = &ctmp_;  double* vtmp = &vtmp_; double* dem = &dem_;
                double* am = &am_; double* bm = &bm_; double* de = &de_; 
                double* ai = &ai_; double* bi = &bi_; double* s = &s_;
                double* eta = & eta_; double* res = & res_; double* ans = &ans_;
                double* gamma = & gamma_; double* MU = & MU_; double* z = &z_[0];
                int im; arma::mat tmp;  arma::mat mattmp;
                arma::vec d = arma::sqrt(arma::diagvec(L));
                for(int i = 0; i<n; i++){
                if(d(i) > 0){
                L(arma::span::all, i) = L(arma::span::all, i)/d(i); 
                L(i, arma::span::all) = L(i, arma::span::all)/d(i);
                ap(i) = ap(i)/d(i); bp(i) = bp(i)/d(i);
                }
                }
                for(int k=0; k<n; k++){
                im = k;  *vtmp = 0;  *dem = 1.0; *s = 0;
                for(int i=k; i<n; i++){
                if(L(i,i) > eps){
                *ctmp = sqrt(R::fmax2(L(i,i), 0 ));
                if(k == 0){*s = 0;}
                if(i > 0 && k > 0){ 
                mattmp = L(i, arma::span(0,k-1))*y(arma::span(0,k-1)); *s = mattmp(0,0);
                }
                *ai = (ap(i)-*s)/(*ctmp);  *bi = (bp(i)-*s)/(*ctmp); 
                *de = R::pnorm(*bi,0,1,1,0) - R::pnorm(*ai,0,1,1,0);
                if(*de<=*dem){ *vtmp = *ctmp;  *dem = *de;  *am = *ai;  *bm = *bi; im = i;}
                }
                }
                if(im > k){
                ap.swap_rows(im,k); bp.swap_rows(im,k);
                L(im, im) = L(k,k);
                if(k > 0){
                L(arma::span::all,arma::span(0,k-1)).swap_rows(im,k);
                }
                if( n-im >= 2){
                L(arma::span(im+1,n-1),arma::span::all).swap_cols(im,k);
                }
                if(im - k >= 2){ 
                tmp = L(arma::span(k+1,im-1),k); L(arma::span(k+1,im-1),k) = L(im,arma::span(k+1,im-1)).t();
                L(im,arma::span(k+1,im-1)) = tmp.t();
                }
                }
                if(n - k >= 2)  L(k,arma::span(k+1,n-1)).zeros();
                if(*vtmp > eps*k){
                L(k,k) = *vtmp; 
                for(int i=k+1; i<n; i++){
                L(i,k) = L(i,k)/(*vtmp);
                L(i, arma::span(k+1,i)) = L(i, arma::span(k+1,i))-L(i,k)*L(arma::span(k+1,i),k).t();
                }
                if(std::abs(*dem) > eps){
                y(k) = (R::dnorm(*am,0,1,0) - R::dnorm(*bm,0,1,0) )/ (*dem);
                }else{
                if(*am < -10) {y(k) = *bm;}
                else if(*bm > 10) {y(k) = *am;}
                else{y(k) = (*am+*bm)/2;}
                }
                }
                else{L(arma::span(k,n-1), k).zeros(); y(k) = 0; }
                }
                for(int ii=0; ii<nrep; ii++) {
                prod = 1.0;
                for(int jj=0; jj<n; jj++) {    
                *MU=0; 
                for(int kk=0; kk<jj; kk++)  *MU += L(jj,kk)*(*(z+kk));
                *eta = Rcpp::stats::pnorm_1((ap(jj)-*MU-mu(jj))/L(jj,jj),0,1,0);
                *gamma = Rcpp::stats::pnorm_1((bp(jj)-*MU-mu(jj))/L(jj,jj),0,1,0);
                prod *= *gamma-*eta; 
                double u = Rf_runif(0, 1);
                *ans=u*(*gamma)+(1-u)*(*eta);
                if(*ans > 1-EPS) *ans = 1-EPS;
                if(*ans < EPS) *ans = EPS ;
                *(z+jj) = Rcpp::stats::qnorm_1(*ans,0,1,0);
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

      likGHK <- function(pars, y, x = NULL, locs, family, corr, effort, nrep = 1000, seed = 123)
     {
        matD <- as.matrix(dist(locs, method = "euclidean", diag = T, upper = T))
        npars <- length(pars)
        nparcorr <- corr$npar.cor
        R <- corr$corr(pars[(npars-nparcorr+1):npars], matD)
        bounds <- family$bounds(y = y, x = x, pars = pars, effort = effort)
        set.seed(seed)
        loglik <- ghkLiko(mu = rep(0, nrow(R)), R = R, lower = bounds$lower, upper = bounds$upper, nrep = nrep)
        if(is.nan(loglik)) loglik = -1e6;  if (loglik == Inf) loglik = 1e6
        return(-loglik)
      }
      
      ########################################################################################
      ########################################################################################
      ###################### FUNCTION RELOADING END ############################################
      ########################################################################################
      
      tmpfun <- function(xtmp) {
        exp(-likGHK(pars = estpar, y = c(y, xtmp), x = rbind(as.matrix(x), pred.x[j,]), 
                    locs = rbind(obs.locs, pred.locs[j,]), family = family, corr = corr,
                    effort = c(sample.effort, pred.effort[j]),nrep = nrep[nnrep], seed = seed) - loglik)
      }
      
      p.m0 <- p.n0 <- tmpfun(m0[j]); mu.m0 <- mu.n0 <- p.m0*m0[j];  mu2.m0 <- mu2.n0 <- p.m0*m0[j]^2
      MM1 <- matrix(0, nrow = 2, ncol = 3); MM2 <- matrix(0, nrow = 2, ncol = 3)
      MM1[1,] <- c(p.m0, mu.m0, mu2.m0);  MM2[1,] <- c(p.n0, mu.n0, mu2.n0)
      
      p.m0 <- tmpfun(m0[j]-1);mu.m0 <- p.m0*(m0[j]-1); mu2.m0 <- p.m0*(m0[j]-1)^2; MM1[2,] <- c(p.m0, mu.m0, mu2.m0)
      
      while( (p.m0 > sqrt(.Machine$double.eps) |  MM1[nrow(MM1), 2] > MM1[nrow(MM1)-1, 2]) & m0[j] > 1)
      {
        p.m0 <- tmpfun(m0[j]-2); mu.m0 <- p.m0*(m0[j]-2);  mu2.m0 <- p.m0*(m0[j]-2)^2
        MM1 <- rbind(MM1, c(p.m0, mu.m0, mu2.m0)); m0[j] <- m0[j]-1
      } 
      
      p.n0 <- tmpfun(n0[j]+1); mu.n0 <- p.n0*(n0[j]+1);  mu2.n0 <- p.n0*(n0[j]+1)^2
      MM2[2, ] <- c(p.n0, mu.n0, mu2.n0)
      
      while( (p.n0 > sqrt(.Machine$double.eps) | MM2[nrow(MM2), 2] > MM2[nrow(MM2)-1, 2]) & n0[j] < max.count)
      {
        p.n0 <- tmpfun(n0[j]+2);  mu.n0 <- p.n0*(n0[j]+2); mu2.n0 <- p.n0*(n0[j]+2)^2
        MM2 <- rbind(MM2, c(p.n0, mu.n0, mu2.n0));  n0[j] <- n0[j]+1
      }
      MM2 <- MM2[-1, ];  MM.all <- rbind(MM1, MM2)
      #### Due to some Monte Carlo Error, put sum to 1 constraint
      
      weight <- 1/sum(MM.all[,1])
      return( c(sum(MM.all[,2])*weight, sum(MM.all[,3]*weight)-(sum(MM.all[,2])*weight)^2) )
    }
    out = sfLapply(1:NPL, par.pred.inner)
    sfStop()
    message("Parallel prediction completed with  ")
    }
  return(out)
}

predGHK.sf(y = simtmp.y, x = cbind(xloc,yloc), obs.locs = cbind(xloc, yloc), 
        pred.locs = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),
        family = binomial_t.gcgc2(df.t = 5), 
        corr = matern.gcgc2(matern.kappa = 0.5, nugget = T), 
        sample.effort = rep(c(11,19),50), pred.effort = c(10,12,12,19), 
        pred.x = matrix(c(0.4,0.4,0.5,0.5,0,0.1,0.9,0.8),4,2,byrow = T),
        nrep = c(100,1000),seed = 124, max.range = 5, max.count = 15, n.cores = 4, cluster.type="SOCK")