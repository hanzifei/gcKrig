selectmodel.gcgc <- function(y, x, obs.locs, corr, effort, max.range, seed, nrep){
  
  x = cbind(1,x); n <- nrow(x)
  matD <- as.matrix(dist(obs.locs, method = "euclidean", diag = T, upper = T))
  
 if(corr == "matern"){
    corr <- matern.gcgc2(nugget = F)
    family <- poisson.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic1 <- llik$AIC; aicc1 <- llik$AICc; bic1 <- llik$BIC
    
    family <- zip.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic2 <- llik$AIC; aicc2 <- llik$AICc; bic2 <- llik$BIC
    
    family <- negbin.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic3 <- llik$AIC; aicc3 <- llik$AICc; bic3 <- llik$BIC
    
    corr <- matern.gcgc2(nugget = T)
    family <- poisson.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic4 <- llik$AIC; aicc4 <- llik$AICc; bic4 <- llik$BIC
    
    family <- zip.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic5 <- llik$AIC; aicc5 <- llik$AICc; bic5 <- llik$BIC
    
    family <- negbin.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic6 <- llik$AIC; aicc6 <- llik$AICc; bic6 <- llik$BIC
}else if(corr == "powerexp"){
    corr <- powerexp.gcgc2(nugget = F)
    family <- poisson.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic1 <- llik$AIC; aicc1 <- llik$AICc; bic1 <- llik$BIC
    
    family <- zip.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic2 <- llik$AIC; aicc2 <- llik$AICc; bic2 <- llik$BIC
    
    family <- negbin.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic3 <- llik$AIC; aicc3 <- llik$AICc; bic3 <- llik$BIC
    
    corr <- matern.gcgc2(nugget = T)
    family <- poisson.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic4 <- llik$AIC; aicc4 <- llik$AICc; bic4 <- llik$BIC
    
    family <- zip.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic5 <- llik$AIC; aicc5 <- llik$AICc; bic5 <- llik$BIC
    
    family <- negbin.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic6 <- llik$AIC; aicc6 <- llik$AICc; bic6 <- llik$BIC
  }else if(corr == "spherical"){
    corr <- spherical.gcgc2(nugget = F)
    family <- poisson.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic1 <- llik$AIC; aicc1 <- llik$AICc; bic1 <- llik$BIC
    
    family <- zip.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic2 <- llik$AIC; aicc2 <- llik$AICc; bic2 <- llik$BIC
    
    family <- negbin.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic3 <- llik$AIC; aicc3 <- llik$AICc; bic3 <- llik$BIC
    
    corr <- spherical.gcgc2(nugget = T)
    family <- poisson.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr,
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic4 <- llik$AIC; aicc4 <- llik$AICc; bic4 <- llik$BIC
    
    family <- zip.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic5 <- llik$AIC; aicc5 <- llik$AICc; bic5 <- llik$BIC
    
    family <- negbin.gcgc2(link = 'log')
    llik <- mleGHK(y = y, x = x, obs.locs = obs.locs, family = family, corr = corr, 
                   effort = effort, nrep = nrep, seed = seed, max.range = max.range)
    aic6 <- llik$AIC; aicc6 <- llik$AICc; bic6 <- llik$BIC
  }
  
  DF <- data.frame(Model = c("Poisson without Nugget", "ZIP without Nugget", "NB2 without Nugget", 
                             "Poisson with Nugget", "ZIP with Nugget", "NB2 with Nugget"),
                   AIC = c(aic1, aic2, aic3, aic4, aic5, aic6),
                   AICc = c(aicc1, aicc2, aicc3, aicc4, aicc5, aicc6),
                   BIC = c(bic1, bic2, bic3, bic4, bic5, bic6))
  return(DF)
}
  
  
  selectmodel.gcgc(y = simtmp.y, x = cbind(xloc,yloc), obs.locs = cbind(xloc,yloc), 
                  corr = 'spherical', effort = 1, max.range = 100, seed = 123, nrep = c(50,500))
   