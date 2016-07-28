library('VGAM')
##### First Create Zero-Inflated Poisson functions (pdf, cdf, quantiles) according to 
##### my OWN parametrization that is similar as NB2
##### Using packages (VGAM) of standard parametrization by definition
##### Size is 1/Overdispersion, which is alphasq
##### Just used for Simulation case: two dimensional region, three beta.

dzip2 <- function(x, size, mu)
{
  lambda <- mu+mu/size
  pi <- 1/(size+1)
  out <- VGAM::dzipois(x = x, lambda = lambda, pstr0 = pi)
  return(out)
}

pzip2 <- function(q, size, mu)
{
  lambda <- mu+mu/size
  pi <- 1/(size+1)
  out <- VGAM::pzipois(q = q, lambda = lambda, pstr0 = pi)
  return(out)
}

qzip2 <- function(p, size, mu)
{
  lambda <- mu+mu/size
  pi <- 1/(size+1)
  out <- VGAM::qzipois(p=p, lambda=lambda, pstr0 = pi)
  return(out)
}