#include <Rmath.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace arma; 
using namespace Rcpp;


// [[Rcpp::export]]
double FHUBNB2 (double m1, double m2, double od1, double od2){
  int N1 = 0; int N2 = 0; int i; int j; double corr;
  double SZ1 = 1/od1; double PB1 = SZ1/(SZ1+m1);
  double SZ2 = 1/od2; double PB2 = SZ2/(SZ2+m2);
  while(R::pnbinom(N1, SZ1, PB1, 1, 0) < 1)   N1 += 1;
  while(R::pnbinom(N2, SZ2, PB2, 1, 0) < 1)   N2 += 1;
  if(N1 > 8000) throw("ERROR");
  if(N2 > 8000) throw("ERROR");
  arma::mat Exytmp(N1, N2);
  for (i=0; i<N1; i++){
    for (j=0; j<N2; j++){
    Exytmp(i,j) = 1-R::fmax2(R::pnbinom(i, SZ1, PB1, 1, 0), 
                             R::pnbinom(j, SZ2, PB2, 1, 0));
    }
  }
  corr = (arma::accu(Exytmp)-m1*m2)/(sqrt(m1*m2*(1+od1*m1)*(1+od2*m2)));
  return corr;
}
  
  
// [[Rcpp::export]]
double FHUBZIP (double m1, double m2, double od1, double od2){
    int N1 = 0; int N2 = 0; int i; int j; double corr;
    double pi1 = od1/(1+od1); double lambda1 = (1+od1)*m1;
    double pi2 = od2/(1+od2); double lambda2 = (1+od2)*m2;
    while(pi1+(1-pi1)*R::ppois(N1, lambda1, 1, 0) < 1)  N1 += 1;
    while(pi2+(1-pi2)*R::ppois(N2, lambda2, 1, 0) < 1)  N2 += 1;
    if(N1 > 8000) throw("ERROR");
    if(N2 > 8000) throw("ERROR");
    arma::mat Exytmp(N1, N2);
    for (i=0; i<N1; i++){
        for (j=0; j<N2; j++){
            Exytmp(i,j) = 1-R::fmax2(pi1+(1-pi1)*R::ppois(i, lambda1, 1, 0), 
                                     pi2+(1-pi2)*R::ppois(j, lambda2, 1, 0));
        }
    }
    corr = (arma::accu(Exytmp)-m1*m2)/(sqrt(m1*m2*(1+od1*m1)*(1+od2*m2)));
    return corr;
}



// [[Rcpp::export]]
double FHUBbinom (double m1, double m2, double n1, double n2){
   int i; int j; double corr;
   double p1 = m1/n1; double p2 = m2/n2;
   if(n1 > 8000) throw("ERROR");
   if(n2 > 8000) throw("ERROR");
   arma::mat Exytmp(n1, n2);
   for (i=0; i<n1; i++){
     for (j=0; j<n2; j++){
       Exytmp(i,j) = 1-R::fmax2(R::pbinom(i, n1, p1, 1, 0), 
                                R::pbinom(j, n2, p2, 1, 0));
    }
  }
   corr = (arma::accu(Exytmp)-m1*m2)/(sqrt(m1*m2*(1-p1)*(1-p2)));
   return corr;
}


// [[Rcpp::export]]
double FHUBZIPNB2 (double zipmu, double nbmu, double zipod, double nbod){
  int N1 = 0; int N2 = 0; int i; int j; double corr;
  double pi = zipod/(1+zipod); double lambda = (1+zipod)*zipmu;
  double SZ = 1/nbod; double PB = SZ/(SZ+nbmu);
  while(pi+(1-pi)*R::ppois(N1, lambda, 1, 0) < 1)  N1 += 1;
  while(R::pnbinom(N2, SZ, PB, 1, 0) < 1)   N2 += 1;
  if(N1 > 8000) throw("ERROR");
  if(N2 > 8000) throw("ERROR");
  arma::mat Exytmp(N1, N2);
  
  for (i=0; i<N1; i++){
    for (j=0; j<N2; j++){
      Exytmp(i,j) = 1-R::fmax2( pi+(1-pi)*R::ppois(i, lambda, 1, 0), 
                                R::pnbinom(j, SZ, PB, 1, 0));
    }
  }
  corr = (arma::accu(Exytmp)-zipmu*nbmu)/(sqrt(zipmu*nbmu*(1+zipod*zipmu)*(1+nbod*nbmu)));
  return(corr);
}


// [[Rcpp::export]]
double FHUBZIPbinomial (double zipmu, double bmu, double zipod, double bn){
  int N1 = 0; int i; int j; double corr;
  double pi = zipod/(1+zipod); double lambda = (1+zipod)*zipmu;
  double p = bmu/bn;
  while(pi+(1-pi)*R::ppois(N1, lambda, 1, 0) < 1)  N1 += 1;
  if(N1 > 8000) throw("ERROR");
  if(bn > 8000) throw("ERROR");
  arma::mat Exytmp(N1, bn);
  
  for (i=0; i<N1; i++){
    for (j=0; j<bn; j++){
      Exytmp(i,j) = 1-R::fmax2( pi+(1-pi)*R::ppois(i, lambda, 1, 0), 
                                R::pbinom(j, bn, p, 1, 0));
    }
  }
  corr = (arma::accu(Exytmp)-zipmu*bmu)/(sqrt(zipmu*bmu*(1+zipod*zipmu)*(1-p)));
  return(corr);
}





// [[Rcpp::export]]
double FHUBNB2binomial (double nbmu, double bmu, double nbod, double bn){
  int N1 = 0; int i; int j; double corr;
  double SZ = 1/nbod; double PB = SZ/(SZ+nbmu);
  double p = bmu/bn;
  while(R::pnbinom(N1, SZ, PB, 1, 0) < 1)   N1 += 1;
  if(N1 > 8000) throw("ERROR");
  if(bn > 8000) throw("ERROR");
  arma::mat Exytmp(N1, bn);
  for (i=0; i<N1; i++){
    for (j=0; j<bn; j++){
      Exytmp(i,j) = 1-R::fmax2( R::pnbinom(i, SZ, PB, 1, 0), 
                                R::pbinom(j, bn, p, 1, 0));
    }
  }
  corr = (arma::accu(Exytmp)-nbmu*bmu)/(sqrt(nbmu*bmu*(1+nbod*nbmu)*(1-p)));
  return(corr);
}

