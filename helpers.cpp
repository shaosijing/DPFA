// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace Rcpp;


// C_kn = array(NA, dim= c(K, numTime, numSample))
//   for (n in 1:numSample){
//     for (k in 1:K){
//       if (ZZip_3D[k,1,n] == 0){
//         C_kn[k,1,n] = 0
//       } else if (ZZip_3D[k,1,n] ==1){
//         C_kn[k,1,n] = rztpois(1,bias_0[k])
//       }
//     }
//   }
//   
//   for (n in 1:numSample){
//     for (t in 2:numTime){
//       for (k in 1:K){
//         if (ZZip_3D[k,t,n]==0){
//           C_kn[k,t,n] = 0
//         }
//         else if (ZZip_3D[k,t,n]==1){
//           C_kn[k,t,n] = rztpois(1, (Phi%*%(W_3D[,t-1,n]*ZZip_3D[,t-1,n]))[k]+bias_0[k])}
// #  if (is.na(C_kn[k,t,n]) == TRUE) print(c(n, t))
//       }
//     }
//   }


// [[Rcpp::export]]
int rztpois_single(double lambda) {
  if (lambda < 0 || !R_FINITE(lambda)) return R_NaN;
  
  /* limiting case as lambda approaches zero is point mass at one */
  if (lambda == 0) return 1.0;
  return Rf_qpois(Rf_runif(exp(-lambda), 1), lambda, 1, 0);
}

// [[Rcpp::export]]
arma::vec rztpois_cpp(unsigned int n, double lambda){
  arma::vec d;
  d.zeros(n);
  for (unsigned int i = 0; i < n; ++i) {
    d(i) = rztpois_single(lambda);
  }
  return d;
}
// [[Rcpp::export]]
arma::Cube<int> calcC_kn(arma::Cube<int> ZZip_3D, arma::vec bias_0, arma::Cube<int> W_3D, arma::mat Phi) {
  Environment pkg = Environment::namespace_env("actuar");
  Function f = pkg["rztpois"];
  
  int K = ZZip_3D.n_rows;
  int nTime = ZZip_3D.n_cols;
  int nSample = ZZip_3D.n_slices;
  arma::Cube<int> C_kn(K, nTime, nSample);
  for (int n = 0; n < nSample; ++n) {
    for (int k = 0; k < K; ++k) {
      if (ZZip_3D(k, 0, n) == 0) {
        C_kn(k, 0, n) = 0;
      } else if (ZZip_3D(k, 0, n) ==1){
        C_kn(k, 0, n) = rztpois_single(bias_0(k));
      }
    }
  }

  for (int n = 0; n < nSample; ++n) {
    for (int t = 1; t < nTime; ++t) {
      for (int k = 0; k < K; ++k) {
        if (ZZip_3D(k, t, n) == 0) {
          C_kn(k, t, n) = 0;
        } else if (ZZip_3D(k, 0, n) ==1){
          auto w1 = W_3D.slice(n);
          auto z2 = ZZip_3D.slice(n);
          arma::mat ret = Phi * (w1.col(t-1) % z2.col(t-1));
          C_kn(k, t, n) = rztpois_single(ret(k) * bias_0(k));
        }
      }
    }
  }
  return C_kn;
}


// [[Rcpp::export]]
arma::mat calcTheta(arma::vec rk, arma::mat ZZip, arma::mat x_kn, double p0) {
  int N = x_kn.n_rows;
  int M = x_kn.n_cols;
  arma::mat Theta(N, M);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      Theta(i,j) = Rf_rgamma(rk(i) * ZZip(i,j) + x_kn(i,j), 1/p0);
    }
  }
  return Theta;
}

// [[Rcpp::export]]
arma::mat calcW(arma::vec sk, arma::mat ZZip, arma::mat C_k1n, double p0) {
  int N = ZZip.n_rows;
  int M = ZZip.n_cols;
  arma::mat W(N, M);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      W(i,j) = Rf_rgamma(sk(i) * ZZip(i,j) + C_k1n(i,j), 1/p0);
    }
  }
  return W;
}