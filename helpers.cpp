// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace Rcpp;


// [[Rcpp::export]]
int rztpois_single(double lambda) {
  if (lambda < 0 || !R_FINITE(lambda)) return R_NaN;
  
  /* limiting case as lambda approaches zero is point mass at one */
  if (lambda < 0.000001) return 1.0;
  return Rf_qpois(Rf_runif(exp(-lambda), 1), lambda, 1, 0);
}

// [[Rcpp::export]]
arma::rowvec rztpois_cpp(unsigned int n, double lambda){
  arma::rowvec d;
  d.zeros(n);
  for (unsigned int i = 0; i < n; ++i) {
    d(i) = rztpois_single(lambda);
  }
  return d;
}
// [[Rcpp::export]]
arma::Cube<int> calcC_kn(arma::Cube<int> ZZip_3D, arma::vec bias_0, arma::Cube<int> W_3D, arma::mat Phi) {
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
        } else if (ZZip_3D(k, 0, n) == 1){
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

// [[Rcpp::export]]
List mult_cpp(arma::mat X_mtn, arma::mat Psi, arma::mat Theta, arma::mat ZZip) {
  int M = Psi.n_rows;
  int K = Psi.n_cols;
  int NT = Theta.n_cols;
  arma::Cube<double> lambda_mknt(M, K, NT);
  arma::Cube<double> lambda_mknt_hat(M, K, NT);
  arma::Cube<int> latent_count(M, K, NT);

  for (int n = 0; n < NT; ++n) {
    for (int m = 0; m < M; ++m) {
      for (int k = 0; k < K; ++k) {
        lambda_mknt(m, k, n) = Psi(m, k) * Theta(k, n) * ZZip(k, n);
      }
    }
  }
  
  for (int n = 0; n < NT; ++n) {
    for (int m = 0; m < M; ++m) {
      for (int k = 0; k < K; ++k) {
        auto s = sum(lambda_mknt.slice(n).col(k));
        if (s != 0)
          lambda_mknt_hat(m, k, n) = lambda_mknt(m, k, n) / s;
        else
          lambda_mknt_hat(m, k, n) = 0;
      }
    }
  }

  
  for (int n = 0; n < NT; ++n) {
    for (int m = 0; m < M; ++m) {
      if (sum(lambda_mknt_hat.slice(n).row(m)) !=0 ) {
        arma::rowvec probs = lambda_mknt_hat.slice(n).row(m);
        double s = arma::sum(probs);
        probs = probs / s;
        arma::irowvec v(K);
        Rf_rmultinom(X_mtn(m, n), probs.begin(), K, v.begin());
        latent_count.slice(n).row(m) = v;
      } else {
        latent_count.slice(n).row(m).fill(0);
      }
    }
  }

  arma:: mat X_mk2(M, K);
  for (int n = 0; n < NT; ++n) {
    X_mk2 = X_mk2 + latent_count.slice(n);
  }

  arma:: mat X_kn2(K, NT);
  for (int k = 0; k < K; ++k) {
    for (int n = 0; n < NT; ++n) {
      for (int m = 0; m < M; ++m) {
        X_kn2(k, n) = X_kn2(k, n) + latent_count(m, k, n);
      }
    }
  }


  return Rcpp::List::create(Rcpp::Named("X_mk2") = X_mk2, Rcpp::Named("X_kn2") = X_kn2);  
}

// [[Rcpp::export]]
arma::rowvec crt_cpp(arma::mat X_kn, arma::vec rk) {
  int K = rk.n_elem;
  int nTotal = X_kn.n_cols;
  arma::rowvec L(K);

  for (int k = 0; k < K; k++) {
    for (int n = 0; n < nTotal; ++n) {
      for (int i = 0; i < X_kn(k, n); ++i) {
				L[k] += (double)( Rf_runif(0, 1) <= rk[k]/( rk[k] + i ) );
      }
    }
  }

  return L;
}