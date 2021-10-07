library(R.utils)

mult_rnd4 = function(x_mtn,psi, theta, ZZip){
  M = nrow(psi)
  K = ncol(psi)
  NT = dim(theta)[2]
 # Time = dim(theta)[3]
  
  lambda_mknt = array(NA, dim= c(M,K,NT))
  for (n in 1:NT){
      for (m in 1: M){
        for (k in 1:K){
          lambda_mknt[m,k,n] = psi[m,k] * (theta[k,n]*ZZip[k,n]) # removed for PFA
        }
      }
  }
  
  lambda_mknt_hat = array(NA, dim= c(M,K,NT))
  for (n in 1:NT){
      for (m in 1:M){
        for (k in 1:K){
          if (sum(lambda_mknt[1:M,k,n]) !=0){
            lambda_mknt_hat[m,k,n] = lambda_mknt[m,k,n]/sum(lambda_mknt[1:M,k,n])
          }
          else {lambda_mknt_hat[m,k,n] = 0}
        }
      }
  }
  
 # lambda_mn = matrix(NA, M,NT)
 # for (m in 1:M){
  #  for (n in 1:NT){
  #    lambda_mn[m,n] = sum(lambda_mknt[m,1:K,n])
  #  }
  #}
  
  latent_count = array(NA, dim=c(M,K,NT))
  for (n in 1:NT){
    for (m in 1:M){
        if (sum(lambda_mknt_hat[m,,n]) !=0){
          latent_count[m,,n] <- t(rmultinom(1,x_mtn[m,n],lambda_mknt_hat[m,,n]))
        }
        else {latent_count[m,,n] <- 0}
      }
    }
  
  X_mk2 = matrix(0,M,K)
  for(n in 1:NT){
    X_mk2 = X_mk2 + latent_count[,,n]
  }
  
  X_kn2 = matrix(0,K,NT)
  for(m in 1:M){
    X_kn2 = X_kn2 + latent_count[m,,]
  }
  
  res=list(X_mk2, X_kn2)
  return(res)
}

  
 