crt <- function(x_kn,rk){
  K = length(rk)
  numTotal = ncol(x_kn)
  Lk = rep(0,K)
  
  for(k in 1:K){
    for(n in 1:numTotal){
      for(i in 1:x_kn[k,n]){
        # print(Lk)
        Lk[k] = Lk[k] + (runif(1) <= rk[k]/(rk[k]+(i-1)))*1
      }
    }
  }
  Lk
}