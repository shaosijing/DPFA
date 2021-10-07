sample_Z_timespan_try3 = function(x_kn , p0, rk, Phi,W, sk, p1,C_k1n, Pi_k, numSample,ZZip){
  numTime = ncol(x_kn)/numSample
  K = nrow(x_kn)
  
  
  BPL = function(lambda){
    1-exp(-lambda)
  }
  
  lambda= Phi %*% W #lambda= Phi %*% W
  Pi = BPL(lambda)
  
  lix = (C_k1n == 0 & x_kn ==0)
  ind = which(C_k1n == 0 & x_kn ==0,arr.ind=T)
  rix = ind[,1];cix=ind[,2]
  
  p_1 = Pi[lix]*(( 1 - p0 )^rk[rix]) * ((1-p1)^sk[rix]) 
  p_0 = 1 - Pi[rix]
  ZZip = matrix(1, nrow = K, ncol = ncol(x_kn))
  #ZZip[lix] =  rbern(length(lix),.1)#((p_1/( p_1 + p_0 ) ) > runif(length(rix)))*1
  ZZip[lix] =  ((p_1/( p_1 + p_0 ) ) > runif(length(rix)))*1
  
  
  lix2 = (C_k1n > 0 | x_kn > 0) 
  ind2 = which(C_k1n > 0 | x_kn > 0,arr.ind=T)
  rix2 = ind2[,1];cix2=ind2[,2]
  ZZip[lix2] = 1
  res = list(ZZip)
  return(res)
}