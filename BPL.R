BPL<-function(lambda,K){
  prob = 1-exp(-lambda)
  out = rbern(rep(1,K),prob)
  return(matrix(out, nrow = K))
}
