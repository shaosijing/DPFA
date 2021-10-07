update_Z <- function(z0,ZZip,time){
  K = nrow(ZZip)
  numTotal = ncol(ZZip)
  t1_ix = seq(1,numTotal,time)
  ZZip2 = cbind(matrix(0,K,1),ZZip[,1:(ncol(ZZip)-1)])
  ZZip2[,t1_ix] = z0
  ZZip2
}
