dynamicPFA = function(dat, K, numTime, niter, burnin){
  M = dim(dat)[1]
  numTotal = dim(dat)[2]
  numSample = numTotal/numTime
  timeSelect = 1:numTime
  timeSpan = c(diff(timeSelect),1)
  # hyperparameters
  hyp = .1;
  alpha_psi = hyp;
  p0 = 0.5;
  p1 = 0.5;
  rk_a = 100;
  rk_b = 1/rk_a;
  alpha_phi = .01
  p1 = .5
  sk_a = 1e5
  sk_b=1/sk_a
  a0 = 1; b0 = 1;
  rk = rep(1,K)
  sk = rgamma(K,sk_a)*sk_b
  p0 = 0.5
  p1 = 0.5
  bias_0 = rep(1,K)

  # initialization

  pca=prcomp(t(Xmtot))
  Psi = abs(pca$rotation[,1:3])

  K1 = K
  K2 = K
  Phi = matrix(rgamma(K1*K2,shape=alpha_phi),K1,K2) #randg(alpha_psi,P,K)
  for (i in 1:K2){
    Phi[,i] = Phi[,i]/sum(Phi[,i])
  }
  Theta<-array(dim=c(K,numTotal))
  for (n in 1:numTotal){
    for (k in 1:K){
      Theta[k,n] = rgamma(1,1,0.5)
    }
  }


z0 = matrix(1,K,numSample)

ZZip = (matrix(runif(K*numTotal),K,numTotal)<1)*1

ZZip2 = update_Z(z0,ZZip,numTime)
W = ZZip2
Pi_k <- BPL(Phi %*% (W *ZZip))

rett <- list()
#fits <- matrix(NA,niter,2)

for(b in 1:niter){
  out = mult_cpp(Xmtot,Psi,Theta, ZZip)
  x_pk = out[[1]]
  x_kn = out[[2]]

  W_time = sweep(W,2,timeSpan,'/')


  W_3D = matrix_to_array(W_time, K, numSample,numTime)
  ZZip_3D = matrix_to_array(ZZip,K, numSample,numTime)
  C_kn <- calcC_kn(ZZip_3D, bias_0, W_3D, Phi)
  C_kn = matrix(C_kn, nrow=K)

  out2 = mult_cpp(C_kn,Phi,W, ZZip)
  C_kk1 = out2[[1]]
  C_k1n = out2[[2]]

  res=sample_Z(x_kn , p0, rk, Phi,W_time, sk, p1,C_k1n,Pi_k, numSample,ZZip)
  ZZip = res[[1]]

  Psi = matrix(NA,dim(x_pk)[1],dim(x_pk)[2])
  for (i in 1:dim(x_pk)[2]){
    Psi[,i] <-  rdirichlet(1,(alpha_psi+x_pk)[,i])
  }

  # chinese restaurant table distribution
  Lk = crt_cpp(x_kn,rk)

  sumbpi = rowSums(ZZip) * log(1-p0)
  rk = rgamma(K, rk_a + Lk)/( rk_b - sumbpi); # from code

  Theta = calcTheta(rk, ZZip, x_kn, p0)

  Phi = matrix(NA,dim(x_pk)[2],dim(x_pk)[2])
  for (i in 1:dim(x_pk)[2]){
    Phi[,i] <-  rdirichlet(1,(alpha_psi+C_kk1)[,i])
  }

  Lk = crt_cpp(C_k1n,sk)
  sumbpi = rowSums(ZZip) * log(1-p0)
  sk = rgamma(K, sk_a + Lk)/( 1/sk_b - sumbpi); # from code
  W = calcW(sk, ZZip, C_k1n, 0.5)

  # Pi_k = rbeta(K,shape1 = a0 + rowSums(ZZip),shape2=b0 + numTotal*numTime - rowSums(ZZip));
  # Pi_k = matrix(rbeta(K*numTotal,a0,b0),K,numTotal)
  #theta = matrix(rgamma(rk*ZZip+ x_kn),dim(x_kn)) * p0;# from code
  #theta = matrix(rgamma(matrix(rk,K,numTotal) + x_kn,1),dim(x_kn)) # from paper

  # nz = (Xmtot[,1:numTime]!=0)*1
  # lambda = Psi%*%Theta[,1:numTime]
  # RMSE = sum((Xmtot[,1:numTime][nz==1]-lambda[nz==1])^2)/numTotal
  # lambda_avg = lambda/colSums(lambda)
  # negLL = -  sum(Xmtot[nz==1]*log(lambda_avg[nz==1]))/sum(Xmtot[nz==1])
  rett[[b]] <- list(Theta=Theta,Lk=Lk,Psi=Psi,rk=rk, Phi = Phi, Xmtot = Xmtot, sumbpi = sumbpi, W = W, x_kn = x_kn, ZZip = ZZip,
                    C_kn = C_kn, C_k1n = C_k1n, C_kk1=C_kk1)
}

  psi.array <- array(NA,c(dim(Psi),(niter-burnin)))

  for(i in 1:(niter-burnin)){
    psi.array[,,i] <- rett[[i+burnin]]$Psi
  }


  sum.psi <- matrix(0,dim(Psi)[1],dim(Psi)[2])

  for(i in 1:(niter-burnin)){
    sum.psi <- sum.psi + psi.array[,,i]
  }
  psi_est = round(sum.psi/(niter-burnin),2)


  Phi.array <- array(NA,c(dim(Phi),(niter-burnin)))
  for(i in 1:(niter-burnin)){
    Phi.array[,,i] <- rett[[burnin+i]]$Phi
  }
  sum.Phi <- matrix(0,dim(Phi)[1],dim(Phi)[2])

  for(i in 1:(niter-burnin)){
    sum.Phi <- sum.Phi + Phi.array[,,i]
  }

  phi_est = round(sum.Phi/(niter-burnin),3)


  Theta.array <- array(NA,c(dim(Theta),(niter-burnin)))
  for(i in 1:(niter-burnin)){
    Theta.array[,,i] <- rett[[burnin+i]]$Theta
  }
  sum.Theta <- matrix(0,dim(Theta)[1],dim(Theta)[2])

  for(i in 1:(niter-burnin)){
    sum.Theta <- sum.Theta + Theta.array[,,i]
  }

  theta_est = round(sum.Theta/(niter-burnin),3)

  W.array <- array(NA,c(dim(rett[[niter]]$W),(niter-burnin)))
  for(i in 1:(niter-burnin)){
    W.array[,,i] <- rett[[burnin+i]]$W
  }
  mean.W <- matrix(0,dim(W)[1],dim(W)[2])

  for(i in 1:(niter-burnin)){
    mean.W <- mean.W + W.array[,,i]
  }

  W_est = round(mean.W/(niter-burnin),3)

  ZZip.array <- array(NA,c(dim(rett[[niter]]$ZZip),(niter-burnin)))
  for(i in 1:(niter-burnin)){
    ZZip.array[,,i] <- rett[[burnin+i]]$ZZip
  }
  mean.ZZip <- matrix(0,dim(ZZip)[1],dim(ZZip)[2])

  for(i in 1:(niter-burnin)){
    mean.ZZip <- mean.ZZip + ZZip.array[,,i]
  }

  ZZip_est = round(mean.ZZip/(niter-burnin),3)
  ZZip_rowsums = rowSums(ZZip_est)

  return<- list(psi_est,phi_est,theta_est,W_est,ZZip_est,ZZip_rowsums)
}

