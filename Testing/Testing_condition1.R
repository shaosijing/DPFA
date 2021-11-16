library(extraDistr)
library(Morpho)
library(actuar)
library(Rcpp)
library(RcppArmadillo)
setwd("H:/Shared drives/SLAM Lab/ema_text/DFA R/dynamic_pfa/CRC/DPFA/")

sourceCpp('helpers.cpp')
#source("crt.R")
source("update_Z.R")
#source("mult.R")
source("sample_Z.R")
source("fit_DIC.R")
setwd("H:/Shared drives/SLAM Lab/ema_text/DFA R/dynamic_pfa/CRC/DPFA/Testing/")

dataset = readRDS("true_para_data1.rds")
psi_tru = dataset[[1]]
theta_tru = dataset[[2]]
Xmtot = dataset[[3]]
phi_tru = dataset[[4]]
ZZip = dataset[[5]]
bias_0 = dataset[[6]]

M=15
K=3
numSample=100
numTime = 10
numTotal = (numTime) * numSample
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

BPL<-function(lambda){
  return(1-exp(-lambda))
}

norm = function(x){
  sum1 = sum(x)
  if (sum1 != 0){
    x = x/sum1
  }
  return(x)
}

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

rk = c(1,1,1)
z0 = matrix(1,K,numSample)
sk = rgamma(K,sk_a)*sk_b

ZZip = (matrix(runif(K*numTotal),K,numTotal)<1)*1

ZZip2 = update_Z(z0,ZZip,numTime)
W = ZZip2
Pi_k <- BPL(Phi %*% (W *ZZip))

niter=5000
rett <- list()
fits <- matrix(NA,niter,2)

t1=proc.time()
for(b in 1:niter){
  out = mult_cpp(Xmtot,Psi,Theta, ZZip)
  x_pk = out[[1]]
  x_kn = out[[2]]

  W_time = sweep(W,2,timeSpan,'/')

  matrix_to_array = function(matrix){
    K = 3
    numSample = 100
    numTime = 10
    output = array(NA, dim = c(K, numTime, numSample))
    for (i in 1:numSample){
      output[,,i]= matrix[,(10*i-9):(10*i)]
    }
    return(output)
  }
  W_3D = matrix_to_array(W_time)
  ZZip_3D = matrix_to_array(ZZip)
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
  # fits[b,] <- c(RMSE,negLL)

  #print(b)
}
t2 = proc.time() - t1
t2

fS <- 2000
psi.array <- array(NA,c(dim(Psi),(niter-fS)))

for(i in 1:(niter-fS)){
  psi.array[,,i] <- rett[[i+fS]]$Psi
}


sum.psi <- matrix(0,dim(Psi)[1],dim(Psi)[2])

for(i in 1:(niter-fS)){
  sum.psi <- sum.psi + psi.array[,,i]
}
psi_est = round(sum.psi/(niter-fS),2)


Phi.array <- array(NA,c(dim(Phi),(niter-fS)))
for(i in 1:(niter-fS)){
  Phi.array[,,i] <- rett[[fS+i]]$Phi
}
sum.Phi <- matrix(0,dim(Phi)[1],dim(Phi)[2])

for(i in 1:(niter-fS)){
  sum.Phi <- sum.Phi + Phi.array[,,i]
}

phi_est = round(sum.Phi/(niter-fS),3)


Theta.array <- array(NA,c(dim(Theta),(niter-fS)))
for(i in 1:(niter-fS)){
  Theta.array[,,i] <- rett[[fS+i]]$Theta
}
sum.Theta <- matrix(0,dim(Theta)[1],dim(Theta)[2])

for(i in 1:(niter-fS)){
  sum.Theta <- sum.Theta + Theta.array[,,i]
}

theta_est = round(sum.Theta/(niter-fS),3)

W.array <- array(NA,c(dim(rett[[niter]]$W),(niter-fS)))
for(i in 1:(niter-fS)){
  W.array[,,i] <- rett[[fS+i]]$W
}
mean.W <- matrix(0,dim(W)[1],dim(W)[2])

for(i in 1:(niter-fS)){
  mean.W <- mean.W + W.array[,,i]
}

W_est = round(mean.W/(niter-fS),3)

ZZip.array <- array(NA,c(dim(rett[[niter]]$ZZip),(niter-fS)))
for(i in 1:(niter-fS)){
  ZZip.array[,,i] <- rett[[fS+i]]$ZZip
}
mean.ZZip <- matrix(0,dim(ZZip)[1],dim(ZZip)[2])

for(i in 1:(niter-fS)){
  mean.ZZip <- mean.ZZip + ZZip.array[,,i]
}

ZZip_est = round(mean.ZZip/(niter-fS),3)
ZZip_rowsums = rowSums(ZZip_est)

return<- list(psi_est,phi_est,theta_est,W_est,ZZip_est,ZZip_rowsums)
round(psi_est,2)
phi_est
# psi_1_1= matrix(NA, nrow = niter,ncol = 1)
# for (i in 1:niter){
#   psi_1_1[i,1] = rett[[i]][["Psi"]][1,1]
# }
# psi_1_2= matrix(NA, nrow = niter,ncol = 1)
# for (i in 1:niter){
#   psi_1_2[i,1] = rett[[i]][["Psi"]][1,2]
# }
# psi_1_3= matrix(NA, nrow = niter,ncol = 1)
# for (i in 1:niter){
#   psi_1_3[i,1] = rett[[i]][["Psi"]][1,3]
# }
# psi_2_1= matrix(NA, nrow = niter,ncol = 1)
# for (i in 1:niter){
#   psi_2_1[i,1] = rett[[i]][["Psi"]][2,1]
# }
# psi_2_2= matrix(NA, nrow = niter,ncol = 1)
# for (i in 1:niter){
#   psi_2_2[i,1] = rett[[i]][["Psi"]][2,2]
# }
# psi_2_3= matrix(NA, nrow = niter,ncol = 1)
# for (i in 1:niter){
#   psi_2_3[i,1] = rett[[i]][["Psi"]][2,3]
# }
# psi13_first_row = cbind(psi_1_1,psi_1_2,psi_1_3)
# psi13_second_row = cbind(psi_2_1,psi_2_2,psi_2_3)
# par(mfrow=c(2,1))
# matplot(psi13_first_row)

fit_DIC(rett, K, N, numTime, M)
matplot(psi13_second_row)
