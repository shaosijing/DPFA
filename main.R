# generate data, with true K = 3
library(extraDistr)
library(Morpho)
library(actuar)
M=15
K=3
numSample=100
numTime = 10
numTotal = (numTime) * numSample
rk_tru = rep(1,K)
sk_tru = rep(1,K)
p0 = 0.5
p1 = 0.5
phi_tru<-diag(1,K)

psi_tru <- matrix(c(0.19,0.005,0.005,
                    0.19,0.005,0.005,
                    0.19,0.005,0.005,
                    0.19,0.005,0.005,
                    0.19,0.005,0.005,
                    0.005,0.19,0.005,
                    0.005,0.19,0.005,
                    0.005,0.19,0.005,
                    0.005,0.19,0.005,
                    0.005,0.19,0.005,
                    0.005,0.005,0.19,
                    0.005,0.005,0.19,
                    0.005,00.005,0.19,
                    0.005,00.005,0.19,
                    0.005,0.005,0.19),nrow = M, byrow = T)
theta_tru<-array(dim=c(K,numTotal))
for (n in 1:numTotal){
  for (k in 1:K){
    theta_tru[k,n] = rgamma(1,rk_tru[k],0.5)
  }
}
W_tru<-array(dim=c(K,numTotal))
for (n in 1:numTotal){
  for (k in 1:K){
    W_tru[k,n] = rgamma(1,rk_tru[k],0.5)
  }
}
p0 = 0.5
p1 = 0.5
bias_0 = rep(1,K)
BPL<-function(lambda,K){
  prob = 1-exp(-lambda)
  out = rbern(rep(1,K),prob)
  return(matrix(out, nrow = K))
}

#sample the initial ZZip for each person
ZZip_1 = matrix(NA, ncol = numSample, nrow = K)
for (i in 1:numSample){
  ZZip_1[,i] = BPL(bias_0,K)
}

#sample each person's ZZip (each person has numTime time points)
ZZip = array(NA, dim =c(K, numTime, numSample))
for (n in 1:numSample){
  for (t in 1: (numTime-1)){
    ZZip[,1,n] = ZZip_1[,n]
    ZZip[,t+1,n] = BPL(phi_tru %*% (W_tru[,t])*ZZip[,1,n]+bias_0,K)
  }
}

lambda = array(NA, dim = c(M,numSample*numTime))
count = 1
for (n in 1:numSample){
  for (t in 1:numTime){
    lambda[,count] = psi_tru %*% (theta_tru[,count] * ZZip[,t,n])
    count = count + 1
  }
}

Xmtot = matrix(NA,M,numTotal)
count = 1
#for(m in 1:M){
for(j in 1:numSample){
  for(t in 1:numTime){
    Xmtot[,count] = rpois(M, lambda[,count])
    count = count + 1
    #    }
  }
}

ZZip = matrix(ZZip, nrow = K)

dataset<- list(psi_tru, theta_tru, Xmtot,phi_tru,ZZip,bias_0)

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
Rep <-1
resfolder <- expand.grid(K,numSample,M,numTime)
A<-expand.grid(K,numSample,M,numTime,Rep)
args<-as.numeric(commandArgs(trailingOnly = TRUE))
K<-A[args,1]
numSample<-A[args,2]
M<-A[args,3]
numTime <-A[args,4]
Rep<-A[args,5]
seed <- Rep+100
source("crt.R")
source("update_Z.R")
source("mult.R")
#out2 = mult_rnd3(C_kn,Phi,W, ZZip)
source("sample_Z.R")
#source("mult_rnd_latent.R")

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


#test with case of 100 participants


#Xmtot[,(numSample*numTime+1):(numSample*(numTime+1))] = rowMeans(Xmtot,na.rm=T)

# initialization

Psi = psi_tru 
#Psi = matrix(NA,M,K)
#for (i in 1:K){
#  Psi[,i] <-  rdirichlet(1,(alpha_psi)[,i])
#}
rk_tru = c(1,1,1)
Theta=theta_tru
Phi = phi_tru

rk = rk_tru
z0 = matrix(1,K,numSample)
sk = rgamma(K,sk_a)*sk_b


#ZZip = t(apply(H_full,1,rbind))
#ZZip = (matrix(runif(K*numTotal),K,numTotal)<.1)*1

ZZip2 = update_Z(z0,ZZip,numTime)
W = ZZip2 # update later for sparse storage
#Pi_k = matrix(rbeta(K,a0,b0),K,1)
#lbd_hat<-rgamma(3,1,1)
#Pi_k<-BPL(lbd_hat)
Pi_k <- BPL(Phi %*% (W *ZZip))



p1m = matrix(.5,1,numTotal)#rep(.5,numTotal)#rep(p1/(p1+(1-p1)*timeSpan),numSample); # seems wrong

# just do 5000 samples to start
niter=5000
rett <- list()
fits <- matrix(NA,niter,2)

t1=proc.time()
for(b in 1:niter){
  
  out = mult(Xmtot,Psi,Theta, ZZip)
  x_pk = out[[1]]
  x_kn = out[[2]]
  #print(x_pk)
  #print(b)
  W_time = sweep(W,2,timeSpan,'/')
  
  #C_kn = matrix(NA,K,ncol(ZZip))
  #for(i in 1:ncol(ZZip)){
  #  for(k in 1:K){
  #    if (ZZip[k,i] >0){
  #      C_kn[k,i] = rztpois(1,Phi[k,] %*% (W_time[,i]*ZZip[,i]) + bias_0[k])
  #    } 
  #    else if (ZZip[k,i]==0){
  #      C_kn[k,i] = 0
  #    }
  #  }
  #}
  
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
  
  C_kn = array(NA, dim= c(K, numTime, numSample))
  for (n in 1:numSample){
    for (k in 1:K){
      if (ZZip_3D[k,1,n] == 0){
        C_kn[k,1,n] = 0
      } else if (ZZip_3D[k,1,n] ==1){
        C_kn[k,1,n] = rztpois(1,bias_0[k])
      }
    }
  }
  
  for (n in 1:numSample){
    for (t in 2:numTime){
      for (k in 1:K){
        if (ZZip_3D[k,t,n]==0){
          C_kn[k,t,n] = 0
        }
        else if (ZZip_3D[k,t,n]==1){
          C_kn[k,t,n] = rztpois(1, (Phi%*%(W_3D[,t-1,n]*ZZip_3D[,t-1,n]))[k]+bias_0[k])}
        #  if (is.na(C_kn[k,t,n]) == TRUE) print(c(n, t))
      }
    }
  }
  
  
  C_kn = matrix(C_kn, nrow=K)
  
  out2 = mult(C_kn,Phi,W, ZZip)
  C_kk1 = out2[[1]]
  C_k1n = out2[[2]]
  
  
  res=sample_Z(x_kn , p0, rk, Phi,W_time, sk, p1,C_k1n,Pi_k, numSample,ZZip)
  ZZip = res[[1]]
  # ZZip[is.na(ZZip)] <- 0
  ### z0 = res[[2]]
  #ZZip = matrix(1,nrow(x_kn),ncol(x_kn));
  Psi = matrix(NA,dim(x_pk)[1],dim(x_pk)[2])
  for (i in 1:dim(x_pk)[2]){
    Psi[,i] <-  rdirichlet(1,(alpha_psi+x_pk)[,i])
  }
  
  
  
  # chinese restaurant table distribution
  Lk = crt(x_kn,rk)
  
  
  
  sumbpi = rowSums(ZZip) * log(1-p0)
  rk = rgamma(K, rk_a + Lk)/( rk_b - sumbpi); # from code
  #rk = rgamma(3,1+Lk,1-log(.5)) # from paper
  #rk=rgamma(3, 1+Lk, 1-sumbpi)
  
  
  Theta = matrix(NA,dim(x_kn)[1],dim(x_kn)[2])
  for(i in 1:dim(x_kn)[1]){
    for(j in 1:dim(x_kn)[2]){
      #theta[i,j] = rgamma(1,shape=rk[i]*ZZip[i,j]+ x_kn[i,j])*p0
      #theta[i,j] = rgamma(1,shape=rk[i]*ZZip[i,j]+ x_kn[i,j])*p0
      Theta[i,j] = rgamma(1,rk[i]*ZZip[i,j]+ x_kn[i,j], p0)
    }
  }
  
  
  # Theta<-apply(Theta, 2, norm)
  
  # for (i in 1:100){
  #  theta[,11*i] <- norm(theta[,11*i])
  # }
  
  Phi = matrix(NA,dim(x_pk)[2],dim(x_pk)[2])
  for (i in 1:dim(x_pk)[2]){
    Phi[,i] <-  rdirichlet(1,(alpha_psi+C_kk1)[,i])
  }
  
  Lk = crt(C_k1n,sk)
  #sumbpi = ZZip2%*%log(1-t(p1m)); # !!!!!! changed from original for p1m
  sumbpi = rowSums(ZZip) * log(1-p0)
  sk = rgamma(K, sk_a + Lk)/( 1/sk_b - sumbpi); # from code
  
  W = matrix(NA,nrow(ZZip),ncol(ZZip))
  for(i in 1:nrow(ZZip)){
    for(j in 1:ncol(ZZip)){
      W[i,j] = rgamma(1,sk[i]*ZZip[i,j] + C_k1n[i,j],0.5) # adjust dimension
    }
  }
  #W<-apply(W, 2, norm)
  
  #Pi_k = rbeta(nrow(z0),a0 + rowSums(z0),b0 + numSample - rowSums(z0))
  
  # tn_ix = seq(numTime+1,numTotal,numTime+1)
  
  # mat_mult = Psi%*%theta[,tn_ix]
  #  XPred = matrix(NA,nrow(mat_mult),ncol(mat_mult))
  #  for(i in 1:nrow(mat_mult)){
  #    for(j in 1:ncol(mat_mult)){
  #      XPred[i,j] = rpois(1,mat_mult[i,j])
  #    }
  #  }
  #  Xmtot[,tn_ix] = XPred
  
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
  
  print(b)
}
t2 = proc.time() - t1
t2

#head(fits)

#plot(fits[,2],type="l")

#rett[[10]]$Xmtot

#View(rett[[30]]$Xmtot)


psi.array <- array(NA,c(dim(Psi),(niter-2000)))

for(i in 1:(niter-2000)){
  psi.array[,,i] <- rett[[i+2000]]$Psi
}


mean.psi <- matrix(0,dim(Psi)[1],dim(Psi)[2])

for(i in 1:(niter-2000)){
  mean.psi <- mean.psi + psi.array[,,i]
}
psi_est = round(mean.psi/(niter-2000),2)


Phi.array <- array(NA,c(dim(Phi),(niter-2000)))
for(i in 1:(niter-2000)){
  Phi.array[,,i] <- rett[[2000+i]]$Phi
}
mean.Phi <- matrix(0,dim(Phi)[1],dim(Phi)[2])

for(i in 1:(niter-2000)){
  mean.Phi <- mean.Phi + Phi.array[,,i]
}

phi_est = round(mean.Phi/(niter-2000),3)


Theta.array <- array(NA,c(dim(Theta),(niter-2000)))
for(i in 1:(niter-2000)){
  Theta.array[,,i] <- rett[[2000+i]]$Theta
}
mean.Theta <- matrix(0,dim(Theta)[1],dim(Theta)[2])

for(i in 1:(niter-2000)){
  mean.Theta <- mean.Theta + Theta.array[,,i]
}

theta_est = round(mean.Theta/(niter-2000),3)

W.array <- array(NA,c(dim(rett[[5000]]$W),(niter-2000)))
for(i in 1:(niter-2000)){
  W.array[,,i] <- rett[[2000+i]]$W
}
mean.W <- matrix(0,dim(W)[1],dim(W)[2])

for(i in 1:(niter-2000)){
  mean.W <- mean.W + W.array[,,i]
}

W_est = round(mean.W/(niter-2000),3)

ZZip.array <- array(NA,c(dim(rett[[5000]]$ZZip),(niter-2000)))
for(i in 1:(niter-2000)){
  ZZip.array[,,i] <- rett[[2000+i]]$ZZip
}
mean.ZZip <- matrix(0,dim(ZZip)[1],dim(ZZip)[2])

for(i in 1:(niter-2000)){
  mean.ZZip <- mean.ZZip + ZZip.array[,,i]
}

ZZip_est = round(mean.ZZip/(niter-2000),3)
ZZip_rowsums = rowSums(ZZip_est)

return<- list(ret,psi_est,phi_est,theta_est,W_est,ZZip_est,ZZip_rowsums)
setwd("/afs/crc.nd.edu/user/s/sshao2/Private/DPFA/RDS_output/")
saveRDS(return,"job14.rds")
