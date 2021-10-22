fit_DIC <- function(rett, K, N, numTime, M){
  theta_sample_keep_SJ<-lapply(rett,'[[', 1)
  psi_sample_keep_SJ <-lapply(rett,'[[', 3)
  theta_sample_average_SJ<-apply(simplify2array(theta_sample_keep_SJ), 1:2, mean)
  psi_sample_average_SJ<-apply(simplify2array(psi_sample_keep_SJ), 1:2, mean)
  
  
  niter_SJ <- length(theta_sample)
 #  k_SJ = 3
 #  n_SJ = 100
 #  Time_SJ = 10
 #  m_SJ = 15
  m_SJ = rep(m_SJ,n_SJ*Time_SJ)
  z_sample_DIC_SJ = Xmtot
  text_token_SJ = matrix(Xmtot,ncol = 1, byrow = F)
  
  a = vector()
  for(i in 1:n_SJ){
    for (t in 1:Time_SJ){
      a=c(a,paste0('id',i,'time',t))
    }
  }
  text_token_SJ = cbind(rep(a, each = 15),text_token_SJ)
  
 
  
  z_sample_DIC_SJ = NULL
  count = 1
  for(i in 1:n_SJ){
    for (t in 1:Time_SJ){
      while (count <n_SJ*Time_SJ){
      w_i_SJ <- text_token_SJ[text_token_SJ[,1]  == paste0('id',i,'time',t),]
      theta_i_t_SJ <- theta_sample_average_SJ[,count]
      count= count+1
      
      }
    for(j in 1:(m_SJ[count])){
      psi_j_SJ <- t(psi_sample_average_SJ)[,j]
      z_sample_DIC_SJ <- c(z_sample_DIC_SJ,log(sum(theta_i_t_SJ*psi_j_SJ)))
      
   }
  
    }
  }
  d_inbar <- sum(z_sample_DIC_SJ)
  
  z_sample_keep_DIC = list()
  count = 1
  for (s in 1:niter){
    z_sample_keep_DIC[[s]] = list()
    for(i in 1:n_SJ){
      z_sample_keep_DIC[[s]][[i]] = list()
      for (t in 1:Time_SJ){
      while (count <n_SJ*Time_SJ){
        w_i_SJ <- text_token_SJ[text_token_SJ[,1]  == paste0('id',i,'time',t),]
        theta_i_t_SJ <- theta_sample_keep_SJ[[s]][,count]
        count= count+1
      }
      z_sample_keep_DIC[[s]][[i]][[t]] = list()
      for(j in 1:(m_SJ[count])){
        psi_j_SJ <- t(psi_sample_keep_SJ[[s]])[,j]
        z_sample_keep_DIC[[s]][[i]][[t]][j]<- log(sum(theta_i_t_SJ*psi_j_SJ))
      }
      sum(z_sample_keep_DIC[[s]][[i]][[t]])
    }
  }
  }
  
  count = 1
  z_sample = list()
  for (i in 1:n_SJ){
    for (t in 1:Time_SJ){
      while (count <n_SJ*Time_SJ){
     z_sample[[count]] =lapply(z_sample_keep_DIC[[i]][[t]], function(x) sum(unlist(x)))
     count= count+1
     print(count)
      }
    }
  }
  
  d_outbar=mean(unlist(z_sample))
  p_d = 2*(d_inbar - d_outbar )
  
  DIC <- (-2)*d_inbar + 2*p_d
  
  DIC
}