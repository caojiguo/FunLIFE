# calculate the FLM estimates

GetFLMEst=function(respon_Y,cov_W,cov_B,pena_mat){
cov_Z=covZ_mat=covZY_mat=NULL
for (i in 1:N){
cov_Z[,,i]=cbind(cov_W[,,i],cov_B[,,i])
covZ_mat=covZ_mat+t(cov_Z[,,i])%*%cov_Z[,,i]
covZY_mat=covZY_mat+t(cov_Z[,,i])%*%respon_Y[,i]
}

flm_theta=ginv(covZ_mat+pena_mat)%*%covZY_mat

flm_alpha=flm_theta[1:dim(X)[2]]
flm_gamma=flm_theta[-(1:dim(X)[2]))]

beta_gamma=matrix(flm_gamma,nrow=dim(X)[2],byrow=T)
flm_beta=beta_gamma%*%beta_estbasis

out=list(flm_beta,flm_gamma,flm_theta)
return(out)
}


