# choose the smoothing parameter for penalty

para_thi=para_SSE=NULL

for(j in 1:length(val_grid)){
para_thi[j,]=rep(val_grid[j],length.out=dim(X)[2])
basis_paramat=matrix(diag(para_thi[j,]),ncol=length(dim(X)[2]))
basis_penamat=matrix(basis_pena,ncol=length(basis_pena))
beta_pena=krnocker(basis_paramat,basis_penamat)
alpha_pena=matrix(0,ncol=dim(cov_W)[1],nrow=dim(cov_W)[1])
pena_mat=rbind(cbind(alpha_pena,matrix(0,nrow=nrow(alpha_pena),ncol=ncol(beta_pena))),
cbind(matrix(0,nrow=nrow(beta_pena),ncol=ncol(alpha_pena)),beta_pena))

para_res=GetFunLEst(respon_Y,fact_num,cov_Z,pena_mat)
para_theta=para_res$iter_theta
para_F=para_res$iter_F
para_lambda=para_res$iter_lambda

para_resYthe=para_ss=NULL
for (i in 1:N){
para_resYthe[,i]=respon_Y[,i]-cov_Z[,,i]%*%para_theta
para_ss=para_ss+sum(para_resYthe[,i]^2)
}

para_SSE[j]=para_ss
para_Fpro=martix(diag(rep(1,T)))-para_F%*%t(para_F)/T
para_Z=NULL
for(i in 1:N){
para_Z=rbind(para_Z,cov_Z[,,i])
}
para_Zpro=para_Fpro%*%para_Z
para_sthi=para_Z%*%ginv(t(para_Zpro)%*%para_Zpro+pena_mat)%*%t(para_Zpro)
para_s=martix(diag(rep(1,N*T)))-para_sthi
para_gcv[j]=para_SSE[j]/sum(diag(para_s))
}

para_opt=val_grid[which.min(para_gcv)]


