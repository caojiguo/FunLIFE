# calculate the Panel data model estimates

pdm_beta=mean(funl_beta)
for(i in 1:N){
for(j in 1:T){
for(p in 1:dim(cov_indX([2]){
cov_X[i,p,j]=mean(cov_indX[j,p,,i])
}}}

error=matrix(rnrom(N*T,0,1),ncol=T)
for(i in 1:N){
for(j in 1:T){
respon_Y[i,j]=t(cov_W[j,,i])%*%alpha+t(cov_X[j,,i])%*%pdm_beta+t(lambda[i,])%*%F[,j]
+error[i,j]
}
}

GetPDMEst=function(respon_Y,cov_W,cov_X){

init_Vval=eigen(respon_Y %*% t(respon_Y))$val
init_Vmat=matrix(diag(init_Vval[1:fact_num]),ncol=fact_num)
init_Fmat=eigen(respon_Y %*% t(respon_Y))$vec
init_F=matrix(init_Fmat[,1:fact_num])
init_lambda=t(t(init_F)%*%respon_Y/N)

resi_Y=NULL
covZ_mat=covZY_mat=NULL
for (i in 1:N){
resi_Y[,i]=respon_Y[,i]-init_F%*%matrix(init_lambda[i,])
covZ_mat=covZ_mat+t(cov_Z[,,i])%*%cov_Z[,,i]
covZY_mat=covZY_mat+t(cov_Z[,,i])%*%resi_Y[,i]
}

init_theta=ginv(covZ_mat)%*%covZY_mat

iter_theta=init_theta
iter_num=1

repeat{

iterep_theta=iter_theta
iter_resY=NULL
for(i in 1:N){
iter_resYF[,i]=respon_Y[,i]-cov_Z[,,i]%*%iter_theta
}

iter_Vval=eigen(iter_YF %*% t(resi_YF))$val
iter_Vmat=matrix(diag(iter_Vval[1:fact_num]),ncol=fact_num)
iter_Fmat=eigen(iter_resYF %*% t(iter_resYF)/(NT))$vec
iter_F=matrix(iter_Fmat[,1:fact_num])
iter_lambda=t(t(iter_YF)%*%iter_resYF/N)

iter_resYthe=NULL
itcov_Z=cov_Z
itcovZ_mat=covZ_mat
itcovZY_mat=NULL
for (i in 1:N){
iter_resYthe[,i]=respon_Y[,i]-iter_F%*%matrix(iter_lambda[i,])
itcovZY_mat=itcovZY_mat+t(itcov_Z[,,i])%*%iter_resYthe[,i]
}

iter_theta=ginv(itcovZ_mat+pena_mat)%*%itcovZY_mat

iter_num=iter_num+1
itnor_theta=norm(iter_theta-iterep_theta,2)
if(iter_num>500 | itnor_theta<0.001 )  break

}
 out=list(iter_theta,iter_F,iter_lambda)
 return(out)
}

pdm_res=GetPDMEst(respon_Y,fact_num,cov_W,cov_X)
pdm_thetaest=pdm_res$iter_theta
pdm_Fest=pdm_res$iter_F
pdm_lambdaest=pdm_res$iter_lambda

pdm_alphaest=funl_theta[1:dim(cov_X)[2]]
pdm_betaest=funl_theta[-(1:dim(cov_X)[2]))]

for(i in 1:N){
for(j in 1:T){
pdm_resest=t(cov_W[j,,i])%*%pdm_alphaest+t(cov_X[j,,i])%*%pdm_betaest+
t(lambdaest[i,])%*%pdm_Fest[,j]

pdm_resrm=sqrt(sum(apply((pdm_resest-respon_Y)^2,1,mean)))

