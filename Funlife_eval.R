# Funlife model and FLM model estimation evaluation

funl_res=GetFunLEst(respon_Y,fact_num,cov_Z,pena_mat)
funl_theta=funl_res$iter_theta
funl_F=funl_res$iter_F
funl_lambda=funl_res$iter_lambda

funl_alpha=funl_theta[1:dim(X)[2]]
funl_gamma=funl_theta[-(1:dim(X)[2]))]

funlbeta_gamma=matrix(funl_gamma,nrow=dim(X)[2],byrow=T)
funl_beta=funlbeta_gamma%*%beta_estbasis

flm_res=GetFLMEst(respon_Y,cov_W,cov_B,pena_mat)
flm_gamma=flm_res$flm_gamma

flmbeta_gamma=matrix(flm_gamma,nrow=dim(X)[2],byrow=T)
flm_beta=flmbeta_gamma%*%beta_estbasis

funl_beta=flm_beta=funl_srm=flm_srm=funl_ssd=flm_ssd=NULL
funl_Fsrm=funl_lambdasrm=NULL
for(sn in 1:sim_num){
funl_res=GetFunLEst(respon_Y,fact_num,cov_Z,pena_mat)
funl_theta=funl_res$iter_theta
funl_gamma=funl_theta[-(1:dim(X)[2]))]
funlbeta_gamma=matrix(funl_gamma,nrow=dim(X)[2],byrow=T)
funl_beta[sn,]=funlbeta_gamma%*%beta_estbasis

funl_Fest=funl_res$iter_F
funl_lambdaest=funl_res$iter_lambda

funl_Fsrm[i]=apply((funl_Fest-funl_F)^2,2,mean)
funl_lambdasrm[i]=apply((funl_lambdaest-funl_lambda)^2),2,mean)

flm_gamma=GetFLMEst(respon_Y,cov_W,cov_B,pena_mat)$flm_gamma
flmbeta_gamma=matrix(flm_gamma,nrow=dim(X)[2],byrow=T)
flm_beta[sn,]=flmbeta_gamma%*%beta_estbasis

funl_srm[i]=sum((funl_beta[i,]-beta)^2)/length(s_grid)
flm_srm[i]=sum((flm_beta[i,]-beta)^2)/length(s_grid)
}
funl_mest=apply(funl_beta,2,mean)
flm_mest=apply(flm_beta,2,mean)

for(s in 1:sim_num){
funl_ssd[i]=sum((funl_beta[i,]-funl_mest)^2)/length(s_grid)
flm_ssd[i]=sum((flm_beta[i,]-flm_mest)^2)/length(s_grid)
}

funl_rm=sqrt(mean(funl_srm))
flm_rm=sqrt(mean(fm_srm))
funl_sd=sqrt(mean(funl_ssd))
flm_sd=sqrt(mean(flm_ssd))
Funl_Frm=sqrt(funl_Fsrm)
Funl_lambdarm=sqrt(funl_lambdasrm)




