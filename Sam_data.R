# generate the data sample

s_grid=seq(0,1,length.out = 500)

beta_1=2+3*s+exp(2*s)
beta_2=5+3*sin(2*pi*s)+2*cos(2*pi*s)
beta=rbind(beta_1,beta_2)

alpha=c(1,0.5)
sigma=matrix(c(1,0,0,1),ncol=2)
lambda=matrix(mvrnorm(N*fact_num,c(0,0),sigma),nrow=N)
F=matrix(mvrnorm(T*fact_num,c(0,0),0.5*sigma),nrow=T)

cov_Ws1=matrix(rexp(N*T,2),ncol=T)
cov_Ws2=matrix(runif(N*T,0,1),ncol=T)

delta_1=matrix(runif(N*T,-1,1),ncol=T)
delta_2=matrix(rnorm(N*T,0,1),ncol=T)

cov_Xs1=cov_Xb1=cov_Xs2=cov_Xb2=cov_W1=NULL
for(i in 1:N){ 
for(j in 1:T){
cov_Xs1[j,,i]=1+c1*t(lambda[i,])%*%F[,j]+delta_1[i,j]*s_grid
cov_Xs2[j,,i]=c2*t(lambda[i,])%*%F[j,]+delta_2[i,j]*sin(2*pi*s_grid)
cov_Xb1[i,j]=t(cov_Xs1[j,,i])%*%beta_1/length(s_grid)
cov_Xb2[i,j]=t(cov_Xs2[j,,i])%*%beta_2/length(s_grid)
cov_W1[i,j]=cov_Ws1+c1*t(lambda[i,])%*%F[,j]
}
}

for(i in 1:N){
cov_X[,,i]=cbind(cov_Xb1[i,],cov_Xb2[i,])
cov_W[,,i]=cbind(cov_Ws1[i,],cov_W2[i,])
}

error=matrix(rnrom(N*T,0,1),ncol=T)
for(i in 1:N){
for(j in 1:T){
respon_Y[i,j]=t(cov_W[j,,i])%*%alpha+t(cov_X[j,,i])%*%c(1,1)+t(lambda[i,])%*%F[,j]
+error[i,j]
}
}

data=respon_Y


