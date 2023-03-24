
#Libraries:
library(MASS)
library(dr)
library(car)
library(olsrr)

#Simulation Study:

#Data 1:
set.seed(1234)
n=500
p=70
Sigma=matrix(rep(0,p^2),nrow=p,ncol=p)
for(i in 1:p){
  for(j in 1:p){
    Sigma[i,j]=0.8^abs(i-j)
  }
}
X=mvrnorm(n,rep(10,p),Sigma)
y=(X[,1]+X[,2]+X[,3])/sqrt(3)+0.1*rnorm(n)
df=data.frame(y,X)
x=rep(0,n-1)
sir=dr(formula = y~ X, data = df, slice.function = dr.slices.arc,
       nslices = 10, chi2approx = "wood", numdir = 8, method = "sir")
dr.basis(sir)
org_dir1=c(rep(1,3),rep(0,p-3))
p1=X%*%org_dir1
sqrt(summary(lm(X%*%dr.basis(sir)[,1]~p1))$r.squared)


#standard Error of Estimate
est=rep(0,200)
for(i in 1:200){
  set.seed(round(rexp(1,500)*10000000))
  n=500
  p=70
  Sigma=matrix(rep(0,p^2),nrow=p,ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      Sigma[i,j]=0.8^abs(i-j)
    }
  }
  X=mvrnorm(n,rep(10,p),Sigma)
  y=(X[,1]+X[,2]+X[,3])/sqrt(3)+0.1*rnorm(n)
  df=data.frame(y,X)
  x=rep(0,n-1)
  sir=dr(formula = y~ X, data = df, slice.function = dr.slices.arc,
         nslices = 10, chi2approx = "wood", numdir = 4, method = "sir")
  dr.basis(sir)
  org_dir1=c(rep(1,3),rep(0,p-3))
  p1=X%*%org_dir1
  est[i]=sqrt(summary(lm(X%*%dr.basis(sir)[,1]~p1))$r.squared)
}
sd(est)

#Choice of effective dimension Size(k):
for(i in 1:p){
  k_eff=i
  if(mean(sir$evalues[(i+1):p])<=qchisq(0.95,p-i,10-i-1)/(n*(p-i))){
    break
  }
}
k_eff
plot(1:p,cumsum(sir$evalues/sum(sir$evalues)))


#Number of Slices
for(h in 2:n){
  sir1=dr(formula = y~ X, data = df, slice.function = dr.slices.arc,
          nslices = h, chi2approx = "wood", numdir = 4, method = "sir")
  x[h]=sqrt(summary(lm(X%*%dr.basis(sir1)[,1]~p1))$r.squared)
}
plot(x,main="Model 1",xlab="Number of Slices",ylab="R_squared(e.d.r 1)")








#Data 2:
set.seed(1234)
n=500
p=70
Sigma=matrix(rep(0,p^2),nrow=p,ncol=p)
for(i in 1:p){
  for(j in 1:p){
    Sigma[i,j]=0.8^abs(i-j)
  }
}
X=mvrnorm(n,rep(10,p),Sigma)
y=1+exp((X[,1]+X[,2]+X[,3])/sqrt(3))+rnorm(n)
df=data.frame(y,X)
x=rep(0,n-1)
sir=dr(formula = y~ X, data = df, slice.function = dr.slices.arc,
       nslices = 10, chi2approx = "wood", numdir = 4, method = "sir")
dr.basis(sir)
org_dir1=c(rep(1,3),rep(0,p-3))
p1=X%*%org_dir1
sqrt(summary(lm(X%*%dr.basis(sir)[,1]~p1))$r.squared)

#Standard Error of Estimate
est=rep(0,200)
for(i in 1:200){
  set.seed(round(rexp(1,500)*10000000))
  n=500
  p=70
  Sigma=matrix(rep(0,p^2),nrow=p,ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      Sigma[i,j]=0.8^abs(i-j)
    }
  }
  X=mvrnorm(n,rep(10,p),Sigma)
  y=1+exp((X[,1]+X[,2]+X[,3])/sqrt(3))+rnorm(n)
  df=data.frame(y,X)
  x=rep(0,n-1)
  sir=dr(formula = y~ X, data = df, slice.function = dr.slices.arc,
         nslices = 10, chi2approx = "wood", numdir = 4, method = "sir")
  dr.basis(sir)
  org_dir1=c(rep(1,3),rep(0,p-3))
  p1=X%*%org_dir1
  est[i]=sqrt(summary(lm(X%*%dr.basis(sir)[,1]~p1))$r.squared)
}
sd(est)

#Choice of effective dimension Size(k):
for(i in 1:p){
  k_eff=i
  if(mean(sir$evalues[(i+1):p])<=qchisq(0.95,p-i,10-i-1)/(n*(p-i))){
    break
  }
}
k_eff
plot(1:p,cumsum(sir$evalues/sum(sir$evalues)))


#Number of Slices
for(h in 2:n){
  sir1=dr(formula = y~ X, data = df, slice.function = dr.slices.arc,
          nslices = h, chi2approx = "wood", numdir = 4, method = "sir")
  x[h]=sqrt(summary(lm(X%*%dr.basis(sir1)[,1]~p1))$r.squared)
}
plot(x,main="Model 2",xlab="Number of Slices",ylab="R_squared(e.d.r 1)")







#Data 3:
set.seed(1234)
n=500
p=70
Sigma=matrix(rep(0,p^2),nrow=p,ncol=p)
for(i in 1:p){
  for(j in 1:p){
    Sigma[i,j]=0.8^abs(i-j)
  }
}
X=mvrnorm(n,rep(2,p),Sigma)
y=(X[,1]+X[,2]+X[,3])/(0.5+(X[,4]+X[,5]+1.5)^2)+0.1*rnorm(n)
df=data.frame(y,X)
x1=rep(0,n-1)
x2=rep(0,n-1)
sir=dr(formula = y~ X, data = df, slice.function = dr.slices.arc,
       nslices = 10, chi2approx = "wood", numdir = 4, method = "sir")
dr.basis(sir)
org_dir1=c(rep(1,3),rep(0,p-3))
org_dir2=c(rep(0,3),1,1,rep(0,p-5))
p1=X%*%org_dir1
p2=X%*%org_dir2
sqrt(summary(lm(X%*%dr.basis(sir)[,1]~p1+p2))$r.squared)

#Standard Error of Estimate
est=rep(0,200)
for(i in 1:200){
  set.seed(round(rexp(1,500)*10000000))
  n=500
  p=70
  Sigma=matrix(rep(0,p^2),nrow=p,ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      Sigma[i,j]=0.8^abs(i-j)
    }
  }
  X=mvrnorm(n,rep(2,p),Sigma)
  y=(X[,1]+X[,2]+X[,3])/(0.5+(X[,4]+X[,5]+1.5)^2)+0.1*rnorm(n)
  df=data.frame(y,X)
  x1=rep(0,n-1)
  x2=rep(0,n-1)
  sir=dr(formula = y~ X, data = df, slice.function = dr.slices.arc,
         nslices = 10, chi2approx = "wood", numdir = 4, method = "sir")
  dr.basis(sir)
  org_dir1=c(rep(1,3),rep(0,p-3))
  org_dir2=c(rep(0,3),1,1,rep(0,p-5))
  p1=X%*%org_dir1
  p2=X%*%org_dir2
  est[i]=sqrt(summary(lm(X%*%dr.basis(sir)[,1]~p1+p2))$r.squared)
}
sd(est)

#Choice of effective dimension Size(k):
for(i in 1:p){
  k_eff=i
  if(mean(sir$evalues[(i+1):p])<=qchisq(0.95,p-i,10-i-1)/(n*(p-i))){
    break
  }
}
k_eff
plot(1:p,cumsum(sir$evalues/sum(sir$evalues)))

#Number of Slices
for(h in 2:n){
  sir1=dr(formula = y~ X, data = df, slice.function = dr.slices.arc,
          nslices = h, chi2approx = "wood", numdir = 4, method = "sir")
  x1[h]=sqrt(summary(lm(X%*%dr.basis(sir1)[,1]~p1+p2))$r.squared)
  x2[h]=sqrt(summary(lm(X%*%dr.basis(sir1)[,2]~p1+p2))$r.squared)
}
plot(x1,main="Model 3",xlab="Number of Slices",ylab="R_squared(e.d.r 1)")





#Data 4:
set.seed(1234)
n=500
p=70
Sigma=matrix(rep(0,p^2),nrow=p,ncol=p)
for(i in 1:p){
  for(j in 1:p){
    Sigma[i,j]=0.8^abs(i-j)
  }
}
X=mvrnorm(n,rep(15,p),Sigma)
y=X[,1]*(X[,1]+X[,2])+rnorm(n)
df=data.frame(y,X)
x1=rep(0,n-1)
x2=rep(0,n-1)
sir=dr(formula = y~ X, data = df, slice.function = dr.slices.arc,
       nslices = 10, chi2approx = "wood", numdir = 4, method = "sir")
dr.basis(sir)
org_dir1=c(1,rep(0,p-1))
org_dir2=c(1,1,rep(0,p-2))
p1=X%*%org_dir1
p2=X%*%org_dir2
sqrt(summary(lm(X%*%dr.basis(sir)[,1]~p1+p2))$r.squared)

#Standard Error of Estimate
est=rep(0,200)
for(i in 1:200){
  set.seed(round(rexp(1,500)*10000000))
  n=500
  p=70
  Sigma=matrix(rep(0,p^2),nrow=p,ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      Sigma[i,j]=0.8^abs(i-j)
    }
  }
  X=mvrnorm(n,rep(15,p),Sigma)
  y=X[,1]*(X[,1]+X[,2])+rnorm(n)
  df=data.frame(y,X)
  x1=rep(0,n-1)
  x2=rep(0,n-1)
  sir=dr(formula = y~ X, data = df, slice.function = dr.slices.arc,
         nslices = 10, chi2approx = "wood", numdir = 4, method = "sir")
  dr.basis(sir)
  org_dir1=c(1,rep(0,p-1))
  org_dir2=c(1,1,rep(0,p-2))
  p1=X%*%org_dir1
  p2=X%*%org_dir2
  est[i]=sqrt(summary(lm(X%*%dr.basis(sir)[,1]~p1+p2))$r.squared)
}
sd(est)

#Choice of effective dimension Size(k):
for(i in 1:p){
  k_eff=i
  if(mean(sir$evalues[(i+1):p])<=qchisq(0.95,p-i,10-i-1)/(n*(p-i))){
    break
  }
}
k_eff
plot(1:p,cumsum(sir$evalues/sum(sir$evalues)))


#Number of Slices
for(h in 2:n){
  sir1=dr(formula = y~ X, data = df, slice.function = dr.slices.arc,
          nslices = h, chi2approx = "wood", numdir = 4, method = "sir")
  x1[h]=sqrt(summary(lm(X%*%dr.basis(sir1)[,1]~p1+p2))$r.squared)
  x2[h]=sqrt(summary(lm(X%*%dr.basis(sir1)[,2]~p1+p2))$r.squared)
}
plot(x1,main="Model 4",xlab="Number of Slices",ylab="R_squared(e.d.r 1)")









#PCA VS SIR

#Data
set.seed(1234)
n=1000
x1=rnorm(n)
x2=rnorm(n)
X=cbind(x1,x2)
y=sin(0.7*x1-0.7*x2)+0.1*rnorm(n)
actual_dir=c(1,-1)
pca=prcomp(X,scale=T,center=T)
cumsum(pca$sdev/sum(pca$sdev))
pca$rotation[,1]
sir=dr(formula = y ~ X, slice.function = dr.slices.arc,
       nslices = 10, chi2approx = "wood", method = "sir")
dr.basis(sir)
plot(x1,x2)
abline(h=0,v=0)
par(mfrow=c(1,2))
plot(X%*%pca$rotation[,1],y,main="PCA",xlab="Projection along 1st direction",ylab="y")
plot(X%*%dr.basis(sir)[,1],y,main="SIR",xlab="Projection along 1st direction",ylab="y")
theta_pca=as.numeric()
theta_sir=as.numeric()
for(i in 1:1000){
  set.seed(round(rexp(1,500)*10000000))
  n=1000
  x1=rnorm(n)
  x2=rnorm(n)
  X=cbind(x1,x2)
  y=sin(0.7*x1-0.7*x2)+0.1*rnorm(n)
  actual_dir=c(1,-1)/sqrt(2)
  pca=prcomp(X,scale=T,center=T)
  pca_dir=pca$rotation[,1]/ sqrt(sum(pca$rotation[,1]^2))
  sir=dr(formula = y ~ X, slice.function = dr.slices.arc,
         nslices = 10, chi2approx = "wood", method = "sir")
  sir_dir=dr.basis(sir)[,1]/ sqrt(sum(dr.basis(sir)[,1]^2))
  theta_pca[i]=acos(sum(actual_dir*pca_dir))
  theta_sir[i]=acos(sum(actual_dir*sir_dir))
}
par(mfrow=c(1,2))
hist(theta_sir,main="SIR",xlab="Angle(in Radians)",ylab="Count")
hist(theta_pca,main="PCA",xlab="Angle(in Radians)",ylab="Count")
par(mfrow=c(1,1))













#Analysis:

#Only Training Set Model
proj=function(t)
{
data=read.csv("NGDP.csv")
forecast=t(as.numeric(data$NGDP2))
fgdp=matrix(forecast,,nrow=40,ncol=10,byrow=T)
agdp=as.array(unique(data$NGDP1))
#Stationary Transformation
par(mfrow=c(5,2))
for(i in 1:ncol(fgdp)){
  # plot(fgdp[,i])
}
par(mfrow=c(1,1))
for(i in 1:nrow(fgdp)){
  fgdp[i,]=fgdp[i,]-agdp[i]
}
par(mfrow=c(5,2))
for(i in 1:ncol(fgdp)){
  # plot(fgdp[,i])
}
par(mfrow=c(1,1))

agdp=diff(agdp)
fgdp=fgdp[-40,]#2006 Q1-2015Q3

#SIR
sir=dr(formula =agdp~fgdp, slice.function = dr.slices.arc,
       nslices = 5,numdir=10, chi2approx = "wood", method = "sir")
dr.basis(sir)


#Choice of effective dimension Size(k):
p=ncol(fgdp)
n=nrow(fgdp)
for(i in 1:p){
  k_eff=i
  if(mean(sir$evalues[(i+1):p])<=qchisq(0.95,p-i,10-i-1)/(n*(p-i))){
    break
  }
}
k_eff
plot(1:p,cumsum(sir$evalues/sum(sir$evalues)))


proj_gdp=fgdp%*%dr.basis(sir)[,1:t]
mod=lm(agdp~0+proj_gdp)
summary(mod)
#plot(mod)

#PCA
pca=prcomp(fgdp,scale=T,center=T)
plot(1:p,cumsum(pca$sdev/sum(pca$sdev)))
proj_gdp_pca=fgdp%*%pca$rotation[,1:t]
mod_pca=lm(agdp~0+proj_gdp_pca)
summary(mod_pca)
#plot(mod_pca)


#Benchmark Model

data1=read.csv("NGDP.csv")
forecast=t(as.numeric(data$NGDP2))
fgdp1=matrix(forecast,,nrow=40,ncol=10,byrow=T)
agdp1=as.array(unique(data$NGDP1))
benchmark_forecast=as.numeric()
for(i in 1:nrow(fgdp1)){
  benchmark_forecast[i]=mean(fgdp1[i,])
}
benchmark_forecast=diff(benchmark_forecast)
benchmark_res=diff(agdp1)-benchmark_forecast

RMSE_benchmark=sqrt(mean(benchmark_res^2))
RMSE_model_sir=sqrt(mean(mod$residuals^2))
RMSE_model_pca=sqrt(mean(mod_pca$residuals^2))
out=c(RMSE_model_sir/RMSE_benchmark,RMSE_model_pca/RMSE_benchmark)
return(out)
}

sir=as.numeric()
pc=as.numeric()
for(i in 1:10){
  sir[i]=proj(i)[1]
  pc[i]=proj(i)[2]
}
  

plot(1:10,sir,type="l",col="BLUE",xlab="Number of Directions",ylab="RMSE Ratio")
lines(pc,col="RED",type="l")
legend("topright", legend=c("SIR", "PCA"), fill = c("blue","red"))




























#One Step forecast Model: Including one pt every iteration in training set
proj=function(t)
{
data=read.csv("NGDP.csv")
forecast=t(as.numeric(data$NGDP2))
fgdp=matrix(forecast,nrow=40,ncol=10,byrow=T)
agdp=as.array(unique(data$NGDP1))
#Stationary Transformation
par(mfrow=c(5,2))
for(i in 1:ncol(fgdp)){
  # plot(fgdp[,i])
}
par(mfrow=c(1,1))

for(i in 1:nrow(fgdp)){
  fgdp[i,]=fgdp[i,]-agdp[i]
}
par(mfrow=c(5,2))
for(i in 1:ncol(fgdp)){
  # plot(fgdp[,i])
}
par(mfrow=c(1,1))

agdp=diff(agdp)
fgdp=fgdp[-40,]#2006 Q1-2015Q3

#Train-Test Split
fgdp_train=fgdp[1:20,]#Pseudo differenced forecasts from 10 individuals from Q1,2006-Q4,2010
agdp_train=agdp[1:20]#Actual differneced GDP values from Q1,2006-Q4,2010
fgdp_test=fgdp[21:39,]#Pseudo differenced forecasts from 10 individuals from Q1,2011-Q3,2015
agdp_test=agdp[21:39]#Pseudo differenced GDP values from Q1,2011-Q3,2015

forecast_residuals_sir=as.numeric()
forecast_residuals_pca=as.numeric()
forecast_residuals_benchmark=as.numeric()

for(i in 21:39)
{
#SIR
sir=dr(formula =agdp_train~fgdp_train, slice.function = dr.slices.arc,
       nslices = 5,numdir=10, chi2approx = "wood", method = "sir")
proj_gdp=fgdp_train%*%dr.basis(sir)[,1:t]
mod=lm(agdp_train~0+proj_gdp)
summary(mod)

#PCA
pca=prcomp(fgdp_train,scale=T,center=T)
proj_gdp_pca=fgdp_train%*%pca$rotation[,1:t]
mod_pca=lm(agdp_train~0+proj_gdp_pca)
summary(mod_pca)


#Performance on Test Set

#SIR
proj_gdp_test_sir=fgdp[i,]%*%dr.basis(sir)[,1:t]
predict_sir=proj_gdp_test_sir%*%mod$coefficients
forecast_residuals_sir[i-20]=agdp[i]-predict_sir


#PCA
proj_gdp_test_pca=fgdp[i,]%*%pca$rotation[,1:t]
predict_pca=proj_gdp_test_pca%*%mod_pca$coefficients
forecast_residuals_pca[i-20]=agdp[i]-predict_pca


#Updation of Training Set including another time point
if(i<39)
{fgdp_train=rbind(fgdp_train,fgdp[i,])
agdp_train=c(agdp_train,agdp[i])}else{
    break
  }
}

#Benchmark Model 1


data1=read.csv("NGDP.csv")
forecast1=t(as.numeric(data1$NGDP2))
fgdp1=matrix(forecast1,nrow=40,ncol=10,byrow=T)
agdp1=as.array(unique(data1$NGDP1))
benchmark_forecast1=as.numeric()
for(i in 1:nrow(fgdp1)){
  benchmark_forecast1[i]=mean(fgdp1[i,])
}
benchmark_forecast1=diff(benchmark_forecast1)
benchmark_res1=diff(agdp1)-benchmark_forecast1
forecast_residuals_benchmark1=benchmark_res1[21:39]


#Benchmark Model 2:
data1=read.csv("NGDP.csv")
forecast1=t(as.numeric(data1$NGDP2))
fgdp1=matrix(forecast1,nrow=40,ncol=10,byrow=T)
agdp1=as.array(unique(data1$NGDP1))
agdp_train1=agdp1[2:21]
fgdp_train1=fgdp1[1:20,]
fgdp_test1=fgdp1[21:40,]
benchmark_forecast2=as.numeric()
for(k in 21:40){
sMAPE_inter=matrix(rep(0,nrow(fgdp_train1)*ncol(fgdp_train1)),nrow=nrow(fgdp_train1),ncol=ncol(fgdp_train1))
for(i in 1:nrow(fgdp_train1)){
  for(j in 1:ncol(fgdp_train1)){
    sMAPE_inter[i,j]=abs(agdp_train1[i]-fgdp_train1[i,j])/abs(agdp_train1[i]+fgdp_train1[i,j])
  }
}
sMAPE=as.numeric()
wts=as.numeric()
for(j in 1:ncol(fgdp_train1)){
  sMAPE[j]=mean(sMAPE_inter[,j])
}
wts=(1/sMAPE)/sum(1/sMAPE)
forecast2=weighted.mean(fgdp1[k,],wts)
benchmark_forecast2=c(benchmark_forecast2,forecast2)

#Updation of Training Set including another time point
if(k<=39)
{fgdp_train1=rbind(fgdp_train1,fgdp1[k,])
agdp_train1=c(agdp_train1,agdp1[k+1])}else{
  break
}
}
benchmark_forecast2=diff(benchmark_forecast2)
forecast_residuals_benchmark2=diff(agdp1)[21:39]-benchmark_forecast2



#Comparison:
RMAE_sir=mean(abs(forecast_residuals_sir))
RMAE_pca=mean(abs(forecast_residuals_pca))
RMAE_benchmark1=mean(abs(forecast_residuals_benchmark1))
RMAE_benchmark2=mean(abs(forecast_residuals_benchmark2))
RMSE_sir=sqrt(mean(forecast_residuals_sir^2))
RMSE_pca=sqrt(mean(forecast_residuals_pca^2))
RMSE_benchmark1=sqrt(mean(forecast_residuals_benchmark1^2))
RMSE_benchmark2=sqrt(mean(forecast_residuals_benchmark2^2))
out=c(RMAE_sir,RMAE_pca,RMAE_benchmark1,RMAE_benchmark2,RMSE_sir,RMSE_pca,RMSE_benchmark1,RMSE_benchmark2)
return(out)

}


RMAE_sir=as.numeric()
RMAE_pc=as.numeric()
RMAE_bench1=as.numeric()
RMAE_bench2=as.numeric()
RMSE_sir=as.numeric()
RMSE_pc=as.numeric()
RMSE_bench1=as.numeric()
RMSE_bench2=as.numeric()
for(i in 1:10){
  RMAE_sir[i]=proj(i)[1]
  RMAE_pc[i]=proj(i)[2]
  RMAE_bench1[i]=proj(i)[3]
  RMAE_bench2[i]=proj(i)[4]
  RMSE_sir[i]=proj(i)[5]
  RMSE_pc[i]=proj(i)[6]
  RMSE_bench1[i]=proj(i)[7]
  RMSE_bench2[i]=proj(i)[8]
}


plot(1:10,RMAE_sir,type="l",col="BLUE",ylim=c(0,200),main="Dynamic Model",xlab="Number of Directions",ylab="MAE")
lines(RMAE_pc,col="RED",type="l")
lines(RMAE_bench1,col="YELLOW",type="l")
lines(RMAE_bench2,col="BLACK",type="l")
legend("topright", legend=c("SIR", "PCA","SA","WA"), fill = c("blue","red","yellow","black"))



plot(1:10,RMSE_sir,type="l",col="BLUE",ylim=c(0,300),main="Dynamic Model",xlab="Number of Directions",ylab="RMSE")
lines(RMSE_pc,col="RED",type="l")
lines(RMSE_bench1,col="YELLOW",type="l")
lines(RMSE_bench2,col="BLACK",type="l")
legend("topright", legend=c("SIR", "PCA","SA","WA"), fill = c("blue","red","yellow","black"))



















#Static Model: No updation of training Set

proj=function(t)
{
  data=read.csv("NGDP.csv")
  forecast=t(as.numeric(data$NGDP2))
  fgdp=matrix(forecast,nrow=40,ncol=10,byrow=T)
  agdp=as.array(unique(data$NGDP1))
  plot(agdp,xlab="Quarters(2006-2015)",ylab="NGDP",main="Nominal GDP US(2006-2015)")
  #Stationary Transformation
  par(mfrow=c(5,2))
  for(i in 1:ncol(fgdp)){
    # plot(fgdp[,i])
  }
  par(mfrow=c(1,1))
  
  for(i in 1:nrow(fgdp)){
    fgdp[i,]=fgdp[i,]-agdp[i]
  }
  par(mfrow=c(5,2))
  for(i in 1:ncol(fgdp)){
    # plot(fgdp[,i])
  }
  par(mfrow=c(1,1))
  
  agdp=diff(agdp)
  fgdp=fgdp[-40,]#2006 Q1-2015Q3
  par(mfrow=c(1,2))
  plot(agdp,ylab="Differenced NGDP")
  plot(fgdp[,1],ylab="Pseudo Differenced Forecast")
  par(mfrow=c(1,1))
  #Train-Test Split
  fgdp_train=fgdp[1:20,]#Pseudo differenced forecasts from 10 individuals from Q1,2006-Q4,2010
  agdp_train=agdp[1:20]#Actual differneced GDP values from Q1,2006-Q4,2010
  fgdp_test=fgdp[21:39,]#Pseudo differenced forecasts from 10 individuals from Q1,2011-Q3,2015
  agdp_test=agdp[21:39]#Pseudo differenced GDP values from Q1,2011-Q3,2015
  
  forecast_residuals_sir=as.numeric()
  forecast_residuals_pca=as.numeric()
  forecast_residuals_benchmark=as.numeric()
  
    #SIR
    sir=dr(formula =agdp_train~fgdp_train, slice.function = dr.slices.arc,
           nslices = 5,numdir=10, chi2approx = "wood", method = "sir")
    proj_gdp=fgdp_train%*%dr.basis(sir)[,1:t]
    mod=lm(agdp_train~0+proj_gdp)
    summary(mod)
    
    #PCA
    pca=prcomp(fgdp_train,scale=T,center=T)
    proj_gdp_pca=fgdp_train%*%pca$rotation[,1:t]
    mod_pca=lm(agdp_train~0+proj_gdp_pca)
    summary(mod_pca)
    
    
    #Performance on Test Set
    
    #SIR
    proj_gdp_test_sir=fgdp_test%*%dr.basis(sir)[,1:t]
    predict_sir=proj_gdp_test_sir%*%mod$coefficients
    forecast_residuals_sir=t(agdp_test)-t(as.matrix(predict_sir))
    
    
    #PCA
    proj_gdp_test_pca=fgdp_test%*%pca$rotation[,1:t]
    predict_pca=proj_gdp_test_pca%*%mod_pca$coefficients
    forecast_residuals_pca=t(agdp_test)-t(as.matrix(predict_pca))
   
    
  #Benchmark Model 1
  
  
  data1=read.csv("NGDP.csv")
  forecast1=t(as.numeric(data1$NGDP2))
  fgdp1=matrix(forecast1,nrow=40,ncol=10,byrow=T)
  agdp1=as.array(unique(data1$NGDP1))
  benchmark_forecast1=as.numeric()
  for(i in 1:nrow(fgdp1)){
    benchmark_forecast1[i]=mean(fgdp1[i,])
  }
  benchmark_forecast1=diff(benchmark_forecast1)
  benchmark_res1=diff(agdp1)-benchmark_forecast1
  forecast_residuals_benchmark1=benchmark_res1[21:39]
  
  
  
  #Benchmark Model 2:
  
  data1=read.csv("NGDP.csv")
  forecast1=t(as.numeric(data1$NGDP2))
  fgdp1=matrix(forecast1,nrow=40,ncol=10,byrow=T)
  agdp1=as.array(unique(data1$NGDP1))
  agdp_train1=agdp1[2:21]
  fgdp_train1=fgdp1[1:20,]
  fgdp_test1=fgdp1[21:40,]
  sMAPE_inter=matrix(rep(0,nrow(fgdp_train1)*ncol(fgdp_train1)),nrow=nrow(fgdp_train1),ncol=ncol(fgdp_train1))
  for(i in 1:nrow(fgdp_train)){
    for(j in 1:ncol(fgdp_train)){
      sMAPE_inter[i,j]=abs(agdp_train1[i]-fgdp_train1[i,j])/abs(agdp_train1[i]+fgdp_train1[i,j])
    }
  }
  sMAPE=as.numeric()
  wts=as.numeric()
    for(j in 1:ncol(fgdp_train1)){
      sMAPE[j]=mean(sMAPE_inter[,j])
    }
  wts=(1/sMAPE)/sum(1/sMAPE)
  benchmark_forecast2=as.numeric()
  for(i in 1:nrow(fgdp_test1)){
    benchmark_forecast2[i]=weighted.mean(fgdp_test1[i,],wts)
  }
  benchmark_forecast2=diff(benchmark_forecast2)
  forecast_residuals_benchmark2=diff(agdp1)[21:39]-benchmark_forecast2


  
  #Comparison:
  RMAE_sir=mean(abs(forecast_residuals_sir))
  RMAE_pca=mean(abs(forecast_residuals_pca))
  RMAE_benchmark1=mean(abs(forecast_residuals_benchmark1))
  RMAE_benchmark2=mean(abs(forecast_residuals_benchmark2))
  RMSE_sir=sqrt(mean(forecast_residuals_sir^2))
  RMSE_pca=sqrt(mean(forecast_residuals_pca^2))
  RMSE_benchmark1=sqrt(mean(forecast_residuals_benchmark1^2))
  RMSE_benchmark2=sqrt(mean(forecast_residuals_benchmark2^2))
  out=c(RMAE_sir,RMAE_pca,RMAE_benchmark1,RMAE_benchmark2,RMSE_sir,RMSE_pca,RMSE_benchmark1,RMSE_benchmark2)
  return(out)
  
}

RMAE_sir=as.numeric()
RMAE_pc=as.numeric()
RMAE_bench1=as.numeric()
RMAE_bench2=as.numeric()
RMSE_sir=as.numeric()
RMSE_pc=as.numeric()
RMSE_bench1=as.numeric()
RMSE_bench2=as.numeric()
for(i in 1:10){
  RMAE_sir[i]=proj(i)[1]
  RMAE_pc[i]=proj(i)[2]
  RMAE_bench1[i]=proj(i)[3]
  RMAE_bench2[i]=proj(i)[4]
  RMSE_sir[i]=proj(i)[5]
  RMSE_pc[i]=proj(i)[6]
  RMSE_bench1[i]=proj(i)[7]
  RMSE_bench2[i]=proj(i)[8]
}


plot(1:10,RMAE_sir,type="l",col="BLUE",ylim=c(0,200),main="Static Model",xlab="Number of Directions",ylab="MAE")
lines(RMAE_pc,col="RED",type="l")
lines(RMAE_bench1,col="YELLOW",type="l")
lines(RMAE_bench2,col="BLACK",type="l")
legend("topright", legend=c("SIR", "PCA","SA","WA"), fill = c("blue","red","yellow","black"))



plot(1:10,RMSE_sir,type="l",col="BLUE",ylim=c(0,300),main="Static Model",xlab="Number of Directions",ylab="RMSE")
lines(RMSE_pc,col="RED",type="l")
lines(RMSE_bench1,col="YELLOW",type="l")
lines(RMSE_bench2,col="BLACK",type="l")
legend("topright", legend=c("SIR", "PCA","SA","WA"), fill = c("blue","red","yellow","black"))

  
  
  
  








#SIR search for better model: Predicting 2011,2012


#Data Preparation:
plodata=read.csv("NGDP.csv")
forecast=t(as.numeric(data$NGDP2))
fgdp=matrix(forecast,nrow=40,ncol=10,byrow=T)
agdp=as.array(unique(data$NGDP1))
#Stationary Transformation
par(mfrow=c(5,2))
for(i in 1:ncol(fgdp)){
  # plot(fgdp[,i])
}
par(mfrow=c(1,1))

for(i in 1:nrow(fgdp)){
  fgdp[i,]=fgdp[i,]-agdp[i]
}
par(mfrow=c(5,2))
for(i in 1:ncol(fgdp)){
  # plot(fgdp[,i])
}
par(mfrow=c(1,1))

agdp=diff(agdp)
fgdp=fgdp[-40,]#2006 Q1-2015Q3

#Train-Test Split
fgdp_train=fgdp[1:20,]#Pseudo differenced forecasts from 10 individuals from Q1,2006-Q4,2010
agdp_train=agdp[1:20]#Actual differneced GDP values from Q1,2006-Q4,2010
fgdp_test=fgdp[21:28,]#Pseudo differenced forecasts from 10 individuals from Q1,2011-Q4,2012
agdp_test=agdp[21:28]#Pseudo differenced GDP values from Q1,2011-Q4,2012

forecast_residuals_sir=as.numeric()
forecast_residuals_pca=as.numeric()
forecast_residuals_benchmark=as.numeric()


#Baseline Model:

#SIR
sir=dr(formula =agdp_train~fgdp_train, slice.function = dr.slices.arc,
       nslices = 6,numdir=10, chi2approx = "wood", method = "sir")
proj_gdp=fgdp_train%*%dr.basis(sir)[,1:2]
mod=lm(agdp_train~0+proj_gdp)
summary(mod)

#PCA
pca=prcomp(fgdp_train,scale=T,center=T)
proj_gdp_pca=fgdp_train%*%pca$rotation[,1:2]
mod_pca=lm(agdp_train~0+proj_gdp_pca)
summary(mod_pca)


#Performance on Test Set

#SIR
proj_gdp_test_sir=fgdp_test%*%dr.basis(sir)[,1:2]
predict_sir=proj_gdp_test_sir%*%mod$coefficients
forecast_residuals_sir=t(agdp_test)-t(as.matrix(predict_sir))
RMAE_sir=mean(abs(forecast_residuals_sir))

#PCA
proj_gdp_test_pca=fgdp_test%*%pca$rotation[,1:2]
predict_pca=proj_gdp_test_pca%*%mod_pca$coefficients
forecast_residuals_pca=t(agdp_test)-t(as.matrix(predict_pca))
RMAE_pca=mean(abs(forecast_residuals_pca))

RMAE_sir
RMAE_pca



#Model Development:

#SIR
sir=dr(formula =agdp_train~fgdp_train, slice.function = dr.slices.arc,
       nslices = 4,numdir=10, chi2approx = "wood", method = "sir")
proj_gdp=fgdp_train%*%dr.basis(sir)[,1:4]
f1=proj_gdp[,1]
f2=proj_gdp[,2]
f3=proj_gdp[,3]
f4=proj_gdp[,4]
df=data.frame(agdp_train,f1,f2,f3,f4)
#Model 1:

mod=lm(agdp_train~.,data=df)
summary(mod)
ols_regress(agdp_train~.,data=df)
ols_step_forward_p(mod)

#Model 2:

mod1=lm(agdp_train~0+f1+f2,data=df)
summary(mod1)
ols_regress(mod1)

#Model 3(Boxcox) :
seq=seq(-5,5,by=0.01)
pr=as.numeric()
for(i in 1:length(seq)){
bc_agdp=bcnPower(agdp_train,lambda=seq[i],gamma=abs(min(agdp_train)))
mod2=lm(bc_agdp~f1)
sum=ols_regress(mod2)
pr[i]=sum$prsq
}
plot(seq,pr,type="l",xlab="Choices of Lambda",ylab="Predicted R_2")
l=seq[which.max(pr)]
bc_agdp=bcnPower(agdp_train,lambda=l,gamma=abs(min(agdp_train)))
mod2=lm(bc_agdp~f1)
ols_regress(mod2)   

#Performance on Test Set

#SIR
proj_gdp_test_sir=fgdp_test%*%dr.basis(sir)[,1]
predict_sir=cbind(rep(1,nrow(proj_gdp_test_sir)),proj_gdp_test_sir)%*%mod2$coefficients
#predict_sir=(proj_gdp_test_sir)%*%mod1$coefficients
forecast_residuals_sir=t(agdp_test)-t(as.matrix(predict_sir))
RMAE_sir=mean(abs(forecast_residuals_sir))

RMAE_sir

#Modification of Model 3(Boxcox):
ols_plot_resid_qq(mod2)
ols_plot_resid_fit(mod2)
ols_plot_dffits(mod2)
ols_plot_dfbetas(mod2)
bc_agdp[c(10,12,13)]
f1[c(10,12,13)]
x=c(10,12,13)
mod3=lm(bc_agdp[-x]~f1[-x])
ols_regress(mod3)


#Performance on Test Set(Modified Model after Residual Analysis)

#SIR
proj_gdp_test_sir=fgdp_test%*%dr.basis(sir)[,1]
predict_sir=cbind(rep(1,nrow(proj_gdp_test_sir)),proj_gdp_test_sir)%*%mod3$coefficients
#predict_sir=(proj_gdp_test_sir)%*%mod1$coefficients
forecast_residuals_sir=t(agdp_test)-t(as.matrix(predict_sir))
RMAE_sir=mean(abs(forecast_residuals_sir))

RMAE_sir

#PCA
pca=prcomp(fgdp_train,scale=T,center=T)
proj_gdp_pca=fgdp_train%*%pca$rotation
d1=proj_gdp[,1]
d2=proj_gdp[,2]
d3=proj_gdp[,3]
d4=proj_gdp[,4]
df=data.frame(agdp_train,d1,d2,d3,d4)

#Model 1:

mod=lm(agdp_train~.,data=df)
summary(mod)
ols_regress(agdp_train~.,data=df)
ols_step_forward_p(mod)

#Model 2:

mod1=lm(agdp_train~d1+d2,data=df)
summary(mod1)
ols_regress(mod1)

#Model 3(Boxcox):
seq=seq(-5,5,by=0.01)
pr=as.numeric()
for(i in 1:length(seq)){
  bc_agdp=bcnPower(agdp_train,lambda=seq[i],gamma=abs(min(agdp_train)))
  mod2=lm(bc_agdp~d1)
  sum=ols_regress(mod2)
  pr[i]=sum$mae
}
plot(seq,pr,type="l")
l=seq[which.max(pr)]
bc_agdp=bcnPower(agdp_train,lambda=l,gamma=abs(min(agdp_train)))
mod2=lm(bc_agdp~d1)
ols_regress(mod2)   


#Performance on Test Set

#PCA
proj_gdp_test_pca=fgdp_test%*%pca$rotation[,1]
predict_pca=cbind(rep(1,nrow(proj_gdp_test_pca)),proj_gdp_test_pca)%*%mod2$coefficients
forecast_residuals_pca=t(agdp_test)-t(as.matrix(predict_pca))
RMAE_pca=mean(abs(forecast_residuals_pca))

RMAE_pca










