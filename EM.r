#=====================================
#R code for chapter 7
#=====================================

#====================================================
#Example: Multinomial observations
x=c(125,18,20,34) #Original observations for the four cells
M=10 #Number of EM iterations
resvec=loglik=rep(NA,M+1)
resvec[1]=0.5 #Initial value

for(i in 2:(M+1)){
	z12=x[1]*resvec[i-1]/(2+resvec[i-1])
	nexttheta=(z12+x[4])/(z12+sum(x[2:4]))
	resvec[i]=nexttheta
	loglik[i]=dmultinom(x,size=sum(x),prob=c(0.5+nexttheta/4,
		(1-nexttheta)/4,(1-nexttheta)/4,nexttheta/4),log=T)
}

par(mfrow=c(1,2))
plot(0:M,resvec,type="o",main="Convergence of estimated theta",xlab="Iter",ylab="Value")
plot(0:M,loglik,type="o",main="Log-likelihood of observed data",xlab="Iter",ylab="Value")

#====================================================
#Example: Censored observations
x=c(0.920,1.382,2.552,0.920,1.333,2.497,0.957,0.509,0.578,0.939,3.519,0.299,0.607,1.679,3.491,3.565)
#Assume also that there are 4 censored observations; we only know that they are >4.
m=length(x); n=4; cc=4 #c is an internal command so we use cc instead
M=20 #Number of EM iterations
resvec=loglik=rep(NA,M+1)
resvec[1]=1/mean(x) #Initial value

for(i in 2:(M+1)){
	nexttheta=(m+n)/(sum(x)+(n*cc*resvec[i-1]+n)/resvec[i-1])
	resvec[i]=nexttheta
	loglik[i]=n*log(1-pexp(cc,rate=nexttheta))+sum(dexp(x,rate=nexttheta,log=T))
}

plot(0:M,resvec,type="o",main="Convergence of estimated theta",xlab="Iter",ylab="Value")
plot(0:M,loglik,type="o",main="Log-likelihood of observed data",xlab="Iter",ylab="Value")
#Data was actually generated based on a true theta of 0.5

#====================================================
#Example: Mixture distributions

#Generate data
set.seed(160)
y1=rnorm(50,mean=4,sd=1); y2=rnorm(50,mean=10,sd=2);
component=sample(1:2,size=50,replace=T,prob=c(0.7,0.3))
x=c(y1[which(component==1)],y2[which(component==2)]) #Note how we pick the correct element from y1/y2
x=sample(x) #Shuffle them
hist(x,breaks=10)

M=30 #Number of EM iterations
resmat=matrix(NA,nrow=M+1,ncol=5) #Columns are mu1,sig1,mu2,sig2,p
resmat[1,]=c(5,2,11,2,0.5) #Initial value
loglik=rep(NA,M+1)

for(i in 2:(M+1)){
	mu10=resmat[i-1,1]; sig10=resmat[i-1,2]
	mu20=resmat[i-1,3]; sig20=resmat[i-1,4]; p0=resmat[i-1,5]
	q0vec=p0*dnorm(x,mu10,sig10)/(p0*dnorm(x,mu10,sig10)+(1-p0)*dnorm(x,mu20,sig20))
	
	mu11=sum(q0vec*x)/sum(q0vec); sig11=sqrt(sum(q0vec*(x-mu11)^2)/sum(q0vec))
	mu21=sum((1-q0vec)*x)/sum((1-q0vec)); sig21=sqrt(sum((1-q0vec)*(x-mu21)^2)/sum((1-q0vec)))
	p1=mean(q0vec)
	resmat[i,]=c(mu11,sig11,mu21,sig21,p1)
	loglik[i]=sum(log(p1*dnorm(x,mu11,sig11)+(1-p1)*dnorm(x,mu21,sig21)))
}

par(mfrow=c(2,3))
plot(0:M,resmat[,1],type="o",main="mu1",ylim=c(3,6),xlab="Iter",ylab="Value"); abline(h=4,lty=2)
plot(0:M,resmat[,2],type="o",main="sig1",ylim=c(0,3),xlab="Iter",ylab="Value"); abline(h=1,lty=2)
plot(0:M,resmat[,3],type="o",main="mu2",ylim=c(9,12),xlab="Iter",ylab="Value"); abline(h=10,lty=2)
plot(0:M,resmat[,4],type="o",main="sig2",ylim=c(0,3),xlab="Iter",ylab="Value"); abline(h=2,lty=2)
plot(0:M,resmat[,5],type="o",main="p",ylim=c(0,1),xlab="Iter",ylab="Value"); abline(h=0.7,lty=2)
plot(0:M,loglik,type="o",main="Log-likelihood of observed data",xlab="Iter",ylab="Value")

#====================================================
#Standard error: Multinomial observations
#Method 1
library(numDeriv)
loglikfun=function(theta,dat){
	return(dmultinom(dat,size=sum(dat),prob=c(0.5+theta/4,(1-theta)/4,(1-theta)/4,theta/4),log=T))
}
dat=c(125,18,20,34)
theta=0.62682149762889672
GG=grad(loglikfun,theta,dat=dat) #Gradient, should be 0 at mle
HH=hessian(loglikfun,theta,dat=dat)
GG; sqrt(1/-HH) #SE

#Method 2
xx=c(rep(1,125),rep(2,18),rep(3,20),rep(4,34))
B=1000 #Number of bootstrap samples
M=20 #Number of EM iterations
resvec=rep(NA,B)
set.seed(123)
for(i in 1:B){
	theta=0.5 #Initial value
	thisx=sample(xx,replace=T) #Sample from original data
	y=table(thisx) #Cell counts for this bootstrap sample
	for(j in 1:M){
		z12=y[1]*theta/(2+theta)
		theta=(z12+y[4])/(z12+sum(y[2:4]))
	}
	resvec[i]=theta
	if(i%%100==0) print(i)
}
hist(resvec,main="Sampling dist. of theta based on 1000 bootstrap samples",xlab="theta",breaks=20)
sd(resvec)

#====================================================
#Standard error: Censored observations
#Method 1
library(numDeriv)
loglikfun=function(theta,dat,cc,n){
	return(n*log(1-pexp(cc,rate=theta))+sum(dexp(dat,rate=theta,log=T)))
}
dat=c(0.920,1.382,2.552,0.920,1.333,2.497,0.957,0.509,0.578,0.939,3.519,0.299,0.607,1.679,3.491,3.565)
n=4; cc=4; theta=0.38326107265
GG=grad(loglikfun,theta,dat=dat,cc=cc,n=n) #Gradient, should be 0 at mle
HH=hessian(loglikfun,theta,dat=dat,cc=cc,n=n)
GG; sqrt(1/-HH) #SE

#Method 2
B=1000 #Number of bootstrap samples
M=30 #Number of EM iterations
resvec=rep(NA,B)
set.seed(123)
for(i in 1:B){
	theta=0.5 #Initial value
	thisind=sample(1:20,replace=T) #Sample the indices
	cc=4; n=sum(thisind>=17); m=20-n #n is number of censored obs
	thisx=dat[thisind[thisind<=16]] #Observed data in bootstrap sample
	for(j in 1:M){
		nexttheta=(m+n)/(sum(thisx)+(n*cc*theta+n)/theta)
		theta=nexttheta
	}
	resvec[i]=theta
	if(i%%100==0) print(i)
}
hist(resvec,main="Sampling dist. of theta based on 1000 bootstrap samples",xlab="theta",breaks=20)
sd(resvec)

#====================================================
#Standard error: Mixture distributions
#Method 1
library(numDeriv)
loglikfun=function(theta,dat){ #Note that theta is a vector parameter (cannot use multiple arguments)
	mu1=theta[1]; sig1=theta[2]; mu2=theta[3]; sig2=theta[4]; p=theta[5]
	return(sum(log(p*dnorm(dat,mu1,sig1)+(1-p)*dnorm(dat,mu2,sig2))))
}
#Generate data set
set.seed(160)
y1=rnorm(50,mean=4,sd=1); y2=rnorm(50,mean=10,sd=2);
component=sample(1:2,size=50,replace=T,prob=c(0.7,0.3))
x=c(y1[which(component==1)],y2[which(component==2)]) #Note how we pick the correct element from y1/y2
dat=sample(x)

thetamle=c(4.09986063703316006,1.06565631064516997,10.24527945167525189,
		1.64816000488422931,0.63159485412309035)
GG=grad(loglikfun,thetamle,dat=dat) #Gradient, should be 0 at mle
HH=hessian(loglikfun,thetamle,dat=dat)
GG; solve(-HH) #Var-cov matrix
sqrt(diag(solve(-HH))) #SE of parameters

#Method 2
B=1000 #Number of bootstrap samples
M=50 #Number of EM iterations
resmat=matrix(NA,nrow=B,ncol=5)
set.seed(123)
for(i in 1:B){
	theta=thetamle #Initial value
	thisx=sample(dat,replace=T)
	for(j in 1:M){
		mu10=theta[1]; sig10=theta[2]
		mu20=theta[3]; sig20=theta[4]; p0=theta[5]
		q0vec=p0*dnorm(thisx,mu10,sig10)/(p0*dnorm(thisx,mu10,sig10)+(1-p0)*dnorm(thisx,mu20,sig20))
		
		mu11=sum(q0vec*thisx)/sum(q0vec); sig11=sqrt(sum(q0vec*(thisx-mu11)^2)/sum(q0vec))
		mu21=sum((1-q0vec)*thisx)/sum((1-q0vec)); sig21=sqrt(sum((1-q0vec)*(thisx-mu21)^2)/sum((1-q0vec)))
		p1=mean(q0vec)
		theta=c(mu11,sig11,mu21,sig21,p1)
	}
	resmat[i,]=theta
	if(i%%100==0) print(i)
}

par(mfrow=c(2,3))
hist(resmat[,1],main="Sampling dist.: mu1",xlab="mu1",breaks=20)
hist(resmat[,2],main="Sampling dist.: sig1",xlab="sig1",breaks=20)
hist(resmat[,3],main="Sampling dist.: mu2",xlab="mu2",breaks=20)
hist(resmat[,4],main="Sampling dist.: sig2",xlab="sig2",breaks=20)
hist(resmat[,5],main="Sampling dist.: p",xlab="p",breaks=20)
apply(resmat,2,sd) #Calculate sd columnwise
