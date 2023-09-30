#=====================================
#R code for chapter 6
#=====================================

#====================================================
#Example: Distribution of sample variance
x=c(55,69,46,51,68,57,48,56,55,52,44,73,55,46,47) #Original observations
varx=var(x) #Sample variance of x
print(varx)
B=1000 #Number of bootstrap samples
bootsam=matrix(NA,nrow=B,ncol=length(x)) #Matrix storing the bootstrap samples
#Note: it is typically not necessary to store the bootstrap samples - can demand a lot of
#memory if sample size is large

set.seed(1000) #Make sure the results are reproducible by setting a seed for the random number generator
#Note: If you are using R with version prior to 3.6.0, you will get a different random sequence when using
#the function "sample" (even when the same seed is used). This is because a patch has been made to the
#random number generator in version 3.6.0.

resvec=rep(NA,B) #Vector storing the bootstrap estimates of the sample variance
for (i in 1:B){
	thissam=sample(x=x,size=length(x),replace=TRUE) #Sample from the x's
	bootsam[i,]=thissam
	resvec[i]=var(thissam)
	if(i%%100==0) print(i) #Print every 100th run
}

hist(resvec,main="Sampling dist. of S-squared based on 1000 bootstrap samples",xlab="S-squared")

mean(resvec) #mean value of the bootstrap estimates
mean(resvec-varx) #bias
mean((resvec-mean(resvec))^2) #variance
mean((resvec-varx)^2) #mse
sort(resvec)[B*0.025] #lower limit of 95% bootstrap CI (percentile)
sort(resvec)[B*0.975] #upper limit of 95% bootstrap CI (percentile)
2*varx-sort(resvec)[B*0.975] #lower limit of 95% bootstrap CI (pivot)
2*varx-sort(resvec)[B*0.025] #upper limit of 95% bootstrap CI (pivot)

#====================================================
#Example: P(X>Y)
mu1=30; sig1=10
mu2=35; sig2=10
trueprob=1-pnorm(0,mu1-mu2,sqrt(sig1^2+sig2^2))
trueprob #True prob that X>Y

estimator=function(x,y){#Define a function that estimates P(X>Y) given data
	count=0
	for(i in 1:length(y)){
		count=count+sum(x>y[i])
	}
	return(count/length(x)/length(y))
}

set.seed(2000) #Make sure the results are reproducible by setting a seed for the random number generator

x=rnorm(50,mu1,sig1); y=rnorm(50,mu2,sig2) #Our original samples
B=1000
resvec=rep(NA,B)
bootsamx=matrix(NA,nrow=B,ncol=length(x))
bootsamy=matrix(NA,nrow=B,ncol=length(y))
theta=estimator(x,y)

for(i in 1:B){
	thissamx=sample(x=x,size=length(x),replace=TRUE) #Sample from the x's
	thissamy=sample(x=y,size=length(y),replace=TRUE) #Sample from the y's
	bootsamx[i,]=thissamx
	bootsamy[i,]=thissamy
	resvec[i]=estimator(thissamx,thissamy)
	if(i%%100==0) print(i) #Print every 100th run
}

hist(resvec,main="Sampling dist. of P(X>Y) based on 1000 bootstrap samples",xlab="Probability",breaks=20)

mean(resvec) #mean value of the bootstrap estimates
mean(resvec-theta) #bias
mean((resvec-mean(resvec))^2) #variance
mean((resvec-theta)^2) #mse
sort(resvec)[B*0.025] #lower limit of 95% bootstrap CI (percentile)
sort(resvec)[B*0.975] #upper limit of 95% bootstrap CI (percentile)
2*theta-sort(resvec)[B*0.975] #lower limit of 95% bootstrap CI (pivot)
2*theta-sort(resvec)[B*0.025] #upper limit of 95% bootstrap CI (pivot)
