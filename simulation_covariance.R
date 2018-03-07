
require("rjags")
require("mvtnorm")

# need both logit and logit^{-1} functions
logit=function(x)log(x/(1-x))
invlogit=function(x){1/(1+exp(-x))}

# set up model specification
n <- 500	    # total number of trials 

# establish the hyperparameters
sig1 <- .5		# std. dev. of single-trial BOLD responses
sig2 <- 1			# std. dev. of item memory strength (logit scale)
rho <- .6			# cor b/n brain activation and memory strength

# set up hyper variance-covariance matrix Sigma
sigma <- matrix(c(sig1^2,         # element [1,1]
                  sig1*sig2*rho,  # element [1,2]
                  sig1*sig2*rho,  # element [2,1]
                  sig2^2          # element [2,2]
),2,2,byrow=TRUE)

# set up hyper mean vector phi
phi <- c(2,0)

# simulate single-trial delta and theta matrix
DeltaTheta <- rmvnorm(n,phi,sigma)

# generate observed variable nodes
ts <- seq(0,4,1)   # scan times
sig <- .5			     # the std. dev. of BOLD responses

# declare some storage objects
N=matrix(NA,n,length(ts))
B=numeric(n)

# loop over trials
for(i in 1:n){
# N is a normal deviate with mean controlled by delta
N[i,]=rnorm(length(ts),DeltaTheta[i,1]*ts,sig) 
# B is a Bernoulli deviate with probability controlled by theta
B[i]=rbinom(1,1,invlogit(DeltaTheta[i,2]))
}

# combine the generated data into a list to pass to JAGS
dat=list('n'=n, 
         'B'=B, 
         'N'=N, 
         'ts'=ts, 
         'Nt'=length(ts), 
         'sig'=sig, 
         'I0'=diag(2), 
         'n0'=2,
         'phi0'=rep(0,2), 
         's0'=diag(2))

############################################################# Sampling from the model

# locate the JAGS code, pass variables, setup sampler
jags <- jags.model('model_hierarchical.txt',
                data = dat,
                n.chains = 4,
                n.adapt = 1000)

# continue adapting the sampler to optimize efficiency
adapt(jags, 1000, end.adaptation=TRUE);

# continue sampling to ensure convergence
update(jags, 1000)

# draw final samples, and monitor important variables
out=jags.samples(jags,
             c('phi', 'Sigma', 'DeltaTheta'),
             1000)

################################################### parameter recovery

# calculate the mean of the posteriors
pms=apply(out$DeltaTheta,c(1,2),mean)
# delta is the first column, theta is the second column
delta=pms[,1]
theta=pms[,2]

lwd=3
cex.big=2.5
cex.tiny=1.5
col=rgb(0,0,1,.25)

mar1=c(.1,.1,.1,.1)

pdf("st_pars.pdf",9,4)
layout(matrix(c(1,2,3,4,5,5),2,3,byrow=TRUE),widths=c(1,7,7),heights=c(7,1))
par(cex=1)

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(.5,.5,"Posterior Mean",srt=90,cex=cex.big)

cor(DeltaTheta[,1],delta)
par(mar=c(3,3,1,1))
lim=range(c(DeltaTheta[,1],delta))
plot(DeltaTheta[,1],delta,xlab="",ylab="",xlim=lim,ylim=lim,pch=16,col=col)
abline(0,1)
text(.7,3.2,expression(delta),cex=cex.big)
text(3.1,.6,"R = 0.99",cex=cex.tiny) #insert correlation coefficient from line 102

cor(DeltaTheta[,2],theta)
par(mar=c(3,3,1,1))
lim=range(c(DeltaTheta[,2],theta))
plot(DeltaTheta[,2],theta,xlab="",ylab="",xlim=lim,ylim=lim,pch=16,col=col)
abline(0,1)
text(-2.6,2.6,expression(theta),cex=cex.big)
text(2.2,-2.6,"R = 0.71",cex=cex.tiny) # insert correlation coefficient from line 110

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(.5,.5,"True Model Parameter",cex=cex.big)

dev.off()


################################################### hyperparameter parameter recovery

require(vioplot)

lwd=3
cex.big=2.5
cex.tiny=1.5
col=rgb(0,0,1,.25)
cex.pts=3

mar1=c(.1,.1,.1,.1)

pdf("hyper_pars.pdf",9,4)
layout(matrix(c(1,2,3,4,5,5),2,3,byrow=TRUE),widths=c(1,7,7),heights=c(7,1))
par(cex=1)

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(.5,.5,"Estimate",srt=90,cex=cex.big)

par(mar=c(3,3,1,1))
xlim=c(.5,2.5)
ylim=c(-.5,2.5)
plot(NA,xlab="",ylab="",xlim=xlim,ylim=ylim,xaxt="n")
vioplot(as.numeric(out$phi[1,,]),at=1,col=col,add=TRUE)
vioplot(as.numeric(out$phi[2,,]),at=2,col=col,add=TRUE)
points(1,phi[1],pch=4,cex=cex.pts,lwd=4)
points(2,phi[2],pch=4,cex=cex.pts,lwd=4)
axis(1,at=c(1,2),labels=c(expression(phi[1]),expression(phi[2])))

par(mar=c(3,3,1,1))
xlim=c(.5,3.5)
ylim=c(.4,2.2)
plot(NA,xlab="",ylab="",xlim=xlim,ylim=ylim,xaxt="n")
sig11=as.numeric(out$Sigma[1,1,,])
sig12=as.numeric(out$Sigma[1,2,,])
sig22=as.numeric(out$Sigma[2,2,,])
esig1=sqrt(sig11)
esig2=sqrt(sig22)
erho=sig12/(sqrt(sig11)*sqrt(sig22))
vioplot(esig1,at=1,col=col,add=TRUE)
vioplot(esig2,at=2,col=col,add=TRUE)
vioplot(erho,at=3,col=col,add=TRUE)
points(1,sig1,pch=4,cex=cex.pts,lwd=4)
points(2,sig2,pch=4,cex=cex.pts,lwd=4)
points(3,rho,pch=4,cex=cex.pts,lwd=4)
axis(1,at=c(1,2,3),labels=c(expression(sigma[1]),expression(sigma[2]),expression(rho)),cex=3)

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")
text(.5,.5,"True Model Parameter",cex=cex.big)

dev.off()


