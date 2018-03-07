
require("rjags")
require("mvtnorm")

# need both logit and logit^{-1} functions
logit=function(x)log(x/(1-x))
invlogit=function(x){1/(1+exp(-x))}

# set up model specification
n <- 500	    # total number of trials 

# establish the hyperparameters for delta
sig1 <- .5		# std. dev. of single-trial BOLD responses, ROI 1
sig2 <- .6		# std. dev. of single-trial BOLD responses, ROI 2
rho <- .4 		# cor b/n brain activations

# set up hyper variance-covariance matrix Sigma
sigma <- matrix(c(sig1^2,         # element [1,1]
                  sig1*sig2*rho,  # element [1,2]
                  sig1*sig2*rho,  # element [2,1]
                  sig2^2          # element [2,2]
),2,2,byrow=TRUE)

# set up hyper mean vector phi
phi <- c(1.5,2)
  
# simulate single-trial delta matrix
Delta <- rmvnorm(n,phi,sigma)

# generate observed variable nodes
ts <- seq(0,4,1)   # scan times
sig <- .5			     # the std. dev. of BOLD responses

Nroi <- 2     # total number of ROIs

# declare some storage objects
N=array(NA,c(n,length(ts),Nroi))
B=numeric(n)
theta=numeric(n)

# set up regression parameters
beta <- c(.5,.3) # two ROIs

# loop over trials
for(i in 1:n){
  for(k in 1:Nroi){
  # N is a normal deviate with mean controlled by delta
  N[i,,k]=rnorm(length(ts),Delta[i,k]*ts,sig) 
  }
  # theta[i] is the single-trial behavioral parameter
  theta[i] <- Delta[i,]%*%beta
  # B is a Bernoulli deviate with probability controlled by theta
  B[i]=rbinom(1,1,invlogit(theta[i]))
}

dat = list('n'=n,
           'Nroi'=Nroi,
           'B'=B, 
           'N'=N, 
           'ts'=ts, 
           'Nt'=length(ts), 
           'sig'=sig, 
           'I0'=diag(2), 
           'n0'=2,
           'phi0'=rep(0,2), 
           's0'=diag(2))

#############################################################

# locate the JAGS code, pass variables, setup sampler
jags <- jags.model('model_directional.txt',
                data = dat,
                n.chains = 4,
                n.adapt = 1000)

# continue adapting the sampler to optimize efficiency
adapt(jags, 1000, end.adaptation=TRUE);

# continue sampling to ensure convergence
update(jags, 1000)

# draw final samples, and monitor important variables
out=jags.samples(jags,
             c('phi', 'Sigma', 'beta'),
             1000)

############################################################# parameter recovery & potting

cex.text=2
mar1=c(.1,.1,.1,.1)
mar2=c(3,3,1,1)
breaks=40

pdf("beta.pdf",10,5)
layout(matrix(c(1,2,3,4,5,6),2,3,byrow=TRUE),widths=c(1,6,6),heights=c(6,1))
par(cex=1.3)

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
text(.5,.5,"Density",srt=90,cex=cex.text)

xs=seq(-5,5,.01)
par(mar=mar2)
k=1
hist(out$beta[k,,],prob=TRUE,main="",xlab="",ylab="",breaks=breaks)
abline(v=beta[k],col="red",lwd=4)
lines(xs,dnorm(xs,0,1/.001^2),lty=3,lwd=3)

par(mar=mar2)
k=2
hist(out$beta[k,,],prob=TRUE,main="",xlab="",ylab="",breaks=breaks)
abline(v=beta[k],col="red",lwd=4)
lines(xs,dnorm(xs,0,1/.001^2),lty=3,lwd=3)

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n",xaxt="n",yaxt="n")

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
text(.5,.5,expression(beta[1]),cex=cex.text)

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
text(.5,.5,expression(beta[2]),cex=cex.text)

dev.off()
