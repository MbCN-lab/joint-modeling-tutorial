# Clear your workspace & set your working directory
rm(list=ls())

# Load required packages and modules
require("rjags")
load.module("wiener") 

# Behavioral and neural data files should be in the same folder.
load("standardizedBOLD.Rdata")
load("block1.Rdata")

# Recode data
temp.rt=block1.data$rt
temp.resp=block1.data$resp
temp.rt[temp.resp==2]<-temp.rt[temp.resp==2]*-1
rt=temp.rt

# Data
TR = 2
lenS = nrow(block1.data$design.matrix)
onset = block1.data$design.matrix[,"onset"]
N = standardized.mean.data[,"block1.v1"] # Block 1, V1

dat = list(N = N, lenN = length(N), TR = TR, t=rt, n.trials=length(rt),
           onset = onset, lenS = length(onset),
           a1 = 6, a2 = 16, b1 = 1, b2 = 1, c = 1/6)

# HRF model 1: Standard double-gamma
model.double.gamma.wiener = "
model{
  # Likelihood
  for (i in 1:lenN) {
    N[i] ~ dnorm(muN[i], sigma)
    Npred[i] ~ dnorm(muN[i], sigma)
    muN[i] = beta0 + inprod(delta[], X[i, ])
  }

  # Define a design matrix
  # HRF model: standard double gamma
  # The model statement is much simpler because the standard double gamma does not require normalization of the curve.
  # h(t) = beta * ( t^(a1-1) * b1^(a1) * exp(-b1*t) / exp(loggam(a1)) - c * t^(a2-1) * b2^(a2) * exp(-b2*t) / exp(loggam(a2)))
  for (i in 1:lenS){
    for (j in 1:lenN){
      temp[j,i] = (j-1) * TR - onset[i]
      Xt[j,i] = ifelse(temp[j,i] >= 0, temp[j,i], 0)
      X[j,i] = (Xt[j,i]^(a1-1) * (b1)^(a1) * exp(-b1*Xt[j,i]) / exp(loggam(a1))) - c * (Xt[j,i]^(a2-1) * (b2)^(a2) * exp(-b2*Xt[j,i]) / exp(loggam(a2)))
    }
  }
  # About JAGS scalar functions:
  # 1. Note that JAGS does not provide a gamma function defined in a linear scale, but only in a log scale.
  # 2. We can also use pow(x, z) to express x^z.

# Prior
sigma ~ dgamma(.001, .001)
beta0 ~ dnorm(0, 0.001)
for (j in 1:lenS){
delta[j] ~ dnorm(0, 0.001)
}
alpha ~ dunif(0.0001,10) 
tau ~ dunif(0,.04)
beta <- 0.5

for (i in 1:n.trials) {
t[i] ~ dwiener(alpha,tau,beta,(delta[2*i]-delta[2*i-1]))
}
}
"

dat.dg = list(N = N, lenN = length(N), TR = TR,
           onset = onset, lenS = length(onset),
           a1 = 6, a2 = 16, b1 = 1, b2 = 1, c = 1/6)

# HRF model 1: Standard double-gamma
model.double.gamma = "
model{
# Likelihood
for (i in 1:lenN) {
N[i] ~ dnorm(muN[i], sigma)
Npred[i] ~ dnorm(muN[i], sigma)
muN[i] = beta0 + inprod(delta[], X[i, ])
}

# Define a design matrix
# HRF model: standard double gamma
# The model statement is much simpler because the standard double gamma does not require normalization of the curve.
# h(t) = beta * ( t^(a1-1) * b1^(a1) * exp(-b1*t) / exp(loggam(a1)) - c * t^(a2-1) * b2^(a2) * exp(-b2*t) / exp(loggam(a2)))
for (i in 1:lenS){
for (j in 1:lenN){
temp[j,i] = (j-1) * TR - onset[i]
Xt[j,i] = ifelse(temp[j,i] >= 0, temp[j,i], 0)
X[j,i] = (Xt[j,i]^(a1-1) * (b1)^(a1) * exp(-b1*Xt[j,i]) / exp(loggam(a1))) - c * (Xt[j,i]^(a2-1) * (b2)^(a2) * exp(-b2*Xt[j,i]) / exp(loggam(a2)))
}
}
# About JAGS scalar functions:
# 1. Note that JAGS does not provide a gamma function defined in a linear scale, but only in a log scale.
# 2. We can also use pow(x, z) to express x^z.

# Prior
sigma ~ dgamma(.001, .001)
beta0 ~ dnorm(0, 0.001)
for (j in 1:lenS){
delta[j] ~ dnorm(0, 0.001)
}
}
"

dat.wiener = list(t=rt, n.trials=length(rt))

# HRF model 1: Standard double-gamma

model.wiener = "
model{
alpha ~ dunif(0.0001,10) 
tau ~ dunif(0,.04)
beta <- 0.5
delta ~ dnorm(0, 0.001)

for (i in 1:n.trials) {
t[i] ~ dwiener(alpha,tau,beta,delta)
}
}
"

### Directed Joint Model 

# Initialization
model.dgw = jags.model(textConnection(model.double.gamma.wiener), data = dat,
                   n.chains = 3, n.adapt = 2000)

# Burn-in
update(model.dgw, n.iter = 4000, progress.bar = "text") 

# Posterior sampling
double.gamma.wiener.out = coda.samples(model = model.dgw, 
                   variable.names = c("beta0", "delta", "sigma", "alpha", "tau", "Npred"), 
                   n.iter = 6000) 

### Double Gamma Only

# Initialization
model.dg = jags.model(textConnection(model.double.gamma), data = dat.dg,
                        n.chains = 3, n.adapt = 2000)

# Burn-in
update(model.dg, n.iter = 4000, progress.bar = "text") 

# Posterior sampling
double.gamma.out = coda.samples(model = model.dg,
                   variable.names = c("beta0", "delta", "sigma", "Npred"), 
                   n.iter = 6000) 

### Wiener Only

# Initialization
model.w = jags.model(textConnection(model.wiener), data = dat.wiener,
                        n.chains = 3, n.adapt = 2000)

# Burn-in
update(model.w, n.iter = 4000, progress.bar = "text") 

# Posterior sampling
wiener.out = coda.samples(model = model.w, 
                   variable.names = c("delta", "alpha", "tau"), 
                   n.iter = 6000) 

############################################################# BOLD Recovery (Figure 12)

# Define HRF function and extract posterior samples 

hrf = function(t, beta){
  beta *( (t^(6-1) * 1^(6) * exp(-1*t)) / gamma(6) - (1/6) * (t^(16-1) * 1^(16) * exp(-1*t)) / gamma(16) )
}

post.samples = function(out){
  n.chains = length(out)
  nrow = nrow(out[[1]])
  ncol = ncol(out[[1]])
  post = array(NA, c(nrow, ncol, n.chains))
  for (i in 1:n.chains){
    post[,,i] = out[[i]]
  }
  colnames(post) = colnames(out[[1]])
  return(post)
}

post = post.samples(double.gamma.wiener.out)

ts = seq(0, 242, 0.01)
onset = block1.data$design.matrix[,1]

var.names = colnames(double.gamma.wiener.out[[1]])
# Posterior mean
beta0 = mean(post[,"beta0",])
betas.mean = apply(post[,124:163,], 2, function(x) mean(x))
betas.CI95 = apply(post[,124:163,], 2, function(x) quantile(x, probs = c(.025, .975)) )

contrast = as.vector(t(block1.data$stim.matrix))

# Design matrix
temp = matrix(NA, nrow = length(ts), ncol = 40)
X = matrix(NA, nrow = length(ts), ncol = 40)
for (i in 1:40){
  temp[,i] = seq(0 - onset[i], 242 - onset[i], 0.01)
  temp[which(temp[,i] < 0),i] = 0
  X[,i] = hrf(temp[,i], betas.mean[i])
}

# Convolved HRF
Xbeta = rowSums(X) + beta0
X.CI95 = apply(post[,1:121,], 2, function(x) quantile(x, probs = c(.025, .975)) )

# Plot BOLD data & model prediction
par(mfrow = c(1,1), mar = c(4,5,1,2))
ts.discrete = seq(0, 240, by = 2)
plot(ts.discrete, standardized.mean.data[,1], ylim = c(-5, 5), pch = 16,
     xlab = "Time (seconds)", ylab = "BOLD response",
     cex.lab = 2, cex.axis = 1.5)
lines(ts, Xbeta, col = "red", lwd = 2)
matlines(ts.discrete, t(X.CI95), col = "red", lty = 3)
#legend(inset = 0.01, "bottomleft", bty = "n",cex = 0.75,
#       lty = c(1,3), lwd = c(2,1), col = "red",
#       c("Model prediction (based on posterior mean of betas)",
#         "95% credible interval of posterior predictive distribution"))

############################################################# Neural Activation vs Behavioral Variables (Figure 13)

# Posterior estimates for delta -- Directed joint model -- replace with your estimates if run independently 
delta2<-c(4.149,9.582,11.286,11.605,-7.324,11.944,8.203,7.965,8.350,10.826,
          6.581,0.597,8.663,1.919,1.440,5.330,6.098,5.277,5.098,5.017) 
delta1<-c(5.720,9.002,12.353,10.594,4.964,-11.670,5.748,4.048,3.862,6.041,
          9.059,2.062,7.062,4.532,2.881,9.885,7.652,10.170,9.192,3.588)

# Remove outlier and nonresponse (for illustrative purposes -- may not be necessary)

comp.rt<-abs(rt) ; comp.rt<- comp.rt[-c(5,6)] 
comp.resp<-temp.resp ; comp.resp<- comp.resp[-c(5,6)] 
dd<-delta2-delta1 ; dd<-dd[-c(5,6)] # calculate difference 

# Plot 
lwd=3
cex.big=2.5
cex.tiny=1.5
col=rgb(0,0,1,.25)

par(mar=c(5,6,2,2))
plot(comp.rt,dd,xlab="Response Time", pch=c(16,3)[as.numeric(comp.resp)], 
     ylab= expression(delta[i]), cex.lab=cex.big)
abline(h=0,lty=2)
legend(2.30,4.9, pch=c(16,3), c("Respond 1st Stimulus", "Respond 2nd Stimulus"), bty='n', 
       cex=1.2)


############################################################# PSD Comparison (Figure 14)

temp1 = post.samples(double.gamma.wiener.out)
temp2 = post.samples(wiener.out)
temp3 = post.samples(double.gamma.out)

sum.dgw = summary(double.gamma.wiener.out)
sum.dg = summary(double.gamma.out)
sum.w = summary(wiener.out)

post.sd.dgw = sum.dgw$statistics[,2]
post.sd.dg = sum.dg$statistics[,2]
post.sd.w = sum.w$statistics[,2]

# Figure 14
layout(matrix(c(1,2,
                1,3), nrow = 2, byrow = T))

par(mar = c(4,5,2.5,2))
plot(post.sd.dgw[which(names(post.sd.dgw) == "delta[1]"):which(names(post.sd.dgw) == "delta[40]")],
     post.sd.dg[which(names(post.sd.dg) == "delta[1]"):which(names(post.sd.dg) == "delta[40]")],
     xlab = "Joint model", ylab = "Neural only",
     cex.lab = 2, cex.axis = 1.5, cex.main = 2.3,
     xlim = c(2, 5), ylim = c(2, 5), pch = 16, col = rgb(0,0,0,.3), cex = 1.2,
     main = expression("Standard deviation: " * delta))
curve(1 * x, lty = 3, col = "gray", add = T)

par(mar = c(4, 5, 1, 2))
hist(temp1[,which(colnames(temp1) == "sigma"),], breaks = seq(0, 4, 0.1),
     col = rgb(1,0,0,0.3), cex.lab = 2, cex.axis = 1.5, cex.main = 2.3,
     xlab = expression(sigma), main = "", #main = expression(sigma),
     freq = F)
hist(temp3[,which(colnames(temp3) == "sigma"),], breaks = seq(0, 4, 0.1),
     col = rgb(0,0,1,0.3),
     freq = F, add = T)
legend(inset = 0.01, "topright", cex = 1.3, pt.cex = 1.5, bty = "n",
       col = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)), pch = 15,
       c("Joint model", "Neural only"))

hist(temp1[,which(colnames(temp1) == "beta0"),], breaks = seq(-2, 0.5, 0.1),
     col = rgb(1,0,0,0.3), cex.lab = 2, cex.axis = 1.5, cex.main = 2.3,
     xlab = expression(beta[0]), main = "", #main = expression(beta[0]),
     freq = F, ylim = c(0, 2.5))
hist(temp3[,which(colnames(temp3) == "beta0"),], breaks = seq(-2, 0.5, 0.1),
     col = rgb(0,0,1,0.3),
     freq = F, add = T)

