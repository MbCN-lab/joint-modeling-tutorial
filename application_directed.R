rm(list=ls())

hrf.conv = function(t, beta, onset, prmt = c(6, 16, 1, 1, 1/6)){
  a1 = prmt[1]
  a2 = prmt[2]
  b1 = prmt[3]
  b2 = prmt[4]
  c = prmt[5]
  
  t = t - onset
  t[t<0] = 0
  
  return(beta * ( ((t^(a1-1) * b1^a1 * exp(-b1 * t) )/gamma(a1) )
                 - c * (((t^(a2 - 1) * b2^a2 * exp(-b2 * t) )/gamma(a2))) ))
}

# Load required packages and modules
require("rjags")
load.module("wiener")

###############################
# Load and extract behavioral log
setwd("~/Dropbox/LabProject/JointModelingTutorial/201712_newdata")
load("application_dataset.Rdata")

# Recode data
rt[resp == 0] = rt[resp == 0] * -1

# Data
TR = 2

lenS = length(onset) # Total number of stiuli presented in the block

dat = list(N = N, lenN = length(N), TR = TR,
           t = rt, n.trials = length(rt),
           onset = onset, lenS = lenS,
           a1 = 6, a2 = 16, b1 = 1, b2 = 1, c = 1/6)


model.double.gamma.wiener = "
model{
  # Likelihood
  for (i in 1:lenN) {
  N[i] ~ dnorm(muN[i], inv.sigma.sq)
  Npred[i] ~ dnorm(muN[i], inv.sigma.sq)
  muN[i] = beta0 + inprod(beta[], X[i, ])
  }

  # Define a design matrix
  for (i in 1:lenS){
    for (j in 1:lenN){
      temp[j,i] = (j-1) * TR - onset[i]
      Xt[j,i] = ifelse(temp[j,i] >= 0, temp[j,i], 0)
      X[j,i] = (Xt[j,i]^(a1-1) * (b1)^(a1) * exp(-b1*Xt[j,i]) / exp(loggam(a1))) - c * (Xt[j,i]^(a2-1) * (b2)^(a2) * exp(-b2*Xt[j,i]) / exp(loggam(a2)))
    }
  }

  for (i in 1:n.trials){
    xi[i] = beta[2*i] - beta[2*i-1]
    t[i] ~ dwiener(alpha, tau, omega, xi[i])
  }

  # Prior
  # The neural submodel
  inv.sigma.sq ~ dgamma(.001, .001)
  sigma.sq = 1/inv.sigma.sq
  beta0 ~ dnorm(0, 0.001)
  for (j in 1:lenS){
    beta[j] ~ dnorm(0, 0.001)
  }
  # The behavioral submodel
  alpha ~ dunif(0.0001, 10)
  tau ~ dunif(0, 0.04)
  omega = 0.5
}
"


# Initialization
model.dgw = jags.model(textConnection(model.double.gamma.wiener), data = dat,
                      n.chains = 3, n.adapt = 2000)

# Burn-in
update(model.dgw, n.iter = 4000, progress.bar = "text") 

# Posterior sampling
dgw.out = coda.samples(model = model.dgw,
                       variable.names = c("beta0", "beta", "sigma.sq", "Npred",
                                          "alpha", "tau", "xi"), 
                       n.iter = 6000) 

out = summary(dgw.out)
save(dgw.out, file = "dgwout.Rdata")
###########################
# Extract estimates for plotting
idx = list()
idx$Npred = seq(which(rownames(out$quantiles) == "Npred[1]"),
                which(rownames(out$quantiles) == paste("Npred[",length(N),"]", sep=""), 1))
idx$beta = seq(which(rownames(out$statistics) == "beta[1]"),
               which(rownames(out$statistics) == paste("beta[",lenS,"]", sep=""), 1))
idx$beta0 = which(rownames(out$statistics) == "beta0")
idx$sigma.sq = which(rownames(out$statistics) == "sigma.sq")
idx$alpha = which(rownames(out$statistics) == "alpha")
idx$tau = which(rownames(out$statistics) == "tau")
idx$xi = seq(which(rownames(out$statistics) == "xi[1]"),
             which(rownames(out$statistics) == paste("xi[",lenS/2,"]", sep=""), 1))

# Recover BOLD responses from beta
ts = seq(0, (length(N)-1) * TR, TR)
X.estim = mapply(hrf.conv, out$statistics[idx$beta, 1], onset = onset, MoreArgs = list(t = ts))
###########################
# Plot
## (1) BOLD recovery
ylim.BOLD = c(floor(min(min(rowSums(X.estim)), min(N), min(out$quantiles[idx$Npred, c(1, 5)]))),
              ceiling(max(max(rowSums(X.estim)), max(N), max(out$quantiles[idx$Npred, c(1, 5)]))))

layout(matrix(c(2,1,0,3), 2, 2, byrow = T),
       widths = c(1, 9), heights = c(8, 1))
par(mar = c(1.8, 1.5, 0.2, 1.5))
plot(N, pch = 16, cex = 1.2, cex.axis = 2, ylim = ylim.BOLD)
lines(rowSums(X.estim) + out$statistics[idx$beta0,1], col = "red", lwd = 2)
matlines(out$quantiles[idx$Npred, c(1, 5)], lty = 3, col = "red")

par(mar = c(0,0,0,0))
plot(0, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white")
text(0, "BOLD Response", srt = 90, cex = 3)

par(mar = c(0,0,0,0))
plot(0, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white")
text(0, "Time (TR)", cex = 3)

## (2) RT vs xi
layout(matrix(c(2,1,0,3), 2, 2, byrow = T),
       widths = c(1, 9), heights = c(8, 1))
par(mar = c(1.5, 1.5, 0.2, 1.5))
plot(abs(rt), out$statistics[idx$xi,1], col = "white", cex.axis = 2.1,
     xlim = c(0.4, 1.2),
     #xlim = c(0, round(max(abs(rt)), 1)),
     ylim = c(-6, 6))
     #ylim = c(floor(min(out$statistics[idx$xi,1])), ceiling(max(out$statistics[idx$xi,1]))))
points(abs(rt[resp == 0]), out$statistics[idx$xi,1][resp == 0], pch = 16, cex = 2)
points(abs(rt[resp == 1]), out$statistics[idx$xi,1][resp == 1], pch = 3, cex = 2)
abline(h = 0, lty = 3)
legend(inset = 0.025, "topright", bty = "n", cex = 2,
       pch = c(16, 3), pt.cex = 2,
       c("Respond 1st Stimulus", "Respond 2nd Stimulus"))

par(mar = c(0,0,0,0))
plot(0, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white")
text(0, expression(xi[i]), srt = 90, cex = 3)

plot(0, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white")
text(0, "Response time", cex = 3)

## (3) contrast vs beta (not included in the text)
beta.mean = out$statistics[idx$beta,1]
temp = as.data.frame(cbind(stim, beta.mean))
beta.cond.mean = aggregate(beta.mean ~ stim, data = temp, FUN = mean)

layout(matrix(c(2,1,0,3), 2, 2, byrow = T),
       widths = c(1, 6), heights = c(6, 1))
par(mar = c(1, 1, 0.2, 0.2))
plot(stim, beta.mean, pch = 16, col = rgb(0,0,0,0.3), cex = 2,
     xlab = "Contrast level", ylab = "Beta")
points(beta.cond.mean[,1], beta.cond.mean[,2], pch = 4, col = "red", lwd = 5)
legend(inset = 0.01, "bottomright", cex = 1.5,
       pch = c(16, 4), col = c(rgb(0,0,0,0.3), "red"), pt.cex = c(2, 1), lwd = c(0, 5), lty = NA,
       title = "Note: Not included in the main text",
       c("Single-trial beta", "Condition mean"))


par(mar = c(0,0,0,0))
plot(0, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white")
text(0, expression(beta[i]), srt = 90, cex = 3)

plot(0, bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white")
text(0, "Contrast level", cex = 3)