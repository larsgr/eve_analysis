library(mvtnorm)
library(tidyverse)
library(ape)

library(evemodel)



# Simulate data without covariances, just define between-species and within-species variance
# The phylogeny allows the two-theta model to fit the species mean with two additional variables 
# (theta2 and alpha).
# The one-theta model can use the alpha to fit any coincidental covariance


tree <- read.tree(file = "data/salmonData/salmonidsPhylo.newick")

n = 4 # n per species
N = Ntip(tree) # number of species

bVar = 1 # between species var
wVar = 1 # within species var

K = 500 # number of "genes"


# generate a multivariate distribution
iExp = rep(1:N,each=n)
sigma = diag(rep(bVar,N))[iExp,iExp] + diag(rep(wVar,N*n))
dat <- rmvnorm(K,sigma = sigma)

colnames(dat) <- rep(tree$tip.label,each=n)
rownames(dat) <- paste0("gene",1:K)

## Simulate with covariance


simEVE <- function(theta, sigma2, alpha, beta, tree, index.expand, n)
{
  # Define the expected species mean for the root and the evolutionary variance for the root
  expected.species.mean.root <- theta
  evol.var.root <- sigma2 / (2 * alpha)
  
  
  expression.var <- calcExpVarOUconst(tree = tree, theta = theta, alpha = alpha, sigma2 = sigma2)
  covar.matrix <- calcCovMatOU(tree, alphas = alpha, evol.variance = expression.var$evol.variance)
  
  expanded.matrix <- expandECovMatrix(expected.mean = expression.var$expected.mean,covar.matrix,
                                      sigma2, alpha, beta, index.expand)
  
  # Simulate values
  rmvnorm(n = n, mean = expanded.matrix$expected.mean, sigma = expanded.matrix$cov.matr)
}

dat2 <- simEVE(theta=0,
               alpha=42,
               sigma2=2*bVar*alpha,
               beta=wVar/bVar,
               tree=tree,
               index.expand=iExp,
               n=K)

colnames(dat2) <- rep(tree$tip.label,each=n)
rownames(dat2) <- paste0("gene",1:K)

t2res2 <-
  twoThetaTest(tree = tree,
               gene.data=dat2,
               shiftSpecies = c("Salp","Okis","Omyk","Ssal"))

# t2res <-
#   twoThetaTest(tree = tree,
#                gene.data=dat,
#                shiftSpecies = c("Salp","Okis","Omyk","Ssal"))


#save(t2res, file="t2res.RData")
load("t2res.RData")

plot(sort(t2res$LRT),ylab="LRT", type="l", lwd=2)
lines(sort(t2res2$LRT),ylab="LRT", lwd=2, lty=2)
lines(qchisq(p = seq(0,1,length.out = K), df = 1)*2,col="red")
lines(qchisq(p = seq(0,1,length.out = K), df = 2),col="blue")
legend("topleft",lty=1,col=c("red","blue","black"), legend=c("chi-squared df = 1","chi-squared df = 2","LRT simulated"))


