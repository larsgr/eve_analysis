library(mvtnorm)
library(tidyverse)
library(ape)

source("scripts/dvdt.R")
source("scripts/EVEcore.R")
source("scripts/twoThetaTest.R")


# Simulate data without covariances, just define between-species and within-species variance
# Use the phylogeny so that the 


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



t2res <-
  twoThetaTest(tree = tree,
               gene.data=dat,
               shiftSpecies = c("Salp","Okis","Omyk","Ssal"))


plot(sort(t2res$LRT),ylab="LRT")
lines(qchisq(p = seq(0,1,length.out = K), df = 1),col="red")
lines(qchisq(p = seq(0,1,length.out = K), df = 2),col="blue")
legend("topleft",lty=1,col=c("red","blue"), legend=c("chi-squared df = 1","chi-squared df = 2"))
