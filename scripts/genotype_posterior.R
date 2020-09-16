# this script implements a toy model to
# calculate the posterior probability 
# of genotype calls for a single sample
# given read data
library(tidyverse)

# likelihood function for biallelic variation modified from equation 2 in Li 2011 
# to have a single error rate for all bases and ploidy = 2
	# g - genotype in number of alternate alleles (e.g. 0,1,2)
	# l - # reference alleles
	# k - total number of alleles sampled, a.k.a coverage of the site (e.g. # alt alleles = k - l )
	# e - base call error rate, e.g. 0.01 (Q20), 0.001 (Q30)...

	# returns the probability of the data given g and e. 
lg <- function(g,l,k,e){
	(((2-g)*e + g*(1-e))^l*((2-g)*(1-e) + g*e)^(k-l))/(2^k)
}

# prior probability distribution based on expected AFS for a single sample
	# p = theta, or expected number of pairwise differences
	# in humans this is 0.001, in drosophila it's more like 0.01
prio <- function(p){
	c((1 - (p + p/2)),p,p/2)
}

# posterior probability of genotype given data, prior and error rate
	# gg = g from lg
	# ll = l from lg
	# kk = k from lg
	# ee = e from lg
	# pp = p from prio
plg <- function(gg,ll,kk,ee,pp){

	lg(gg,ll,kk,ee)*prio(pp)/sum(lg(gg,ll,kk,ee)*prio(pp))

}


n <- 15
posteriors <- c()
for(i in 0:n){

	posteriors <- rbind(posteriors,
						plg(c(0,1,2),i,n,0.001,0.01)
					)

}

plot(0:n,posteriors[,1],type="b",lwd=2,lty=3,xlab="# alt alleles",ylab="probability of genotype")
lines(0:n,posteriors[,2],type="b",lwd=2,lty=3,col="gray")
lines(0:n,posteriors[,3],type="b",lwd=2,lty=3,col="blue")


rbinom(n=1000,prob=0.5,size=15) %>% table() %>% plot()

(dbinom(x=0:n,size=15,prob=0.5)*1000) %>% points(x=0:n,y=.)
