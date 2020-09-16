# this script implements a toy model to
# calculate the posterior probability 
# of genotype calls for a single sample
# given read data
library(tidyverse)

# likelihood function from Li 2011
lg <- function(g,l,k,e){
	(((2-g)*e + g*(1-e))^l*((2-g)*(1-e) + g*e)^(k-l))/(2^k)
}

# prior based on expected AFS for a single sample
prio <- function(p){
	c((1 - (p + p/2)),p,p/2)
}

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
