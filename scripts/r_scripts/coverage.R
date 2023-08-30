library(tidyverse)

# focal region is: NC_000020.11:29400000-34400000

# read in short read coverage table
sc <- read.table("../../results/03_AlignmentAndCoverage/short_read_coverage/coverage_1kb.bed.gz") %>%
	mutate(V4=as.numeric(V4),V5=as.numeric(V5))

# filter down to focal region	
sc <- sc %>% filter(V1=="NC_000020.11" & V2 >= 29400000 & V3 <= 34400000)

# plot full range of coverage in focal region
plot(sc[,2],sc[,4],pch=20,cex=.2)

# zoom in in y-axis to see a little better
plot(sc[,2],sc[,4],pch=20,cex=.2,ylim=c(0,1000))
abline(h=c(120,250),col="red")

# plot a histogram of coverage
hist(sc[,4],breaks=1000,xlim=c(0,1000))
abline(v=c(120,250),col="red")

# based on eye-balling the graph, avoid variant calling using short reads in windows where mean coverage < 120 or >250


# have a look at long read coverage as well. 
# read in long read coverage table
lc <- read.table("../../results/03_AlignmentAndCoverage/long_read_coverage/coverage_1kb.bed.gz") %>%
	mutate(V4=as.numeric(V4),V5=as.numeric(V5))

# filter down to focal region	
lc <- lc %>% filter(V1=="NC_000020.11" & V2 >= 29400000 & V3 <= 34400000)

# plot full range of coverage in focal region
plot(lc[,2],lc[,4],pch=20,cex=.2)

# zoom in in y-axis to see a little better
plot(lc[,2],lc[,4],pch=20,cex=.2,ylim=c(0,450))
abline(h=c(170,270),col="red")

# plot a histogram of coverage
hist(lc[,4],breaks=1000,xlim=c(0,500))
abline(v=c(170,270),col="red")


# plot long read and short read coverage against each other
plot(sc[,4],lc[,4],pch=20,cex=.2,xlim=c(0,400),ylim=c(0,400))


# plot together
par(mfrow=c(2,1),mar=c(0,0,0,0),oma=c(3,3,1,1))

plot(sc[,2],sc[,4],pch=20,cex=.2,ylim=c(0,450))
abline(h=c(120,250),col="red")

plot(lc[,2],lc[,4],pch=20,cex=.2,ylim=c(0,450))
abline(h=c(170,270),col="red")




