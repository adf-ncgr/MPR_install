#!/bin/env Rscript
library(Biobase)
library(MPR)
library(optparse)
option_list = list(
    make_option(c("--het"), type="double", default=0.5, dest="het",
            help="probability of heterozygous genotype", metavar="0 <= x <= 1"),
    make_option(c("--err"), type="double", default=0.01, dest="err",
            help="probability of genotype error", metavar="0 <= x <= 1"),
    make_option(c("--recomb"), type="integer", default=100, dest="recomb",
            help="number of observations from last recombination until another recombination event is allowed", metavar="positive integer")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

`makeMissingDataEmissionFUN` <-
function (errorRate = 0.01) 
{
    E <- log(errorRate)
    E2 <- log(1 - errorRate)
    E3 <- log(0.5)
    #not sure what I was trying to accomplish with this!?
    #E3 <- log(0.75)
    #h = hidden state, x = observed genotype, n => see line 13 R/hmm.vitFUN.rils.R; looks like this is tracking the length of the run of a specific genotype?
    function(h, x, n) {
        #this was h before; this configuration seems to do well with missing data, but I am no longer getting het calls at all.
        #if (x == 0 || x == 3) 
            #return(n * E3)
        #this config produces what appears to be a very similar if not identical result to the original configuration; too many het calls, esp in the case of many missing data.
        #if (h == 0 || h == 3) 
            #return(n * E3)
#produces way too many hets, now starting to see that the n * E3 is basically trying to account for the likelihood of seeing a stretch of consecutive homozygous genotypes of the same parental allele if the state is really equally likely to produce either
#        if (h == 0 || h == 3) 
#            return(E3)
	#having trouble staying in het state, even when setting het prob to 1.0!
        #if (h == 3) 
        if (h == 3 && (x == 1 || x == 2)) 
            return(n * E3)
        else if (x == 0) 
            return(E3)
        else 
            return(ifelse(h == x, E2, E))
    }
}
myTransitionFun <- function(fromState, toState, physicalDistance, lastRecombNobs) {
    d <- physicalDistance/(100000000000*100)
    p <- (1 - exp(-2 * d))
    p <- p/(1 + p)
    #write(c("lastRecombNobs = ", lastRecombNobs), stderr())
    if (lastRecombNobs < opt$recomb) p <- 0
    #p=0.000000001
    ifelse(fromState == toState, 1 - p, p)
}
`hmm.vitFUN.suppressQuickRecomb` <-
function (geno, position, geno.probability, transitionFUN = phy2get.haldane.rils, 
    emissionFUN = makeEmissionFUN(errorRate = 0.01), ...) 
{
    n.obs <- length(geno)
    n.state <- length(geno.probability)
    #psi.con tracks the number of consecutive observations since the last recombination
    psi <- delta <- psi.con <- matrix(0, nrow = n.state, ncol = n.obs)
    n.con <- geno.cr <- numeric(n.obs)
    geno.dis <- abs(diff(position))
    n.con[1] <- 1
    g <- geno[1]
    for (i in 2:n.obs) {
        n.con[i] <- ifelse(geno[i] == g, n.con[i - 1] + 1, 1)
        g <- geno[i]
    }
    for (i in 1:n.state) {
        delta[i, 1] <- log(geno.probability[i]) + emissionFUN(i, geno[1], n.con[1])
        psi.con[i,1] <- 1
    }
    preProb <- numeric(n.state)
    for (t in 2:n.obs) {
        for (j in 1:n.state) {
            for (i in 1:n.state) preProb[i] <- delta[i, t - 1] + 
                log(transitionFUN(i, j, geno.dis[t - 1], psi.con[j,t-1]))
            psi[j, t] <- which.max(preProb)
            if (psi[j,t] == j) {
                psi.con[j,t] = psi.con[j,t-1] + 1
            }
            else {
                psi.con[1:n.state,t] = 1
            }
            delta[j, t] <- max(preProb) + emissionFUN(j, geno[t], 
                n.con[t])
        }
    }
    geno.cr[n.obs] <- which.max(delta[, n.obs])
    for (t in seq(n.obs - 1, 1, by = -1)) geno.cr[t] <- psi[geno.cr[t + 
        1], t + 1]
    geno.cr
}
#debug(hmm.vitFUN.rils)
#trace("hmm.vitFUN.rils", quote(if (T) {browser()}), at=12, print=T)
#trace("hmm.vitFUN.rils", quote(if (!(which.max(preProb))) {browser()}), at=12, print=F)
t<-read.table("stdin", header=T, check.names=F)
O.pos<-t[,2]
for (i in 3:ncol(t)) {
#for (i in 3:4) {
    O<-t[,i]
    #increasing error rate to try to prevent miscalls from affecting overall genotype blocks didn't seem to prevent switching states frequently, in fact it seemed to exacerbate it
    #O.cr <- hmm.vitFUN.rils(geno=O,position=O.pos,geno.probability=c(0.25, 0.25,0.50),transitionFUN =phy2get.haldane.rils, emissionFUN = makeMissingDataEmissionFUN(errorRate = 0.1))
    #O.cr <- hmm.vitFUN.rils(geno=O,position=O.pos,geno.probability=c(0.25, 0.25,0.50),transitionFUN =myTransitionFun, emissionFUN = makeMissingDataEmissionFUN(errorRate = 0.05))
    #this was working well except for the double-recombinant issue
    O.cr <- hmm.vitFUN.suppressQuickRecomb(geno=O,position=O.pos,geno.probability=c(opt$het/2, opt$het/2, opt$het),transitionFUN =myTransitionFun, emissionFUN = makeMissingDataEmissionFUN(errorRate = opt$err))
    #write.table(O.cr)
    t[,i]<-O.cr
}
write.table(t, file="", sep="\t", quote=F, row.names = FALSE, col.names=TRUE)
