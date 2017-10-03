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
    #cost of getting it wrong
    E <- log(errorRate)
    #cost of getting it right
    E2 <- log(1 - errorRate)
    #cost of making successive HOM calls of the same kind when you are in HET region
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
myTransitionFun <- function(fromState, toState, physicalDistance, lastRecombDist) {
    d <- physicalDistance/(100000000000*100)
    p <- (1 - exp(-2 * d))
    p <- p/(1 + p)
    #write(c("lastRecombNobs = ", lastRecombNobs), stderr())
    if (lastRecombDist < opt$recomb) p <- 0
    #p=0.000000001
    ifelse(fromState == toState, 1 - p, p)
}
#equivalent to original version packaged with MPR, but arguments for being called in suppression context
myTransitionFun2 <- function(fromState, toState, physicalDistance, lastRecombNobs) {
    d <- physicalDistance/(._rice_phy2get_factor_ * 100)
    #d <- physicalDistance/(100000000000*100)
    p <- (1 - exp(-2 * d))
    p <- p/(1 + p)
    #write(c("lastRecombNobs = ", lastRecombNobs), stderr())
    #if (lastRecombNobs < opt$recomb) p <- 0
    #p=0.000000001
    ifelse(fromState == toState, 1 - p, p)
}
#reinstitutes the key difference line proving that the behavior is caused by line
#if (lastRecombNobs < opt$recomb) p <- 0
myTransitionFun3 <- function(fromState, toState, physicalDistance, lastRecombNobs) {
    d <- physicalDistance/(._rice_phy2get_factor_ * 100)
    #d <- physicalDistance/(100000000000*100)
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
    psi.dist <- matrix(0, nrow = n.state, ncol = n.obs)
    n.con <- geno.cr <- numeric(n.obs)
    n.dist <- numeric(n.obs)
    geno.dis <- abs(diff(position))
    n.con[1] <- 1
    n.dist[1] <- 1
    g <- geno[1]
    for (i in 2:n.obs) {
        n.con[i] <- ifelse(geno[i] == g, n.con[i - 1] + 1, 1)
        #augment consecutive obs counts used in emission prob with distances for transition prob
        n.dist[i] <- ifelse(geno[i] == g, n.dist[i - 1] + geno.dis, 1)
        g <- geno[i]
    }
    for (i in 1:n.state) {
        delta[i, 1] <- log(geno.probability[i]) + emissionFUN(i, geno[1], n.con[1])
        psi.con[i,1] <- 1
        psi.dist[i,1] <- 1
    }
    preProb <- numeric(n.state)
    for (t in 2:n.obs) {
        for (j in 1:n.state) {
            for (i in 1:n.state) preProb[i] <- delta[i, t - 1] + 
                #adf: extra arg is for myTransitionFun
                log(transitionFUN(i, j, geno.dis[t - 1], psi.dist[j,t-1]))
                #log(transitionFUN(i, j, geno.dis[t - 1]))
            psi[j, t] <- which.max(preProb)
            if (psi[j,t] == j) {
                psi.con[j,t] = psi.con[j,t-1] + 1
                psi.dist[j,t] = psi.dist[j,t-1] + geno.dis[t-1]
            }
            else {
                psi.con[1:n.state,t] = 1
                psi.dist[1:n.state,t] = 1
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
    #adf: this was the version we started with when examining the recombination pileup issue
    #O.cr <- hmm.vitFUN.suppressQuickRecomb(geno=O,position=O.pos,geno.probability=c(opt$het/2, opt$het/2, opt$het),transitionFUN =myTransitionFun, emissionFUN = makeMissingDataEmissionFUN(errorRate = opt$err))
    #adf: this was the earlier version we tested and found worked much better for the recombnation pileup issue
    #O.cr <- hmm.vitFUN.rils(geno=O,position=O.pos,geno.probability=c(opt$het/2, opt$het/2, opt$het),emissionFUN = makeMissingDataEmissionFUN(errorRate = opt$err))
    #adf: this users suppressQuickRecomb but not the special transition function: myTransitionFun; it produced results equivalent to not using suppressQuickRecomb at all.
    #O.cr <- hmm.vitFUN.suppressQuickRecomb(geno=O,position=O.pos,geno.probability=c(opt$het/2, opt$het/2, opt$het),emissionFUN = makeMissingDataEmissionFUN(errorRate = opt$err))
    #adf: as expected, this produces output identical to "theoldway" (because it is the old way)
    #O.cr <- hmm.vitFUN.suppressQuickRecomb(geno=O,position=O.pos,geno.probability=c(opt$het/2, opt$het/2, opt$het),transitionFUN =myTransitionFun, emissionFUN = makeMissingDataEmissionFUN(errorRate = opt$err))
    O.cr <- hmm.vitFUN.suppressQuickRecomb(geno=O,position=O.pos,geno.probability=c(opt$het/2, opt$het/2, opt$het),transitionFUN =myTransitionFun, emissionFUN = makeMissingDataEmissionFUN(errorRate = opt$err))
    #write.table(O.cr)
    t[,i]<-O.cr
}
write.table(t, file="", sep="\t", quote=F, row.names = FALSE, col.names=TRUE)
