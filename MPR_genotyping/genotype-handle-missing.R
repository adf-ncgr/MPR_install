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
    #h = hidden state, x = observed genotype, n => see line 13 R/hmm.vitFUN.rils.R; looks like this is tracking the length of the run of a specific genotype?
    function(h, x, n) {
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
    #if (lastRecombDist < opt$recomb) p <- 0
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
    geno.dist <- abs(diff(position))
    #track number of consecutive observations of same genotype, which is
    #used in estimating probabilities for hidden het state
    n.con[1] <- 1
    g <- geno[1]
    for (i in 2:n.obs) {
        n.con[i] <- ifelse(geno[i] == g || g == 0, n.con[i - 1] + 1, 1)
        g <- geno[i]
    }
    #state probabilities of starting point of path
    #delta is probability of most likely path to this point leading to the given state for that point in the observation set
    #psi gives the previous state on the most likely path leading to the given state
    for (i in 1:n.state) {
        delta[i, 1] <- log(geno.probability[i]) + emissionFUN(i, geno[1], n.con[1])
        psi.con[i,1] <- 1
        #we want to start out free to change if needed
        psi.dist[i,1] <- opt$recomb
    }
    #probabilities of getting to j state from any previous state, including transition costs
    preProb <- numeric(n.state)
    for (t in 2:n.obs) {
        for (j in 1:n.state) {
            for (i in 1:n.state) {
                preProb[i] <- delta[i, t - 1] + 
                log(transitionFUN(i, j, geno.dist[t - 1], psi.dist[j,t-1]))
                #print(sprintf("i=%i,delta=%f,trans=%f",i, delta[i,t-1],log(transitionFUN(i, j, geno.dist[t - 1], psi.dist[j,t-1]))));
                #print(sprintf("i=%i,preProb=%f",i, preProb[i]));
            }
            psi[j, t] <- which.max(preProb)
            if (psi[j,t] == j) {
                psi.con[j,t] = psi.con[j,t-1] + 1
                psi.dist[j,t] = psi.dist[j,t-1] + geno.dist[t-1]
                #psi.dist[j,t] = psi.dist[j,t-1] + geno.dist[t]
            }
            else {
                #psi.con[1:n.state,t] = 1
                #psi.dist[1:n.state,t] = 1
                psi.con[j,t] = 1
                psi.dist[j,t] = 1
            }
            delta[j, t] <- max(preProb) + emissionFUN(j, geno[t], 
                n.con[t])
            #print(sprintf("t=%i,j=%i,delta=%f,psi=%i,psi.dist=%i",t,j,delta[j,t],psi[j,t],psi.dist[j,t]));
        }
    }
    geno.cr[n.obs] <- which.max(delta[, n.obs])
    for (t in seq(n.obs - 1, 1, by = -1)) geno.cr[t] <- psi[geno.cr[t + 
        1], t + 1]

    i<-1
    last_g <- geno.cr[i];
    while (i <= n.obs) {
       if (geno.cr[i] != last_g) {
            j <- i
            new_g <- geno.cr[i]
            while (geno.cr[j] == new_g && j < n.obs) {j <- j+1}
            if (position[j] - position[i] < opt$recomb) {
                geno.cr[i:j] <- last_g
            }
            else {
                last_g = new_g;
            }
            i <- j
       }
       i<-i+1
    }
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
    O.cr <- hmm.vitFUN.suppressQuickRecomb(geno=O,position=O.pos,geno.probability=c(opt$het/2, opt$het/2, opt$het),transitionFUN =myTransitionFun, emissionFUN = makeMissingDataEmissionFUN(errorRate = opt$err))
    #write.table(O.cr)
    t[,i]<-O.cr
}
cat("#VERSION postprocessing_recombination_initial-2-ga552139-8;", " recomb=",opt$recomb, ", het=", opt$het, ", error=", opt$err, "\n");
write.table(t, file="", sep="\t", quote=F, row.names = FALSE, col.names=TRUE)
