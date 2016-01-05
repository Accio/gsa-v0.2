library(limma)
library(lattice)
library(ribiosPlot)

figfile <- function(x) file.path("figures", x)
outfile <- function(x) file.path("output", x)

set.seed(1887)

y <- matrix(rnorm(1000*6),1000,6)
design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))

## index1: genuinely significantly up-regulated
index1 <- 1:20
y[index1,4:6] <- y[index1,4:6]+1
## index2: not DEs
index2 <- 21:40

camera(y, index1, design)
camera(y, index2, design)
camera(y, list(set1=index1,set2=index2), design) ## FDR 0.003 versus 0.156 - set1 is significant

## however the FDR is quite sensitive to number of genesets that are tested. In fact if a long list of randomly selected gene sets is tested along with index1, the FDR value of index1 will drop (see below), though it is still ranked the highest

gsN <- 100
indexes <- lapply(1:gsN, function(x) sample(21:1000, 20))
names(indexes) <- sprintf("RandomGeneSet%d", 1:gsN)
head(fdrDep <- camera(y, c(list(set1=index1), indexes), design))

## is this caused only by the BH (global?) FDR correction? The answer is no: local false discovery rate estimated for index1 is 0.12, 
library(fdrtool)
fdrDepFDR <- fdrtool(fdrDep$PValue, statistic="pvalue")
head(fdrDepFDR$qval)

## can we draw a curve showing how BH FDR values change simply with the increased lengths of gene lists?
fdrAfterAppendRandom <- function(matrix, index, design, nRandom=10) {
    if(nRandom>0) {
        indLen <- length(index)
        restInd <- setdiff(1:nrow(matrix), index)
        rl <- lapply(1:nRandom, function(x) sample(restInd, indLen, replace=FALSE))
        names(rl) <- sprintf("Random%d", 1:nRandom)
        index <- c(list(index=index), rl)
    } else {
        index <- list(index=index)
    }
    cameraRes <- camera(matrix, index=index, design=design, sort=FALSE)
    if(nRandom==0) {
        p <- fdr <- cameraRes["index", "PValue"] ## in case of single gene set no FDR correction is made
        pRank <- fdrRank <- 1
    } else {
        p <- cameraRes["index", "PValue"]
        fdr <- cameraRes["index", "FDR"]
        pRank <- sum(cameraRes$PValue<=p)
        fdrRank <- sum(cameraRes$FDR<=fdr)
    }
    return(list(p=p, pRank=pRank, FDR=fdr, FDRrank=fdrRank))
}

(nRandoms <- as.integer(c(0, 1,5,10^seq(1,4,0.125), 5*10^4,1*10^5,2*10^5)))
library(multicore)
randomFDRs <- mclapply(nRandoms, function(x)
    fdrAfterAppendRandom(y, index1, design, nRandom=x))
randomFDRresult <- data.frame(nRandom=nRandoms,
                              p=sapply(randomFDRs, function(x) x$p),
                              pRank=sapply(randomFDRs, function(x) x$pRank),
                              FDR=sapply(randomFDRs, function(x) x$FDR),
                              FDRrank=sapply(randomFDRs, function(x) x$FDRrank))
logSteps <- c(1,2,5,10,50,100,200, 500,1000,
              2000,5000,10000, 20000, 50000, 1*10^5, 2*10^5)
xyplot(FDR~nRandom, data=randomFDRresult,
       abline=list(h=c(log2(1), log2(0.10), log2(randomFDRresult$p[1])),
           col=c(rep("black", 2), "red"), lwd=c(rep(1,2), 1.5), lty=c(rep(2,2), 1)),
       main="FDR of a genuine DE-geneset (camera)", 
       xlab="Number of random gene sets tested", ylab="FDR (Benjamini-Hochberg)",
       scales=list(tck=c(1,0), alternating=1L,
           y=list(log=2, at=c(1E-4,1E-3,1E-2,0.02, 0.05, 0.10, 0.25, 0.5, 1)),
           x=list(at=logSteps,
               log=TRUE), rot=45), type="b", fill="darkgray", pch=21)
ipdf(figfile("FDR-numberOfRandomGeneSets.pdf"), width=6L, height=6L)

xyplot(FDRrank~nRandom, data=randomFDRresult,
       abline=list(a=0,b=1, lty=2),
       main="FDR rank of a genuine DE-geneset (camera)", 
       xlab="Number of random gene sets tested", ylab="FDR rank (Benjamini-Hochberg)",
       scales=list(tck=c(1,0), alternating=1L,
           y=list(log=10, at=logSteps),
           x=list(at=c(1,10,50,100,200, 500,1000,
                      2000,5000,10000, 20000, 50000, 1*10^5, 2*10^5),
               log=TRUE), rot=45), type="b", fill="darkgray", pch=21)
ipdf(figfile("FDRrank-numberOfRandomGeneSets.pdf"), width=6L, height=6L)

xyplot(pRank~nRandom, data=randomFDRresult,
       abline=list(a=0,b=1, lty=2),
       main="pValue rank of a genuine DE-geneset (camera)", 
       xlab="Number of random gene sets tested", ylab="pValue rank",
       scales=list(tck=c(1,0), alternating=1L,
           y=list(log=10, at=logSteps),
           x=list(at=c(1,10,50,100,200, 500,1000,
                      2000,5000,10000, 20000, 50000, 1*10^5, 2*10^5),
               log=TRUE), rot=45), type="b", fill="darkgray", pch=21)
ipdf(figfile("pRank-numberOfRandomGeneSets.pdf"), width=6L, height=6L)
