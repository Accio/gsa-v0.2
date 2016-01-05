library(limma)
library(lattice)
library(ribiosPlot)
library(ribiosIO)

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

## is this caused only by the BH (global?) FDR correction? The answer is no: local false discovery rate estimated for index1 is 0.309, 
library(fdrtool)
fdrDepFDR <- fdrtool(fdrDep$PValue, statistic="pvalue")
head(fdrDepFDR$qval)

## can we draw a curve showing how BH FDR values change simply with the increased lengths of gene lists?
fdrAfterAppendRandom <- function(matrix, index, design, nRandom=10, use.ranks=FALSE) {
    if(nRandom>0) {
        indLen <- length(index)
        restInd <- setdiff(1:nrow(matrix), index)
        rl <- lapply(1:nRandom, function(x) sample(restInd, indLen, replace=FALSE))
        names(rl) <- sprintf("Random%d", 1:nRandom)
        index <- c(list(index=index), rl)
    } else {
        index <- list(index=index)
    }
    cameraRes <- camera(matrix, index=index, design=design, sort=FALSE, use.ranks=use.ranks)
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

(nRandoms <- as.integer(c(0, 1,2, 5,10^seq(1,4,0.125), 5*10^4,1*10^5,2*10^5)))
library(multicore)
getFDRresult <- function(matrix, index, design, nRandoms, use.ranks=FALSE) {
    randomFDRs <- mclapply(nRandoms, function(x)
        fdrAfterAppendRandom(matrix, index, design, nRandom=x, use.ranks=use.ranks))
    randomFDRresult <- data.frame(nRandom=nRandoms,
                                  p=sapply(randomFDRs, function(x) x$p),
                                  pRank=sapply(randomFDRs, function(x) x$pRank),
                                  FDR=sapply(randomFDRs, function(x) x$FDR),
                                  FDRrank=sapply(randomFDRs, function(x) x$FDRrank))
    return(randomFDRresult)
}
randomFDRresult <- getFDRresult(y, index1, design, nRandoms)



logSteps <- function(x) return(c(1,2,5,10,20, 50,100,200, 500,1000,
                                 2000,5000,10000, 20000, 50000, 1*10^5, 2*10^5))
getNbyFDR <- function(fdrResult, threshold=.10) {
    isNearest <- with(fdrResult, which.min(abs(FDR-threshold)))
    if(fdrResult$FDR[isNearest]==threshold) {
        return(fdrResult$nRandom[isNearest])
    } else if (fdrResult$FDR[isNearest]<threshold)  {
        a <- with(fdrResult, nRandom[isNearest+1])
        m <- with(fdrResult, FDR[isNearest+1])
        b <- with(fdrResult, nRandom[isNearest])
        j <- with(fdrResult, FDR[isNearest])

    } else {
        a <- with(fdrResult, nRandom[isNearest])
        m <- with(fdrResult, FDR[isNearest])
        b <- with(fdrResult, nRandom[isNearest-1])
        j <- with(fdrResult, FDR[isNearest-1])
    }
    x <- (b*(m-threshold)+a*(threshold-j))/(m-j)
    x <- as.integer(x)
    if(length(x)==0) return(0L)
    return(x)
}

plotFDR <- function(fdrResult,main="FDR of a genuine DE-geneset (camera)") {
    n <- getNbyFDR(fdrResult, thr=0.1)
    res <- xyplot(FDR~nRandom, data=fdrResult,
                  abline=list(h=c(log2(1), log2(0.10), log2(fdrResult$p[1])),v=log10(n),
                      col=c(rep("black", 2), "red"), lwd=c(rep(1,2), 1.5),
                      lty=c(rep(2,2), 1)),
                  main=main, 
                  xlab="Number of random gene sets tested", ylab="FDR (Benjamini-Hochberg)",
                  key=list(space="top", lines=list(col=c("red")),
                      text=list("pValue")),
                  scales=list(tck=c(1,0), alternating=1L,
                      y=list(log=2, at=c(1E-4,1E-3,1E-2,0.02, 0.05, 0.10, 0.25, 0.5, 1)),
                      x=list(at=logSteps(),
                          log=TRUE, rot=45)), type="b", fill="darkgray", pch=21)
    print(res)
    return(invisible(res))

}
plotFDR(randomFDRresult)
ipdf(figfile("FDR-numberOfRandomGeneSets.pdf"), width=6L, height=6L)

plotFDRrank <- function(fdrResult,main="FDR rank of a genuine DE-geneset (camera)") {
    res <- xyplot(FDRrank~nRandom, data=fdrResult,
           abline=list(a=0,b=1, lty=2),
           main=main,
           xlab="Number of random gene sets tested", ylab="FDR rank (Benjamini-Hochberg)",
           scales=list(tck=c(1,0), alternating=1L,
               y=list(log=10, at=logSteps()),
               x=list(at=c(1,10,50,100,200, 500,1000,
                          2000,5000,10000, 20000, 50000, 1*10^5, 2*10^5),
                   log=TRUE, rot=45)), type="b", fill="darkgray", pch=21)
    print(res)
    return(invisible(res))

}
plotFDRrank(randomFDRresult)
ipdf(figfile("FDRrank-numberOfRandomGeneSets.pdf"), width=6L, height=6L)

plotPrank <- function(fdrResult,main="pValue rank of a genuine DE-geneset (camera)") {
    res <- xyplot(pRank~nRandom, data=fdrResult,
           abline=list(a=0,b=1, lty=2),
           main=main, 
           xlab="Number of random gene sets tested", ylab="pValue rank",
           scales=list(tck=c(1,0), alternating=1L,
               y=list(log=10, at=logSteps()),
               x=list(at=c(1,10,50,100,200, 500,1000,
                          2000,5000,10000, 20000, 50000, 1*10^5, 2*10^5),
                   log=TRUE, rot=45)), type="b", fill="darkgray", pch=21)
        print(res)
    return(invisible(res))
}
plotPrank(randomFDRresult)
ipdf(figfile("pRank-numberOfRandomGeneSets.pdf"), width=6L, height=6L)

## export FDR simulation results
writeMatrix(randomFDRresult, outfile("randomFDRresult.txt"))

## is this because of the relative moderate change (1-sigma?)
## 2sigma
y2 <- y
y2[index1, 4:6] <- y2[index1, 4:6]+1 ## 2-sigma changes
## fdr gets better: 0.0009 
head(fdrDep2sigma <- camera(y2, c(list(set1=index1), indexes), design))
fdrDepFDR2sigma <- fdrtool(fdrDep2sigma$PValue, statistic="pvalue")
head(fdrDepFDR2sigma$qval) ## local fdr 0.005, much better!

generateFDRreport <- function(matrix, index, design, nRandoms, use.ranks=FALSE, suffix="") {
    randomFDRresult <- getFDRresult(matrix, index, design, nRandoms, use.ranks=use.ranks)
    plotFDR(randomFDRresult,
            main=sprintf("FDR of a genuine DE-geneset (%s)", suffix))
    ipdf(figfile(sprintf("FDR-numberOfRandomGeneSets-%s.pdf", suffix)),
         width=6L, height=6L)
    plotFDRrank(randomFDRresult,
                main=sprintf("FDR rank of a genuine DE-geneset (%s)", suffix))
    ipdf(figfile(sprintf("FDRrank-numberOfRandomGeneSets-%s.pdf", suffix)),
         width=6L, height=6L)
    plotPrank(randomFDRresult,
              main=sprintf("pValue rank of a genuine DE-geneset (%s)", suffix))
    ipdf(figfile(sprintf("pRank-numberOfRandomGeneSets-%s.pdf", suffix)),
         width=6L, height=6L)
    
    writeMatrix(randomFDRresult,
                outfile(sprintf("randomFDRresult-%s.txt", suffix)))
    return(invisible(randomFDRresult))

}
## nRandomsSmall <- c(0,1,2,5,10)
system.time(randomFDRresult.2sigma <- generateFDRreport(y2, index1, design, nRandoms, suffix="2sigma"))

system.time(randomFDRresult.2sigma.ranks <- generateFDRreport(y2, index1, design, nRandoms, suffix="2sigmaUsingRanks", use.ranks=TRUE))

## note: t-test takes about 32s, and wilcoxon-test takes about 106s (note by David: likely to be improved by BioQC)

## 3sigma
y3 <- y
y3[index1, 4:6] <- y3[index1, 4:6]+2 ## 3-sigma changes

randomFDRresult.3sigma <- generateFDRreport(y3, index1, design, nRandoms, suffix="3sigma")

## half of genes with 1 sigma, the rest not changing
y1half <- y
y1half[index1[11:20], 4:6] <- y1half[index1[11:20], 4:6]-1 ## 1-sigma changes

randomFDRresult.sigmaHalf <- generateFDRreport(y1half, index1, design, nRandoms, suffix="sigmaHalf")

## half of genes with 2 sigma, the rest not changing
y2half <- y
y2half[index1[1:10], 4:6] <- y2half[index1[1:10], 4:6]+1 ## 2-sigma changes

randomFDRresult.2sigmaHalf <- generateFDRreport(y2half, index1, design, nRandoms, suffix="2sigmaHalf")

## half of genes with 3 sigma, the rest not changing
y3half <- y
y3half[index1[1:10], 4:6] <- y3half[index1[1:10], 4:6]+2 ## 3-sigma changes

randomFDRresult.3sigmaHalf <- generateFDRreport(y3half, index1, design, nRandoms, suffix="3sigmaHalf")


## if we set the threshold of FDR to 0.1, what is the maximal gene sets that can be simultaneously tested in order to not to miss the genuine-DE geneset?


effects <- c("sigmaHalf", "sigma", "2sigmaHalf","2sigma", "3sigmaHalf", "3sigma")
thrs <- c(0.05, 0.10, 0.25)
fdrThres <- data.frame(effect=rep(effects, each=length(thrs)),
                       FDRthreshold=thrs,
                       N=c(sapply(thrs, function(x) getNbyFDR(randomFDRresult.sigmaHalf,x)),
                           sapply(thrs, function(x) getNbyFDR(randomFDRresult,x)),
                           sapply(thrs, function(x) getNbyFDR(randomFDRresult.2sigmaHalf,x)),
                           sapply(thrs, function(x) getNbyFDR(randomFDRresult.2sigma,x)),
                           sapply(thrs, function(x) getNbyFDR(randomFDRresult.3sigmaHalf,x)),
                           sapply(thrs, function(x) getNbyFDR(randomFDRresult.3sigma,x))))

writeMatrix(fdrThres, outfile("effect-FDRthreshold-N-table.txt"), row.names=FALSE)

## alternatively: given a number of gene sets (say 3000), what is the FDR value of the genuine-DE geneset if there is one and only one such gene set?


