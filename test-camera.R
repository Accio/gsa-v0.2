library(limma)

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
        fdr <- cameraRes["index", "PValue"] ## in case of single gene set no FDR correction is made
        fdrRank <- 1
    } else {
        fdr <- cameraRes["index", "FDR"]
        fdrRank <- sum(cameraRes$FDR<=fdr)
    }
    return(list(FDR=fdr, FDRrank=fdrRank))
}

(nRandoms <- as.integer(c(0, 1,5,10,20,50,10^seq(2,6,0.125))))
library(multicore)
randomFDRs <- mclapply(nRandoms, function(x)
    fdrAfterAppendRandom(y, index1, design, nRandom=x))
