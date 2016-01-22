library(MASS)
library(testthat)
library(limma)

## from covariance matrix to correlation matrix (column-wise)
myCorr <- function(mat) {
    matVar <- var(mat)
    sigma <- diag(diag(matVar)^(-.5), nrow=ncol(mat), ncol=ncol(mat))
    return(sigma %*% matVar %*% sigma)
}

myMatrix <- mvrnorm(1000, mu=rep(0,2), Sigma=matrix(c(10,2,3,10), nrow=2, byrow=TRUE))
expect_equal(cor(myMatrix), myCorr(myMatrix))

simulateTwoGroupData <- function(Ngenes=1000, Nsamples=c(3,3),
                                 indGeneSet=1:20, 
                                 deltaMean=1, cor=0, seed=1887) {
    set.seed(seed)
    stopifnot(length(Nsamples)==2)
    mat <- matrix(rnorm(Ngenes*sum(Nsamples)), nrow=Ngenes)
    rownames(mat) <- paste("gene", 1:Ngenes, sep="_")
    NgeneSet <- length(indGeneSet)
    isGroup2 <- 1:ncol(mat) > Nsamples[1]
    if(all(cor==0)) {
        mat[indGeneSet,isGroup2] <-     mat[indGeneSet,isGroup2]+deltaMean
    } else {
        if(length(cor)==1)
            cor <- rep(cor, 2)
        sigmaGroup1 <- matrix(cor[1], nrow=NgeneSet, ncol=NgeneSet)
        diag(sigmaGroup1) <- 1
        sigmaGroup2 <- matrix(cor[2], nrow=NgeneSet, ncol=NgeneSet)
        diag(sigmaGroup2) <- 1
        mvrnormGroup1 <- t(mvrnorm(n=Nsamples[1],
                                   mu=rep(0, NgeneSet),
                                   Sigma=sigmaGroup1))
        mvrnormGroup2 <- t(mvrnorm(n=Nsamples[2],
                                   mu=rep(deltaMean, NgeneSet),
                                   Sigma=sigmaGroup2))
        mat[indGeneSet,] <- cbind(mvrnormGroup1, mvrnormGroup2)
    }
    return(mat)
}

test <- simulateTwoGroupData(Nsamples=c(3,3), cor=c(0.25, 0.25))
testDesign <- matrix(c(rep(1,6), rep(0,3), rep(1,3)), nrow=6, byrow=FALSE)
camera(test, index=1:20, design=testDesign, contrast=2L)

myTests <- lapply(1:1000, function(x)
    simulateTwoGroupData(Nsamples=c(3,3), deltaMean=1, cor=c(0.5, 0.5), seed=x))
myCamera <- sapply(myTests, function(x)
    camera(x, index=1:20, design=testDesign, contrast=2L)$PValue)
summary(myCamera>=0.05)

myNegTests <- lapply(1:1000, function(x)
    simulateTwoGroupData(Nsamples=c(3,3), deltaMean=0, cor=c(0.5, 0.5), seed=x))
myNegCamera <- sapply(myNegTests, function(x)
    camera(x, index=1:20, design=testDesign, contrast=2L)$PValue)
summary(myNegCamera>=0.05)

ps <- seq(0,1,0.01)
fpr <- sapply(ps, function(x) mean(myNegCamera<=x))
tpr <- sapply(ps, function(x) mean(myCamera<=x))
plot(tpr~fpr, type="l")
abline(0,1,lwd=2)

