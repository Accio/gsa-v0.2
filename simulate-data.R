library(MASS)
library(testthat)
library(limma)
library(lattice)
library(ribiosExpression)
library(latticeExtra)

library(ribiosPlot)
figfile <- function(x) file.path("figures", x)
##------------------------------------------------------------##
## Education: correlation and covariance matrix
## from covariance matrix to correlation matrix (column-wise)
##------------------------------------------------------------##
myCorr <- function(mat) {
    matVar <- var(mat)
    sigma <- diag(diag(matVar)^(-.5), nrow=ncol(mat), ncol=ncol(mat))
    return(sigma %*% matVar %*% sigma)
}

myMatrix <- mvrnorm(1000, mu=rep(0,2), Sigma=matrix(c(10,2,3,10), nrow=2, byrow=TRUE))
expect_equal(cor(myMatrix), myCorr(myMatrix))

validateTwoGroupExprsSimulator <- function(object) {
    if(length(object@nGenes)<=1) {
        warning("nGenes length less than 1")
        return(FALSE)
    }
    if(length(object@nSamples)!=2) {
        warning("nSamples must be of two values")
        return(FALSE)
    }
    if(length(object@tpGeneSetInd)<1) {
        warning("tpGeneSetInd length less than 1")
        return(FALSE)
    }
    if(length(object@tpGeneSetCor)!=2) {
        warning("tpGeneSetCor length must equal 2")
        return(FALSE)
    }
    return(TRUE)
}
setClass("TwoGroupExprsSimulator",
         representation(nGenes="integer",
                        nSamples="integer",
                        tpGeneSetInd="integer",
                        deltaMean="numeric",
                        tpGeneSetCor="numeric",
                        randomSeed="integer",
                        matrix="matrix"),
         prototype=list(nGenes=12000L,
             nSamples=c(3L,3L),
             tpGeneSetInd=1:20,
             deltaMean=1,
             tpGeneSetCor=c(0,0),
             randomSeed=1887L,
             matrix=matrix()),
         validity=validateTwoGroupExprsSimulator)

setGeneric("nGenes", function(object) standardGeneric("nGenes"))
setGeneric("nSamples", function(object) standardGeneric("nSamples"))
setGeneric("tpGeneSetInd", function(object) standardGeneric("tpGeneSetInd"))
setGeneric("deltaMean", function(object) standardGeneric("deltaMean"))
setGeneric("tpGeneSetCor", function(object) standardGeneric("tpGeneSetCor"))
setGeneric("randomSeed", function(object) standardGeneric("randomSeed"))

setGeneric("nGenes<-", function(object,value) standardGeneric("nGenes<-"))
setGeneric("nSamples<-", function(object,value) standardGeneric("nSamples<-"))
setGeneric("tpGeneSetInd<-", function(object,value) standardGeneric("tpGeneSetInd<-"))
setGeneric("deltaMean<-", function(object,value) standardGeneric("deltaMean<-"))
setGeneric("tpGeneSetCor<-", function(object,value) standardGeneric("tpGeneSetCor<-"))
setGeneric("randomSeed<-", function(object,value) standardGeneric("randomSeed<-"))

setMethod("nGenes", "TwoGroupExprsSimulator", function(object) object@nGenes)
setMethod("nSamples", "TwoGroupExprsSimulator", function(object) object@nSamples)
setMethod("tpGeneSetInd", "TwoGroupExprsSimulator", function(object) object@tpGeneSetInd)
setMethod("deltaMean", "TwoGroupExprsSimulator", function(object) object@deltaMean)
setMethod("tpGeneSetCor", "TwoGroupExprsSimulator", function(object) object@tpGeneSetCor)
setMethod("randomSeed", "TwoGroupExprsSimulator", function(object) object@randomSeed)

setMethod("nGenes<-",
          c("TwoGroupExprsSimulator", "numeric"), function(object,value) {
              object@nGenes<-as.integer(value)
              return(object)
          })
setMethod("nSamples<-",
          c("TwoGroupExprsSimulator", "numeric"), function(object, value) {
              if(length(value)>2)
                  stop("value must be of length 2")
              if(length(value)==1) {
                  value <- c(value,value)
              }
              value <- as.integer(value)
              object@nSamples <- value
              return(object)
          })
setMethod("tpGeneSetInd<-",
          c("TwoGroupExprsSimulator", "numeric"), function(object,value) {
              object@tpGeneSetInd <- as.integer(value)
              return(object)
          })
setMethod("deltaMean<-",
          c("TwoGroupExprsSimulator", "numeric"), function(object,value) {
              object@deltaMean <- value
              return(object)
          })
setMethod("tpGeneSetCor<-",
          c("TwoGroupExprsSimulator", "numeric"), function(object,value) {
              if(length(value)>2)
                  stop("value must be of length 2")
              if(length(value)==1) {
                  value <- c(value,value)
              }
              object@tpGeneSetCor <- value
              return(object)
          })
setMethod("randomSeed<-",
          c("TwoGroupExprsSimulator", "numeric"), function(object,value) {
              object@randomSeed <- as.integer(value)
              return(object)
          })
setMethod("designMatrix", "TwoGroupExprsSimulator", function(object) {
              matrix(c(rep(1, sum(nSamples(object))),
                       rep(c(0,1), nSamples(object))),
                     ncol=2, byrow=FALSE)
          })
setMethod("contrastMatrix", "TwoGroupExprsSimulator", function(object) {
              matrix(c(0,1), byrow=FALSE, ncol=1)
          })
as.matrix.TwoGroupExprsSimulator <- function(x) return(x@matrix)
setAs(from="TwoGroupExprsSimulator", to="matrix", function(from) return(from@matrix))

simulateTwoGroupData <- function(tgSim) {
    set.seed(randomSeed(tgSim))
    Ngenes <- nGenes(tgSim)
    Nsamples <- nSamples(tgSim)
    mat <- matrix(rnorm(Ngenes*sum(Nsamples)), nrow=Ngenes)
    delta <- deltaMean(tgSim)

    tpGeneInd <- tpGeneSetInd(tgSim)
    tpGeneCount <- length(tpGeneInd)
    tpCor <- tpGeneSetCor(tgSim)
    isGroup2 <- 1:ncol(mat) > Nsamples[1]
    if(all(tpCor==0)) {
        mat[tpGeneInd,isGroup2] <-     mat[tpGeneInd,isGroup2]+delta
    } else {
        sigmaGroup1 <- matrix(tpCor[1], nrow=tpGeneCount, ncol=tpGeneCount)
        diag(sigmaGroup1) <- 1
        sigmaGroup2 <- matrix(tpCor[2], nrow=tpGeneCount, ncol=tpGeneCount)
        diag(sigmaGroup2) <- 1
        mvrnormGroup1 <- t(mvrnorm(n=Nsamples[1],
                                   mu=rep(0, tpGeneCount),
                                   Sigma=sigmaGroup1))
        mvrnormGroup2 <- t(mvrnorm(n=Nsamples[2],
                                   mu=rep(delta, tpGeneCount),
                                   Sigma=sigmaGroup2))
        mat[tpGeneInd,] <- cbind(mvrnormGroup1, mvrnormGroup2)
    }

    rownames(mat) <- paste("Gene", 1:Ngenes, sep="")
    colnames(mat) <- paste("Group", rep(c(1,2), Nsamples), ".Sample", c(1:Nsamples[1], 1:Nsamples[2]),sep="")
    return(mat)
}

simulateTwoGroupExprs <- function(nGenes=12000,
                                  nSamples=c(3,3),
                                  tpGeneSetInd=1:20,
                                  deltaMean=1.0,
                                  tpGeneSetCor=0,
                                  randomSeed=1887) {
    obj <- new("TwoGroupExprsSimulator")
    nGenes(obj) <- nGenes
    nSamples(obj) <- nSamples
    tpGeneSetInd(obj) <- tpGeneSetInd
    deltaMean(obj) <- deltaMean
    tpGeneSetCor(obj) <- tpGeneSetCor
    randomSeed(obj) <- randomSeed
    obj@matrix <- simulateTwoGroupData(obj)
    return(obj)
}



pFuncCamera <- function(tgSim, index) {
    res <- camera(as.matrix(tgSim),
                  index=index,
                  design=designMatrix(tgSim),
                  contrast=contrastMatrix(tgSim),
                  sort=FALSE)
    return(res$PValue)
}

setClass("BenchmarkDataset",
         representation("simulators"="list",
                        "genesets"="list"))
newBenchmarkDataset <- function(simulators, genesets) {
    new("BenchmarkDataset", simulators=simulators, genesets=genesets)
}

setGeneric("simulators", function(object) standardGeneric("simulators"))
setGeneric("genesets", function(object) standardGeneric("genesets"))
setMethod("simulators", "BenchmarkDataset", function(object) object@simulators)
setMethod("genesets", "BenchmarkDataset", function(object) object@genesets)

## pFunc: a function accepts an TwoGroupExprsSimulator object and a list of gene set indices
## expected return: p-values of the gene sets in the same order as a numeric vector

setClass("Benchmarker",
         representation("pFunc"="function"),
         contains="BenchmarkDataset")

setGeneric("pFunc", function(object) standardGeneric("pFunc"))
setGeneric("pFunc<-", function(object,value) standardGeneric("pFunc<-"))
setGeneric("pValues", function(object) standardGeneric("pValues"))
setMethod("pFunc", "Benchmarker", function(object) return(object@pFunc));
setMethod("pFunc<-", c("Benchmarker", "function"), function(object,value) {
              object@pFunc <- value
              return(object)
          })


generateBenchmarkData <- function(tgSim, ngsr=99, B=100) {
    geneSetInd <- tpGeneSetInd(tgSim)
    geneSets <- c(list(tpGeneSet=geneSetInd),
                  lapply(1:99, function(x) sample(setdiff(1:nGenes(tgSim), geneSetInd), length(geneSetInd))))
    names(geneSets) <- c("truePositive", paste("GeneSet", 2:length(geneSets)))
    mySims <- lapply(1:B, function(x)
        simulateTwoGroupExprs(tpGeneSetCor=tpGeneSetCor(tgSim), randomSeed=x))
    res <- newBenchmarkDataset(mySims, geneSets)
    return(res)
}

newBenchmarker <- function(tgSim, ngsr=99, B=100, pFunc) {
    bd <- generateBenchmarkData(tgSim, ngsr=ngsr, B=B)
    ber <- as(bd, "Benchmarker")
    pFunc(ber) <- pFunc
    return(ber)
}

setMethod("pValues", "Benchmarker", function(object) {
              lapply(simulators(object), function(sim) {
                         do.call(pFunc(object),
                                 list(sim, index=genesets(object)))
                     })
          })



validBenchmarkResult <- function(object) {
    if(ncol(object@ROC)!=3) {
        warning("'ROC' must be of three columns")
        return(FALSE)
    }
    if(!identical(colnames(object@ROC),c("pThr", "FPR", "TPR"))) {
        warning("'ROC' must contain 3 columns: pThr, FPR, TPR")
        return(FALSE)
    }
    return(TRUE)
}
setClass("BenchmarkResult",
         representation(ROC="data.frame",
                        AUC="numeric",
                        ranks="numeric"),
         validity=validBenchmarkResult)


setGeneric("benchmark", function(object) standardGeneric("benchmark"))
setGeneric("ROC", function(object) standardGeneric("ROC"))
setGeneric("AUC", function(object) standardGeneric("AUC"))
setGeneric("ranks", function(object) standardGeneric("ranks"))
setGeneric("avgRank", function(object) standardGeneric("avgRank"))
setGeneric("minRank", function(object) standardGeneric("minRank"))
setGeneric("maxRank", function(object) standardGeneric("maxRank"))
setGeneric("rankStat", function(object) standardGeneric("rankStat"))

setMethod("ROC", "BenchmarkResult", function(object) object@ROC)
setMethod("AUC", "BenchmarkResult", function(object) object@AUC)
setMethod("ranks", "BenchmarkResult", function(object) object@ranks)
setMethod("avgRank", "BenchmarkResult", function(object) mean(object@ranks, na.rm=TRUE))
setMethod("minRank", "BenchmarkResult", function(object) min(object@ranks))
setMethod("maxRank", "BenchmarkResult", function(object) max(object@ranks))
setMethod("rankStat", "BenchmarkResult", function(object) paste("min=",minRank(object),
                                                                "/avg=", avgRank(object),
                                                                "/max=", maxRank(object), sep=""))


newBenchmarkResult <- function(ROC, AUC, ranks) {
    return(new("BenchmarkResult",
               ROC=ROC, AUC=AUC, ranks=ranks))
}
setMethod("benchmark", "Benchmarker", function(object) {
              pval <- pValues(object)
              tpVals <- sapply(pval, "[[", 1L)
        
              pStep <- 0.001
              ps <- seq(0, 1, pStep)
              fpr <- sapply(ps, function(p0) mean(sapply(pval, function(pp) mean(pp[-1]<p0))))
              tpr <- sapply(ps, function(p0) mean(tpVals<=p0))
              roc <- data.frame(pThr=ps,
                                FPR=fpr,
                                TPR=tpr)
              auc <- sum(tpr)*pStep
              ranks <- sapply(pval, function(x) rank(x)[1])
              res <- new("BenchmarkResult", ROC=roc, AUC=auc, ranks=ranks)
              return(res)
          })

drawAt <- function(at, pThr, x, y, ...) {
    isAt <- which.min(abs(pThr-at))
    panel.points(x[isAt], y[isAt], ...)
}
panel.benchmarkResult <- function(x,y,pThr,...) {
    panel.xyplot(x,y,...)
    panel.abline(v=c(0.01, 0.05), lty=2, col="darkgray")
    drawAt(0.05, pThr, x, y, col="black", pch=24, fill="orange", cex=1.3)
    drawAt(0.01, pThr, x, y, col="black", pch=24, fill="red", cex=1.3)

}
xyplot.BenchmarkResult <- function(x, data,...) {
    roc <- ROC(x)
    xyplot(TPR~FPR, type="l", data=roc, abline=c(0,1),
           xlab="False positive rate", ylab="True positiev rate",
           panel=panel.benchmarkResult, pThr=roc$pThr,
           key=list(space="top", columns=2,
               text=list(labels=c("nominal p=0.05", "nominal p=0.01")),
               points=list(fill=c("orange", "red"), pch=24, cex=1.5, col="black")),
           scales=list(tck=c(1,0), alternating=1L),
           ...)
}
tt <- simulateTwoGroupExprs()
cameraBenchmarker <- newBenchmarker(tt, pFunc=pFuncCamera)
cameraBenchmarkResult <- benchmark(cameraBenchmarker)
(roc1 <- xyplot(cameraBenchmarkResult))

varParData <- function(varPar=c("nGenes", "nSamples", "tpGeneSetInd", "deltaMean", "tpGeneSetCor"),
                       varParList, ...) {
    varPar <- match.arg(varPar)
    benchSimParas <- lapply(varParList, function(x) {
                                 basic <- list()
                                 basic[[varPar]] <- x
                                 return(c(basic, list(...)))
                             })
    benchData <- lapply(benchSimParas, function(x) {
                            do.call(simulateTwoGroupExprs, x)
                    })
    return(benchData)
}


varParPerformance <- function(varPar=c("nGenes", "nSamples", "tpGeneSetInd", "deltaMean", "tpGeneSetCor"),
                              varParList, pFunc, ...) {
    varData <- varParData(varPar=varPar,
                          varParList=varParList, ...)
    performance <- lapply(varData, function(x) {
                              benchmarker <- newBenchmarker(x, pFunc=pFunc)
                              benchmarkResult <- benchmark(benchmarker)
                              return(benchmarkResult)
                          })
    return(performance)
}

colRamp <- function(cols, n) colorRampPalette(cols)(n)

plotROC <- function(performanceList,
                    cols=c("black", "red"),
                    main, key.title, values) {
    rocPlots <- lapply(performanceList, xyplot)
    rocCols <- colRamp(cols, length(performanceList))
    rocComb <- update(rocPlots[[1]], col=rocCols[1])
    if(length(rocPlots)>1) {
        for(i in 2:length(rocPlots)) {
            rocComb <- rocComb+update(rocPlots[[i]], col=rocCols[[i]])
        }
    }
    plot(update(rocComb, main=main,
            key=list(space="right", cex=0.85,title=key.title,
                text=list(labels=sprintf("%1.1f", values))
              , lines=list(col=rocCols))))
}

panel.AUC <- function(x,y,col,...) {
    panel.xyplot(x,y,pch=21, fill=col,cex=1.25, col="black",...)
    for(i in 1:(length(x)-1)) {
        panel.segments(x[i], y[i],
                      x[i+1], y[i+1], col=col[i])
    }
}
plotAUC <- function(performanceList,
                    cols=c("black", "red"),
                    main, key.title, values,
                    scales=list(tck=c(1,0), alternating=1L, x=list(at=values)),
                    ...) {
    rocCols <- colRamp(cols, length(performanceList))
    auc <- sapply(performanceList, AUC)
    xyplot(auc ~ values, col=rocCols, panel=panel.AUC,
           scales=scales,
           xlab=key.title, ylab="AUC", main=main,...)
}

panel.ranks <- function(x,y,col,...) {
    panel.dotplot(x,y,col=col,...)
    yMean <- tapply(y, x, mean)
    xUniq <- 1:ribiosUtils::ulen(x)
    ucol <- unique(col)
    barWidth <- 0.2
    for(i in 1:(length(yMean)-1)) {
        panel.segments(xUniq[i], yMean[i],
                       xUniq[i+1], yMean[i+1], col=ucol[i])
    }
    for(i in 1:length(yMean)) {
        panel.segments(xUniq[i]-barWidth, yMean[i],
                      xUniq[i]+barWidth, yMean[i], lwd=2,
                       col=ucol[i])
    }
}
plotRanks <- function(performanceList,
                    cols=c("black", "red"),
                    main, key.title, values) {
    rocCols <- colRamp(cols, length(performanceList))
    ranks <- lapply(performanceList, ranks)
    df <- data.frame(Rank=unlist(ranks),
                     Value=rep(values, sapply(ranks, length)))
    pointCols <- rep(rocCols, sapply(ranks, length))
    dotplot(Rank ~ Value, col=pointCols, data=df, horizontal=FALSE,
            panel=panel.ranks,
            scales=list(tck=c(1,0), alternating=1L, x=list(at=seq(along=values), labels=sapply(values, "[[", 1L))),
            xlab=key.title, ylab="Ranks", main=main)
}

##----------------------------------------##
## benchmarks
##----------------------------------------##

## test
tt1 <- simulateTwoGroupExprs(nSamples=c(3,3), tpGeneSetCor=0.1)
tt2 <- simulateTwoGroupExprs(nSamples=c(15,15), tpGeneSetCor=0.1)

tt1Benchmark <- newBenchmarker(tt1, pFunc=pFuncCamera)
tt2Benchmark <- newBenchmarker(tt2, pFunc=pFuncCamera)

tt1BenchmarkRes <- benchmark(tt1Benchmark)
tt2BenchmarkRes <- benchmark(tt2Benchmark)

plotROC(list(list(tt1=tt1BenchmarkRes, tt2=tt2BenchmarkRes)), key.title="Sample size", values=c(3,15))

## corrlation
testCors <- c(seq(0, 0.5, 0.1), seq(0.6, 1, 0.2))
cameraCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=testCors, pFunc=pFuncCamera)
plotROC(cameraCorPerformance, main="Camera's performance", key.title="Correlation", values=testCors)
ipdf(figfile("camera-ROC-cor.pdf"), width=6.5, height=6L)
plotAUC(cameraCorPerformance, main="Camera's performance", key.title="Correlation", values=testCors)
ipdf(figfile("camera-AUC-cor.pdf"), width=6, height=6L)
plotRanks(cameraCorPerformance, main="Camera's performance", key.title="Correlation", values=testCors)
ipdf(figfile("camera-Ranks-cor.pdf"), width=6, height=6L)


## test sample sizes: under the simulation scenario the sample number does not make a difference
testSampleSizes <- lapply(c(seq(2,5), seq(6,12, 2), seq(15,30,5)), function(x) rep(x,2))
cameraSamplePerformance <- varParPerformance(varPar="nSample", varParList=testSampleSizes, pFunc=pFuncCamera, tpGeneSetCor=0.1)
plotROC(cameraSamplePerformance, main="Camera's performance", key.title="nSample", values=sapply(testSampleSizes, "[[", 1L))
ipdf(figfile("camera-ROC-sampleSize.pdf"), width=6.5, height=6L)

## test gene size
testGeneSetSize <- c(2,5,10,15,20,30,50,100,200,500,1000)
cameraGSsizePerformance <- varParPerformance(varPar="tpGeneSetInd",
                                             varParList=lapply(testGeneSetSize, function(x) 1:x),
                                             pFunc=pFuncCamera, tpGeneSetCor=0)
gsSizeCol <- c(royalbluered(3)[1],"lightgray", royalbluered(3)[3])
plotROC(cameraGSsizePerformance, cols=gsSizeCol,
        main="Camera's performance", key.title="Gene set size", values=testGeneSetSize)
ipdf(figfile("camera-ROC-geneSetSize.pdf"), width=6.5, height=6L)
plotAUC(cameraGSsizePerformance, main="Camera's performance", key.title="Gene set size", values=testGeneSetSize,
        cols=gsSizeCol,
        scales=list(tck=c(1,0), alternating=1L, x=list(at=testGeneSetSize, log=2)))
ipdf(figfile("camera-AUC-geneSetSize.pdf"), width=6, height=6L)
plotRanks(cameraGSsizePerformance, main="Camera's performance", key.title="Gene set size", values=testGeneSetSize,
          cols=gsSizeCol)
ipdf(figfile("camera-Ranks-geneSetSize.pdf"), width=5, height=5)

## test effect size
testDeltaMean <- c(seq(-3, -0.5, 0.5), 0, seq(0.5, 3, 0.5))
cameraESPerformance <- varParPerformance(varPar="deltaMean",
                                         varParList=testDeltaMean, 
                                         pFunc=pFuncCamera, tpGeneSetCor=0.1)
plotROC(cameraESPerformance, 
        main="Camera's performance", key.title="Effect size", values=testDeltaMean)
plotRanks(cameraESPerformance, main="Camera's performance", key.title="Effect size", values=testDeltaMean)
