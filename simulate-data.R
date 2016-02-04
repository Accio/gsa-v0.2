library(ribiosPlot)
library(ribiosIO)
figfile <- function(x) file.path("figures", x)

source("simulate-data-funcs.R")

##----------------------------------------##
## test
##----------------------------------------##
tt <- newTwoGroupExprsSimulator(tpGeneSetCor=0.1)

cameraBenchmarker <- newBenchmarker(tt, pFunc=pFuncCamera)
cameraBenchmarkResult <- benchmark(cameraBenchmarker)
(roc1 <- xyplot(cameraBenchmarkResult))

tt2 <- cloneTwoGroupExprsSimulator(tt, randomSeed=1008)
expect_equal(randomSeed(tt2), 1008)

## test
set.seed(1887); twoGroupMvrnorm1887 <- twoGroupMvrnorm(20, c(3,3), 1, 0.1)
stopifnot(all(tt@matrix[1:20,]-twoGroupMvrnorm1887==0))

## mutate byground
ttMut <- mutateBg(tt, bgDgeInd=30:50, bgDgeDeltaMean=rnorm(20),
                  bgCorInd=40:60,
                  bgCorCluster=gl(3,7),
                  bgCorSigma=.95)

ttMut2 <- randomlyMutateBg(tt, bgDgePerc=0.2, bgCorPerc=0.2)
ttMut3 <- randomlyMutateBg(tt, bgDgePerc=0.2, bgCorPerc=0.2, bgDgeCorPerc=0.02)

##----------------------------------------##
## benchmarks
##----------------------------------------##
outfile <- function(x) file.path("output", x)
cameraBenchmarkFile <- outfile("cameraBenchmark.RData")
if(!loadFile(cameraBenchmarkFile)) {
    ## correlation
    testCors <- c(seq(0, 0.5, 0.1), seq(0.6, 1, 0.2))
    cameraCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=testCors, pFunc=pFuncCamera)

    testSampleSizes <- lapply(c(seq(2,5), seq(6,12, 2), seq(15,30,5)), function(x) rep(x,2))
    testSampleSizeHalf <- sapply(testSampleSizes, "[[", 1L)
    cameraSamplePerformance <- varParPerformance(varPar="nSample",
                                                 varParList=testSampleSizes, pFunc=pFuncCamera, tpGeneSetCor=0.1)

    testGeneSetSize <- c(2,5,10,15,20,30,50,100,200,300,400, 500,1000)
    cameraGSsizePerformance <- varParPerformance(varPar="tpGeneSetInd",
                                             varParList=lapply(testGeneSetSize, function(x) 1:x),
                                             pFunc=pFuncCamera, tpGeneSetCor=0)

    testDeltaMean <- c(seq(-3, -0.5, 0.5), 0, seq(0.5, 3, 0.5))
    cameraESPerformance <- varParPerformance(varPar="deltaMean",
                                             varParList=testDeltaMean, 
                                             pFunc=pFuncCamera, tpGeneSetCor=0.1)

    testOppDeltaMeanVal <- seq(0.1, 0.5, 0.1)

    testOppDeltaMean <- lapply(testOppDeltaMeanVal,
                               function(x) {
                                   negs <- rep(-1.5, as.integer(20*x))
                                   pos <- rep(1.5, 20-length(negs))
                                   return(c(negs, pos))
                               })
    cameraOppESPerformance <- varParPerformance(varPar="deltaMean",
                                                varParList=testOppDeltaMean, 
                                                pFunc=pFuncCamera, tpGeneSetCor=0.1)
    
    testPartDeltaMeanVal <- seq(0.1, 0.9, 0.1)
    testPartDeltaMean <- lapply(testPartDeltaMeanVal,
                                function(x) {
                                    noneff <- rep(0, as.integer(20*x))
                                    pos <- rep(1.5, 20-length(noneff))
                                    return(c(noneff, pos))
                               })
    cameraPartESPerformance <- varParPerformance(varPar="deltaMean",
                                                varParList=testPartDeltaMean, 
                                                pFunc=pFuncCamera, tpGeneSetCor=0.1)

    save(testCors, cameraCorPerformance,
         testSampleSizes,testSampleSizeHalf, cameraSamplePerformance,
         testGeneSetSize,cameraGSsizePerformance,
         testDeltaMean,cameraESPerformance,
         testOppDeltaMeanVal, testOppDeltaMean,cameraOppESPerformance,
         testPartDeltaMeanVal, testPartDeltaMean,cameraPartESPerformance,
         file=cameraBenchmarkFile)
}

plotROC(cameraCorPerformance, main="Camera's performance", key.title="Correlation", values=testCors)
ipdf(figfile("camera-ROC-cor.pdf"), width=6.5, height=6L)
plotAUC(cameraCorPerformance, main="Camera's performance", key.title="Correlation", values=testCors)
ipdf(figfile("camera-AUC-cor.pdf"), width=6, height=6L)
plotRanks(cameraCorPerformance, main="Camera's performance", key.title="Correlation", values=testCors)
ipdf(figfile("camera-Ranks-cor.pdf"), width=6, height=6L)


## test sample sizes: under the simulation scenario the sample number does not make a difference
plotROC(cameraSamplePerformance, main="Camera's performance", key.title="nSample", values=testSampleSizeHalf)
ipdf(figfile("camera-ROC-sampleSize.pdf"), width=6.5, height=6L)
plotAUC(cameraSamplePerformance, main="Camera's performance", key.title="nSample", values=testSampleSizeHalf)
ipdf(figfile("camera-AUC-sampleSize.pdf"), width=6, height=6L)
plotRanks(cameraSamplePerformance, main="Camera's performance", key.title="nSample", values=testSampleSizeHalf)
ipdf(figfile("camera-Ranks-sampleSize.pdf"), width=6, height=6L)

## test gene set size
plotROC(cameraGSsizePerformance, 
        main="Camera's performance", key.title="Gene set size", values=testGeneSetSize)
ipdf(figfile("camera-ROC-geneSetSize.pdf"), width=6.5, height=6L)
plotAUC(cameraGSsizePerformance, main="Camera's performance", key.title="Gene set size", values=testGeneSetSize,
        scales=list(tck=c(1,0), alternating=1L, x=list(at=testGeneSetSize, log=2)))
ipdf(figfile("camera-AUC-geneSetSize.pdf"), width=6, height=6L)
plotRanks(cameraGSsizePerformance, main="Camera's performance", key.title="Gene set size", values=testGeneSetSize)
ipdf(figfile("camera-Ranks-geneSetSize.pdf"), width=5, height=5)


## test effect size
plotROC(cameraESPerformance, 
        main="Camera's performance", key.title="Effect size", values=testDeltaMean)
esCol <- c(royalbluered(3)[1],"lightgray", royalbluered(3)[3])
plotROC(cameraESPerformance, main="Camera's performance", key.title="Effect size", values=testDeltaMean, cols=esCol)
ipdf(figfile("camera-ROC-effectSize.pdf"), width=6.5, height=6L)
plotAUC(cameraESPerformance, main="Camera's performance", key.title="Effect size", values=testDeltaMean, cols=esCol,
        ylim=c(0.45,1.05))
ipdf(figfile("camera-AUC-effectSize.pdf"), width=6, height=6)
plotRanks(cameraESPerformance, main="Camera's performance", key.title="Effect size", values=testDeltaMean, cols=esCol)
ipdf(figfile("camera-Ranks-effectSize.pdf"), width=5, height=5)

## test opposite effect size

plotROC(cameraOppESPerformance, 
        main="Camera's performance", key.title="Opp. effect", values=testOppDeltaMeanVal)
ipdf(figfile("camera-ROC-OppositeEffectSize.pdf"), width=6.5, height=6L)
plotAUC(cameraOppESPerformance, main="Camera's performance", key.title="Opp. effect", values=testOppDeltaMeanVal, 
        ylim=c(0.45,1.05))
ipdf(figfile("camera-AUC-OppositeEffectSize.pdf"), width=6, height=6)
plotRanks(cameraOppESPerformance, main="Camera's performance", key.title="Opp. effect", values=testOppDeltaMeanVal)
ipdf(figfile("camera-Ranks-OppositeEffectSize.pdf"), width=5, height=5)

## test partsite effect size

plotROC(cameraPartESPerformance, 
        main="Camera's performance", key.title="Part. effect", values=testPartDeltaMeanVal)
ipdf(figfile("camera-ROC-PartEffectSize.pdf"), width=6.5, height=6L)
plotAUC(cameraPartESPerformance, main="Camera's performance", key.title="Part. effect", values=testPartDeltaMeanVal, 
        ylim=c(0.45,1.05))
ipdf(figfile("camera-AUC-PartEffectSize.pdf"), width=6, height=6)
plotRanks(cameraPartESPerformance, main="Camera's performance", key.title="Part. effect", values=testPartDeltaMeanVal)
ipdf(figfile("camera-Ranks-PartEffectSize.pdf"), width=5, height=5)


##----------------------------------------##
## Other methods
##----------------------------------------##
bioqcTstat <- newBenchmarker(tt, pFunc=pFuncBioQCtStat)
bioqcTstatResult <- benchmark(bioqcTstat)
(rocBioQCTstat <- xyplot(bioqcTstatResult))

bioqcTstatCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=testCors, pFunc=pFuncBioQCtStat)

plotROC(bioqcTstatCorPerformance, key.title="Cor", values=testCors)
ipdf(figfile("BioQC.log10P.Tstat-ROC-cor.pdf"), width=6.5, height=6L)
xyplot(sapply(bioqcTstatCorPerformance, AUC)~sapply(cameraCorPerformance, AUC), type="b", abline=c(0,1))
ipdf(figfile("ROC-correlation-BioQC.log10P.Tstat-Camera.pdf"), width=6.5, height=6L)



## camera ~ 14s
(cameraSpeed <- system.time(testCameraPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncCamera)))

## BioQC ~ 14s
(bioqcSpeed <- system.time(testBioQCTstatPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncBioQCtStat)))

## cameraRank ~ 75s
(cameraRankSpeed <- system.time(cameraRankCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncCameraRank)))

## fisher's method ~ 99s
system.time(fisherMethodCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncFisherMethod))

## mroast ~ 126s 
(mroastSpeed <- system.time(mroastCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncMroast)))

## fisher: 131s
(fisherExactSpeed <- system.time(testFisherPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncFisherExact)))

## globaltest 290s
(gtSpeed <- system.time(gtCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncGlobaltest)))


## one-sample z-test of t-statistic: 122s
system.time(tStatZtestCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncTstatZtest))
## two-sample t-test of t-statistic: 144s
(tStatTtestSpeed <- system.time(tStatTtestCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncTstatTtest)))

## WMW of t-statistic: 112s
(tStatWmwSpeed <- system.time(tStatWmwCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncTstatWMW)))

## Chisq^test of t-statistic: 110s
(chisqSpeed <- system.time(chisqCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncChisq)))


## test pFuncLimmaAggregated
## system.time(pFuncLimmaAggregated(tt, list(1:20, 21:40)))

## romer ~ 336s
system.time(romerCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncRomer))

## two-sample t-test of t-statistic
## TODO
## (1) Script
## (2) GSEA
