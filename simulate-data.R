library(ribiosPlot)
library(ribiosIO)
figfile <- function(x) file.path("figures", x)

source("simulate-data-funcs.R")

testSimulator <- function(tgSim, pFunc=pFuncCamera) {
    bench <- newBenchmarker(tgSim, pFunc=pFunc)
    benchRes <- benchmark(bench)
    print(xyplot(benchRes))
    return(invisible(benchRes))
}
##----------------------------------------##
## test
##----------------------------------------##
tt <- newTwoGroupExprsSimulator(tpGeneSetCor=0.1)

cameraBenchmarkResult <- testSimulator(tt, pFuncCamera)
bioqcBenchmarkResult <- testSimulator(tt, pFuncBioQCtStat)
(rocCamera <- xyplot(cameraBenchmarkResult))
(rocBioQC <- xyplot(bioqcBenchmarkResult))

## tt2 <- cloneTwoGroupExprsSimulator(tt, randomSeed=1008)
## expect_equal(randomSeed(tt2), 1008)

## test
## set.seed(1887); twoGroupMvrnorm1887 <- twoGroupMvrnorm(20, c(3,3), 1, 0.1)
## stopifnot(all(tt@matrix[1:20,]-twoGroupMvrnorm1887==0))

## mutate backgroundground

## ttMut <- mutateBgByParams(tt, bgDgeInd=30:50,
##                          bgDgeDeltaMean=rnorm(21),
##                          bgCorInd=40:60,
##                          bgCorCluster=gl(7,3),
##                          bgCorSigma=0.95)

## ttMutRes <- testSimulator(ttMut)
## (rocMut <- xyplot(ttMutRes))

trellisLineCols <- trellis.par.get()$superpose.line$col

## if the background genes have random correlation structures, it's not a big problem for either camera or BioQC
ttMutRandCor <- randomlyMutateBg(tt, bgCorPerc=0.20)
ttMutRandCorCamera <- testSimulator(ttMutRandCor, pFunc=pFuncCamera)
ttMutRandCorBioQC <- testSimulator(ttMutRandCor, pFunc=pFuncBioQCtStat)
update(rocCamera, sub="CAMERA:Blue; BioQC: Red. Solid/Dash: iid/20% random cor", main="Random correlations in background genes") + update(rocBioQC, col="red") + xyplot(ttMutRandCorCamera,col=trellisLineCols[1], lty=2) + xyplot(ttMutRandCorBioQC, col="red", lty=2)
ipdf(figfile("ROC-randomCorrelation.pdf"), width=5L, height=5L)

## what happens if the background genes have increased expression?
ttMutPos <- randomlyMutateBg(tt, bgDgePerc=0.3,
                             bgDgeDeltaMeanFunc=function(n) rnorm(n, mean=1, sd=1))
ttMutPosCamera <- testSimulator(ttMutPos, pFunc=pFuncCamera)
ttMutPosBioQC <- testSimulator(ttMutPos, pFunc=pFuncBioQCtStat)
update(rocCamera, sub="CAMERA:Blue; BioQC: Red. Solid/Dash: iid/20% with logFC~N(1,1)", main="30% DE in background genes") + update(rocBioQC, col="red") + xyplot(ttMutPosCamera, lty=2, col=trellisPar$superpose.line$col[1]) + xyplot(ttMutPosBioQC, col="red", lty=2)
ipdf(figfile("ROC-randomDE.pdf"), width=5L, height=5L)

update(rocCamera, sub="CAMERA. Solid/Dash: iid/30% with logFC~N(1,1)", main="DE in background genes")+xyplot(ttMutPosCamera, lty=2, col=trellisPar$superpose.line$col[1])
ipdf(figfile("ROC-randomDE-cameraOnly.pdf"), width=5L, height=5L)

update(rocBioQC, sub="BioQC. Solid/Dash: iid/30% with logFC~N(1,1)", main="DE in background genes", col="red")+xyplot(ttMutPosBioQC, lty=2, col="red")
ipdf(figfile("ROC-randomDE-BioQConly.pdf"), width=5L, height=5L)

## a mild change?
ttMutMildPos <- randomlyMutateBg(tt, bgDgePerc=0.05,
                             bgDgeDeltaMeanFunc=function(n) rnorm(n, mean=1, sd=1))
ttMutMildPosCamera <- testSimulator(ttMutMildPos, pFunc=pFuncCamera)
xyplot(ttMutMildPosCamera)

## if random gene sets are correlated
testMutCor <- function(tgSim, pFunc=pFuncCamera) {
    bench <- newBenchmarker(tgSim, pFunc=pFunc)
    gsIndList <- bench@genesets[-1]
    gsInd <- unlist(gsIndList)
    gsFac <- factor(rep(seq(along=gsIndList), sapply(gsIndList,length)))
    tgSimMut <- mutateBgByParams(tgSim, bgCorInd=gsInd, bgCorCluster=gsFac, bgCorSigma=0.5)
    bench <- newBenchmarker(tgSimMut, pFunc=pFunc, geneSets=genesets(bench))
    benchRes <- benchmark(bench)
    print(xyplot(benchRes))
    return(invisible(benchRes))
}
(rocMut3 <- testMutCor(tt))
(rocMut3bioqc <- testMutCor(tt, pFunc=pFuncBioQCtStat))


xyplot(rocMut3, main="Strong correlation (0.5) in GS_R", sub="CAMERA:Blue; BioQC: Red. Solid/Dash: iid/50% cor in GS_R", lty=2)+xyplot(rocMut3bioqc, col="red", lty=2)+update(rocBioQC, col="red")+rocCamera
ipdf(figfile("ROC-GS_Rcor.pdf"), width=5L, height=5L)

xyplot(rocMut3, main="Strong correlation (0.5) in GS_R", sub="CAMERA. Solid/Dash: iid/50% cor in GS_R", lty=2)+rocCamera
ipdf(figfile("ROC-GS_Rcor-cameraOnly.pdf"), width=5L, height=5L)

update(rocBioQC, col="red", main="Strong correlation (0.5) in GS_R", sub="BioQC. Solid/Dash: iid/50% cor in GS_R")+xyplot(rocMut3bioqc, col="red", lty=2)
ipdf(figfile("ROC-GS_Rcor-BioQConly.pdf"), width=5L, height=5L)

AUC(rocMut3bioqc) ## AUC=.9
summary(ranks(rocMut3bioqc)<=5) ## in 90% case ranks highest

## as expected, if the expression of GS_R are also changed, the curve will look much worsex
testMutDeltaMean <- function(tgSim, pFunc=pFuncCamera) {
    bench <- newBenchmarker(tgSim, pFunc=pFunc)
    gsIndList <- bench@genesets[-1]
    gsInd <- unlist(gsIndList)
    gsFac <- factor(rep(seq(along=gsIndList), sapply(gsIndList,length)))
    tgSimMut <- mutateBgByParams(tgSim, bgDgeInd=gsInd, bgDgeDeltaMean=0.5)
    bench <- newBenchmarker(tgSimMut, pFunc=pFunc, geneSets=genesets(bench))
    benchRes <- benchmark(bench)
    print(xyplot(benchRes))
    return(invisible(benchRes))
}
(rocMut4 <- testMutDeltaMean(tt))
(rocMut4bioqc <- testMutDeltaMean(tt, pFunc=pFuncBioQCtStat))

xyplot(rocMut4, lty=2)+xyplot(rocMut4bioqc, col="red", lty=2)+update(rocBioQC, col="red")+rocCamera

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


## test pFuncLimmaAggregated ~ 180s
## system.time(pFuncLimmaAggregated(tt, list(1:20, 21:40)))

## romer ~ 336s
system.time(romerCorPerformance <- varParPerformance(varPar="tpGeneSetCor", varParList=0.1, pFunc=pFuncRomer))

## two-sample t-test of t-statistic
## TODO
## (1) Script
## (2) GSEA


## generate data
biosBenchmarkData <- function(nSamples,
                              deltaMean,
                              tpGeneSetCor,
                              tpGeneSetInd,
                              bgDgePerc) {
    obj <- newTwoGroupExprsSimulator(nSamples=nSamples,
                                     deltaMean=deltaMean,
                                     tpGeneSetCor=tpGeneSetCor,
                                     tpGeneSetInd=tpGeneSetInd)
    obj <- randomlyMutateBg(obj, bgDgePerc=bgDgePerc)
    return(obj)
                              
}
biosNs <- 2:6
biosDeltas <- c(0.5, 1, 1.5, 2, 3)
biosRos <- c(0, 0.1, 0.2, 0.5)
biosThetas <- c(5, 10, 20, 50, 100)
biosEpsilon <- c(0.01, 0.05, 0.25)

system.time(test <- biosBenchmarkData(2, 0.5, 0, 3, 0))

