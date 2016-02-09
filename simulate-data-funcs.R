library(MASS)
library(testthat)
library(limma)
library(lattice)
library(ribiosExpression)
library(ribiosUtils)
library(latticeExtra)
library(BioQC)
library(globaltest)


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
