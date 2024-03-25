\name{PVD}
\alias{PVD}
\title{PVD: A fast computational procedure to detect the genetic pleiotropy effect for multiple variants}
\usage{
PVD(G, E)
}
\arguments{
\item{G}{the similarity matrix for genetic data with nrow(G)>5}

\item{E}{the similarity matrix for phenotypic data with nrow(E)=nrow(G)}
}
\value{
the p value for the fast method
}
\description{
PVD: A fast computational procedure to detect the genetic pleiotropy effect for multiple variants
}
\examples{
library(MASS)
n = 100
m1 = 10
m2 = 10
betax = matrix(rnorm(m1*m2,0,1),m1,m2)
x = mvrnorm(n, rep(0,m1), diag(m1))
y = x\%*\%betax + mvrnorm(n, rep(0,m2), diag(m2))
f <- function(x){
  n = nrow(x)
  m = ncol(x)
  M = matrix(NA, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      M[i,j]=(sum(2-abs(x[i,]-x[j,])))/(2*m)
    }
  }
  return(M)
}
G = f(x)
E = f(y)
PVD(G, E)
}
