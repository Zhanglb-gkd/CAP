#' Ep calculates the expectation of tr(A\%*\%B) based on a permutation test
#'
#' @param A a centralized and symmetric matrix with nrow(A)>1
#' @param B a centralized and symmetric matrix with the same dimension as A
#' @return the expectation of tr(A\%*\%B) based on a permutation test
#' @export
Ep <- function(A, B){
  j = sum(diag(A))
  T = sum(diag(B))
  n = nrow(A)
  Expectation = j*T/(n-1)
  return(Expectation)
}


#' Vp calculates the variance of tr(A\%*\%B) based on a permutation test
#'
#' @param A a centralized and symmetric matrix with nrow(A)>3
#' @param B a centralized and symmetric matrix with the same dimension as A
#' @return the variance of tr(A\%*\%B) based on a permutation test
#' @export
Vp <- function(A, B){
  n = nrow(A)
  j = sum(diag(A))
  T = sum(diag(B))
  t2 = sum(diag(A %*% A))
  T2 = sum(diag(B %*% B))
  s2 = sum(diag(A)*diag(A))
  S2 = sum(diag(B)*diag(B))
  Variance =
    (s2*S2*(n-1)^2*(n-2)*(n-3)+((j^2-s2)*(T^2-S2)+2*(t2-s2)*(T2-S2)+4*s2*S2)*(n-1)*(n-2)*(n-3)+(4*(2*s2-t2)*(2*S2-T2)+
                                                                                                  2*(2*s2-j^2)*(2*S2-T^2))*(n-1)*(n-3)+(2*t2-6*s2+j^2)*(2*T2-6*S2+T^2)*(n-1)- (j^2)*(T^2)*(n-2)*(n-3)*n)/(n*((n-1)^2)*(n-2)*(n-3))
  #s2*S2/n + ((j^2-s2)*(T^2-S2)+2*(t2-s2)*(T2-S2)+4*s2*S2)/(n*(n-1)) + (4*(2*s2-t2)*(2*S2-T2) +
  #2*(2*s2-j^2)*(2*S2-T^2))/(n*(n-1)*(n-2)) + (2*t2-6*s2+j^2)*(2*T2-6*S2+T^2)/(n*(n-1)*(n-2)*(n-3)) - (j^2)*(T^2)/((n-1)^2)
  return(Variance)
}


#' Tp calculates the third moment of tr(A\%*\%B) based on a permutation test
#'
#' @param A a centralized and symmetric matrix with nrow(A)>5
#' @param B a centralized and symmetric matrix with the same dimension as A
#' @return the third moment of tr(A\%*\%B) based on a permutation test
#' @export
Tp <- function(A, B){
  n = nrow(A)
  j = sum(diag(A))
  T = sum(diag(B))
  t2 = sum(diag(A %*% A))
  T2 = sum(diag(B %*% B))
  s2 = sum(diag(A)*diag(A))
  S2 = sum(diag(B)*diag(B))
  t3 = sum(diag(A %*% A %*% A))
  T3 = sum(diag(B %*% B %*% B))
  s3 = sum(diag(A)*diag(A)*diag(A))
  S3 = sum(diag(B)*diag(B)*diag(B))
  u = sum(A*A*A)
  U = sum(B*B*B)
  d = sum(t(diag(A)) %*% diag(A %*% A))
  D = sum(t(diag(B)) %*% diag(B %*% B))
  b = sum(t(diag(A)) %*% A %*% diag(A))
  B = sum(t(diag(B)) %*% B %*% diag(B))
  ThirdMoment = ((n^2)*(n+1)*(n^2+15*n-4)*s3*S3+4*(n^4-8*n^3+19*n^2-4*n-16)*u*U+24*(n^2-n-4)*(u*B+b*U)
                 +6*(n^4-8*n^3+21*n^2-6*n-24)*b*B+12*(n^4-n^3-8*n^2+36*n-48)*d*D
                 +12*(n^3-2*n^2+9*n-12)*(j*s2*D+d*T*S2)+3*(n^4-4*n^3-2*n^2+9*n-12)*j*T*s2*S2
                 +24*((n^3-3*n^2-2*n+8)*(d*U+u*D)+(n^3-2*n^2-3*n+12)*(d*B+b*D))
                 +12*(n^2-n+4)*(j*s2*U+u*T*S2)+6*(2*n^3-7*n^2-3*n+12)*(j*s2*B+b*T*S2)
                 -2*n*(n-1)*(n^2-n+4)*((2*u+3*b)*S3+(2*U+3*B)*s3)
                 -3*n*((n-1)^2)*(n+4)*((j*s2+4*d)*S3+(T*S2+4*D)*s3)
                 +2*n*(n-1)*(n-2)*((j^3+6*j*t2+8*t3)*S3+(T^3+6*T*T2+8*T3)*s3)
                 +(j^3)*((n^3-9*(n^2)+23*n-14)*T^3+6*(n-4)*T*T2+8*T3)
                 +6*j*t2*((n-4)*T^3+(n^3-9*n^2+24*n-14)*T*T2+4*(n-3)*T3)
                 +8*t3*(T^3+3*(n-3)*T*T2+(n^3-9*n^2+26*n-22)*T3)
                 -16*((j^3)*U+u*(T^3))-12*(n^2-5*n+8)*(j*t2*U+u*T*T2)-8*(3*n^2-15*n+16)*(t3*U+u*T3)
                 -6*(n^2-5*n+4)*((j^3)*B+b*(T^3))-24*((n^2-5*n+8)*(t3*B+b*T3)+(n^2-5*n+6)*(j*t2*B+b*T*T2))
                 -3*(n-2)*(8*((j^3)*D+d*(T^3))+4*(n^2-5*n+12)*(j*t2*D+d*T*T2)
                           +8*(n^2-5*n+8)*(t3*D+d*T3)+(n^2-5*n+2)*((j^3)*T*S2+(T^3)*j*s2)
                           +2*(n^2-5*n+6)*(j*t2*T*S2+j*s2*T*T2)+16*(t3*T*S2+j*s2*T3)))/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
  return(ThirdMoment)
}


#' PVD: A fast computational procedure to detect the genetic pleiotropy effect for multiple variants
#'
#' @param G the similarity matrix for genetic data with nrow(G)>5
#' @param E the similarity matrix for phenotypic data with nrow(E)=nrow(G)
#' @return the p value for the fast method
#' @export
#' @importFrom PearsonDS ppearsonIII
#' @examples
#' library(MASS)
#' n = 100
#' m1 = 10
#' m2 = 10
#' betax = matrix(rnorm(m1*m2,0,1),m1,m2)
#' x = mvrnorm(n, rep(0,m1), diag(m1))
#' y = x\%*\%betax + mvrnorm(n, rep(0,m2), diag(m2))
#' f <- function(x){
#'   n = nrow(x)
#'   m = ncol(x)
#'   M = matrix(NA, nrow = n, ncol = n)
#'   for(i in 1:n){
#'     for(j in 1:n){
#'       M[i,j]=(sum(2-abs(x[i,]-x[j,])))/(2*m)
#'     }
#'   }
#'   return(M)
#' }
#' G = f(x)
#' E = f(y)
#' PVD(G, E)
PVD <- function(G,E){
  n = nrow(G)
  mr = nrow(E)
  nc = ncol(G)
  mc = ncol(E)
  if((n>5) & (n == mr) & (mr == nc) & (nc == mc)) {
  H = diag(n) - (1/n)*rep(1,n)%*%t(rep(1,n))

  A = H %*% G %*% H
  B = H %*% E %*% H

  eigA = eigen(A)
  valA0 = eigA$values
  kA = sum(Re(valA0)>0.01)
  valA = valA0[1:kA]
  vecA = eigA$vectors[,1:kA]

  eigB = eigen(B)
  valB0 = eigB$values
  kB = sum(Re(valB0)>0.01)
  valB = valB0[1:kB]
  vecB = eigB$vectors[,1:kB]

  X = c(1/16,1/4,1,4)
  P = rep(NA,4)

  for(i in 1:4){
    valA1 = valA^X[i]
    A1 = vecA %*% diag(valA1) %*% t(vecA)
    valB1 = valB^X[i]
    B1 = vecB %*% diag(valB1) %*% t(vecB)

    s = sum(A1*B1)

    x = Ep(A1,B1)
    y = Vp(A1,B1)
    z = Tp(A1,B1)

    a = (z-3*x*y-x^3)/(y^(3/2))
    b = 4/(a^2)

    tt = (s-x)/(sqrt(y))

    if(Re(a)>0){
      P[i] = ppearsonIII(q=Re(tt), shape=Re(b), location=Re(-sqrt(b)), scale=Re(1/sqrt(b)), lower.tail=F)
    }
    else{
      P[i] = ppearsonIII(q=Re(tt), shape=Re(b), location=Re(sqrt(b)), scale=Re(-1/sqrt(b)), lower.tail=F)
    }

  }

  t0 = sum((1/4)*tan((1/2-P)*pi))
  p0 = 0.5 - atan(t0)/pi

  return(p0)
  }
  else{
    print("error:it is necessary that nrow(G) = ncol(G) = nrow(E) = ncol(E) > 5")
  }
}
