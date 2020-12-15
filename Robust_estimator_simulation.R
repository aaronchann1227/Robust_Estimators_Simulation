##################################################################
#                        STAT 527
#                        Fall 2020
##################################################################
#                 Final Project, Robust Regression
#                        Mehrdad Mohammadi
#                        Man Chong Chan
#                        Yilin Zhu
##################################################################
#install.packages("MASS")
#install.packages("robust")
#install.packages("robustbase")
#install.packages("quantreg")
library(MASS)
library(robust)
library(quantreg)
library(robustbase)
##################################################################
# M estimators
# data
set.seed(114396701)
l1.data <- function(n1=100, n2){
  data <- mvrnorm(n=n1,mu=c(0, 0),
                  Sigma=matrix(c(1, .9, .9, 1), ncol=2))
  # generate n2 'bad' cases
  data <- rbind(data, mvrnorm(n=n2,
                              mu=c(1.5, -1.5), Sigma=.2*diag(c(1, 1))))
  data <- data.frame(data)
  names(data) <- c("X", "Y")
  ind <- c(rep(1, n1),rep(2, n2))
  plot(Y ~ X, data, pch=c("x", "o")[ind],
       col=c("black", "red")[ind],
       main=substitute(list(N[1] == n1, N[2] == n2), list(n1=n1, n2=n2)))
  # Regression Lines
r1  = rlm(Y~X, data=data)
r2 = rlm(Y~X, data=data, psi = psi.bisquare)
r3 = lm(Y ~ X, data)
r4 = lm(Y ~ X, data, subset=1:n1)
abline(r1, lwd=2)
abline(r2, lwd=2,col="darkgreen")
abline(r3, lwd=2, col="magenta")
abline(r4, lty=4, lwd=2, col="blue")
legend("topleft", c("Huber's","bisquare","OLS","OLS on good"),
       inset=0.02, lty=c(1, 1,1, 4), lwd=2, col=c("black", "darkgreen","magenta", "blue"),
       cex=.9)}
tiff("M.tiff", height = 10,width = 15,units = 'in',res=200)
par(mfrow=c(2, 2))
# different simulations for N
# function (N1,N2)
l1.data(100, 20)
l1.data(100, 30)
l1.data(100, 75)
l1.data(100, 100)
dev.off()


##################################################################
# s-estimator and MM-estimator

set.seed(114396701)
l1.data <- function(n1=100, n2){
  data <- mvrnorm(n=n1,mu=c(0, 0),
                  Sigma=matrix(c(1, .9, .9, 1), ncol=2))
  # generate n2 'bad' cases
  data <- rbind(data, mvrnorm(n=n2,
                              mu=c(1.5, -1.5), Sigma=.2*diag(c(1, 1))))
  data <- data.frame(data)
  names(data) <- c("X", "Y")
  ind <- c(rep(1, n1),rep(2, n2))
  plot(Y ~ X, data, pch=c("x", "o")[ind],
       col=c("black", "red")[ind],
       main=substitute(list(N[1] == n1, N[2] == n2), list(n1=n1, n2=n2)))
  # Regression Lines
  r1 = rlm(Y~X, data=data)
  abline(r1, lwd=2)
  abline(rlm(Y~X, data=data, psi = psi.bisquare), lwd=2,col="darkgreen")
  s1 = lmRob(Y ~ X,data,estim="initial")
  s2 = rlm(Y ~ X,data,method="MM")
  abline(s1, lwd=2, col="gold4")
  abline(s2, lwd=2,col="red")
  abline(lm(Y ~ X, data), lwd=2, col="magenta")
  abline(lm(Y ~ X, data, subset=1:n1), lty=2, lwd=2, col="blue")
  legend("topleft", c("Huber's","bisquare","S","MM","OLS","OLS on good"),
         inset=0.02, lty=c(1,1,1, 1,1, 4), lwd=2, col=c("darkgreen","gold4","black","red", "magenta", "blue"),
         cex=.9)}
# different simulations for N
# function (N1,N2)
tiff("S.tiff", height = 10,width = 15,units = 'in',res=200)
par(mfrow=c(2, 2))
l1.data(100, 20)
l1.data(100, 30)
l1.data(100, 75)
l1.data(100, 100)
dev.off()
# regressions

##################################################################
#LTS and MLS
set.seed(114396701)
l1.data <- function(n1=100, n2){
  data <- mvrnorm(n=n1,mu=c(0, 0),
                  Sigma=matrix(c(1, .9, .9, 1), ncol=2))
  # generate n2 'bad' cases
  data <- rbind(data, mvrnorm(n=n2,
                              mu=c(1.5, -1.5), Sigma=.2*diag(c(1, 1))))
  data <- data.frame(data)
  names(data) <- c("X", "Y")
  ind <- c(rep(1, n1),rep(2, n2))
  plot(Y ~ X, data, pch=c("x", "o")[ind],
       col=c("black", "red")[ind],
       main=substitute(list(N[1] == n1, N[2] == n2), list(n1=n1, n2=n2)))
  # Regression Lines
  lts1 = lqs(Y ~ X, data, method = "lts")
  lms1 = lqs(Y ~ X, data, method = "lms")
  abline(lts1, lwd=2,col="gold4")
  abline(lms1, lwd=2, col="red")
  abline(lm(Y ~ X, data), lwd=2, col="magenta")
  abline(lm(Y ~ X, data, subset=1:n1), lty=2, lwd=2, col="blue")
  legend("topleft", c("LTS","LMS","OLS","OLS on good"),
         inset=0.02, lty=c(1, 1,1, 4), lwd=2, col=c("gold4","red", "magenta", "blue"),
         cex=.9)}
tiff("LTS.tiff", height = 10,width = 15,units = 'in',res=200)
par(mfrow=c(2, 2))
l1.data(100, 20)
l1.data(100, 30)
l1.data(100, 75)
l1.data(100, 100)
dev.off()

##################################################################
#All in one graph 
set.seed(114396701)
l1.data <- function(n1=100, n2){
  data <- mvrnorm(n=n1,mu=c(0, 0),
                  Sigma=matrix(c(1, .9, .9, 1), ncol=2))
  # generate n2 'bad' cases
  data <- rbind(data, mvrnorm(n=n2,
                              mu=c(1.5, -1.5), Sigma=.2*diag(c(1, 1))))
  data <- data.frame(data)
  names(data) <- c("X", "Y")
  ind <- c(rep(1, n1),rep(2, n2))
  plot(Y ~ X, data, pch=c("x", "o")[ind],
       col=c("black", "red")[ind],
       main=substitute(list(N[1] == n1, N[2] == n2), list(n1=n1, n2=n2)))
  # Regression Lines

  abline(rlm(Y~X, data=data), lwd=2,col="darkorange")
  abline(rlm(Y~X, data=data, psi = psi.bisquare), lwd=2,col="darkgreen")
  abline(lmRob(Y ~ X,data,estim="initial"), lwd=2, col="brown4")
  abline(rlm(Y ~ X,data,method="MM"), lwd=2,col="gray0")
  abline(lqs(Y ~ X, data, method = "lts"), lwd=2,col="gold4")
  abline(lqs(Y ~ X, data, method = "lms"), lwd=2, col="red")
  abline(lm(Y ~ X, data), lwd=2, col="magenta")
  abline(lm(Y ~ X, data, subset=1:n1), lty=2, lwd=2, col="blue")
  abline(rq(Y ~ X, data=data, tau=0.5, ci=FALSE),lty=2, lwd=2, col="green4")
  legend("topleft", c("Huber's","bisquare","S","MM","LTS","LMS","OLS","OLS on good","L1"),
         inset=0.02, lty=c(1,1,1,1,1, 1,1, 4,4), lwd=2, col=c("darkorange","darkgreen","brown4","gray0","gold4","red", "magenta", "blue","green4"),
         cex=.9)}
tiff("All.tiff", height = 10,width = 15,units = 'in',res=200)
par(mfrow=c(1, 2))
l1.data(100, 20)
l1.data(100, 30)
dev.off()


