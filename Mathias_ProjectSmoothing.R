#####################################
#########       SCRIPT   ############
#####################################

# require packages
library(ggplot2)
theme_set(theme_bw())
library(lokern)


# Function and asymptotic confidence interval construction

set.seed(28112022)

# sample size
n <- 50 

# m(x) the known mean the regression
m <- function(x) (sin(2*pi*x^3))^3

# X from uniform distribution
X <- seq(from =1, to = n)/n

# Y obtain from formula
Y <- m(X) + 0.5*rnorm(n, 0,1)

# chosen points to estimate m(x) 
# x <- sort(sample(x = X, size = n/2))

x <- seq(from = min(X), to = max(X), length = n)

# global bandwidth to use
h <- bw.ucv(X)

# kernel to use : epanechnikov
epanKer <- function(u) 0.75*(1-u^2)*(abs(u) <=1)

# Nadaraya-Watson estimate to use
nwregEst <- function(x, X, Y, h, epanKer) sum(Y*epanKer((x-X)/h))/sum(epanKer((x-X)/h))

# NW on x
nwEst <- sapply(X, function(x) nwregEst(x, X, Y, h, epanKer))


# residu for whole data
residu <- Y - sapply(X, function(x) nwregEst(x, X, Y, h, epanKer))

# integral of square of kernel
R_K <- integrate(function(u) epanKer(u)^2, -Inf, Inf)$value

# sigma_hat for chosen x
h2 <- bw.SJ(X)
sigma_hat <- function(x) sum(epanKer((X-x)/h2)*residu^2)/sum(epanKer((X-x)/h2))


# estimation of f(x)
fest <- function(x) (1/50)*sum((1/h2)*epanKer((X-x)/h2))

# lowerbound of CI at 95%
lowb <- function(xi) nwregEst(xi, X, Y, h, epanKer) - 1.96*(sqrt(1/(n*h))*sqrt(R_K * sigma_hat(xi) /fest(xi)))

lowerbound <- rep(NA, 50)

for (j in 1:50){
  lowerbound[j] <- lowb(x[j])
}


# upperband of CI at 95%

upb <- function(xi) nwregEst(xi, X, Y, h, epanKer) + 1.96*(sqrt(1/(n*h))*sqrt(R_K * sigma_hat(xi) /fest(xi)))


upperBound <- rep(NA, 50)
for (j in 1:50){
  upperBound[j] <- upb(x[j])
}




# Plot of X, Y, m(x) and CI
plot(X, Y, pch = "+")
lines(X, m(X), type = "l", col = 1)
lines(X, nwEst, type = "l", col = 2)
lines(X, lowerbound, col = 3)
lines(x, upperBound, col = 4)
rug(X)
title(main = "Plot of X, Y, m(x) estimate and CI")
legend("topright", legend = c("ci_upb", "ci_lowb", "NW_Est", "m(x)"), col = c(4,3,2,1), lty = c(1,1,1,1), cex = 0.75)


# chosen reference points
x1 <- 0.18; x2 <- 0.39; x3 <- 0.52; x4 <- 0.68; x5 <- 0.81; x6 <- 0.92
m1 <- m(0.18) ; m2 <- m(0.39); m3 <- m(0.52); m4 <- m(0.68) ; m5 <- m(0.81); m6 <- m(0.92)



# Confidence interval construction from monte carlo
# asymptotic approach

ci_construct <- function(x_cho, m_cho, n, level){
  
  X50 <- seq(from =1, to = n)/n
  h_to_use <- bw.ucv(X50)
  Y50 <- m(X50) + 0.5*rnorm(n, 0,1)
  
  
  # integral of square of kernel
  R_K <- integrate(function(u) epanKer(u)^2, -Inf, Inf)$value
  
  # sigma_hat for chosen x
  h2 <- bw.SJ(X50)
  sigma_hat <- function(x) sum(epanKer((X50-x)/h2)*residu^2)/sum(epanKer((X50-x)/h2))
  
  
  # estimation of f(x)
  fest <- function(x) (1/n)*sum((1/h2)*epanKer((X50-x)/h2))
  
  coverage <- rep(NA, 1000)
  # mesti <- rep(NA, 1000)
  
  for (j in 1:1000){
    
    # resampling from X50
    Xsa <- sample(X50, replace = TRUE)
    
    
    Ysa <- m(Xsa) + 0.5*rnorm(n, 0,1)
    
    mhat <- nwregEst(x_cho, Xsa, Ysa, h_to_use, epanKer)
    
    # mesti[j] <- mhat
    
    # residu for whole data
    residu <- Ysa - sapply(Xsa, function(x) nwregEst(x, Xsa, Ysa, h_to_use, epanKer))
    
    lb <- mhat - level*(sqrt(1/(n*h_to_use))*sqrt(R_K * sigma_hat(x_cho) /fest(x_cho)))
    
    ub <- mhat + level*(sqrt(1/(n*h_to_use))*sqrt(R_K * sigma_hat(x_cho) /fest(x_cho)))
    
    coverage[j] <- (lb < m_cho & m_cho < ub)*1
  }
  
  return(round(mean(coverage),2))
  
}


# quantile based approach

ci2_construct <- function(x_cho, m_cho, n, level){
  
  X50 <- seq(from =1, to = n)/n
  h_to_use <- bw.ucv(X50)
  Y50 <- m(X50) + 0.5*rnorm(n, 0,1)
  
  # mhat0 <- nwregEst(x_cho, X50, Y50, h_to_use, epanKer)
  
  
  # 
  # # integral of square of kernel
  # R_K <- integrate(function(u) epanKer(u)^2, -Inf, Inf)$value
  # 
  # # sigma_hat for chosen x
  # h2 <- bw.SJ(X50)
  # sigma_hat <- function(x) sum(epanKer((X50-x)/h2)*residu^2)/sum(epanKer((X50-x)/h2))
  # 
  # 
  # # estimation of f(x)
  # fest <- function(x) (1/n)*sum((1/h2)*epanKer((X50-x)/h2))
  
  coverage <- rep(NA, 500)
  
  
  for (i in 1:500){
    mesti <- rep(NA, 100)
    for (j in 1:100){
      
      # resampling from X50
      Xsa <- sample(X50, replace = TRUE)
      
      Ysa <- m(Xsa) + 0.5*rnorm(n, 0,1)
      
      mesti[j] <- nwregEst(x_cho, Xsa, Ysa, h_to_use, epanKer)
      
    }
    
    # residu for whole data
    # residu <- Ysa - sapply(Xsa, function(x) nwregEst(x, Xsa, Ysa, h_to_use, epanKer))
    
    lb <- quantile(mesti, level)
    
    ub <- quantile(mesti, 1-level)
    
    coverage[i] <- (lb <= m_cho & m_cho <= ub)*1
    
  }
  
  return(round(mean(coverage),2))
  
}