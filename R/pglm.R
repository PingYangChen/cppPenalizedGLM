n <- 200
d <- 10
trueb0 <- .1
trueact <- cbind(
  c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0),
  c(0, 0, 1, 1, 1, 1, 1, 0, 0, 0),
  c(0, 0, 0, 0, 1, 1, 1, 1, 1, 0)
)
trueb <- matrix(runif(d*3, -1, 1)*10, d, 3)
for (m in 1:3) { trueb[which(trueact[,m] == 0),m] <- 0 }

x <- matrix(rnorm(n*d), n, d)

family <- "multinomial" #  "gaussian" #  #"binomial" #
if (family == "gaussian") {
  y <- trueb0[1] + x %*% trueb[,1] + rnorm(n, 0, 0.5)
} else if (family == "binomial") {
  eta <- trueb0[1] + x %*% trueb[,1]; yp <- exp(eta)/(1 + exp(eta))
  y <- rbinom(n, 1, yp)
} else if (family == "multinomial") {
  eta <- trueb0 + x %*% trueb
  yp <- exp(eta)/matrix(rowSums(exp(eta)), n, 3)
  y <- matrix(0, n, 3)
  for (i in 1:n) {
    y[i,] <- rmultinom(1, 1, yp[i,])
  }
}


lambda = 0.05
alpha = 0.5
maxit = 100
lambdaLength = 100
minLambdaRatio = 1e-3
tol = 1e-4

penalizedGLM <- function(
  y, x, lambdaLength = 100, minLambdaRatio = 1e-3,
  alpha = 0.5, maxit = 100)
{

  xm <- colMeans(x)
  xs <- sapply(1:d, function(j) sd(x[,j]))
  xs[which(xs == 0)] <- 1
  x <- (x - matrix(xm, n, d, byrow = TRUE))/xs

  lambdaMax <- calcMaxLambda(family, x, y, alpha)
  lambdaVec <- seq(lambdaMax, lambdaMax*minLambdaRatio, length = lambdaLength)

  if (family == "gaussian") {
    mdl <- gaussian_elastic(y, x, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)
  } else if (family == "binomial") {
    mdl <- binomial_elastic(y, x, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)
  } else if (family == "multinomial") {
    mdl <- multinomial_elastic(y, x, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)
  }


}



pglmEigen <- function(y, x) {

}

pglmArma <- function(y, x) {

}
