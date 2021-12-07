#'
#'
#'
calcMaxLambda <- function(family, x, y, alpha) {
  if (family == "gaussian") {
    lambdaMax <- max(t(x) %*% y)/(n*alpha)
  } else if (family == "binomial") {
    lambdaMax <- max( t(x) %*% (y - mean(y)*(1 - mean(y))) )/(n*alpha)
  } else if (family == "multinomial") {
    lambdaMax <- max(sapply(1:ncol(y), function(m) {
      max( t(x) %*% (y[,m] - mean(y[,m])*(1 - mean(y[,m]))) )/(n*alpha)
    }))
  }
  return(lambdaMax)
}

#'
gaussian_elastic <- function(
  y, x, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)
{
  #
  n <- nrow(x)
  d <- ncol(x)
  b0 <- 0
  b <- rep(0, d)
  lambdaLength <- length(lambdaVec)
  #
  b0Vec <- matrix(0, 1, lambdaLength)
  bMat <- matrix(0, d, lambdaLength)
  sy <- sum(y)
  for (lm in 1:lambdaLength) {
    lambda <- lambdaVec[lm]
    for (i in 1:maxit) {
      b0 <- (sy - sum(x %*% b))/n
      for (j in 1:d) {
        xb <- x %*% b
        yj <- b0 + xb - x[,j]*b[j]
        b[j] <- shooting((t(x[,j]) %*% (y - yj))/n, lambda*alpha)/(1 + lambda*(1 - alpha))
      }
      # Compute the mean absolute difference
      dval <- (abs(b0Vec[lm] - b0) + sum(abs(bMat[,lm] - b)))/(d + 1)
      # Update coefficients
      b0Vec[lm] <- b0
      bMat[,lm] <- b
      # Check stopping criterion
      print(c(lm, m, i, dval))
      if (dval < tol) { break }
    }
  }
  return(list(
    #conv = conv,
    b0 = b0Vec,
    b = bMat,
    lambda = lambdaVec
  ))
}

#'
binomial_elastic <- function(
  y, x, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)
{
  #
  n <- nrow(x)
  d <- ncol(x)
  b0 <- 0
  b <- rep(0, d)
  lambdaLength <- length(lambdaVec)
  #
  b0Vec <- matrix(0, 1, lambdaLength)
  bMat <- matrix(0, d, lambdaLength)
  for (lm in 1:lambdaLength) {
    lambda <- lambdaVec[lm]
    for (i in 1:maxit) {
      wz <- calcWZ(y, x, b0, b)
      b0 <- (sum(wz$w*(wz$z - x %*% b)))/n
      for (j in 1:d) {
        xb <- x %*% b
        zj <- b0 + xb - x[,j]*b[j]
        b[j] <- shooting((t(wz$w*x[,j]) %*% (wz$z - zj))/n, lambda*alpha)/(1 + lambda*(1 - alpha))
      }
      # Compute the mean absolute difference
      dval <- (abs(b0Vec[lm] - b0) + sum(abs(bMat[,lm] - b)))/(d + 1)
      # Update coefficients
      b0Vec[lm] <- b0
      bMat[,lm] <- b
      # Check stopping criterion
      print(c(lm, m, i, dval))
      if (dval < tol) { break }
    }
  }
  return(list(
    #conv = conv,
    b0 = b0Vec,
    b = bMat,
    lambda = lambdaVec
  ))
}

#'
multinomial_elastic <- function(
  y, x, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)
{
  #
  n <- nrow(x)
  d <- ncol(x)
  yd <- ncol(y)
  b0 <- rep(0, yd)
  b <- matrix(0, d, yd)
  lambdaLength <- length(lambdaVec)
  #
  b0Array <- matrix(0, yd, lambdaLength)
  bArray <- array(0, c(d, yd, lambdaLength))
  for (lm in 1:lambdaLength) {
    lambda <- lambdaVec[lm]
    for (m in 1:ncol(y)) {
      for (i in 1:maxit) {
        wz <- calcWZ(y[,m], x, b0[m], b[,m])
        b0[m] <- (sum(wz$w*(wz$z - x %*% b[,m])))/n
        for (j in 1:d) {
          wz <- calcWZ(y[,m], x, b0[m], b[,m])
          xb <- x %*% b[,m]
          zj <- b0[m] + xb - x[,j]*b[j,m]
          b[j,m] <- shooting((t(wz$w*x[,j]) %*% (wz$z - zj))/n, lambda*alpha)/(1 + lambda*(1 - alpha))
        }
        # Compute the mean absolute difference
        dval <- (sum(abs(b0Array[,lm] - b0)) + sum(abs(bArray[,,lm] - b)))/(d*yd + yd)
        # Update coefficients
        b0Array[,lm] <- b0
        bArray[,,lm] <- b
        # Check stopping criterion
        print(c(lm, m, i, dval))
        if (dval < tol) { break }
      }
    }
  }
  return(list(
    #conv = conv,
    b0 = b0Array,
    b = bArray,
    lambda = lambdaVec
  ))
}

prob_value <- function(x, b0, b) {
  val <- 1/(exp(-(b0 + x %*% b)))
  val[which(val < 1e-5)] <- 1e-5
  val[which(val > (1 - 1e-5))] <- 1 - 1e-5
  return(val)
}

calcWZ <- function(y, x, b0, b) {
  pv <- prob_value(x, b0, b)
  wv <- pv*(1 - pv)
  return(list(
    w = wv,
    z = b0 + x %*% b + (y - pv)/wv
  ))
}

#'
gaussian_predict <- function(object, x) {

}

#'
binomial_predict <- function(object, x) {
  eta <- b0 + x %*% b
  yp <- exp(eta)/(1 + exp(eta))
  as.integer(yp > 0.5)
}

#'
multinomial_predict <- function(object, x) {
  eta <- b0 + x %*% b
  yp <- exp(eta)/matrix(rowSums(exp(eta)), n, 3)
  yd <- yp/matrix(sapply(1:n, function(i) max(yp[i,])), n, 3)
  yd*(yd == 1) - y
}

shooting <- function(z, r) {
  if (z > r) {
    return(z - r)
  } else if (z < -r) {
    return(z + r)
  } else {
    return(0)
  }
}
