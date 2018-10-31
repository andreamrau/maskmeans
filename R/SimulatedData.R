#' Simulate multi-view data
#' 
#' Generate simulated multi-view data to illustrate the agglomerative and splitting versions of the
#' multi-view K-means algorithm. 
#' 
#' @param beta Multiplicative factor controlling the spread of values around the origin in the first (\code{beta}), 
#' second (\code{beta}), and sixth (\code{beta_V6_multiplier} x \code{beta}) views. Defaults to 4.
#' @param n  Number of observations per cluster for simulations. Defaults to 100.
#' @param K Number of clusters in the first view for simulations. As this view is made up
#' of a single cluster at the origin surrounded by evenly spaced clusters in a circular pattern around
#' it, this number should be odd. Defaults to 7.
#' @param sigma Variance of noise to be added to views 2, 4, 5, and 6.
#' Defaults to 1.5.
#' @param beta_V6_multiplier Multiplicative parameter for the spread of values around the origin in view 6 
#' (\code{beta_V6_multiplier} x \code{beta}). Defaults to 1.5.
#' @param ... Additional optional parameters. 
#' @return
#' \item{data }{Multi-view simulated data}
#' \item{labels }{Matrix of dimension \code{n} x \code{v}, where \code{n} is the number
#' of observations and \code{v} the number of views, representing the true labels used
#' to generate the data}
#' @export
#'
#' @examples
#' set.seed(12345)
#' sim_1 <- mv_simulate(type = "D1")
#' sim_2 <- mv_simulate(type = "D2")
#' sim_3 <- mv_simulate(type = "D3")
#' sim_4 <- mv_simulate(type = "D4")
#' sim_5 <- mv_simulate(type = "D5")
#' sim_6a <- mv_simulate(type = "D6")
#' sim_6b <- mv_simulate(type = "D6", delta=7, n=200, K=5, sigma=0.5)

mv_simulate <- function(beta=4, n=100, K=7, sigma=1.5, beta_V6_multiplier=1.5, ...) {
  
  ## Parse ellipsis function
  providedArgs <- list(...)
  arg.user <- list(type="D6")
  arg.user[names(providedArgs)] <- providedArgs
  type <- arg.user$type

  delta <- beta
  sigma_V6 <- sigma
  delta_V6 <- beta_V6_multiplier * beta
  
  if(type == "D1") {
    sim <- simuD1()
  } else if(type == "D2") {
    sim <- simuD2()
  } else if(type == "D3") {
    sim <- simuD3()
  } else if(type == "D4") {
    sim <- simuD4()
  } else if (type == "D5") {
    sim <- simuD5()
  } else if (type == "D6") {
    if(! K %% 2) stop("Simulation D6 is only supported for an odd number of clusters.")
    if(delta < 0) stop("delta should be a nonnegative multiplicative factor.")
    if(n < 0 | n %% 1) stop("n should be a nonnegative number of observations.")
    if(K < 0 | K %% 1) stop("K should be a nonnegative number of clusters")
    if(sigma < 0) stop("sigma should be a nonnegative variance for the added noise.")
    sim <- simuD6(delta = delta, n=n, K=(K-1)/2, sigma=sigma, delta_V6 = delta_V6,
                  sigma_V6 = sigma_V6)
  } else  stop("Only types D1, D2, D3, D4, D5, and D6 are currently supported.")
  return(list(data = sim$data, labels=sim$lab))
}


#----------------------------------------------------------------------------------
## Unexported simulation functions

simuD1 <- function() {
  n = 1000
  lab <- NULL
  # 4 classes centrées sur le losange
  data1 = matrix(rnorm(2 * n), nrow = n)
  data1[1:250, ] = data1[1:250, ] + matrix(rep(c(0, 6), 250), nrow = 250, byrow =
                                             T)
  data1[251:500, ] = data1[251:500, ] + matrix(rep(c(0, -6), 250), nrow =
                                                 250, byrow = T)
  data1[501:750, ] = data1[501:750, ] + matrix(rep(c(6, 0), 250), nrow =
                                                 250, byrow = T)
  data1[751:1000, ] = data1[751:1000, ] + matrix(rep(c(-6, 0), 250), nrow =
                                                   250, byrow = T)
  lab <- c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250))
  # 4 classes centrées sur un losange plus resserre
  data2 = matrix(rnorm(2000), nrow = 1000)
  data2[1:250, ] = data2[1:250, ] + matrix(rep(c(0, 3), 250), nrow = 250, byrow =
                                             T)
  data2[251:500, ] = data2[251:500, ] + matrix(rep(c(3, 0), 250), nrow =
                                                 250, byrow = T)
  data2[501:750, ] = data2[501:750, ] + matrix(rep(c(0, -3), 250), nrow =
                                                 250, byrow = T)
  data2[751:1000, ] = data2[751:1000, ] + matrix(rep(c(-3, 0), 250), nrow =
                                                   250, byrow = T)
  lab <- cbind(lab, c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250)))
  # 4 classes mais variation que sur y
  data3 = matrix(rnorm(2000), nrow = 1000)
  data3[1:250, ] = data3[1:250, ] + matrix(rep(c(0, -10), 250), nrow = 250, byrow =
                                             T)
  data3[251:500, ] = data3[251:500, ] + matrix(rep(c(0, 10), 250), nrow =
                                                 250, byrow = T)
  data3[501:750, ] = data3[501:750, ] + matrix(rep(c(0, 4), 250), nrow =
                                                 250, byrow = T)
  data3[751:1000, ] = data3[751:1000, ] + matrix(rep(c(0, -4), 250), nrow =
                                                   250, byrow = T)
  lab <- cbind(lab, c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250)))
  # que du bruit N(0,1)
  data4 = 4 * matrix(rnorm(2000), nrow = 1000)
  lab <- cbind(lab, rep(1, 1000))
  
  colnames(lab) <- c("V1", "V2", "V3", "V4")
  D1 = cbind(data1, data2, data3, data4)
  colnames(D1) <-
    c("V1.1", "V1.2", "V2.1", "V2.2", "V3.1", "V3.2", "V4.1", "V4.2")
  
  return(list(data = D1, lab = lab))
}

simuD2 <- function() {
  # plus de variance et plus resserre par rapport à D1
  data1 = matrix(rnorm(2000), nrow = 1000)
  data1[1:250, ] = 2 * data1[1:250, ] + matrix(rep(c(0, 5), 250), nrow =
                                                 250, byrow = T)
  data1[251:500, ] = 2 * data1[251:500, ] + matrix(rep(c(5, 0), 250), nrow =
                                                     250, byrow = T)
  data1[501:750, ] = 2 * data1[501:750, ] + matrix(rep(c(0, -5), 250), nrow =
                                                     250, byrow = T)
  data1[751:1000, ] = 2 * data1[751:1000, ] + matrix(rep(c(-5, 0), 250), nrow =
                                                       250, byrow = T)
  lab <- c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250))
  # plus resserre par rapport à D1
  data2 = matrix(rnorm(2000), nrow = 1000)
  data2[1:250, ] = data2[1:250, ] + matrix(rep(c(0, 2), 250), nrow = 250, byrow =
                                             T)
  data2[251:500, ] = data2[251:500, ] + matrix(rep(c(2, 0), 250), nrow =
                                                 250, byrow = T)
  data2[501:750, ] = data2[501:750, ] + matrix(rep(c(0, -2), 250), nrow =
                                                 250, byrow = T)
  data2[751:1000, ] = data2[751:1000, ] + matrix(rep(c(-2, 0), 250), nrow =
                                                   250, byrow = T)
  lab <- cbind(lab, c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250)))
  ## Plus resserré par rapport à D1
  data3 = matrix(rnorm(2000), nrow = 1000)
  data3[1:250, ] = 0.5 * data3[1:250, ] + matrix(rep(c(0, 3), 250), nrow =
                                                   250, byrow = T)
  data3[251:500, ] = 0.5 * data3[251:500, ] + matrix(rep(c(0, 4), 250), nrow =
                                                       250, byrow = T)
  data3[501:750, ] = 0.5 * data3[501:750, ] + matrix(rep(c(0, -4), 250), nrow =
                                                       250, byrow = T)
  data3[751:1000, ] = 0.5 * data3[751:1000, ] + matrix(rep(c(0, -3), 250), nrow =
                                                         250, byrow = T)
  lab <- cbind(lab, c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250)))
  # variable de bruit N(0,5)
  data4 = 25 * matrix(rnorm(2000), nrow = 1000)
  lab <- cbind(lab, rep(1, 1000))
  
  colnames(lab) <- c("V1", "V2", "V3", "V4")
  D2 = cbind(data1, data2, data3, data4)
  colnames(D2) <-
    c("V1.1", "V1.2", "V2.1", "V2.2", "V3.1", "V3.2", "V4.1", "V4.2")
  return(list(data = D2, lab = lab))
}

simuD3 <- function() {
  data1 = matrix(rnorm(2000), nrow = 1000)
  data1[1:250,] = 2 * data1[1:250,] + matrix(rep(c(0, 5), 250), nrow =
                                               250, byrow = T)
  data1[251:500,] = 2 * data1[251:500,] + matrix(rep(c(5, 0), 250), nrow =
                                                   250, byrow = T)
  data1[501:750,] = 2 * data1[501:750,] + matrix(rep(c(0,-5), 250), nrow =
                                                   250, byrow = T)
  data1[751:1000,] = 2 * data1[751:1000,] + matrix(rep(c(-5, 0), 250), nrow =
                                                     250, byrow = T)
  lab <- c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250))
  
  data2 = matrix(rnorm(2000), nrow = 1000)
  data2[1:250,] = data2[1:250,] + matrix(rep(c(0, 2), 250), nrow = 250, byrow =
                                           T)
  data2[251:500,] = data2[251:500,] + matrix(rep(c(2, 0), 250), nrow =
                                               250, byrow = T)
  data2[501:750,] = data2[501:750,] + matrix(rep(c(0,-2), 250), nrow =
                                               250, byrow = T)
  data2[751:1000,] = data2[751:1000,] + matrix(rep(c(-2, 0), 250), nrow =
                                                 250, byrow = T)
  lab <-
    cbind(lab, c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250)))
  
  data3 = matrix(rnorm(2000), nrow = 1000)
  data3[1:250,] = 0.5 * data3[1:250,] + matrix(rep(c(0, 3), 250), nrow =
                                                 250, byrow = T)
  data3[251:500,] = 0.5 * data3[251:500,] + matrix(rep(c(0, 4), 250), nrow =
                                                     250, byrow = T)
  data3[501:750,] = 0.5 * data3[501:750,] + matrix(rep(c(0,-4), 250), nrow =
                                                     250, byrow = T)
  data3[751:1000,] = 0.5 * data3[751:1000,] + matrix(rep(c(0,-3), 250), nrow =
                                                       250, byrow = T)
  lab <-
    cbind(lab, c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250)))
  
  data4 = 0.1 * matrix(rnorm(2000), nrow = 1000)
  
  lab <- cbind(lab, rep(1, 1000))
  
  colnames(lab) <- c("V1", "V2", "V3", "V4")
  D3 = cbind(data1, data2, data3, data4)
  colnames(D3) <-
    c("V1.1", "V1.2", "V2.1", "V2.2", "V3.1", "V3.2", "V4.1", "V4.2")
  return(list(data = D3, lab = lab))
}

simuD4 <- function() {
  data1 = matrix(rnorm(2000), nrow = 1000)
  data1[1:250,] = data1[1:250,] + matrix(rep(c(0, 5), 250), nrow = 250, byrow =
                                           T)
  data1[251:500,] = data1[251:500,] + matrix(rep(c(5, 0), 250), nrow =
                                               250, byrow = T)
  data1[501:750,] = data1[501:750,] + matrix(rep(c(0,-5), 250), nrow =
                                               250, byrow = T)
  data1[751:1000,] = data1[751:1000,] + matrix(rep(c(-5, 0), 250), nrow =
                                                 250, byrow = T)
  lab <- c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250))
  
  data2 = matrix(rnorm(2000), nrow = 1000)
  data2[1:250,] = data2[1:250,] + matrix(rep(c(0, 2), 250), nrow = 250, byrow =
                                           T)
  data2[251:500,] = data2[251:500,] + matrix(rep(c(2, 0), 250), nrow =
                                               250, byrow = T)
  data2[501:750,] = data2[501:750,] + matrix(rep(c(0,-2), 250), nrow =
                                               250, byrow = T)
  data2[751:1000,] = data2[751:1000,] + matrix(rep(c(-2, 0), 250), nrow =
                                                 250, byrow = T)
  lab <-
    cbind(lab, c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250)))
  
  data3 = matrix(rnorm(2000), nrow = 1000)
  data3[1:250,] = 0.5 * data3[1:250,] + matrix(rep(c(0, 3), 250), nrow =
                                                 250, byrow = T)
  data3[251:500,] = 0.5 * data3[251:500,] + matrix(rep(c(0, 4), 250), nrow =
                                                     250, byrow = T)
  data3[501:750,] = 0.5 * data3[501:750,] + matrix(rep(c(0,-4), 250), nrow =
                                                     250, byrow = T)
  data3[751:1000,] = 0.5 * data3[751:1000,] + matrix(rep(c(0,-3), 250), nrow =
                                                       250, byrow = T)
  lab <-
    cbind(lab, c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250)))
  
  data4 = matrix(rnorm(2000), nrow = 1000)
  data4[1:125,] = 0.5 * data4[1:125,] + matrix(rep(c(0, 5), 125), nrow =
                                                 125, byrow = T)
  data4[126:250,] = 0.5 * data4[126:250,] + matrix(rep(c(0,-5), 125), nrow =
                                                     125, byrow = T)
  data4[251:375,] = 0.5 * data4[251:375,] + matrix(rep(c(0, 5), 125), nrow =
                                                     125, byrow = T)
  data4[376:500,] = 0.5 * data4[376:500,] + matrix(rep(c(0,-5), 125), nrow =
                                                     125, byrow = T)
  data4[501:625,] = 0.5 * data4[501:625,] + matrix(rep(c(5, 0), 125), nrow =
                                                     125, byrow = T)
  data4[626:750,] = 0.5 * data4[626:750,] + matrix(rep(c(-5, 0), 125), nrow =
                                                     125, byrow = T)
  data4[751:875,] = 0.5 * data4[751:875,] + matrix(rep(c(5, 0), 125), nrow =
                                                     125, byrow = T)
  data4[876:1000,] = 0.5 * data4[876:1000,] + matrix(rep(c(-5, 0), 125), nrow =
                                                       125, byrow = T)
  lab <- cbind(lab, sort(rep(seq(1:8), 125)))
  colnames(lab) <- c("V1", "V2", "V3", "V4")
  D4 = cbind(data1, data2, data3, data4)
  colnames(D4) <-
    c("V1.1", "V1.2", "V2.1", "V2.2", "V3.1", "V3.2", "V4.1", "V4.2")
  return(list(data = D4, lab = lab))
}

simuD5 <- function() {
  data1 = matrix(rnorm(4000), nrow = 2000)
  data1[1:250,] = data1[1:250,] + matrix(rep(c(0, 6), 250), nrow = 250, byrow =
                                           T)
  data1[251:500,] = data1[251:500,] + matrix(rep(c(0,-6), 250), nrow =
                                               250, byrow = T)
  data1[501:750,] = data1[501:750,] + matrix(rep(c(6, 0), 250), nrow =
                                               250, byrow = T)
  data1[751:1000,] = data1[751:1000,] + matrix(rep(c(-6, 0), 250), nrow =
                                                 250, byrow = T)
  data1[1251:1500,] = data1[1251:1500,] + matrix(rep(c(-6, 6), 250), nrow =
                                                   250, byrow = T)
  data1[1001:1250,] = data1[1001:1250,] + matrix(rep(c(-6,-6), 250), nrow =
                                                   250, byrow = T)
  data1[1501:1750,] = data1[1501:1750,] + matrix(rep(c(6, 6), 250), nrow =
                                                   250, byrow = T)
  data1[1751:2000,] = data1[1751:2000,] + matrix(rep(c(6,-6), 250), nrow =
                                                   250, byrow = T)
  lab <- sort(rep(seq(1:8), 250))
  data2 = matrix(rnorm(4000), nrow = 2000)
  data2[1:250,] = data2[1:250,] + matrix(rep(c(0, 3), 250), nrow = 250, byrow =
                                           T)
  data2[251:500,] = data2[251:500,] + matrix(rep(c(3, 0), 250), nrow =
                                               250, byrow = T)
  data2[501:750,] = data2[501:750,] + matrix(rep(c(0,-3), 250), nrow =
                                               250, byrow = T)
  data2[751:1000,] = data2[751:1000,] + matrix(rep(c(-3, 0), 250), nrow =
                                                 250, byrow = T)
  data2[1251:1500,] = data2[1251:1500,] + matrix(rep(c(-3, 3), 250), nrow =
                                                   250, byrow = T)
  data2[1001:1250,] = data2[1001:1250,] + matrix(rep(c(-3,-3), 250), nrow =
                                                   250, byrow = T)
  data2[1501:1750,] = data2[1501:1750,] + matrix(rep(c(3, 3), 250), nrow =
                                                   250, byrow = T)
  data2[1751:2000,] = data2[1751:2000,] + matrix(rep(c(3,-3), 250), nrow =
                                                   250, byrow = T)
  lab <- cbind(lab, sort(rep(seq(1:8), 250)))
  data3 = matrix(rnorm(4000), nrow = 2000)
  data3[1:250,] = data3[1:250,] + matrix(rep(c(0,-6), 250), nrow = 250, byrow =
                                           T)
  data3[251:500,] = data3[251:500,] + matrix(rep(c(0, 6), 250), nrow =
                                               250, byrow = T)
  data3[501:750,] = data3[501:750,] + matrix(rep(c(0, 3), 250), nrow =
                                               250, byrow = T)
  data3[751:1000,] = data3[751:1000,] + matrix(rep(c(0,-3), 250), nrow =
                                                 250, byrow = T)
  data3[1251:1500,] = data3[1251:1500,] + matrix(rep(c(-3, 0), 250), nrow =
                                                   250, byrow = T)
  data3[1001:1250,] = data3[1001:1250,] + matrix(rep(c(-6, 0), 250), nrow =
                                                   250, byrow = T)
  data3[1501:1750,] = data3[1501:1750,] + matrix(rep(c(3, 0), 250), nrow =
                                                   250, byrow = T)
  data3[1751:2000,] = data3[1751:2000,] + matrix(rep(c(6, 0), 250), nrow =
                                                   250, byrow = T)
  lab <- cbind(lab, sort(rep(seq(1:8), 250)))
  data4 = 4 * matrix(rnorm(4000), nrow = 2000)
  lab <- cbind(lab, rep(1, 2000))
  colnames(lab) <- c("V1", "V2", "V3", "V4")
  D5 = cbind(data1, data2, data3, data4)
  colnames(D5) <-
    c("V1.1", "V1.2", "V2.1", "V2.2", "V3.1", "V3.2", "V4.1", "V4.2")
  return(list(data = D5, lab = lab))
}


simuD6 <- function(delta, n, K, sigma, sigma_V6, delta_V6) {
  data1 = matrix(rnorm(2 * n * (2 * K + 1)), ncol = 2)
  colnames(data1) <- c("V1.1", "V1.2")
  mvg = ((1:(2 * K + 1)) * n) - n + 1
  mvd = (1:(2 * K + 1)) * n
  theta = (2 * pi * (1:(2 * K))) / (2 * K)
  for (k in 1:(2 * K)) {
    data1[mvg[k]:mvd[k],] = data1[mvg[k]:mvd[k],] + 
      matrix(rep(delta * c(cos(theta[k]), sin(theta[k])), n), ncol = 2, byrow = T)
  }
  lab <- sort(rep(seq(1:(2 * K + 1)), n))
  D6 <- data1
  
  data2 = matrix(rnorm(2 * n * (2 * K + 1), mean = 0, sd = sigma), ncol = 2)
  colnames(data2) <- c("V2.1", "V2.2")
  for (k in 1:K) {
    thetabis = mean(c(theta[(2 * k) - 1] , theta[2 * k]))
    data2[mvg[(2 * k) - 1]:mvd[(2 * k) - 1],] = 
      data2[mvg[(2 * k) - 1]:mvd[(2 * k) - 1],] + matrix(rep(delta * c(cos(thetabis), sin(thetabis)), n), 
                                                         ncol = 2, byrow = T)
    data2[mvg[2 * k]:mvd[2 * k], ] =  data2[mvg[2 * k]:mvd[2 * k], ] + 
      matrix(rep(delta * c(cos(thetabis), sin(thetabis)), n), ncol = 2, byrow = T)
  }
  lab <- cbind(lab, c(sort(rep(seq(
    1:K
  ), 2 * n)), rep(K + 1, n)))
  D6 = cbind(D6, data2)
  
  perm = sample(n * (2 * K + 1))
  lab <- cbind(lab, lab[perm, 1])
  data3 = data1[perm, ]
  colnames(data3) <- c("V3.1", "V3.2")
  D6 = cbind(D6, data3)
  
  I1 = which(data1[, 1] < -delta / 2)
  I2 = which(data1[, 1] > delta / 2)
  V.4 = rnorm(n * (2 * K + 1), mean = 0, sd = sigma)
  V.4[I1] = V.4[I1] - delta / 2
  V.4[I2] = V.4[I2] + delta / 2
  lab4 = rep(3, n * (2 * K + 1))
  lab4[I1] = 1
  lab4[I2] = 2
  lab <- cbind(lab, lab4)
  D6 = cbind(D6, V.4)
  
  V.5 = V.4
  V.5[I1] = V.4[I1] - delta / 2
  V.5[I2] = V.4[I2] + delta / 2
  lab <- cbind(lab, lab4)
  D6 = cbind(D6, V.5)
  
  V6 = matrix(rnorm(2 * n * (2 * K + 1), sd = sigma_V6), ncol = 2)
  colnames(V6) <- c("V6.1", "V6.2")
  lab6 = rep(4, n * (2 * K + 1))
  a = sample(2 * K)[1:3]
  for (i in 1:length(a)) {
    I = which(lab[, 1] == a[i])
     V6[I, ] = V6[I, ] + matrix(rep(delta_V6 * c(cos(theta[a[i]]), sin(theta[a[i]])), n), ncol =
                                  2, byrow = TRUE)
 #   V6[I, ] = V6[I, ] + matrix(rep(delta * c(cos(theta[a[i]]), sin(theta[a[i]])), n), ncol =
  #                               2, byrow = TRUE)
    lab6[I] = i
  }
  lab <- cbind(lab, lab6)
  D6 = cbind(D6, V6)

  return(list(data = data.frame(D6), lab = data.frame(lab)))
}
