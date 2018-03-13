#############################################
##      JEUX DE DONNEES UTILISEES POUR LES TESTS
#############################################


####   Jeu de donnees D1
########
#' Simulate multi-view data
#' @param type One of \lbrace\code{D1}, \code{D2}, \code{D3}, \code{D4}, \code{D5}\rbrace specifying the
#' type of simulated data to generate
#' @return
#' @export
#'
#' @examples

mvsimulate <- function(type = "D1") {
  if(type == "D1") {
    n <- 1000
    lab <- NULL
    # 4 classes centrées sur le losange
    data1 <- matrix(rnorm(2 * n), nrow = n)
    data1[1:250, ] <- data1[1:250, ] + matrix(rep(c(0, 6), 250), nrow = 250, byrow = TRUE)
    data1[251:500, ] <- data1[251:500, ] + matrix(rep(c(0, -6), 250), nrow = 250, byrow = TRUE)
    data1[501:750, ] <- data1[501:750, ] + matrix(rep(c(6, 0), 250), nrow =  250, byrow = TRUE)
    data1[751:1000, ] <- data1[751:1000, ] + matrix(rep(c(-6, 0), 250), nrow = 250, byrow = TRUE)
    lab <- c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250))
    # 4 classes centrées sur un losange plus resserre
    data2 <- matrix(rnorm(2000), nrow = 1000)
    data2[1:250, ] <- data2[1:250, ] + matrix(rep(c(0, 3), 250), nrow = 250, byrow = TRUE)
    data2[251:500, ] <- data2[251:500, ] + matrix(rep(c(3, 0), 250), nrow = 250, byrow = TRUE)
    data2[501:750, ] <- data2[501:750, ] + matrix(rep(c(0, -3), 250), nrow = 250, byrow = TRUE)
    data2[751:1000, ] <- data2[751:1000, ] + matrix(rep(c(-3, 0), 250), nrow = 250, byrow = TRUE)
    lab <- cbind(lab, c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250)))
    # 4 classes mais variation que sur y
    data3 <- matrix(rnorm(2000), nrow = 1000)
    data3[1:250, ] <- data3[1:250, ] + matrix(rep(c(0, -10), 250), nrow = 250, byrow =TRUE)
    data3[251:500, ] <- data3[251:500, ] + matrix(rep(c(0, 10), 250), nrow = 250, byrow = TRUE)
    data3[501:750, ] <- data3[501:750, ] + matrix(rep(c(0, 4), 250), nrow = 250, byrow = TRUE)
    data3[751:1000, ] <- data3[751:1000, ] + matrix(rep(c(0, -4), 250), nrow = 250, byrow = TRUE)
    lab <- cbind(lab, c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250)))
    # que du bruit N(0,1)
    data4 <- 4 * matrix(rnorm(2000), nrow = 1000)
    lab <- cbind(lab, rep(1, 1000))
    colnames(lab) <- c("V1", "V2", "V3", "V4")
    D1 <- cbind(data1, data2, data3, data4)
    colnames(D1) <- c("V1.1", "V1.2", "V2.1", "V2.2", "V3.1", "V3.2", "V4.1", "V4.2")
    final_D <- D1
    final_lab <- lab
  } else if(type == "D2") {
    # plus de variance et plus resserre par rapport à D1
    data1 <- matrix(rnorm(2000), nrow = 1000)
    data1[1:250, ] <- 2 * data1[1:250, ] + matrix(rep(c(0, 5), 250), nrow = 250, byrow = TRUE)
    data1[251:500, ] <- 2 * data1[251:500, ] + matrix(rep(c(5, 0), 250), nrow = 250, byrow = TRUE)
    data1[501:750, ] <- 2 * data1[501:750, ] + matrix(rep(c(0, -5), 250), nrow = 250, byrow = TRUE)
    data1[751:1000, ] <- 2 * data1[751:1000, ] + matrix(rep(c(-5, 0), 250), nrow = 250, byrow = TRUE)
    lab <- c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250))
    # plus resserre par rapport à D1
    data2 <- matrix(rnorm(2000), nrow = 1000)
    data2[1:250, ] <- data2[1:250, ] + matrix(rep(c(0, 2), 250), nrow = 250, byrow =TRUE)
    data2[251:500, ] <- data2[251:500, ] + matrix(rep(c(2, 0), 250), nrow = 250, byrow = TRUE)
    data2[501:750, ] <- data2[501:750, ] + matrix(rep(c(0, -2), 250), nrow = 250, byrow = TRUE)
    data2[751:1000, ] <- data2[751:1000, ] + matrix(rep(c(-2, 0), 250), nrow = 250, byrow = TRUE)
    lab <- cbind(lab, c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250)))
    ## Plus resserré par rapport à D1
    data3 <- matrix(rnorm(2000), nrow = 1000)
    data3[1:250, ] <- 0.5 * data3[1:250, ] + matrix(rep(c(0, 3), 250), nrow = 250, byrow = TRUE)
    data3[251:500, ] <- 0.5 * data3[251:500, ] + matrix(rep(c(0, 4), 250), nrow = 250, byrow = TRUE)
    data3[501:750, ] <- 0.5 * data3[501:750, ] + matrix(rep(c(0, -4), 250), nrow = 250, byrow = TRUE)
    data3[751:1000, ] <- 0.5 * data3[751:1000, ] + matrix(rep(c(0, -3), 250), nrow = 250, byrow = TRUE)
    lab <- cbind(lab, c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250)))
    # variable de bruit N(0,5)
    data4 <- 25 * matrix(rnorm(2000), nrow = 1000)
    lab <- cbind(lab, rep(1, 1000))
    
    colnames(lab) <- c("V1", "V2", "V3", "V4")
    D2 <- cbind(data1, data2, data3, data4)
    colnames(D2) <-
      c("V1.1", "V1.2", "V2.1", "V2.2", "V3.1", "V3.2", "V4.1", "V4.2")
    final_D <- D1
    final_lab <- lab
  } else if(type == "D3") {
    
  } else if(type == "D4") {
  
  } else if (type == "D5") {
    
  } else error("Only types D1, D2, D3, and D4 are currently supported.")
  
  return(list(data = final_D, labels=final_lab))
}




####    Jeu de donnees D3
##########
#' Title
#'
#' @return
#' @export
#'
#' @examples
simuD3 <- function() {
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
  
  data4 = 0.1 * matrix(rnorm(2000), nrow = 1000)
  
  lab <- cbind(lab, rep(1, 1000))
  
  colnames(lab) <- c("V1", "V2", "V3", "V4")
  D3 = cbind(data1, data2, data3, data4)
  colnames(D3) <-
    c("V1.1", "V1.2", "V2.1", "V2.2", "V3.1", "V3.2", "V4.1", "V4.2")
  return(list(data = D3, lab = lab))
}

####    Jeu de donnees D4
##########
#' Title
#'
#' @return
#' @export
#'
#' @examples
simuD4 <- function() {
  data1 = matrix(rnorm(2000), nrow = 1000)
  data1[1:250, ] = data1[1:250, ] + matrix(rep(c(0, 5), 250), nrow = 250, byrow =
                                             T)
  data1[251:500, ] = data1[251:500, ] + matrix(rep(c(5, 0), 250), nrow =
                                                 250, byrow = T)
  data1[501:750, ] = data1[501:750, ] + matrix(rep(c(0, -5), 250), nrow =
                                                 250, byrow = T)
  data1[751:1000, ] = data1[751:1000, ] + matrix(rep(c(-5, 0), 250), nrow =
                                                   250, byrow = T)
  lab <- c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250))
  
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
  
  data4 = matrix(rnorm(2000), nrow = 1000)
  data4[1:125, ] = 0.5 * data4[1:125, ] + matrix(rep(c(0, 5), 125), nrow =
                                                   125, byrow = T)
  data4[126:250, ] = 0.5 * data4[126:250, ] + matrix(rep(c(0, -5), 125), nrow =
                                                       125, byrow = T)
  data4[251:375, ] = 0.5 * data4[251:375, ] + matrix(rep(c(0, 5), 125), nrow =
                                                       125, byrow = T)
  data4[376:500, ] = 0.5 * data4[376:500, ] + matrix(rep(c(0, -5), 125), nrow =
                                                       125, byrow = T)
  data4[501:625, ] = 0.5 * data4[501:625, ] + matrix(rep(c(5, 0), 125), nrow =
                                                       125, byrow = T)
  data4[626:750, ] = 0.5 * data4[626:750, ] + matrix(rep(c(-5, 0), 125), nrow =
                                                       125, byrow = T)
  data4[751:875, ] = 0.5 * data4[751:875, ] + matrix(rep(c(5, 0), 125), nrow =
                                                       125, byrow = T)
  data4[876:1000, ] = 0.5 * data4[876:1000, ] + matrix(rep(c(-5, 0), 125), nrow =
                                                         125, byrow = T)
  lab <- cbind(lab, sort(rep(seq(1:8), 125)))
  colnames(lab) <- c("V1", "V2", "V3", "V4")
  D4 = cbind(data1, data2, data3, data4)
  colnames(D4) <-
    c("V1.1", "V1.2", "V2.1", "V2.2", "V3.1", "V3.2", "V4.1", "V4.2")
  return(list(data = D4, lab = lab))
}


####    Jeu de donnees D5
##########
#' Title
#'
#' @return
#' @export
#'
#' @examples
simuD5 <- function() {
  data1 = matrix(rnorm(4000), nrow = 2000)
  data1[1:250, ] = data1[1:250, ] + matrix(rep(c(0, 6), 250), nrow = 250, byrow =
                                             T)
  data1[251:500, ] = data1[251:500, ] + matrix(rep(c(0, -6), 250), nrow =
                                                 250, byrow = T)
  data1[501:750, ] = data1[501:750, ] + matrix(rep(c(6, 0), 250), nrow =
                                                 250, byrow = T)
  data1[751:1000, ] = data1[751:1000, ] + matrix(rep(c(-6, 0), 250), nrow =
                                                   250, byrow = T)
  data1[1251:1500, ] = data1[1251:1500, ] + matrix(rep(c(-6, 6), 250), nrow =
                                                     250, byrow = T)
  data1[1001:1250, ] = data1[1001:1250, ] + matrix(rep(c(-6, -6), 250), nrow =
                                                     250, byrow = T)
  data1[1501:1750, ] = data1[1501:1750, ] + matrix(rep(c(6, 6), 250), nrow =
                                                     250, byrow = T)
  data1[1751:2000, ] = data1[1751:2000, ] + matrix(rep(c(6, -6), 250), nrow =
                                                     250, byrow = T)
  lab <- sort(rep(seq(1:8), 250))
  data2 = matrix(rnorm(4000), nrow = 2000)
  data2[1:250, ] = data2[1:250, ] + matrix(rep(c(0, 3), 250), nrow = 250, byrow =
                                             T)
  data2[251:500, ] = data2[251:500, ] + matrix(rep(c(3, 0), 250), nrow =
                                                 250, byrow = T)
  data2[501:750, ] = data2[501:750, ] + matrix(rep(c(0, -3), 250), nrow =
                                                 250, byrow = T)
  data2[751:1000, ] = data2[751:1000, ] + matrix(rep(c(-3, 0), 250), nrow =
                                                   250, byrow = T)
  data2[1251:1500, ] = data2[1251:1500, ] + matrix(rep(c(-3, 3), 250), nrow =
                                                     250, byrow = T)
  data2[1001:1250, ] = data2[1001:1250, ] + matrix(rep(c(-3, -3), 250), nrow =
                                                     250, byrow = T)
  data2[1501:1750, ] = data2[1501:1750, ] + matrix(rep(c(3, 3), 250), nrow =
                                                     250, byrow = T)
  data2[1751:2000, ] = data2[1751:2000, ] + matrix(rep(c(3, -3), 250), nrow =
                                                     250, byrow = T)
  lab <- cbind(lab, sort(rep(seq(1:8), 250)))
  data3 = matrix(rnorm(4000), nrow = 2000)
  data3[1:250, ] = data3[1:250, ] + matrix(rep(c(0, -6), 250), nrow = 250, byrow =
                                             T)
  data3[251:500, ] = data3[251:500, ] + matrix(rep(c(0, 6), 250), nrow =
                                                 250, byrow = T)
  data3[501:750, ] = data3[501:750, ] + matrix(rep(c(0, 3), 250), nrow =
                                                 250, byrow = T)
  data3[751:1000, ] = data3[751:1000, ] + matrix(rep(c(0, -3), 250), nrow =
                                                   250, byrow = T)
  data3[1251:1500, ] = data3[1251:1500, ] + matrix(rep(c(-3, 0), 250), nrow =
                                                     250, byrow = T)
  data3[1001:1250, ] = data3[1001:1250, ] + matrix(rep(c(-6, 0), 250), nrow =
                                                     250, byrow = T)
  data3[1501:1750, ] = data3[1501:1750, ] + matrix(rep(c(3, 0), 250), nrow =
                                                     250, byrow = T)
  data3[1751:2000, ] = data3[1751:2000, ] + matrix(rep(c(6, 0), 250), nrow =
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
