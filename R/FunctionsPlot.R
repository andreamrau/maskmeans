#################################################
##      FONCTIONS DE PLOT
#################################################

#' Title
#'
#' @param data 
#' @param lab 
#'
#' @return
#' @export
#'
#' @examples
PlotSimuAntoine <- function(data, lab) {
  g1 <- ggplot(data, aes(x = V1.1, y = V1.2)) +
    geom_point(aes(colour = factor(lab[, 1])))
  g2 <- ggplot(data, aes(x = V2.1, y = V2.2)) +
    geom_point(aes(colour = factor(lab[, 2])))
  g3 <- ggplot(data, aes(x = V3.1, y = V3.2)) +
    geom_point(aes(colour = factor(lab[, 3])))
  g4 <- ggplot(data, aes(x = V4.1, y = V4.2)) +
    geom_point(aes(colour = factor(lab[, 4])))
  grid.arrange(g1, g2, g3, g4, ncol = 2, nrow = 2)
}

#' Title
#'
#' @param data 
#' @param classif 
#'
#' @return
#' @export
#'
#' @examples
PlotClassif <- function(data, classif) {
  don <- data.frame(classif = as.factor(classif), data = data)
  ggpairs(don,
          mapping = aes(color = classif),
          columns = colnames(don)[-1])
}

################################
# exploitation graphique de la sortie de hmv1
################################
#' Title
#'
#' @param Res 
#' @param dendro 
#' @param heights 
#' @param weights 
#' @param CRIT 
#'
#' @return
#' @export
#'
#' @examples
PlotResHMV <- function(Res,
                       dendro = TRUE,
                       heights = TRUE,
                       weights = TRUE,
                       CRIT = TRUE) {
  if (heights == TRUE) {
    # plot des hauteurs de l'arbre
    g1 <-
      ggplot(data.frame(x = as.factor(seq(
        2, length(Res$height)
      )), y = rev(Res$height)[-1]), aes(x = x, y = y)) +
      geom_bar(stat = "identity") +
      ylab("Height") +
      xlab("Nb of clusters")
    print(g1)
  }
  
  if (dendro == TRUE) {
    # dendrogramme
    g2 <- ggdendrogram(Res$R, rotate = FALSE, size = 2)
    print(g2)
  }
  
  if (weights == TRUE) {
    # graphe de l'évolution des poids
    df <- data.frame(cbind(seq(1, ncol(Res$weights)), t(Res$weights)))
    rownames(df) <- NULL
    colnames(df) <- c("time", as.character(seq(1:nrow(Res$weights))))
    ## !!!
    Test <- melt(df, "time")
    g3 <- ggplot(Test, aes(time, value)) +
      geom_area(aes(fill = variable, group = variable)) +
      ylab("weights")
    print(g3)
    
    g4 <- ggplot(Test, aes(time, value, group = variable, colour = variable)) +
      geom_point() +
      geom_line(aes(lty = variable))
    print(g4)
  }
  
  if (CRIT == TRUE) {
    g5 <-
      ggplot(data.frame(x = seq(2, length(Res$order)), y = rev(Res$CRIT)), aes(x =
                                                                                 x, y = y)) +
      geom_point() +
      geom_line() +
      ylab("CRIT") +
      xlab("Nb of clusters")
    print(g5)
  }
}

#############################################
##   exploitation graphique de la sortie splittingClusters
#############################################
#' Title
#'
#' @param Res 
#'
#' @return
#' @export
#'
#' @examples
PlotResSplitting <- function(Res) {
  # plot de CRIT
  aux1 <-
    data.frame(K = apply(Res$clustersplithist, 2, max),
               CRIT = Res$CRIT[-1])
  g1 <- ggplot(aux1, aes(x = K, y = CRIT)) +
    geom_point() +
    geom_line() +
    xlab("Nb of clusters") +
    ylab("Criterion value")
  print(g1)
  
  # plot des poids
  # graphe de l'évolution des poids
  df <- data.frame(cbind(seq(1, ncol(Res$weights)), t(Res$weights)))
  rownames(df) <- NULL
  colnames(df) <- c("time", as.character(seq(1:nrow(Res$weights))))
  #!!!!!!!
  Test <- melt(df, "time")
  g2 <- ggplot(Test, aes(time, value)) +
    geom_area(aes(fill = variable, group = variable)) +
    ylab("weights")
  print(g2)
  
  g3 <- ggplot(Test, aes(time, value, group = variable, colour = variable)) +
    geom_point() +
    geom_line(aes(lty = variable))
  print(g3)
  
  # plot du splitting   ????
}




#############################################
##   exploitation graphique de la sortie splittingClustersbis
#############################################
#' Title
#'
#' @param Res 
#' @param plotCRIT 
#' @param weights_step 
#' @param weights_clust 
#'
#' @return
#' @export
#'
#' @examples
PlotResSplittingbis <-
  function(Res,
           plotCRIT = TRUE,
           weights_step = FALSE,
           weights_clust = FALSE) {
    # plot de CRIT
    if (plotCRIT == TRUE) {
      aux1 <-
        data.frame(K = apply(Res$clustersplithist, 2, max),
                   CRIT = Res$CRIT[-1])
      g1 <- ggplot(aux1, aes(x = K, y = CRIT)) +
        geom_point() +
        geom_line() +
        xlab("Nb of clusters") +
        ylab("Criterion value")
      print(g1)
    }
    
    # plot des poids
    # graphes de l'évolution des poids
    if (weights_step != FALSE) {
      for (h in 1:length(weights_step)) {
        j <- weights_step[h]
        #!!!!!!!
        A <- melt(Res$weights[[j]])
        A[, 1] <- as.factor(A[, 1])
        A[, 2] <- as.factor(A[, 2])
        colnames(A) = c("clust", "mv", "weights")
        g <- ggplot(A, aes(fill = mv, y = weights, x = clust)) +
          geom_bar(stat = "identity") +
          ggtitle(paste("step ", j - 1, " - split ", Res$ksplit[j]))
        print(g)
      }
    }
    
    if (weights_clust != FALSE) {
      #!!!!!!
      A <- melt(Res$weights[[length(Res$weights)]])
      colnames(A) <- c("clust", "mv", "weights")
      I <- NULL
      for (h in 1:length(weights_clust))
        I <- c(I, which(A$clust == weights_clust[h]))
      A[, 1] <- as.factor(A[, 1])
      A[, 2] <- as.factor(A[, 2])
      g2 <- ggplot(A[I, ], aes(fill = mv, y = weights, x = clust)) +
        geom_bar(stat = "identity")
      print(g2)
    }
    # plot du splitting   ????
  }


###########################
##   boxplot des probapost
###########################
#' Title
#'
#' @param probapost 
#'
#' @return
#' @export
#'
#' @examples
BoxplotProbaPost <- function(probapost) {
  aux <- data.frame(probamax = apply(probapost, 1, max),
                   lab = as.factor(apply(probapost, 1, which.max)))
  p <- ggplot(aux, aes(x = lab, y = probamax)) +
    geom_boxplot()
  print(p)
}

###########################
##    % de probapost > seuil
############################

#' Title
#'
#' @param Res 
#' @param probapost.init 
#' @param seuil 
#'
#' @return
#' @export
#'
#' @examples
plotprobaseuil <- function(Res, probapost.init, seuil = 0.8) {
  aux <- NULL
  for (k in 2:length(Res$CRIT)) {
    aux <-
      c(aux, sum(apply(
        cutreeNewProbapost(Res$R, K = k, probapost.init)$probapost,
        1,
        max
      ) > seuil) * 100 / nrow(probapost.init))
  }
  p <- ggplot(data.frame(k = 2:length(Res$CRIT), aux = aux), aes(k, aux)) +
    geom_line() +
    geom_point() +
    xlab("nb of clusters") +
    ylab("percent of probamax>seuil")
  print(p)
}

###########################
##    fonction clustree pour visualiser les étapes de splitting
##     prend l'objet Res venant de splittingClusters et les données X
############################
#' Title
#'
#' @param Res 
#' @param X 
#'
#' @return
#' @export
#'
#' @examples
clustreebis <- function(Res, X) {
  aux <- Res$clustersplithist
  colnames(aux) <- paste0("K", apply(Res$clustersplithist, 2, max))
  dataPlot <- cbind(X, aux)
  clustree(data.frame(dataPlot), prefix = "K")
}
