##   Splitting de clusters inspirée du multiview.
## Entrées :
#  cluster.init : une classif du  re
#  X : les données (le 1er bloc correspond aux données principales utilisées pour obtenir classif)
# mv : la taille de chaque multiview
# weightsopt : if TRUE, multiview version; if FALSE weights = 1/V


#' Title
#'
#' @param X 
#' @param mv 
#' @param gamma 
#' @param Kmax 
#' @param cluster.init 
#' @param weightsopt 
#' @param testkmeans 
#'
#' @return
#' @export
#'
#' @examples
splittingClusters = function(X,
                             mv,
                             gamma,
                             Kmax,
                             cluster.init,
                             weightsopt = T,
                             testkmeans = T) {
  ref = c(0, cumsum(mv))
  CRIT <- 0    #evolution critere general
  clustersplithist = matrix(c(cluster.init, cluster.init), ncol = 2)  #pour garder l'historique du splitting
  
  # initialisation des centres
  centers <-
    matrix(0, nrow = max(cluster.init), ncol = ncol(X))      # matrice K x V
  for (k in 1:max(cluster.init)) {
    if (length(which(cluster.init == k)) > 1) {
      centers[k, ] <- apply(X[which(cluster.init == k), ], 2, mean)
    } else{
      centers[k, ] <- X[which(cluster.init == k), ]
    }
  }
  rownames(centers) <- 1:max(cluster.init)
  weights.sauv <- NULL
  # calcul des poids initiaux
  if (weightsopt == T) {
    w = weightcalculate(X, mv, centers, cluster.init, gamma)
  } else{
    w = data.frame(weights = rep(1 / length(mv), length(mv)) ,
                   weightsmv = rep(1 / length(mv), sum(mv)))
  }
  weights.sauv <- w$weights
  
  #calcul des variances pondérées par classe
  varclass = rep(0, max(cluster.init))
  for (k in 1:max(cluster.init)) {
    I = which(cluster.init == k)
    for (v in 1:length(mv)) {
      dv = (ref[v] + 1):ref[v + 1]
      varclass[k] = varclass[k] + ((w$weights[v] ^ gamma) * sum((
        as.matrix(X)[I, dv] - matrix(
          rep(centers[k, dv], length(I)),
          nrow = length(I),
          byrow = T
        )
      ) ^ 2))
    }
    rm(I, dv)
  }
  # calcul du critere general
  CRIT <- c(CRIT, sum(varclass))
  
  ksplit.sauv <- NULL
  
  ###   début de la boucle de splitting
  iter = 1
  nbcluster = max(cluster.init)
  while (nbcluster < Kmax) {
    # cluster avec varclass max
    ksplit = which.max(varclass)
    I = which(clustersplithist[, iter] == ksplit)
    ksplit.sauv = c(ksplit.sauv, ksplit)
    
    # Kmeans dans la classe à splitter
    #pb avec les data dans le kmeans  -> la multiplication par les poids déforme trop
    if (testkmeans == T) {
      A <- NULL
      for (i in 1:length(I)) {
        #A<-rbind(A,(w$weightsmv^gamma)*X[I[i],])
        A <- rbind(A, (w$weightsmv ^ (gamma / 2)) * X[I[i], ])
      }
      a = kmeans(A ,
                 centers = 2,
                 nstart = 50,
                 iter.max = 50)   # (w$weightsmv^gamma)*X[I,] -> ne fait pas ce qu'il faut
    } else{
      a = kmeans(X[I, ] ,
                 centers = 2,
                 nstart = 50,
                 iter.max = 50)
    }
    
    J1 = which(a$cluster == 1)
    J2 = which(a$cluster == 2)
    # mise à jour des centres
    if (length(J1) > 1) {
      centers[ksplit, ] = apply(X[I[J1], ], 2, mean)
    } else{
      centers[ksplit, ] = X[I[J1], ]
    }
    if (length(J2) > 1) {
      centers = rbind(centers, apply(X[I[J2], ], 2, mean))
    } else{
      centers = rbind(centers, X[I[J2], ])
    }
    nbcluster = nbcluster + 1
    rownames(centers) = 1:nbcluster
    
    # historique clustering
    if (iter > 1)
      clustersplithist = cbind(clustersplithist, clustersplithist[, iter])
    clustersplithist[I[J1], iter + 1] <- ksplit
    clustersplithist[I[J2], iter + 1] <-
      nbcluster                     #length(varclass) +1
    #mise à jour des poids
    if (weightsopt == T) {
      w = weightcalculate(X, mv, centers, clustersplithist[, iter + 1], gamma)
    } else{
      w = data.frame(
        weights = rep(1 / length(mv), length(mv)) ,
        weightsmv = rep(1 / length(mv), sum(mv))
      )
    }
    weights.sauv = cbind(weights.sauv, w$weights)
    #mise à jour des variances pondérées par cluster
    varclass = rep(0, nbcluster)
    for (k in 1:nbcluster) {
      I = which(clustersplithist[, iter + 1] == k)
      for (v in 1:length(mv)) {
        dv = (ref[v] + 1):ref[v + 1]
        varclass[k] = varclass[k] + ((w$weights[v] ^ gamma) * sum((
          as.matrix(X)[I, dv] - matrix(
            rep(centers[k, dv], length(I)),
            nrow = length(I),
            byrow = T
          )
        ) ^ 2))
      }
      rm(I, dv)
    }
    
    # calcul du critere general
    CRIT <- c(CRIT, sum(varclass))
    
    #
    iter <- iter + 1
  }
  
  return(
    list(
      weights = weights.sauv,
      CRIT = CRIT,
      withinss = varclass,
      clustersplithist = clustersplithist,
      ksplit = ksplit.sauv
    )
  )
}

############################################################################################################
############################################################################################################
##   Splitting de clusters inspirée du multiview  -   poids par cluster
## Entrées :
#     cluster.init : une classif dure
#     X : les données (le 1er bloc correspond aux données principales utilisées pour obtenir classif)
#     mv : la taille de chaque multiview
#     weightsopt : if TRUE, multiview version avec poids par cluster; if FALSE weights = 1/V

#' Title
#'
#' @param X 
#' @param mv 
#' @param gamma 
#' @param Kmax 
#' @param cluster.init 
#'
#' @return
#' @export
#'
#' @examples
splittingClustersbis = function(X, mv, gamma = 1, Kmax, cluster.init) {
  ref = c(0, cumsum(mv))
  CRIT <- 0    #evolution critere general
  clustersplithist = matrix(c(cluster.init, cluster.init), ncol = 2)  #pour garder l'historique du splitting
  
  # initialisation des centres
  centers <-
    matrix(0, nrow = max(cluster.init), ncol = ncol(X))      # matrice K x V
  for (k in 1:max(cluster.init)) {
    if (length(which(cluster.init == k)) > 1) {
      centers[k, ] <- apply(X[which(cluster.init == k), ], 2, mean)
    } else{
      centers[k, ] <- X[which(cluster.init == k), ]
    }
  }
  rownames(centers) <- 1:max(cluster.init)
  weights.sauv <- NULL
  
  # calcul des poids initiaux
  w = weightcalculatebis(X, mv, centers, cluster.init, gamma)
  weights.sauv[[1]] = w$weights
  
  #calcul des variances pondérées par classe
  varclass = rep(0, max(cluster.init))
  for (k in 1:max(cluster.init)) {
    I = which(cluster.init == k)
    for (v in 1:length(mv)) {
      dv = (ref[v] + 1):ref[v + 1]
      varclass[k] = varclass[k] + ((w$weights[k, v] ^ gamma) * sum((
        as.matrix(X)[I, dv] - matrix(
          rep(centers[k, dv], length(I)),
          nrow = length(I),
          byrow = T
        )
      ) ^ 2))
    }
    rm(I, dv)
  }
  # calcul du critere general
  CRIT <- c(CRIT, sum(varclass))
  
  ksplit.sauv <- NULL
  
  ###   début de la boucle de splitting
  iter = 1
  nbcluster = max(cluster.init)
  while (nbcluster < Kmax) {
    # cluster avec varclass max
    ksplit = which.max(varclass)
    I = which(clustersplithist[, iter] == ksplit)
    ksplit.sauv = c(ksplit.sauv, ksplit)
    
    # kmeans sur la classe splittée
    A <- NULL
    for (i in 1:length(I)) {
      A <- rbind(A, (w$weightsmv[k, ] ^ (gamma / 2)) * X[I[i], ])
    }
    a = kmeans(A ,
               centers = 2,
               nstart = 50,
               iter.max = 50)
    
    J1 = which(a$cluster == 1)
    J2 = which(a$cluster == 2)
    # mise à jour des centres
    if (length(J1) > 1) {
      centers[ksplit, ] = apply(X[I[J1], ], 2, mean)
    } else{
      centers[ksplit, ] = X[I[J1], ]
    }
    if (length(J2) > 1) {
      centers = rbind(centers, apply(X[I[J2], ], 2, mean))
    } else{
      centers = rbind(centers, X[I[J2], ])
    }
    nbcluster = nbcluster + 1
    rownames(centers) = 1:nbcluster
    
    # historique clustering
    if (iter > 1)
      clustersplithist = cbind(clustersplithist, clustersplithist[, iter])
    clustersplithist[I[J1], iter + 1] <- ksplit
    clustersplithist[I[J2], iter + 1] <-
      nbcluster                     #length(varclass) +1
    #mise à jour des poids
    w = weightcalculatebis(X, mv, centers, clustersplithist[, iter + 1], gamma)
    
    weights.sauv[[iter + 1]] = w$weights
    #mise à jour des variances pondérées par cluster
    varclass = rep(0, nbcluster)
    for (k in 1:nbcluster) {
      I = which(clustersplithist[, iter + 1] == k)
      for (v in 1:length(mv)) {
        dv = (ref[v] + 1):ref[v + 1]
        varclass[k] = varclass[k] + ((w$weights[k, v] ^ gamma) * sum((
          as.matrix(X)[I, dv] - matrix(
            rep(centers[k, dv], length(I)),
            nrow = length(I),
            byrow = T
          )
        ) ^ 2))
      }
      rm(I, dv)
    }
    
    # calcul du critere general
    CRIT <- c(CRIT, sum(varclass))
    
    #
    iter <- iter + 1
  }
  
  return(
    list(
      weights = weights.sauv,
      CRIT = CRIT,
      withinss = varclass,
      clustersplithist = clustersplithist,
      ksplit = ksplit.sauv
    )
  )
}






##################################
# fonction calcul annexe des poids
#' Title
#'
#' @param X 
#' @param mv 
#' @param centers 
#' @param clustering 
#' @param gamma 
#'
#' @return
#' @export
#'
#' @examples
weightcalculatebis <- function(X, mv, centers, clustering, gamma) {
  ref = c(0, cumsum(mv))
  weights = matrix(0, nrow = nrow(centers), ncol = length(mv))
  weightsmv = matrix(0, nrow = nrow(centers), ncol = sum(mv))
  for (k in 1:nrow(centers)) {
    aux = NULL
    for (v in 1:length(mv)) {
      dv = (ref[v] + 1):ref[v + 1]
      aux = c(aux, sum((
        as.matrix(X)[which(clustering == k), dv] - matrix(
          rep(centers[k, dv], sum(clustering == k)),
          nrow = sum(clustering == k),
          byrow = T
        )
      ) ^ 2))
    }
    if (gamma > 1) {
      aux = aux ^ (1 / (1 - gamma))
      aux = aux / sum(aux)
    }
    if (gamma == 1) {
      aux = aux ^ (-1)
      aux = aux / sum(aux)
    }
    aux1 <- NULL
    for (j in 1:length(mv))
      aux1 <- c(aux1, rep(aux[j], mv[j]))
    weights[k, ] = aux
    weightsmv[k, ] = aux1
  }
  return(list(weights = weights, weightsmv = weightsmv))
}

