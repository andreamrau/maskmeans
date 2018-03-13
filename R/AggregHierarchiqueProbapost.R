##   Agrégation de clusters inspirée du multiview avec poids sur indiv ("probapost").
## Entrées :
#  probapost.init : une classif fuzzy (probapost)  n x K
#  X : les données (le 1er bloc correspond aux données principales utilisées pour obtenir classif)
# mv : la taille de chaque multiview

#' Aggregation of fuzzy clusters based on multi-view inputs
#'
#' @param probapost.init Matrix of available conditional probabilities of cluster membership
#' for \code{X}
#' @inheritParams hmv1
#'
#' @return
#' @export
#'
#' @examples
#' 
hmvprobapost <- function(X, mv, gamma, probapost.init)
{
  weights.sauv <- NULL
  Kmax <- ncol(probapost.init)
  
  # calcul des centres initiaux
  centers.init <- matrix(0, nrow = Kmax, ncol = ncol(X))
  for (k in 1:Kmax) {
    centers.init[k,] <-
      apply(matrix(rep(probapost.init[, k], sum(mv)), nrow = nrow(X)) * X, 2, sum) / 
      sum(probapost.init[, k])
  }
  rownames(centers.init) <- 1:Kmax
  # calcul des poids initiaux
  w <- weightcalculateprobapost(X, mv, centers.init, probapost.init, gamma)
  
  weights.sauv <- w$weights
  
  # pour avoir une structure de hclust
  R <- hclust(dist(centers.init))
  R$method <- "AgregMultiv"
  R$call <- ""
  
  compt <- 0
  aggreg <- NULL
  height <- NULL
  noeud <-
    (-1) * seq(1, Kmax, 1)     # les étiquettes initiales sont -1,-2, ... -Kmax dans hclust
  cluster <- -(1:Kmax)     # on marque les clusters initiaux avec des "-" pour le hclust ensuite
  probapost <- probapost.init
  
  # calcul initial de la matrice des mesures d'aggregation D(Ck,Ck')
  D <- matrix(0, nrow = Kmax, ncol = Kmax)
  for (i in 1:(Kmax - 1)) {
    for (j in (i + 1):Kmax) {
      ni <- sum(probapost.init[, i])
      nj <- sum(probapost.init[, j])
      D[i, j] = D[j, i] = (ni * nj / (ni + nj)) * sum((w$weightsmv ^ gamma) *
                                                        (centers.init[i,] - centers.init[j,]) ^ 2)
    }
  }
  rownames(D) <- (-1) * seq(1, Kmax, 1)
  colnames(D) <- (-1) * seq(1, Kmax, 1)
  
  CRIT <- Crit3(X, mv, gamma, w$weights, probapost, cluster)
  
  while (ncol(D) > 2) {
    # couple ayant valeur hors 0 minimale dans d
    a <- which(D == min(D[D > 0]), arr.ind = TRUE)
    aggreg <- rbind(aggreg, rownames(a))
    compt <- compt + 1
    noeud <-
      c(setdiff(noeud, as.numeric(aggreg[compt,])), compt)  # les étiquettes restantes à aggréger
    height <- c(height, D[aggreg[compt, 1], aggreg[compt, 2]])
    
    #indice des deux classes à aggreger
    indice <-
      c(which(colnames(D) == aggreg[compt, 1]), which(colnames(D) == aggreg[compt, 2]))
    
    #creation d'une nvelle ligne et colonne de fusion
    D <- rbind(D, t(rep(0, ncol(D))))
    D <- cbind(D, rep(0, nrow(D)))
    rownames(D)[nrow(D)] <- compt
    colnames(D)[nrow(D)] <- compt
    
    # On retire les deux lignes et colonnes des classes fusionnées
    D <- D[-sort(indice),-sort(indice)]
    
    # Dans le vecteur cluster, on indique la fusion
    cluster[which(cluster == aggreg[compt, 1])] <- compt
    cluster[which(cluster == aggreg[compt, 2])] <- compt
    
    # probapost
    probapost <- cbind(probapost, rep(0, nrow(probapost)))
    colnames(probapost)[ncol(probapost)] <- compt
    probapost[, ncol(probapost)] <-
      apply(probapost[, indice], 1, sum)
    probapost <- probapost[,-sort(indice)]
    
    # Calcul des centres
    centers <- matrix(0, nrow = length(noeud), ncol = sum(mv))
    for (k in 1:length(noeud)) {
      centers[k,] <-
        apply(matrix(rep(probapost[, k], sum(mv)), nrow = nrow(X)) * X, 2, sum) / sum(probapost[, k])
    }
    rownames(centers) <- noeud
    
    # Calcul des poids
    w <- weightcalculateprobapost(X, mv, centers, probapost, gamma)
    weights.sauv <- cbind(weights.sauv, w$weights)
    
    # Calcul de la matrice D
    for (i in 1:(nrow(D) - 1)) {
      for (j in (i + 1):nrow(D)) {
        ni <- sum(probapost[, i])
        nj <- sum(probapost[, j])
        D[i, j] <- D[j, i] <- (ni * nj / (ni + nj)) * sum((w$weightsmv ^
                                                           gamma) * (centers[i,] - centers[j,]) ^ 2)
      }
    }
    
    CRIT <-
      c(CRIT, Crit3(X, mv, gamma, w$weights, probapost, cluster))
    
  }# fin du while
  
  a <- which(D == min(D[D > 0]), arr.ind = TRUE)
  aggreg <- rbind(aggreg, rownames(a))
  aggreg <- t(apply(aggreg, 1, function(z)
    as.numeric(z)))
  compt <- compt + 1
  height <- c(height, D[1, 2])
  ord <- (-1) * t(aggreg)[which(t(aggreg) < 0)]
  labels <- as.character(sort((-1) * aggreg[which(aggreg < 0)]))
  
  R$merge <- aggreg
  R$height <- height
  R$order <- ord
  R$labels <- labels
  
  return(
    list(
      merge = aggreg,
      height = height,
      order = ord,
      labels = labels,
      R = R,
      weights = weights.sauv,
      CRIT = CRIT
    )
  )
}

#' Calculation of multi-view weights for hard clusters
#'
#' @inheritParams hmv1
#' @param centers Matrix of cluster centers
#' @param probapost Matrix of conditional probabilities of cluster membership for \code{X}
#'
#' @return
#' @export
#'
#' @examples
weightcalculateprobapost <-
  function(X, mv, centers, probapost, gamma) {
    ref <- c(0, cumsum(mv))
    labcluster <- as.numeric(rownames(centers))
    aux <- matrix(0, nrow = length(labcluster), ncol = length(mv))
    for (k in 1:length(labcluster)) {
      # for each cluster k
      for (v in 1:length(mv)) {
        dv <- (ref[v] + 1):ref[v + 1]
        aux[k, v] <- sum(probapost[, k] * rowSums((sweep(
          X[, dv, drop = FALSE], 2, centers[k, dv], FUN = "-")) ^ 2))
      }
    }
    aux1 <- apply(aux, 2, sum)
    
    if (gamma > 1) {
      aux1 <- aux1 ^ (1 / (1 - gamma))
      aux1 <- aux1 / sum(aux1)
    }
    if (gamma == 1) {
      v0 <- which.min(aux1)
      aux1[v0] <- 1
      aux1[-v0] <- 0
    }
    
    weightsmv <- NULL
    for (j in 1:length(mv))
      weightsmv <- c(weightsmv, rep(aux1[j], mv[j]))
    return(list(weights = aux1, weightsmv = weightsmv))
  }


#' Calculate the value of criterion 3
#'
#' @inheritParams hmv1
#' @param weights Blah blah
#' @param probapost Blah blah
#' @param cluster Blah blah
#'
#' @return
#' @export
#'
#' @examples
Crit3 <- function(X, mv, gamma, weights, probapost, cluster) {
  ref <- c(0, cumsum(mv))
  clustername <- unique(cluster)
  nbcluster <- length(clustername)
  centers <- matrix(0, nrow = nbcluster, ncol = ncol(X))
  for (k in 1:nbcluster)
    centers[k,] <-
    apply(matrix(rep(probapost[, k], sum(mv)), nrow = nrow(X)) * X, 2, sum) / sum(probapost[, k])
  
  varclass <- rep(0, nbcluster)
  for (k in 1:nbcluster) {
    for (v in 1:length(mv)) {
      dv <- (ref[v] + 1):ref[v + 1]
      varclass[k] <- varclass[k] + ((weights[v] ^ gamma) * sum(probapost[, k] *
                                                                apply((
                                                                  X[, dv, drop = FALSE] - matrix(
                                                                    rep(centers[k, dv], nrow(X), drop = FALSE),
                                                                    nrow = nrow(X),
                                                                    byrow = T
                                                                  )
                                                                ) ^ 2, 1, sum)))
    }
  }
  return(sum(varclass))
}





#' Title
#'
#' @param R ??
#' @param K ??
#' @param probapost.init  ??
#'
#' @return
#' @export
#'
#' @examples
cutreeNewProbapost <- function(R, K, probapost.init) {
  cluster.init <- apply(probapost.init, 1, which.max)
  cc <- cutree(R, K)
  classif <- cluster.init
  probapost <- matrix(0, nrow = nrow(probapost.init), ncol = K)
  for (k in 1:length(cc))
    classif[which(classif == names(cc)[k])] <- cc[k]
  
  for (k in 1:K) {
    if (length(which(cc == k)) > 1) {
      probapost[, k] = apply(probapost.init[, which(cc == k)], 1, sum)
    } else{
      probapost[, k] = probapost.init[, which(cc == k)]
    }
  }
  return(list(classif = classif, probapost = probapost))
}
