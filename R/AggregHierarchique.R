hmv1 <- function(X, mv, gamma, cluster.init, weightsopt = TRUE) {
  weights.sauv <- NULL
  Kmax <- max(cluster.init)
  # calcul des centres initiaux
  centers.init <- matrix(0, nrow = Kmax, ncol = ncol(X))
  for (k in seq_len(Kmax)) {
    centers.init[k, ] <-
      apply(X[which(cluster.init == k), , drop = FALSE], 2, mean)
    #if (length(which(cluster.init==k))>1){
    #centers.init[k,]<-apply(X[which(cluster.init==k),],2,mean)
    #}else{
    #centers.init[k,]<-X[which(cluster.init==k),]
    #}
  }
  rownames(centers.init) <- seq_len(Kmax)
  # calcul des poids initiaux
  if (weightsopt == TRUE) {
    w <- weightcalculate(X, mv, centers.init, cluster.init, gamma)
  } else{
    w <- data.frame(weights = rep(1 / length(mv), length(mv)) ,
                   weightsmv = rep(1 / length(mv), sum(mv)))
  }
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
  cluster <- -cluster.init      # on marque les clusters initiaux avec des "-" pour le hclust ensuite
  
  # calcul initial de la matrice des mesures d'aggregation D(Ck,Ck')
  D <- matrix(0, nrow = Kmax, ncol = Kmax)
  for (i in 1:(Kmax - 1)) {
    for (j in (i + 1):Kmax) {
      I <- which(cluster.init == i)
      J <- which(cluster.init == j)
      D[i, j] = D[j, i] = (length(I) * length(J) / (length(I) + length(J))) * 
        sum((w$weightsmv ^ gamma) * (centers.init[i, ] - centers.init[j, ]) ^ 2)
    }
  }
  rownames(D) <- (-1) * seq(1, Kmax, 1)
  colnames(D) <- (-1) * seq(1, Kmax, 1)
  
  CRIT <- Crit2(X, mv, gamma, w$weights, cluster.init)
  
  while(ncol(D) > 2) {
    # couple ayant valeur hors 0 minimale dans d
    a <- which(D == min(D[D > 0]), arr.ind = TRUE)
    aggreg <- rbind(aggreg, rownames(a))
    compt <- compt + 1
    noeud <-
      c(setdiff(noeud, as.numeric(aggreg[compt, ])), compt)  # les étiquettes restantes à aggréger
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
    D <- D[-sort(indice), -sort(indice)]
    
    # Dans le vecteur cluster, on indique la fusion
    cluster[which(cluster == aggreg[compt, 1])] <- compt
    cluster[which(cluster == aggreg[compt, 2])] <- compt
    
    # Calcul des centres
    centers <- matrix(0, nrow = length(noeud), ncol = sum(mv))
    for (k in 1:length(noeud)) {
      if (length(which(cluster == noeud[k])) > 1) {
        centers[k, ] <- apply(X[which(cluster == noeud[k]), ], 2, mean)
      } else{
        centers[k, ] <- X[which(cluster == noeud[k]), ]
      }
    }
    rownames(centers) <- noeud
    
    # Calcul des poids
    if (weightsopt == TRUE) {
      w <- weightcalculate(X, mv, centers, cluster, gamma)
    } else{
      w <- data.frame(
        weights <- rep(1 / length(mv), length(mv)) ,
        weightsmv <- rep(1 / length(mv), sum(mv))
      )
    }
    weights.sauv <- cbind(weights.sauv, w$weights)
    
    # Calcul de la matrice D
    for (i in 1:(nrow(D) - 1)) {
      for (j in (i + 1):nrow(D)) {
        I <- which(cluster == rownames(D)[i])
        J <- which(cluster == rownames(D)[j])
        D[i, j] <- D[j, i] <- (length(I) * length(J) / (length(I) + length(J))) * 
          sum((w$weightsmv ^ gamma) * (centers[i, ] - centers[j, ]) ^ 2)
      }
    }
    
    CRIT <- c(CRIT, Crit2(X, mv, gamma, w$weights, cluster))
    
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
weightcalculate <- function(X, mv, centers, cluster, gamma) {
  ref <- c(0, cumsum(mv))
  labcluster <- as.numeric(rownames(centers))
  aux <- matrix(0, nrow = length(labcluster), ncol = length(mv))
  for (k in 1:length(labcluster)) {
    # for each cluster k
    I <- which(cluster == labcluster[k])
    if (length(I) > 0) {
      for (v in 1:length(mv)) {
        dv <- (ref[v] + 1):ref[v + 1]
        aux[k, v] <- sum((as.matrix(X)[I, dv] - matrix(
          rep(centers[k, dv], length(I)), nrow = length(I), byrow = T)) ^ 2)
      }
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


Crit2 <- function(X, mv, gamma, weights, cluster) {
  ref <- c(0, cumsum(mv))
  clustername <- unique(cluster)
  nbcluster <- length(clustername)
  centers <- matrix(0, nrow = nbcluster, ncol = ncol(X))
  for (k in seq_len(nbcluster))
    centers[k, ] <- apply(X[which(cluster == clustername[k]), , drop = FALSE], 2, mean)
  
  varclass <- rep(0, nbcluster)
  for (k in seq_len(nbcluster)) {
    I <- which(cluster == clustername[k])
    for (v in 1:length(mv)) {
      dv <- (ref[v] + 1):ref[v + 1]
      varclass[k] <- varclass[k] + ((weights[v] ^ gamma) * sum((
        X[I, dv, drop = FALSE] - 
          matrix(rep(centers[k, dv], length(I), drop = FALSE),nrow = length(I),byrow = T)) ^ 2))
    }
#    rm(I, dv)
  }
  return(sum(varclass))
}

cutreeNew <- function(R, K, cluster.init) {
  cc <- cutree(R, K)
  classif <- cluster.init
  for (j in 1:length(cc))
    classif[which(classif == names(cc)[j])] <- cc[j]
  return(classif)
}


