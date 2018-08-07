#' Aggregation of hard or fuzzy clusters based on multi-view data
#'
#' @param X Matrix of multi-view data, where the first view corresponds to the 
#' principal data used to obtain the partition or fuzzy clustering in \code{cluster_init}
#' @param mv (Optional unless \code{X} is a matrix.) If \code{X} is a matrix, vector 
#' corresponding to the size of each data view. 
#' The sum of \code{mv} should correspond to the number of columns in \code{X}.
#' @param gamma Parameter that controls the distribution of view weights. Default value is 2. 
#' @param clustering_init Either a vector of available cluster labels (for hard clustering) or a matrix
#' of fuzzy classification labels (For fuzzy clustering)
#' @param use_mv_weights If \code{TRUE}, run algorithm in weighted multi-view mode; if FALSE, the
#' weight for each view is set to be equal. This option is only used for hard clustering.
#' @param verbose If \code{TRUE}, provide verbose output
#'
#' @return
#' \item{merged_clusters }{Matrix providing each pair of merged clusters at each iteration of the algorithm}
#' \item{hclust }{Object of class "hclust" to be used for plotting cluster agglomerations}
#' \item{weights }{Matrix of dimension \code{v} x \code{niterations}, where \code{v} is the number of
#' views and \code{niterations} is the number of successive agglomerations, providing the multi-view weights}
#' \item{criterion }{Value taken on by the agglomerative criterion at each iteration}
#' @export
#'
mv_aggregation <- function(X, mv, clustering_init, gamma=2, use_mv_weights = TRUE, verbose=TRUE) {
  
  if(gamma < 1) stop("gamma must be greater than or equal to 1.")
  mode <- ifelse(is.vector(clustering_init), "hard", "fuzzy")
  cat("Running in mode:", mode, "\n")

  weights_save <- NULL
  cluster_init <- clustering_init
  Kmax <- ifelse(mode == "hard", max(cluster_init), ncol(cluster_init))
  
  if(mode == "hard") {
    # Calculate initial cluster centers
    centers_init <- rowsum(X, group=cluster_init) / as.vector(table(cluster_init))
    ## Calculate initial weights
    if(use_mv_weights) {
      w <- mv_weights(X=X, mv=mv, centers=as.matrix(centers_init), 
                      cluster=cluster_init, gamma=gamma, mode=mode)
    } else{
      w <- data.frame(weights = rep(1 / length(mv), length(mv)),
                      weightsmv = rep(1 / length(mv), sum(mv)))
    }
  }
  
  if(mode == "fuzzy") {
    probapost_init <- probapost <- cluster_init
    # Calculate initial cluster centers
    centers_init <- matrix(0, nrow = Kmax, ncol = ncol(X))
    for (k in 1:Kmax) {
      centers_init[k,] <-
        apply(matrix(rep(probapost_init[, k], sum(mv)), nrow = nrow(X)) * X, 2, sum) / 
        sum(probapost_init[, k])
    }
    rownames(centers_init) <- 1:Kmax
    # Calculate initial weights
    w <- mv_weights(X=X, mv=mv, centers=centers_init, cluster=probapost_init, gamma=gamma, mode=mode)
  }

  weights_save <- w$weights

  # Format as hclust object
  R <- hclust(dist(centers_init))
  R$method <- "AggregMultiv"
  R$call <- ""
  compt <- 0
  aggreg <- height <- NULL
  noeud <-
    (-1) * seq(1, Kmax, 1)     # les étiquettes initiales sont -1,-2, ... -Kmax dans hclust
  if(mode == "hard") cluster <- -cluster_init      # on marque les clusters initiaux avec des "-" pour le hclust ensuite
  if(mode == "fuzzy") cluster <- -(1:Kmax)

  # Initial calculation of aggregation matrices D(Ck,Ck'): Possibly to be optimized? TODO 
  D <- matrix(0, nrow = Kmax, ncol = Kmax)
  for (i in 1:(Kmax - 1)) {
    for (j in (i + 1):Kmax) {
      if(mode == "hard") {
        I <- which(cluster_init == i)
        J <- which(cluster_init == j)
        ni <- length(I)
        nj <- length(J)
      }
      if(mode == "fuzzy") {
        ni <- sum(probapost_init[, i])
        nj <- sum(probapost_init[, j])
      }
      D[i, j] = D[j, i] = (ni * nj / (ni + nj)) * sum((w$weightsmv ^ gamma) *
                                                        (centers_init[i,] - centers_init[j,]) ^ 2)
    }
  }
  rownames(D) <- (-1) * seq(1, Kmax, 1)
  colnames(D) <- (-1) * seq(1, Kmax, 1)
  
  ## Calculate criterion
  CRIT <- ifelse(mode == "hard",
                 criterion(X=X, mv=mv, gamma=gamma, weights=w$weights, cluster=cluster_init, mode=mode),
                 criterion(X=X, mv=mv, gamma=gamma, weights=w$weights, cluster=cluster, 
                           mode=mode, probapost=probapost_init))
  
  while(ncol(D) > 2) {
    # Find pair of clusters with smallest non-zero value in D
    a <- which(D == min(D[D > 0]), arr.ind = TRUE)
    aggreg <- rbind(aggreg, rownames(a))
    compt <- compt + 1
    noeud <-
      c(setdiff(noeud, as.numeric(aggreg[compt, ])), compt)  ## Labels remaining to be aggregated
    height <- c(height, D[aggreg[compt, 1], aggreg[compt, 2]])
    
    # Index of two classes to aggregate
    indice <-
      c(which(colnames(D) == aggreg[compt, 1]), which(colnames(D) == aggreg[compt, 2]))
    if(verbose) {
      cat("  => Aggregating clusters", aggreg[compt, 1], "and", aggreg[compt, 2], "\n")
    }
    
    # Creation of a new row/column for fusion column
    D <- rbind(D, t(rep(0, ncol(D))))
    D <- cbind(D, rep(0, nrow(D)))
    rownames(D)[nrow(D)] <- compt
    colnames(D)[nrow(D)] <- compt
    
    # Delete the rows and columns for the fused clusters
    D <- D[-sort(indice), -sort(indice)]
    
    # In the cluster vector, indicate the fused cluster
    cluster[which(cluster == aggreg[compt, 1])] <- compt
    cluster[which(cluster == aggreg[compt, 2])] <- compt
    
    if(mode == "fuzzy") {
      
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
      w <- mv_weights(X=X, mv=mv, centers=centers, cluster=probapost, gamma=gamma, mode=mode)
    }
  
    if(mode == "hard") {
      # Recalculate the centers
      centers <- rowsum(X, group=cluster) / as.vector(table(cluster))
      ## Make sure the order of the rows does not change
      centers <- centers[match(rownames(centers), noeud),]
      
      # Calcul des poids
      if (use_mv_weights) {
        w <- mv_weights(X=X, mv=mv, centers=as.matrix(centers), 
                        cluster=cluster, gamma=gamma, mode=mode)
      } else{
        w <- data.frame(
          weights <- rep(1 / length(mv), length(mv)) ,
          weightsmv <- rep(1 / length(mv), sum(mv))
        )
      }
    }
    
    weights_save <- cbind(weights_save, w$weights)
    
    # Calculate the D matrix
    for (i in 1:(nrow(D) - 1)) {
      for (j in (i + 1):nrow(D)) {
        if(mode == "hard") {
          I <- which(cluster == rownames(D)[i])
          J <- which(cluster == rownames(D)[j]) 
          ni <- length(I)
          nj <- length(J)
        }
        if(mode == "fuzzy") {
          ni <- sum(probapost[, i])
          nj <- sum(probapost[, j])
        }
        D[i, j] <- D[j, i] <- (ni * nj / (ni + nj)) * 
          sum((w$weightsmv ^ gamma) * (centers[i, ] - centers[j, ]) ^ 2)
      }
    }
    
#    CRIT <- c(CRIT, criterion(X=X, mv=mv, gamma=gamma, weights=w$weights, cluster=cluster))
    CRIT <- c(CRIT, ifelse(mode == "hard",
                   criterion(X=X, mv=mv, gamma=gamma, weights=w$weights, cluster=cluster, mode=mode),
                   criterion(X=X, mv=mv, gamma=gamma, weights=w$weights, cluster=cluster, 
                             mode=mode, probapost=probapost)))
    
  } ## End of while loop
  
  a <- which(D == min(D[D > 0]), arr.ind = TRUE)
  aggreg <- rbind(aggreg, rownames(a))
  aggreg <- t(apply(aggreg, 1, function(z) as.numeric(z)))
  compt <- compt + 1
  height <- c(height, D[1, 2])
  ord <- (-1) * t(aggreg)[which(t(aggreg) < 0)]
  labels <- as.character(sort((-1) * aggreg[which(aggreg < 0)]))
  
  R$merge <- aggreg
  R$height <- height
  R$order <- ord
  R$labels <- labels
  
  res <- list(
    merged_clusters = aggreg,
    #     height = height, order = ord, labels = labels,
    hclust = R,
    weights = weights_save,
    criterion = CRIT
  )
  class(res) <- "maskmeans"
  return(res)
}


## NOT exported: function for calculating multi-view weights

mv_weights <- function(X, mv, centers, cluster, gamma, mode, delta=NULL) {
  ref <- c(0, cumsum(mv))
  labcluster <- as.numeric(rownames(centers))
  if(mode == "fuzzy" & !is.matrix(cluster)) stop("Fuzzy clustering requires a matrix of posterior probabilities.")
  if(mode == "hard" & !is.vector(cluster)) stop("Hard clustering requires a vector of cluster assignments.")
  
  if(mode == "hard") {
    tmp <- centers[match(cluster, rownames(centers)),]
    if(length(unique(cluster)) > 1) {
      aux1 <- as.numeric(rowsum(colSums(rowsum((X - tmp)^2, group=cluster)), 
                                group=rep(LETTERS[1:length(mv)], times=mv)))
    } else {
      aux1 <- as.numeric(rowsum(colSums((X-tmp)^2), group=rep(LETTERS[1:length(mv)], 
                                                              times=mv)))
    }
  }
  if(mode == "fuzzy") {
    aux <- matrix(0, nrow = length(labcluster), ncol = length(mv))
    for (k in 1:length(labcluster)) {
      # for each cluster k
      for (v in 1:length(mv)) {
        dv <- (ref[v] + 1):ref[v + 1]
        delta <- ifelse(is.null(delta), 1, delta)  ## Include delta if provided
        aux[k, v] <- sum(cluster[,k]^delta * rowSums((sweep(
          X[, dv, drop = FALSE], 2, centers[k, dv], FUN = "-")) ^ 2))
      }
    }
    aux1 <- apply(aux, 2, sum)
  }
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

## NOT exported: function to calculate criteria for agglomeration

criterion <- function(X, mv, gamma, weights, cluster, probapost=NULL, mode="hard") {
  if(mode == "fuzzy" & missing(probapost)) stop("Fuzzy clustering requires probapost argument.")
  ref <- c(0, cumsum(mv))
  clustername <- unique(cluster)
  nbcluster <- length(clustername)
  if(mode == "hard") {
    centers <- as.matrix(rowsum(X, group=cluster) / as.vector(table(cluster)))
    tmp <- centers[match(cluster, rownames(centers)),]
    varclass <- as.matrix(rowsum(t(rowsum((X-tmp)^2, group=cluster)), group=rep(LETTERS[1:length(mv)], times=mv)))
    varclass_weighted <- sum(varclass * (weights^gamma) )
  }
  if(mode == "fuzzy") {
    centers <- matrix(0, nrow = nbcluster, ncol = ncol(X))
    for (k in 1:nbcluster)
      centers[k,] <-
        apply(matrix(rep(probapost[, k], sum(mv)), nrow = nrow(X)) * X, 2, sum) / sum(probapost[, k])
    varclass <- rep(0, nbcluster)
    for (k in 1:nbcluster) {
      for (v in 1:length(mv)) {
        dv <- (ref[v] + 1):ref[v + 1]
        varclass[k] <- varclass[k] + ((weights[v] ^ gamma) * 
                                        sum(probapost[, k] * 
                                              rowSums((X[, dv, drop = FALSE] - 
                                                         matrix(rep(centers[k, dv], nrow(X), drop = FALSE),
                                                                nrow = nrow(X), byrow = T)) ^ 2)))
      }
    }
    varclass_weighted <- sum(varclass)
  }
  return(varclass_weighted) 
}

#' Cut an aggregation tree from maskmeans output
#' 
#' Cut an aggregation tree into several groups after running \code{maskmeans_aggregation} or
#' \code{maskmeans(..., type = "aggregation")}.
#'
#' @param obj An object of class \code{"maskmeans"}
#' @param K An integer scalar with the desired number of clusters
#' @param clustering_init Initial partition or posterior probabilities used to cluster the multi-view data
#'
#' @return
#' \item{classif }{Vector of group memberships after cutting tree}
#' \item{probapost }{Matrix of group posterior probabilities after cutting tree}
#' @export
maskmeans_cutree <- function(obj, K, clustering_init) {
  hclust_obj <- obj$hclust
  cluster <- clustering_init
  if(is.vector(cluster)) mode <- "hard"
  if(!is.vector(cluster)) mode <- "fuzzy"
  R <- hclust_obj
  cc <- cutree(R, K)
  if(mode == "fuzzy") {
    probapost_init <- cluster
    classif <- apply(probapost_init, 1, which.max)
  } else {
    classif <- cluster
  }
  for (j in 1:length(cc))
    classif[which(classif == names(cc)[j])] <- cc[j]
  
  if(mode == "hard") return(list(classif = classif))
  if(mode == "fuzzy") {
    probapost <- matrix(0, nrow = nrow(probapost_init), ncol = K)
    for (k in 1:K) {
      if (length(which(cc == k)) > 1) {
        probapost[, k] = apply(probapost_init[, which(cc == k)], 1, sum)
      } else{
        probapost[, k] = probapost_init[, which(cc == k)]
      }
    }
    return(list(classif = classif, probapost = probapost))
  }
}

