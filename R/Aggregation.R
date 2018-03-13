#' Aggregation of hard clusters based on multi-view input
#'
#' @param X Matrix of multi-view data, where the first view corresponds to the 
#' principal data used to obtain the partition in \code{cluster.init}
#' @param mv Vector corresponding to the size of each data view. The sum of \code{mv} should 
#' correspond to the number of columns in \code{X}.
#' @param gamma Parameter that controls the distribution of view weights. Default value is 2. 
#' @param cluster_init Vector of available cluster labels 
#' @param use_mv_weights If \code{TRUE}, run algorithm in weighted multi-view mode; if FALSE, the
#' weight for each view is set to be equal.
#'
#' @return
#' \item{merged_clusters }{Matrix providing each pair of merged clusters at each iteration of the algorithm}
#' \item{hclust }{Object of class "hclust" to be used for plotting cluster agglomerations}
#' \item{weights }{Matrix of dimension \code{v} x \code{niterations}, where \code{v} is the number of
#' views and \code{niterations} is the number of successive agglomerations, providing the multi-view weights}
#' \item{criterion }{Value taken on by the agglomerative criterion at each iteration}
#' @export
#'
#' @examples

mv_aggregation <- function(X, mv, cluster_init, gamma=2, use_mv_weights = TRUE) {
  if(gamma <= 1) stop("gamma must be greater than 1.")
  weights_save <- NULL
  Kmax <- max(cluster_init)
  
  # Calculate initial cluster centers
  centers_init <- rowsum(X, group=cluster_init) / table(cluster_init)
  ## Calculate initial weights
  if(use_mv_weights) {
    w <- mv_weights(X=X, mv=mv, centers=as.matrix(centers_init), 
                    cluster=cluster_init, gamma=gamma)
  } else{
    w <- data.frame(weights = rep(1 / length(mv), length(mv)),
                    weightsmv = rep(1 / length(mv), sum(mv)))
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
  cluster <- -cluster_init      # on marque les clusters initiaux avec des "-" pour le hclust ensuite
  
  # Initial calculation of aggregation matrices D(Ck,Ck'): Possibly to be optimized? TODO 
  D <- matrix(0, nrow = Kmax, ncol = Kmax)
  for (i in 1:(Kmax - 1)) {
    for (j in (i + 1):Kmax) {
      I <- which(cluster.init == i)
      J <- which(cluster.init == j)
      D[i, j] = D[j, i] = (length(I) * length(J) / (length(I) + length(J))) * 
        sum((w$weightsmv ^ gamma) * (centers_init[i, ] - centers_init[j, ]) ^ 2)
    }
  }
  rownames(D) <- (-1) * seq(1, Kmax, 1)
  colnames(D) <- (-1) * seq(1, Kmax, 1)
  ## Calculate criterion 
  CRIT <- criterion(X=X, mv=mv, gamma=gamma, weights=w$weights, cluster=cluster_init)
  
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
    
    # Recalculate the centers
    centers <- rowsum(X, group=cluster) / table(cluster)
    
    # Calcul des poids
    if (use_mv_weights) {
      w <- mv_weights(X=X, mv=mv, centers=as.matrix(centers), 
                      cluster=cluster, gamma=gamma)
    } else{
      w <- data.frame(
        weights <- rep(1 / length(mv), length(mv)) ,
        weightsmv <- rep(1 / length(mv), sum(mv))
      )
    }
    weights_save <- cbind(weights_save, w$weights)
    
    # Calcul de la matrice D
    for (i in 1:(nrow(D) - 1)) {
      for (j in (i + 1):nrow(D)) {
        I <- which(cluster == rownames(D)[i])
        J <- which(cluster == rownames(D)[j])
        D[i, j] <- D[j, i] <- (length(I) * length(J) / (length(I) + length(J))) * 
          sum((w$weightsmv ^ gamma) * (centers[i, ] - centers[j, ]) ^ 2)
      }
    }
    CRIT <- c(CRIT, criterion(X=X, mv=mv, gamma=gamma, weights=w$weights, cluster=cluster))
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
  
  return(
    list(
      merged_clusters = aggreg,
#     height = height, order = ord, labels = labels,
      hclust = R,
      weights = weights_save,
      criterion = CRIT
    )
  )
}


## NOT exported: function for calculating multi-view weights

mv_weights <- function(X, mv, centers, cluster, gamma) {
  ref <- c(0, cumsum(mv))
  labcluster <- as.numeric(rownames(centers))

  tmp <- centers[match(cluster, rownames(centers)),]
  aux1 <- as.numeric(rowsum(colSums(rowsum((X - tmp)^2, group=cluster)), 
                            group=rep(LETTERS[1:length(mv)], times=mv)))
  
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

criterion <- function(X, mv, gamma, weights, cluster) {
  ref <- c(0, cumsum(mv))
  clustername <- unique(cluster)
  nbcluster <- length(clustername)
  centers <- as.matrix(rowsum(X, group=cluster) / table(cluster))
  
  tmp <- centers[match(cluster, rownames(centers)),]
  varclass <- as.matrix(rowsum(t(rowsum((X-tmp)^2, group=cluster)), group=rep(LETTERS[1:length(mv)], times=mv)))
  varclass_weighted <- sum(varclass * (weights^gamma) )
  return(varclass_weighted) 
}

## NOT exported: function to scale and normalize by view size

scaleview <- function(Data, mv){
  X <- scale(Data)
  ref <- cumsum(c(0,mv))
  for (v in 1:length(mv))
    X[,seq((ref[v]+1),ref[v+1])] = X[,seq((ref[v]+1),ref[v+1])] / mv[v]
  return(X)
}

## NOT exported: function to cut tree

cutreeNew <- function(hclust_obj, K, cluster) {
  R <- hclust_obj
  cc <- cutree(R, K)
  classif <- cluster
  for (j in 1:length(cc))
    classif[which(classif == names(cc)[j])] <- cc[j]
  return(classif)
}