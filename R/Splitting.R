#' Splitting of hard or fuzzy clusters based on multi-view data
#'
#' @param X Matrix of multi-view data, where the first view corresponds to the 
#' principal data used to obtain the partition or fuzzy clustering in \code{cluster_init}
#' @param mv (Optional unless \code{X} is a matrix.) If \code{X} is a matrix, vector 
#' corresponding to the size of each data view. 
#' The sum of \code{mv} should correspond to the number of columns in \code{X}.
#' @param Kmax Maximum number of clusters for splitting
#' @param gamma Parameter that controls the distribution of view weights. Default value is 2. 
#' @param clustering_init Either a vector of available cluster labels (for hard clustering) 
#' or a matrix
#' of fuzzy classification labels (For fuzzy clustering)
#' @param use_mv_weights If \code{TRUE}, run algorithm in weighted multi-view mode; 
#' if FALSE, the
#' weight for each view is set to be equal. This option is only used for hard clustering.
#' @param perCluster_mv_weights If \code{TRUE}, use cluster-specific multi-view weights.
#' Otherwise use classic multi-view weights.
#'
#' @return
#' \item{split_clusters }{Matrix providing the history of each cluster splitting at 
#' each iteration of the algorithm}
#' \item{weights }{Matrix of dimension \code{v} x \code{niterations}, where \code{v} 
#' is the number of
#' views and \code{niterations} is the number of successive splits, providing the 
#' multi-view weights}
#' \item{criterion }{Value taken on by the splitting criterion at each iteration}
#' \item{withnss }{The within sum-of-squares for each cluster at the last iteration}
#' \item{ksplit }{Vector identifying which cluster was split at each iteration of the 
#' algorithm}
#' 
#' @export
mv_splitting <- function(X, mv, clustering_init, Kmax, gamma=2, 
                         use_mv_weights = TRUE, 
                         perCluster_mv_weights = TRUE) {
  
  ## TODO fuzzy splitting not added yet
  cluster_init <- clustering_init
  
  if(gamma < 1) stop("gamma must be greater than or equal to 1.")
  mode <- ifelse(is.vector(clustering_init), "hard", "fuzzy")
  if(mode == "fuzzy") stop("Only hard clustering is currently supported.")
  cat("Running in mode:", mode, "\n")
  
  ref <- c(0, cumsum(mv))
  CRIT <- 0    ## Evolution of general criterion
  clustersplithist <- matrix(c(cluster_init, cluster_init), ncol = 2)  ## Keep track of splits
  
  # Calculate initial cluster centers
  centers_init <- rowsum(X, group=cluster_init) / as.numeric(table(cluster_init))

  ## Calculate initial weights
  if(!perCluster_mv_weights) {
    if(use_mv_weights) {
      w <- mv_weights(X=X, mv=mv, centers=as.matrix(centers_init), 
                      cluster=cluster_init, gamma=gamma, mode=mode)
    } else{
      w <- data.frame(weights = rep(1 / length(mv), length(mv)),
                      weightsmv = rep(1 / length(mv), sum(mv)))
    }
    weights_save <- w$weights    
  } else {
    w <- mv_weights_perCluster(X=X, mv=mv, centers=as.matrix(centers_init), 
                               cluster=cluster_init, gamma=gamma,
                               mode = mode)
    weights_save <- c()
    weights_save[[1]] <- w$weights
  }

  
  ## Calculate weighted per-cluster variances
  if(!perCluster_mv_weights) {
    tmp <- centers_init[match(cluster_init, rownames(centers_init)),]
    varclass <- (w$weights ^ gamma) * 
      rowsum(colSums(rowsum((X - tmp)^2, group=cluster_init)), 
             group=rep(LETTERS[1:length(mv)], times=mv))
  } else {
    # TODO
     varclass <- rep(0, max(cluster_init))
     for (k in 1:max(cluster_init)) {
       I <- which(cluster_init == k)
       for (v in 1:length(mv)) {
         dv <- (ref[v] + 1):ref[v + 1]
         varclass[k] <- varclass[k] + ((w$weights[k, v] ^ gamma) * sum((
           as.matrix(X)[I, dv] - matrix(
             rep(centers_init[k, dv], length(I)),
             nrow = length(I),
             byrow = T)) ^ 2))
       }
     }
   }
  
  # Calculate general criterion
  CRIT <- c(CRIT, sum(varclass)) 
  ksplit_save <- NULL
  
  ## Start splitting loop
  iter <- 1
  nbcluster <- max(cluster_init)
  centers <- as.data.frame(centers_init)
  while(nbcluster < Kmax) {
    
    # Cluster with max varclass
    ksplit <- which.max(varclass)
    I <- which(clustersplithist[, iter] == ksplit)
    ksplit_save <- c(ksplit_save, ksplit)
    
    # Kmeans in the cluster to be split
    #pb avec les data dans le kmeans  -> la multiplication par les poids déforme trop
    if(!perCluster_mv_weights) {
      A <- t(w$weightsmv ^ (gamma/2) * t(X[I,])) 
    } else {
      A <- NULL
      for (i in 1:length(I)) {
        A <- rbind(A, (w$weightsmv[k, ] ^ (gamma / 2)) * X[I[i], ])
      }
    }
    a <- kmeans(A, centers = 2, nstart = 50, iter.max = 50)        
    
    J1 <- which(a$cluster == 1)
    J2 <- which(a$cluster == 2)
    
    # Update centers
    centers[c(ksplit, nrow(centers)+1),] <- 
      rowsum(X[I,], group=a$cluster) / as.numeric(table(a$cluster))
    nbcluster <- nbcluster + 1

    # History of clustering
    if(iter > 1) clustersplithist <- cbind(clustersplithist, clustersplithist[, iter])
    clustersplithist[I[J1], iter + 1] <- ksplit
    clustersplithist[I[J2], iter + 1] <- nbcluster    
    #length(varclass) +1
    
    # Update weights
    if(!perCluster_mv_weights) {
      if(use_mv_weights) {
        w <- mv_weights(X=X, mv=mv, centers=as.matrix(centers), 
                        cluster=clustersplithist[, iter + 1], 
                        gamma=gamma, mode=mode)
      } else{
        w <- data.frame(weights = rep(1 / length(mv), length(mv)),
                        weightsmv = rep(1 / length(mv), sum(mv)))
      }
      weights_save <- cbind(weights_save, w$weights)
    } else {
      w <- mv_weights_perCluster(X=X, mv=mv, centers=as.matrix(centers), 
                                 cluster=clustersplithist[, iter + 1], gamma=gamma,
                                 mode = mode)
      weights_save[[iter + 1]] <- w$weights
    }
    
    ## Update weighted variances by cluster
    if(!perCluster_mv_weights) {
      tmp <- centers[match(clustersplithist[, iter + 1], rownames(centers)),]
      varclass <- (w$weights ^ gamma) * 
        rowsum(colSums(rowsum((X - tmp)^2, group=clustersplithist[, iter + 1])), 
               group=rep(LETTERS[1:length(mv)], times=mv))
    } else {
      varclass <- rep(0, nbcluster)
      for (k in 1:nbcluster) {
        I <- which(clustersplithist[, iter + 1] == k)
        for (v in 1:length(mv)) {
          dv = (ref[v] + 1):ref[v + 1]
          varclass[k] <- varclass[k] + ((weights_save[[iter + 1]][k, v] ^ gamma) * sum((
            as.matrix(X)[I, dv] - matrix(rep(as.numeric(centers[k, dv]), length(I)), 
                                         nrow = length(I), byrow = T)) ^ 2))
        }
      }
    }
    
    # Calculate general criterion
    CRIT <- c(CRIT, sum(varclass))
    iter <- iter + 1
  }
  return(
    list(
      weights = weights_save,
      criterion = CRIT,
      withinss = varclass,
      split_clusters = clustersplithist,
      ksplit = ksplit_save
    )
  )
}


## NOT exported: function for calculating per-cluster multi-view weights

mv_weights_perCluster <- function(X, mv, centers, cluster, gamma, mode) {
  ref <- c(0, cumsum(mv))
  if(mode == "fuzzy" & !is.matrix(cluster)) 
    stop("Fuzzy clustering requires a matrix of posterior probabilities.")
  if(mode == "hard" & !is.vector(cluster)) 
    stop("Hard clustering requires a vector of cluster assignments.")
  if(mode == "hard") {
    weights = matrix(0, nrow = nrow(centers), ncol = length(mv))
    weightsmv = matrix(0, nrow = nrow(centers), ncol = sum(mv))
    for (k in 1:nrow(centers)) {
      aux = NULL
      for (v in 1:length(mv)) {
        dv = (ref[v] + 1):ref[v + 1]
        aux = c(aux, sum((
          as.matrix(X)[which(cluster == k), dv] - matrix(
            rep(centers[k, dv], sum(cluster == k)),
            nrow = sum(cluster == k),
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
  }
  if(mode == "fuzzy") {
    stop("Fuzzy clustering not currently supported in splitting algorithm")
  }
  return(list(weights = weights, weightsmv = weightsmv))
}
  