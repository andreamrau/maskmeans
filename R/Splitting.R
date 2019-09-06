#' Splitting of hard or soft clusters based on multi-view data
#'
#' @param X Matrix of multi-view data, where the first view corresponds to the 
#' principal data used to obtain the partition or soft clustering in \code{cluster_init}
#' @param mv (Optional unless \code{X} is a matrix.) If \code{X} is a matrix, vector 
#' corresponding to the size of each data view. 
#' The sum of \code{mv} should correspond to the number of columns in \code{X}.
#' @param Kmax Maximum number of clusters for splitting
#' @param gamma Parameter that controls the distribution of view weights. Default value is 2. 
#' @param clustering_init Either a vector of available cluster labels (for hard clustering) 
#' or a matrix
#' of soft classification labels (for soft clustering)
#' @param use_mv_weights If \code{TRUE}, run algorithm in weighted multi-view mode; 
#' if FALSE, the
#' weight for each view is set to be equal. This option is only used for hard clustering.
#' @param perCluster_mv_weights If \code{TRUE}, use cluster-specific multi-view weights.
#' Otherwise use classic multi-view weights.
#' @param delta Parameter that controls the weights on the soft classifications pi(i,k)
#' @param verbose If \code{TRUE}, provide verbose output
#' @param parallel If \code{FALSE}, no parallelization. If \code{TRUE}, parallel
#' execution using BiocParallel (see next argument \code{BPPARAM}) for the soft splitting algorithm. A note on running
#' in parallel using BiocParallel: it may be advantageous to remove large, unneeded objects
#' from the current R environment before calling the function, as it is possible that R's
#' internal garbage collection will copy these files while running on worker nodes.
#' @param BPPARAM Optional parameter object passed internally to \code{bplapply} when
#' \code{parallel=TRUE}. If not specified, the parameters last registered with \code{register}
#' will be used.
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
#' \item{all_probapost }{List of conditional probabilities for each split for the 
#' soft algorithm}
#' 
#' @export
mv_splitting <- function(X, mv, clustering_init, Kmax, gamma=2, 
                         use_mv_weights = TRUE, delta = 2,
                         perCluster_mv_weights = TRUE,
                         verbose=TRUE, parallel=TRUE, BPPARAM=bpparam()) {
  
  cluster_init <- clustering_init
  
  if(gamma < 1) stop("gamma must be greater than or equal to 1.")
  if(delta <= 1) stop("delta must be strictly greater than 1.")
  mode <- ifelse(is.vector(clustering_init), "hard", "soft")
  cat("Running in mode:", mode, "\n")
  
  ref <- c(0, cumsum(mv))
  CRIT <- c()    ## Evolution of general criterion
  
  if(mode == "hard") { ## Hard mode
    clustersplithist <- matrix(c(cluster_init, cluster_init), ncol = 2)  ## Keep track of splits
    # Calculate initial cluster centers
    centers_init <- rowsum(X, group=cluster_init) / as.numeric(table(cluster_init))    
  } else {  ## Fuzzy mode
    ## Added clustersplithist here
    clustersplithist <- matrix(c(apply(cluster_init, 1, which.max), apply(cluster_init, 1, which.max)), ncol=2) ## Keep track of splits
    Kinit<- ncol(cluster_init)
    # Initialize centers: TODO fix this loop (?)
    centers_init <- matrix(0,nrow=Kinit, ncol=ncol(X))
    for(k in 1:Kinit){
      centers_init[k,]<- colSums(matrix(rep(cluster_init[,k]^delta,sum(mv)),nrow=nrow(X)) * X,2) / 
        sum(cluster_init[,k]^delta)
    }
    rownames(centers_init)<-1:Kinit
  }

  ## Calculate initial weights
  if(!perCluster_mv_weights) { ## Hard mode
    if(mode == "hard") {
      if(use_mv_weights) {
        w <- mv_weights(X=X, mv=mv, centers=as.matrix(centers_init), 
                        cluster=cluster_init, gamma=gamma, mode=mode)
      } else{
        w <- data.frame(weights = rep(1 / length(mv), length(mv)),
                        weightsmv = rep(1 / length(mv), sum(mv)))
      }
    } else { ## Fuzzy mode
      w <- mv_weights(X=X, mv=mv, centers=as.matrix(centers_init),
                      cluster=cluster_init, gamma=gamma, delta=delta, mode=mode)
    }
    weights_save <- w$weights    
  } else {
    if(mode == "hard") { ## Hard mode, per-cluster weights
      w <- mv_weights_perCluster(X=X, mv=mv, centers=as.matrix(centers_init), 
                                 cluster=cluster_init, gamma=gamma, mode = mode)
    } else { ## Fuzzy mode, per-cluster weights
      w <- mv_weights_perCluster(X=X, mv=mv, centers=centers_init, cluster=cluster_init, 
                                 gamma=gamma, delta=delta, mode=mode)
    }
    weights_save <- c()
    weights_save[[1]] <- w$weights
  }
  
  ## Calculate weighted per-cluster variances
  if(!perCluster_mv_weights) {
    if(mode == "hard") { ## Hard mode
      tmp <- centers_init[match(cluster_init, rownames(centers_init)),]
      varclass <- rowSums(rowsum(t(w$weights ^ gamma * 
                                     rowsum(t((X - tmp)^2), group=rep(LETTERS[1:length(mv)], times=mv))), 
                                 group=cluster_init))      
    } else { ## Fuzzy mode:
      varclass <- rep(0, Kinit)
      for(k in 1:Kinit) {
        for(v in 1:length(mv)) {
          dv <- (ref[v]+1):ref[v+1]
          varclass[k] <- varclass[k] + sum(w$weights[v]^gamma * rowSums(cluster_init[,k]^delta * (t(t(X[,dv]) - centers_init[k,dv]))^2))
        }
      }
    }
  } else {
    if(mode == "hard") {  ## Hard mode
      tmp <- centers_init[match(cluster_init, rownames(centers_init)),]
      varclass_tmp <- t(rowsum(t((X - tmp)^2), group=rep(LETTERS[1:length(mv)], times=mv)))
      weights_tmp <- w$weights
      rownames(weights_tmp) <- 1:nrow(weights_tmp)
      weights_tmp <- weights_tmp[match(cluster_init, rownames(weights_tmp)),]
      varclass <- rowSums(rowsum((weights_tmp ^ gamma) * varclass_tmp, group = cluster_init))
    } else {  ## Fuzzy mode:
      varclass <- rep(0, Kinit)
      for (k in 1:Kinit){
        for (v in 1:length(mv)){
          dv <- (ref[v]+1):ref[v+1]
          varclass[k] <- varclass[k] + sum(w$weights[k,v]^gamma * 
                                             rowSums(cluster_init[,k]^delta * (t(t(X[,dv]) - centers_init[k,dv]))^2))
        }
      }
    }
   }
  
  # Calculate general criterion
  CRIT <- c(CRIT, sum(varclass)) 
  ksplit_save <- NULL
  cluster <- cluster_init
  
  ## For fuzzy algorithm, keep track of the probaposts
  if(mode == "soft") {
    all_probapost <- list()
    all_probapost_i <- 1
    all_probapost[[all_probapost_i]] <- cluster_init
    all_probapost_i <- all_probapost_i + 1
  }
  
  ## Start splitting loop
  iter <- 1
  nbcluster <- ifelse(mode == "hard", max(cluster_init), Kinit)
  centers <- as.data.frame(centers_init)
  while(nbcluster < Kmax) {
    
    # Cluster with max varclass
    ksplit <- which.max(varclass)
    ksplit_save <- c(ksplit_save, ksplit)
    
    if(verbose) {
      cat("  => Splitting cluster", ksplit, "\n")
    }
    
    # Kmeans in the cluster to be split
    # pb avec les data dans le kmeans  -> la multiplication par les poids déforme trop
    if(mode == "hard") {
      I <- which(clustersplithist[, iter] == ksplit)
      if(!perCluster_mv_weights) {
        A <- t(w$weightsmv ^ (gamma/2) * t(X[I,])) 
      } else {
        A <- NULL
        for (i in 1:length(I)) {
          A <- rbind(A, (w$weightsmv[ksplit, ] ^ (gamma / 2)) * X[I[i], ])
          #       ### !!!  Error in previous code had this line:
          #       A <- rbind(A, (w$weightsmv[max(nbcluster), ] ^ (gamma / 2)) * X[I[i], ])
        } 
      }
      if(length(I) > 2) {
        a <- kmeans(A, centers = 2, nstart = 50, iter.max = 50)        
        J1 <- which(a$cluster == 1)
        J2 <- which(a$cluster == 2)
      } else {
        J1 <- 1
        J2 <- 2
        a<- list(cluster = c(1,2))
      }
      
      # Update centers
      centers[c(ksplit, nrow(centers)+1),] <- 
        rowsum(X[I,], group=a$cluster) / as.numeric(table(a$cluster))
      nbcluster <- nbcluster + 1
      
      # History of clustering
      if(iter > 1) clustersplithist <- cbind(clustersplithist, clustersplithist[, iter])
      clustersplithist[I[J1], iter + 1] <- as.numeric(ksplit)
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
        varclass <- rowSums(rowsum(t(w$weights ^ gamma * 
                                       rowsum(t((X - tmp)^2), group=rep(LETTERS[1:length(mv)], times=mv))), 
                                   group=clustersplithist[, iter + 1]))
      } else {
        tmp <- centers[match(clustersplithist[, iter + 1], rownames(centers)),]
        varclass_tmp <- t(rowsum(t((X - tmp)^2), group=rep(LETTERS[1:length(mv)], times=mv)))
        weights_tmp <- weights_save[[iter + 1]]
        rownames(weights_tmp) <- 1:nrow(weights_tmp)
        weights_tmp <- weights_tmp[match(clustersplithist[, iter + 1], rownames(weights_tmp)),]
        varclass <- rowSums(rowsum((weights_tmp ^ gamma) * varclass_tmp, group = clustersplithist[, iter + 1]))
      }
    } else {  ## Fuzzy mode
      ## Add clustersplithist here
      I <- which(clustersplithist[,iter] == ksplit)
      
      ## Splitting in annex functions
      A <- .FuzzyKmeansSplit(X=X, mv=mv, gamma=gamma, delta=delta, w=w, ksplit=ksplit,
                            Pik=cluster[,ksplit], ninit=20, niter=20, K=2, parallel=parallel,
                            BPPARAM=BPPARAM)
      ## AR: I think we only want to split the observations from the chosen cluster here
#      A <- .FuzzyKmeansSplit(X=X[I,], mv=mv, gamma=gamma, delta=delta, w=w, ksplit=ksplit,
#                             Pik=cluster[I,ksplit], ninit=20, niter=20, K=2, parallel=parallel,
#                             BPPARAM=BPPARAM)
      
      ## Update centers
      centers[ksplit,] <- A$centersNew[1,]
      centers <- rbind(centers,A$centersNew[2,])
      nbcluster <- nbcluster+1
      rownames(centers) <- 1:nbcluster
      
      ## Update posterior probabilities
      ## AR: I think we only want to split the observations from the chosen cluster here
      cluster[,ksplit] <- A$Pinew[,1]
      cluster <- cbind(cluster, A$Pinew[,2])
      ## AR: Should we divide the conditional probability for observations in unsplit classes between the two new ones?
#      cluster[I,ksplit] <- A$Pinew[,1]
#      cluster <- cbind(cluster, 0)
#      cluster[I,ncol(cluster)] <- A$Pinew[,2]
#      cluster[-I,ksplit] <- cluster[-I,ncol(cluster)] <- cluster[-I,ksplit]/2
      
      ## Keep track of the probapost
      all_probapost[[all_probapost_i]] <- cluster
      all_probapost_i <- all_probapost_i + 1
      
      ## Added update of clustering history here
      if(iter > 1) clustersplithist <- cbind(clustersplithist, clustersplithist[, iter])
      clustersplithist[, iter+1] <- apply(cluster, 1, which.max)
      
      #mise à jour des poids alpha
      if(perCluster_mv_weights) {
        w <- mv_weights_perCluster(X=X, mv=mv, centers=as.matrix(centers), cluster=cluster, 
                                   gamma=gamma, delta=delta, mode=mode)
        weights_save[[iter + 1]] <- w$weights  
        varclass <- rep(0,nbcluster)
        for (k in 1:nbcluster){
          for (v in 1:length(mv)){
            dv <- (ref[v]+1):ref[v+1]
            varclass[k] <- varclass[k] + sum(w$weights[k,v]^gamma * 
                                               rowSums(cluster[,k]^delta * (t(t(X[,dv]) - unlist(centers[k,dv])))^2))
          }
        }
      } else {
        w <- mv_weights(X=X, mv=mv, centers=as.matrix(centers),
                        cluster=cluster, gamma=gamma, delta=delta, mode=mode)
        weights_save <- cbind(weights_save, w$weights)
        ## Update weighted variances per cluster
        varclass <- rep(0,nbcluster)
        for (k in 1:nbcluster){
          for (v in 1:length(mv)){
            dv <- (ref[v]+1):ref[v+1]
            varclass[k] <- varclass[k] + sum(w$weights[v]^gamma * 
              rowSums(cluster[,k]^delta * (t(t(X[,dv]) - unlist(centers[k,dv])))^2))
          }
        }
      }
    }

    # Calculate general criterion
    CRIT <- c(CRIT, sum(varclass))
    iter <- iter + 1
  }
  ret <- list(
    weights = weights_save,
    criterion = CRIT,
    withinss = varclass,
    ksplit = ksplit_save, 
    split_clusters = clustersplithist
  )
  if(mode == "soft") {
    ret$all_probapost <- all_probapost
  }
  return(ret)
}


## NOT exported: function for calculating per-cluster multi-view weights

mv_weights_perCluster <- function(X, mv, centers, cluster, gamma, mode, delta=NULL) {
  ref <- c(0, cumsum(mv))
  if(mode == "soft" & !is.matrix(cluster)) 
    stop("Soft clustering requires a matrix of posterior probabilities.")
  if(mode == "hard" & !is.vector(cluster)) 
    stop("Hard clustering requires a vector of cluster assignments.")
  if(mode == "soft" & is.null(delta))
    stop("Soft clustering requires a value for delta.")
  weights = matrix(0, nrow = nrow(centers), ncol = length(mv))
  weightsmv = matrix(0, nrow = nrow(centers), ncol = sum(mv))
  for (k in 1:nrow(centers)) {
    aux <- NULL
    for (v in 1:length(mv)) {
      dv = (ref[v] + 1):ref[v + 1]
      if(mode == "hard") {  ## Hard mode
        aux = c(aux, sum((
          as.matrix(X)[which(cluster == k), dv] - matrix(
            rep(centers[k, dv], sum(cluster == k)),
            nrow = sum(cluster == k),
            byrow = T)) ^ 2))
      } else {  ## Fuzzy mode
          aux=c(aux, sum((cluster[,k]^delta) * 
                           rowSums((sweep(X[,dv,drop=FALSE],2,
                                          centers[k,dv],FUN="-"))^2)) )
      }
    }
    if (gamma > 1) {
      if(sum(aux) != 0) {
        aux = aux ^ (1 / (1 - gamma))
        
        aux = aux / sum(aux)
      }
    }
    if (gamma == 1) {
      if(sum(aux) != 0) {
        aux = aux ^ (-1)
        aux = aux / sum(aux)
      }
    }
    aux1 <- NULL
    for (j in 1:length(mv))
      aux1 <- c(aux1, rep(aux[j], mv[j]))
    weights[k, ] = aux
    weightsmv[k, ] = aux1
  }
  return(list(weights = weights, weightsmv = weightsmv))
}

## NOT exported: fuzzy K-means with vector of probabilities

###  j'ai du mettre plein d'exceptions à cause du cas delta=1 qui mene à des CRIT=NaN
.FuzzyKmeansSplit <- function(X,mv,gamma,delta,w,ksplit,Pik,ninit=20,
                              niter=20, K=2, parallel=FALSE, BPPARAM=bpparam()) {
  #pour stocker les réponses
  Pinew <- matrix(0,nrow=nrow(X),ncol=K)
  centersNew <- matrix(0,nrow=K,ncol=ncol(X))
  A1 <- .FuzzyKmeansaux(X,w,gamma,delta, ksplit, Pik, K,niter)
  
  ## Not parallel mode
  if(!parallel) {
    for(z in 1:ninit) {   #boucle sur le nombre d'essai du fuzzy Kmeans
      A2 <- .FuzzyKmeansaux(X,w,gamma,delta, ksplit, Pik, K,niter)
      if(!is.na(A1$CRIT)){
        if (!is.na(A2$CRIT)){
          if (A2$CRIT<A1$CRIT){
            A1 <- A2
          }
        }
      } else{
        if(!is.na(A2$CRIT)){
          A1 <- A2
        }
      }
    }  
    Pinew <- A1$Pinew
    centersNew <- A1$centersNew
  } else{
    ## Parallel mode
    tmp <- bplapply(1:ninit, 
                    FUN=function(ii,P_X,P_w,P_gamma,P_delta,P_ksplit,
                                          P_Pik,P_K,P_niter) {
      set.seed(ii)
      res <- .FuzzyKmeansaux(X=P_X,w=P_w,gamma=P_gamma,delta=P_delta,ksplit=P_ksplit,
                             Pik=P_Pik,K=P_K,niter=P_niter)
      return(res)
    }, P_X=X,P_w=w,P_gamma=gamma,P_delta=delta,P_ksplit=ksplit,
       P_Pik=Pik,P_K=K,P_niter=niter,BPPARAM=BPPARAM)
    CRIT_all <- unlist(lapply(tmp, function(x) x$CRIT))
    choose <- which.min(CRIT_all)
    Pinew <- tmp[[choose]]$Pinew
    centersNew <- tmp[[choose]]$centersNew
  }

  return(result=list(Pinew= Pinew, centersNew= centersNew))
}

## NOT exported: fuzzy K-means auxiliary function with vector of probabilities
.FuzzyKmeansaux <- function(X, w,gamma,delta, ksplit, Pik, K,niter){
  # Z_i^{(v)} = alpha_v^{gamma/2} X_i^{(v)}
  if (is.vector(w$weights)) {       # même poids par classe (\alpha_v)
    Z <- data.frame(matrix(rep(w$weightsmv^(gamma/2),nrow(X)),nrow=nrow(X),byrow=T) * X)
  } else{       #poids(\alpha_{k,v})
    Z <- data.frame(matrix(rep(w$weightsmv[ksplit,]^(gamma/2),nrow(X)),nrow=nrow(X),byrow=T) * X)  
  }
  #init de stockage des Pinew et centersNew
  Pinew <- matrix(0,nrow=nrow(X),ncol=K)
  centersNew <- matrix(0,nrow=K,ncol=ncol(X))
  eta <- matrix(0,nrow=K,ncol=ncol(X))
  m <- matrix(0,nrow=nrow(X),ncol=K)   # \|Z_i - \eta_k\|^2
  CRIT <- NULL
  
  #initialisation des centres
  centersNew <- kmeans(X,K)$centers
  if (delta<= 1)  stop("Error delta value")
  
    for(l in 1:niter){
      # Mise à jour des poids
      # for (k in 1:K){  
      #   if (is.vector(w$weights)){       # même poids par classe (\alpha_v)
      #     eta[k,] <- centersNew[k,] * w$weightsmv^(gamma/2)
      #   } else {
      #     eta[k,] <- centersNew[k,] * w$weightsmv[ksplit,]^(gamma/2)  
      #   }
      #   m[,k] <- rowSums((matrix(rep(eta[k,],nrow(Z)),byrow=TRUE, nrow=nrow(Z))-Z)^2)
      # }
      if (is.vector(w$weights)){       # même poids par classe (\alpha_v)
        eta <- t(t(centersNew) * w$weightsmv^(gamma/2))
      } else {
        eta <- t(t(centersNew) * w$weightsmv[ksplit,]^(gamma/2))  
      }
      etamat <- eta[rep(1:2, each=nrow(Z)),]
      m <- matrix(rowSums((etamat - Z[rep(1:nrow(Z), 2),])^2), ncol=K)
      
      Pinew <- (m^(1/(1-delta))) / matrix(rep(apply(m^(1/(1-delta)),1,sum),K),ncol=K)
      Pinew <-  matrix(rep(Pik,K),ncol=K) * Pinew
      
      # Mise à jour des centres
      # for (k in 1:K)
      #   centersNew[k,] <- colSums(matrix(rep((Pinew[,k]^delta),ncol(X)),nrow=nrow(X)) * X) / 
      #     sum(Pinew[,k]^delta)
      centersNew <- rowsum(matrix(as.vector(Pinew)^delta, nrow=K*nrow(Pinew), ncol=ncol(X)) * X[rep(1:nrow(X), K),], 
                           group = rep(1:K, each=nrow(X))) /
        colSums(Pinew^delta)
    } # fin de la boucle ninter
    
    # calcul du critere
    # for (k in 1:K){ 
    #   if (is.vector(w$weights)){       # même poids par classe (\alpha_v)
    #     eta[k,] <- centersNew[k,] * w$weightsmv^(gamma/2)
    #   } else{
    #     eta[k,] <- centersNew[k,] * w$weightsmv[ksplit,]^(gamma/2)  
    #   }
    #   m[,k]=rowSums((matrix(rep(eta[k,],nrow(Z)),byrow=T,nrow=nrow(Z))-Z)^2)
    # }
    if (is.vector(w$weights)){       # même poids par classe (\alpha_v)
      eta <- t(t(centersNew) * w$weightsmv^(gamma/2))
    } else {
      eta <- t(t(centersNew) * w$weightsmv[ksplit,]^(gamma/2))  
    }
    etamat <- eta[rep(1:2, each=nrow(Z)),]
    m <- matrix(rowSums((etamat - Z[rep(1:nrow(Z), 2),])^2), ncol=K)
    CRIT <- sum((Pinew^delta) * m)
  
  return(list(Pinew=Pinew,centersNew=centersNew,CRIT=CRIT))  
}




#' Extract the posterior probabilities from a specified step in the splitting tree
#'
#' @param probapost Final matrix of posterior probabilities
#' @param ksplit Desired level of the splitting tree
#' @param K Desired number of clusters
#'
#' @return Matrix of posterior probabilities.
#' @export
#'
probapost_K <- function(probapost, ksplit, K){
  Kmax <- ncol(probapost)
  Msplit <- .MatrixKsplit(ksplit,Kmax)
  I <- which(apply(Msplit,1,max)==K)    # la ligne de Msplit qui contient la répartition des Kmax classes en Kref classes
  probapostextract <- matrix(0,nrow=nrow(probapost),ncol=K)
  for (k in 1:K)
    probapostextract[,k] <- apply(probapost[,which(Msplit[I,]==k),drop=FALSE],1,sum)
  return(probapostextract)
}


#' Extract the classification by MAP estimator (complete or "smooth")
#'
#' @param probapost Matrix of posterior probabilities
#' @param seuil Threshold, by default 0.80
#'
#' @return Vector of cluster labels. Observations with maximum conditional probabilithy less than \code{threshold} are
#' assigned a cluster label of 0.
#' @export
MAP_smooth <- function(probapost, seuil=0.8){
  if(seuil<1 & 0<=seuil){
    label <- rep(0,nrow(probapost))
    label <- apply(probapost,1,which.max) * (apply(probapost,1,max)>seuil) 
    return(label)
  }
  else{cat("Error for threshold value")}
}



## Probapost helper functions (not exported)

#######################
##   exploitation de l'historique du splitting 
##  la matrice Msplit explique sur chaque ligne la répartition des classes 
##  (1ere ligne le nb classes totales 1 à Kmax jusqu'à sur la dernière ligne la répartition en 1:Kinit classes)
.MatrixKsplit <- function(ksplit,Kmax){
  Kref <- Kmax
  Msplit <- matrix(rep(1:Kref,2),nrow=2,byrow=TRUE)
  iter <- length(ksplit)
  while(iter>0){
    I <- c(Kref,which(Msplit[nrow(Msplit)-1,]==Kref))
    Msplit[nrow(Msplit),I] <- ksplit[iter]
    iter <- iter-1
    Kref <- Kref-1
    Msplit <- rbind(Msplit,Msplit[nrow(Msplit),])
  }
  return(Msplit[-nrow(Msplit),])
}






