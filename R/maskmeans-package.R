#' Multi-view agglomeration/splitting K-means clustering algorithm
#'
#' Brief description of package ...
#'
#' \tabular{ll}{ Package: \tab maskmeans\cr Type: \tab Package\cr Version:
#' \tab 0.0.13\cr Date: \tab 2018-10-31\cr License: \tab GPL (>=3.3.1)\cr LazyLoad:
#' \tab yes\cr }
#'
#' @name maskmeans-package
#' @aliases maskmeans-package
#' @docType package
#' @author Andrea Rau, Cathy Maugis-Rabusseau, Antoine Godichon-Baggioni
#'
#' Maintainer: Andrea Rau <\url{andrea.rau@@inra.fr}>
#' @references
#' To be completed ...
#'
#' @keywords models cluster
#' @importFrom stats hclust cutree prcomp dist rnorm kmeans
#' @import ggplot2
#' @import clustree
#' @import BiocParallel
#' @importFrom viridis scale_color_viridis scale_fill_viridis
#' @importFrom ggdendro ggdendrogram
#' @importFrom tidyr gather_
#' @importFrom cowplot plot_grid
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom capushe capushe
#' @importFrom ggraph scale_edge_color_continuous
NULL

#gridExtra, GGally, reshape2, cluster, 

#' Multi-view agglomeration or splitting K-means clustering algorithm
#' 
#' This is the primary function to run either the agglomeration or splitting version of the multi-view
#' K-means clustering algorithm
#'
#' @param mv_data Multi-view data, either in the form of a concatenated matrix (where the columns of the views have
#' been combined) or in the form of a list, where each element contains one of the views. In the former case, the argument
#' \code{mv} should also be provided to indicate the dimension of each of the views (\code{sum(mv)} should be equal to
#' the number of columns in the matrix). In the latter case, each matrix in the list should have observations (rows) sorted
#' in the same order, as matrices will be combined by columns within the function; note that \code{mv} will be detected 
#' automtically based on the number of columns in each matrix.
#' @param clustering_init Initial hard or fuzzy clustering to be used for aggregating or splitting clusters.
#' @param type Either \code{"splitting"} or \code{"aggregation"}.
#' @param verbose If \code{TRUE}, provide verbose output.
#' @param parallel If \code{FALSE}, no parallelization. If \code{TRUE}, parallel
#' execution using BiocParallel (see next argument \code{BPPARAM}) for the fuzzy splitting algorithm. A note on running
#' in parallel using BiocParallel: it may be advantageous to remove large, unneeded objects
#' from the current R environment before calling the function, as it is possible that R's
#' internal garbage collection will copy these files while running on worker nodes.
#' @param BPPARAM Optional parameter object passed internally to \code{bplapply} when
#' \code{parallel=TRUE}. If not specified, the parameters last registered with \code{register}
#' will be used.
#' @param ... Additional optional parameters. See \code{?mv_aggregation} and \code{?mv_splitting} for more details.
#'
#' @return Output from either \code{mv_aggregation} or \code{mv_splitting}, according to the \code{type} of algorithm
#' specified above.
#' 
#' For the aggregation algorithm, the outputs are as follows:
#' \item{merged_clusters }{Matrix providing each pair of merged clusters at each iteration of the algorithm}
#' \item{hclust }{Object of class "hclust" to be used for plotting cluster aggregations}
#' \item{weights }{Matrix of dimension \code{v} x \code{niterations}, where \code{v} is the number of
#' views and \code{niterations} is the number of successive agglomerations, providing the multi-view weights}
#' \item{criterion }{Value taken on by the agglomerative criterion at each iteration}
#' \item{final_classification }{Final classification of observations using \code{final_K} clusters (obtained by other cutting the
#' true at the appropriate level for the aggregation algorithm)}
#' \item{final_K }{Number of clusters chosen via model selection by the data-drive slope estimation (DDSE) slope heuristics}
#' 
#' 
#' For the splitting algorithm, the outputs are as follows:
#' \item{split_clusters }{Matrix providing the history of each cluster splitting at 
#' each iteration of the algorithm}
#' \item{weights }{If \code{perCluster_mv_weights=FALSE}, a matrix of dimension 
#' \code{v} x \code{niterations}, where \code{v} is the number of
#' views and \code{niterations} is the number of successive agglomerations/splits, providing the multi-view weights. If 
#' \code{perCluster_mv_weights=TRUE}, a list of length \code{niterations}, where each element
#' of the list corresponds to per-cluster multi-view weights in matrices of dimension \code{K} x \code{v}, where \code{K} is the
#' number of clusters at a particular iteration and \code{v} is the number of views}
#' \item{criterion }{Value taken on by the splitting criterion at each iteration}
#' \item{withnss }{The within sum-of-squares for each cluster at the last iteration}
#' \item{ksplit }{Vector identifying which cluster was split at each iteration of the 
#' algorithm}
#' \item{final_classification }{Final classification of observations using \code{final_K} clusters (obtained by fixing the number 
#' of splits in the splitting algorithm)}
#' \item{final_K }{Number of clusters chosen via model selection by the data-drive slope estimation (DDSE) slope heuristics}
#' \item{all_probapost }{List of conditional probabilities for each split for the 
#' fuzzy splitting algorithm}
#' \item{final_probapost }{Matrix of conditional probabilities of cluster membership for
#' each observation for the model with \code{final_K} clusters}
#' 
#' @export
#' @example /inst/examples/maskmeans-package.R
maskmeans <- function(mv_data, clustering_init, type = "splitting", parallel=TRUE, BPPARAM=bpparam(), verbose=TRUE, ...) {
  ## Parse ellipsis function
  providedArgs <- list(...)
  arg.user <- list(gamma=2, use_mv_weights=TRUE, Kmax=NULL, perCluster_mv_weights=FALSE, mv=NULL, parallel=FALSE, BPPARAM=bpparam())
  arg.user[names(providedArgs)] <- providedArgs
  
  ## Format data: X, mv
  X <- mv_data
  if((is.matrix(mv_data) | is.data.frame(mv_data)) & is.null(arg.user$mv))
    stop("If multi-view data are provided as a matrix, the dimension of each view must be specified in mv.")
  if(is.list(mv_data) & !is.data.frame(mv_data)) {
    arg.user$mv <- unlist(lapply(mv_data, ncol))
    ## Sanity check on dimensions
    nr <- unlist(lapply(mv_data, nrow))
    if(sum(diff(nr))) stop("All views must be measured on the same set of observations.")
    ## Sanity check on rownames
    X <- do.call("cbind", mv_data)
    colnames(X) <- unlist(lapply(mv_data, colnames))
    for(i in 2:length(mv_data)) {
      if(is.null(rownames(mv_data[[i]]))) next;
      if(!all.equal(rownames(mv_data[[1]]), rownames(mv_data[[i]])))
        warning(paste0("Rownames of view", i, "do not match the first view. Are you sure that rows in each view are in the same order?"))
    }
    rownames(X) <- rownames(mv_data[[1]])
  }
  if(is.data.frame(mv_data)) {
    X <- as.matrix(mv_data)
    rownames(X) <- rownames(mv_data)
    colnames(X) <- colnames(mv_data)
  }
  
  ## Scale data
  X <- scaleview(X, unlist(arg.user$mv))
  
  ## Feed into appropriate algorithm
  all_probapost <- NULL
  if(type == "aggregation") {
    mv_run <- mv_aggregation(X=X, mv=arg.user$mv, clustering_init=clustering_init, 
                             gamma = arg.user$gamma, use_mv_weights = arg.user$use_mv_weights,
                             verbose=verbose)
    if(length(mv_run$hclust$labels[-1]) < 10) {
      message("DDSE for model selection is only possible if at least 10 cluster merges are performed.\nYou can use maskmeans_cutree() to cut the tree at a specific value of K if desired.")
      final_classification <- NULL
      final_probapost <- NULL
      final_K <- NULL
    } else {
      KDDSE <- suppressWarnings(selectK_aggregation(mv_run, X))
      ct <- maskmeans_cutree(mv_run, K=KDDSE, clustering_init=clustering_init)
      final_classification <- ct$classif
      final_probapost <- ct$probapost 
      final_K <- KDDSE
    }
  } else {
    if(is.null(arg.user$Kmax)) stop("Splitting algorithm requires the user to specify Kmax, the maximum number of clusters")
    mv_run <- mv_splitting(X=X, mv=arg.user$mv, clustering_init=clustering_init,
                           Kmax=arg.user$Kmax, gamma=arg.user$gamma, 
                           use_mv_weights = arg.user$use_mv_weights,
                           perCluster_mv_weights = arg.user$perCluster_mv_weights,
                           verbose=verbose,
                           parallel=arg.user$parallel, 
                           BPPARAM=arg.user$BPPARAM)
    if(is.null(mv_run$split_clusters)) {
      final_classification <- final_probapost <- final_K <- all_probapost <- NA
    } else if(length(apply(mv_run$split_clusters, 2, max)) < 10) {
      message("DDSE for model selection is only possible if at least 10 cluster splits are performed.\nYou can use maskmeans_cutree() to cut the tree at a specific value of K if desired.")
      final_classification <- NULL
      final_probapost <- NULL
      final_K <- NULL
      all_probapost <- mv_run$all_probapost
    } else {
      KDDSE <- suppressWarnings(selectK_splitting(mv_run, X))
      K <- seq(from = length(unique(mv_run$split_clusters[,1])),
               to = arg.user$Kmax,
               by=1)
      final_classification <- mv_run$split_clusters[,which(K == KDDSE)]
      all_probapost <- mv_run$all_probapost
      final_probapost <- mv_run$all_probapost[[which(unlist(lapply(mv_run$all_probapost, ncol)) == 
                                                         KDDSE)]]
      final_K <- KDDSE
    } 
  }
  mv_run[["final_classification"]] <- final_classification
  mv_run[["final_probapost"]] <- final_probapost
  mv_run[["final_K"]] <- final_K
  mv_run[["all_probapost"]] <- all_probapost
  
  class(mv_run) <- "maskmeans"
  return(mv_run)
}

## NOT exported: function to scale and normalize by view size
scaleview <- function(Data, mv){
  X <- scale(Data)
  ref <- cumsum(c(0,mv))
  for (v in 1:length(mv))
    X[,seq((ref[v]+1),ref[v+1])] = X[,seq((ref[v]+1),ref[v+1])] / mv[v]
  return(X)
}


## NOT exported: function to select the number of clusters for aggregation
##   Selection par capushe pour résultat de hv1 (hierarchie dure)   
##     --> il faut revoir ce qui joue le rôle de la dimension
selectK_aggregation <- function(obj, X) {
  M <- as.numeric(obj$hclust$labels[-1])
  M <- cbind(M, matrix(rep(as.numeric(obj$hclust$labels[-1]),2), ncol=2),
             rev(obj$criterion))
  selK <- capushe::capushe(M, n=nrow(data.frame(X)))
  K <- selK@DDSE@model
  return(as.numeric(K))
}


## NOT exported: function to select the number of clusters for splitting
##   Selection par capushe pour résultat de hv1 (splitting dure)  
##     --> il faut revoir ce qui joue le rôle de la dimension
selectK_splitting <- function(obj, X) {
  K <- apply(obj$split_clusters, 2, max)
#  M <- as.numeric(Res$labels[-1])
  M <- cbind(matrix(rep(K,3), ncol=3), obj$criterion[-1])
  selK <- capushe(M, n=nrow(data.frame(X)))
  K <- selK@DDSE@model
  return(as.numeric(K))
}

