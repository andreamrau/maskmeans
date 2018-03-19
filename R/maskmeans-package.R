#' Multi-view agglomeration/splitting K-means clustering algorithm
#'
#' Brief description of package ...
#'
#' \tabular{ll}{ Package: \tab maskmeans\cr Type: \tab Package\cr Version:
#' \tab 0.0.1\cr Date: \tab 2018-03-13\cr License: \tab GPL (>=3)\cr LazyLoad:
#' \tab yes\cr }
#'
#' @name coseq-package
#' @aliases maskmeans-package
#' @docType package
#' @author Andrea Rau, Cathy Maugis-Rabusseau, Antoine Godichon-Baggioni
#'
#' Maintainer: Andrea Rau <\url{andrea.rau@@inra.fr}>
#' @references
#' To be completed ...
#'
#' @keywords models cluster
#' @importFrom stats hclust cutree
NULL

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
#' @param ... Additional optional parameters.
#'
#' @return Output from either \code{mv_aggregation} or \code{mv_splitting}, according to the \code{type} of algorithm
#' specified above.
#' 
#' @export
maskmeans <- function(mv_data, clustering_init, type = "splitting", ...) {
  ## Parse ellipsis function
  providedArgs <- list(...)
  arg.user <- list(gamma=2, use_mv_weights=TRUE, Kmax=NULL, perCluster_mv_weights=TRUE, mv=NULL)
  arg.user[names(providedArgs)] <- providedArgs
  
  ## Format data: X, mv
  X <- mv_data
  if((is.matrix(mv_data) | is.data.frame(mv_data)) & is.null(arg.user$mv))
    stop("If multi-view data are provided as a matrix, the dimension of each view must be specified in mv.")
  if(is.list(mv_data) & !is.data.frame(mv_data)) {
    arg.user$mv <- lapply(mv_data, ncol)
    ## Sanity check on dimensions
    nr <- unlist(lapply(mv_data, nrow))
    if(sum(diff(nr))) stop("All views must be measured on the same set of observations.")
    ## Sanity check on rownames ? TODO
    X <- do.call("cbind", mv_data)
    colnames(X) <- unlist(lapply(mv_data, colnames))
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
  if(type == "aggregation") {
    mv_run <- mv_aggregation(X=X, mv=arg.user$mv, clustering_init=clustering_init, 
                             gamma = arg.user$gamma, use_mv_weights = arg.user$use_mv_weights)
  } else {
    if(is.null(arg.user$Kmax)) stop("Splitting algorithm requires the user to specify Kmax, the maximum number of clusters")
    mv_run <- mv_splitting(X=X, mv=arg.user$mv, clustering_init=clustering_init,
                           Kmax=arg.user$Kmax, gamma=arg.user$gamma, 
                           use_mv_weights = arg.user$use_mv_weights,
                           perCluster_mv_weights = arg.user$perCluster_mv_weights)
  }
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
