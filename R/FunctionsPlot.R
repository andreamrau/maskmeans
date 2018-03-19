#' Overview plot of multi-view data
#' 
#' Global graphical representation of multi-view data.
#'
#' @param mv_data Multi-view data, either in the form of a concatenated matrix (where the columns of the views have
#' been combined) or in the form of a list, where each element contains one of the views. In the former case, the argument
#' \code{mv} should also be provided to indicate the dimension of each of the views (\code{sum(mv)} should be equal to
#' the number of columns in the matrix). In the latter case, each matrix in the list should have observations (rows) sorted
#' in the same order, as matrices will be combined by columns within the function; note that \code{mv} will be detected 
#' automtically based on the number of columns in each matrix.
#' @param scale If \code{TRUE}, data will be scaled so that all variables are normalized and views are standardized according 
#' to their size.
#' @param ... Additional optional parameters. In particular, if points should be colored by cluster membership, a vector
#' \code{labels} providing cluster assigments should be provided.
#'
#' @return A multi-facetted plot, with one facet per view. Univariate views are represented as densities, bivariate views
#' are represented as scatterplots, and multivariate views are represented by plotting scatterplots of the first two 
#' principal components.
#' 
#' @export
mv_plot <- function(mv_data, scale=TRUE, ...) {
  ## Parse ellipsis function
  providedArgs <- list(...)
  arg.user <- list(mv=NULL, labels=NULL)
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
  
  if(is.null(arg.user$labels)) arg.user$labels <- rep(NA, nrow(X))
  
  ## Scale data if desired
  if(scale) X <- scaleview(X, unlist(arg.user$mv))
  
  Xplot <- matrix(nrow=0, ncol=4)
  colnames(Xplot) <- c("x", "y", "view", "labels")
  ref <- c(0, cumsum(arg.user$mv))
  for(i in 1:length(arg.user$mv)) {
    dv <- (ref[i] + 1):ref[i + 1]
    if(arg.user$mv[i] == 1) {
      tmp <- data.frame(x=X[,dv], y=rep(NA, nrow(X)), view=rep(paste0("View_", i), nrow(X)),
                        labels=arg.user$labels)
      Xplot <- rbind(Xplot, tmp)
    } else if(arg.user$mv[i] == 2) {
      tmp <- data.frame(X[,dv], rep(paste0("View_", i), nrow(X)), labels=arg.user$labels)
      colnames(tmp) <- c("x", "y", "view", "labels")
      Xplot <- rbind(Xplot, tmp)
    } else {
      pca1 <- prcomp(X[,dv])
      scores <- data.frame(pca1$x)[,1:2]
      tmp <- cbind(scores, rep(paste0("View_", i), nrow(X)), labels=arg.user$labels)
      colnames(tmp) <- c("x", "y", "view", "labels")
      Xplot <- rbind(Xplot,tmp)
    }
  }
  
  if(sum(is.na(arg.user$labels))) {
    g <- ggplot(Xplot) + 
      geom_point(data = subset(Xplot, !is.na(y)), aes(x=x, y=y), alpha=0.25) + 
      geom_density(data = subset(Xplot, is.na(y)), aes(x=x), alpha=0.25) + 
      facet_wrap(~view, scales="free_y") +
      theme_bw() 
  } else {
    g <- ggplot(Xplot) + 
      geom_point(data = subset(Xplot, !is.na(y)), 
                 aes(x=x, y=y, color=factor(labels), fill=factor(labels)), alpha=0.25) + 
      geom_density(data = subset(Xplot, is.na(y)), aes(x=x, fill=factor(labels), color=factor(labels)), alpha=0.25) + 
      facet_wrap(~view, scales="free_y") +
      guides(color=FALSE, fill=FALSE) + 
      theme_bw() +
      scale_color_viridis(discrete=TRUE) + scale_fill_viridis(discrete=TRUE) 
  }
  print(g)
}


maskmeans_plot <- function(obj, 
                           type = c("dendrogram", "heights", "weights_line", "weights_area", "criterion")) {
  
  ## TODO: Add splitting plot ?
  ## TODO: Add weight plots for per cluster weights
  
  if(class(obj) != "maskmeans") stop("This plot function expects an object of class maskmeans.")

  if("heights" %in% type) { ## The heights of the tree
    if(!"merged_clusters" %in% names(obj)) 
      stop("dendrogram heights plot is only possible for the aggregation algorithm.")
    g1 <- ggplot(data.frame(x = as.factor(seq(2, length(obj$hclust$height))), 
                            y = rev(obj$hclust$height)[-1]), 
                 aes(x=x, y=y)) +
      geom_bar(stat = "identity") +
      ylab("Height") +
      xlab("Nb of clusters")
    print(g1)
  }
  if("dendro" %in% type) { ## Dendrogram
    if(!"merged_clusters" %in% names(obj)) 
      stop("dendrogram heights plot is only possible for the aggregation algorithm.")
    g2 <- ggdendrogram(obj$hclust, rotate=FALSE, size=2)
    print(g2)
  }
  if("criterion" %in% type) { ## Evolution of criterion
    if("merged_clusters" %in% names(obj)) {
      df <- data.frame(x = seq(2, length(obj$hclust$order)), 
                       y = rev(obj$criterion))
    } else {
      df <- data.frame(x = apply(obj$split_clusters, 2, max),
                 y = obj$criterion[-1])
    }
    g5 <- ggplot(df, aes(x = x, y = y)) +
      geom_point() +
      geom_line() +
      ylab("Criterion value") +
      xlab("Number of clusters")
    print(g5)
  }
  if("weights_area" %in% type) { ## Evolution of weights
    if(class(obj$weights) == "list") stop("I'm still working on the plot functions for per-cluster weights in the splitting algorithm.")
    
    df <- data.frame(cbind(seq(1, ncol(obj$weights)), t(obj$weights)))
    rownames(df) <- NULL
    colnames(df) <- c("iteration", as.character(seq(1:nrow(obj$weights))))
    df <- gather(df, key=view, value=value, -iteration)
    g3 <- ggplot(df, aes(iteration, value)) +
      geom_area(aes(fill = view, group = view)) +
      ylab("weights") +
      scale_fill_viridis(discrete=TRUE)
    print(g3)
  }
  if("weights_line" %in% type) {
    if(class(obj$weights) == "list") stop("I'm still working on the plot functions for per-cluster weights in the splitting algorithm.")
  
    df <- data.frame(cbind(seq(1, ncol(obj$weights)), t(obj$weights)))
    rownames(df) <- NULL
    colnames(df) <- c("iteration", as.character(seq(1:nrow(obj$weights))))
    df <- gather(df, key=view, value=value, -iteration)
    g4 <- ggplot(df, aes(iteration, value, group = view, colour = view)) +
      geom_point() +
      geom_line(aes(lty = view)) +
      scale_fill_viridis(discrete=TRUE)
    print(g4)
  }
}

#############################################
##   exploitation graphique de la sortie splittingClustersbis
#############################################
PlotResSplittingbis <-
  function(Res,
           weights_step = FALSE,
           weights_clust = FALSE) {
    
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
##    fonction clustree pour visualiser les étapes de splitting
##     prend l'objet Res venant de splittingClusters et les données X
############################
clustreebis <- function(Res, X) {
  aux <- Res$clustersplithist
  colnames(aux) <- paste0("K", apply(Res$clustersplithist, 2, max))
  dataPlot <- cbind(X, aux)
  clustree(data.frame(dataPlot), prefix = "K")
}

##  Not exported:
probapost_boxplot <- function(probapost) {
  aux <- data.frame(probamax = apply(probapost, 1, max),
                    lab = as.factor(apply(probapost, 1, which.max)))
  p <- ggplot(aux, aes(x = lab, y = probamax)) +
    geom_boxplot()
  print(p)
}

## Not exported:
probapost_threshold <- function(obj, probapost, threshold = 0.8) {
  aux <- NULL
  for (k in 2:length(obj$criteiron)) {
    aux <- c(aux, sum(apply(
      cutreeNewProbapost(obj$hclust, K=k, probapost)$probapost, 1, max) > threshold) * 100 / 
        nrow(probapost))
  }
  p <- ggplot(data.frame(k = 2:length(obj$criterion), aux = aux), aes(k, aux)) +
    geom_line() +
    geom_point() +
    xlab("Number of clusters") +
    ylab("Percent of max posterior prob > threshold")
  print(p)
}


##!!  I have not implemented the following function:
# PlotClassif <- function(data, classif) {
#   don <- data.frame(classif = as.factor(classif), data = data)
#   ggpairs(don,
#           mapping = aes(color = classif),
#           columns = colnames(don)[-1])
# }
