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
    arg.user$mv <- unlist(lapply(mv_data, ncol))
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
      tmp <- data.frame(x=X[,dv], y=rep(NA, nrow(X)), view=rep(paste0("View ", i), nrow(X)),
                        labels=arg.user$labels)
      Xplot <- rbind(Xplot, tmp)
    } else if(arg.user$mv[i] == 2) {
      tmp <- data.frame(X[,dv], rep(paste0("View ", i), nrow(X)), labels=arg.user$labels)
      colnames(tmp) <- c("x", "y", "view", "labels")
      Xplot <- rbind(Xplot, tmp)
    } else {
      pca1 <- prcomp(X[,dv])
      scores <- data.frame(pca1$x)[,1:2]
      tmp <- cbind(scores, rep(paste0("View ", i), nrow(X)), labels=arg.user$labels)
      colnames(tmp) <- c("x", "y", "view", "labels")
      Xplot <- rbind(Xplot,tmp)
    }
  }
  
  if(sum(is.na(arg.user$labels))) {
    g <- ggplot(Xplot) + 
      geom_point(data = Xplot[which(is.na(Xplot$y)==FALSE),], aes_string(x="x", y="y"), alpha=0.25) + 
      geom_density(data = Xplot[which(is.na(Xplot$y)==TRUE),], aes_string(x="x"), alpha=0.25) + 
      facet_wrap("view", scales="free_y") +
      xlab("") + ylab("") +
      labs(caption="Univariate views are represented by density plots, bivariate views by scatterplots,\nand multivariate views by scatterplots of the first two principal components.") +
      theme_bw() 
  } else {
    Xplot$labels <- factor(Xplot$labels)
    g <- ggplot(Xplot) + 
      geom_point(data = Xplot[which(is.na(Xplot$y)==FALSE),], 
                 aes_string(x="x", y="y", color="labels", fill="labels"), alpha=0.25) + 
      geom_density(data = Xplot[which(is.na(Xplot$y)==TRUE),], 
                   aes_string(x="x", fill="labels", color="labels"), alpha=0.25) + 
      facet_wrap("view", scales="free_y") +
      guides(color=FALSE, fill=FALSE) + 
      xlab("") + ylab("")+
      labs(caption="Univariate views are represented by density plots, bivariate views by scatterplots,\nand multivariate views by scatterplots of the first two principal components.") +
      theme_bw() 
#    +  viridis::scale_color_viridis(discrete=TRUE) + viridis::scale_fill_viridis(discrete=TRUE) 
  }
  return(g)
}

#' Plot results of the multi-view aggregation/splitting K-means algorithm
#'
#' Produce a variety of plots after running the multi-view aggregation or splitting
#' version of the K-means algorithm.
#' 
#' @param obj Object of class \code{"maskmeans"} resulting from a call to the \code{maskemans}
#' function
#' @param type Graphic to be produced: one or more of \code{c("dendrogram", "heights", 
#' "weights_line", "weights", "criterion", "tree", "tree_perClusterWeights", "probapost_boxplot")}.
#' @param tree_type Either \code{"final_K"}, \code{"all"}, or a numerical value for the final number of clusters for plots of type \code{"tree"} or \code{"tree_perClusterWeights"}.
#' @param mv_names If desired, a vector of multiview names to be used for the \code{tree_perClusterWeights} plot
#' @param ... Additional optional parameters. 
#'
#' @return List of one or more ggplot2 objects.
#' @export
maskmeans_plot <- function(obj, 
                           type = c("dendrogram", "heights", "weights_line", 
                                    "weights", "criterion", "tree", "tree_perClusterWeights"),
                           tree_type = "final_K",
                           mv_names = NULL,
                           ...) {
  ## Parse ellipsis function
  providedArgs <- list(...)
  arg.user <- list(mv_data=NULL, edge_arrow=TRUE)
  arg.user[names(providedArgs)] <- providedArgs
  
  if(class(obj) != "maskmeans") stop("This plot function expects an object of class maskmeans.")

  g <- vector("list", length=0)

#  if(!is.null(arg.user$mv_data) & !("tree" %in% type))
#    message("Note: you have included mv_data, but this argument is used only when the 'tree' graph is produced for the splitting algorithm.")
  
  if("heights" %in% type) { ## The heights of the tree
    if(!"merged_clusters" %in% names(obj)) {
      message("-- dendrogram heights plot is only possible for the aggregation algorithm.")
    } else {
      g1 <- ggplot(data.frame(x = as.factor(seq(2, length(obj$hclust$height))), 
                              y = rev(obj$hclust$height)[-1]), 
                   aes_string(x="x", y="y")) +
        geom_bar(stat = "identity") +
        ylab("Height") +
        xlab("Nb of clusters") + theme_bw()
      g[["heights"]] <- g1 
    }
  }
  if("dendrogram" %in% type) { ## Dendrogram
    if(!"merged_clusters" %in% names(obj)) {
      message("-- dendrogram plot is only possible for the aggregation algorithm.")
    } else {
      g2 <- ggdendro::ggdendrogram(obj$hclust, rotate=FALSE, size=2)
      g[["dendrogram"]] <- g2 
    }
  }
  if("criterion" %in% type) { ## Evolution of criterion
    if("merged_clusters" %in% names(obj)) {
      df <- data.frame(x = seq(2, length(obj$hclust$order)), 
                       y = rev(obj$criterion))
    } else {
      df <- data.frame(x = apply(obj$split_clusters, 2, max),
                 y = obj$criterion)
    }
    g5 <- ggplot(df, aes_string(x = "x", y = "y")) +
      geom_point() +
      geom_line() +
      ylab("Criterion value") +
      xlab("Number of clusters") + theme_bw()
    g[["criterion"]] <- g5
  }
  if("weights" %in% type) { ## Evolution of weights
    if(class(obj$weights) == "list") {
      tmp <- cbind(do.call("rbind", obj$weights),
                   unlist(lapply(obj$weights, function(x) rep(nrow(x), nrow(x)))),
                   rep(1:length(obj$weights), 
                       times = unlist(lapply(obj$weights, nrow))))
      if(is.null(mv_names)) {
        colnames(tmp) <- c(paste0("View ", 1:(ncol(tmp)-2)), "Cluster", "iteration")
      } else {
        if(length(mv_names) != ncol(tmp)-2) stop("mv_names must be the same length as the number of views")
        colnames(tmp) <- c(mv_names, "Cluster", "iteration")
      }
      h <- Heatmap(tmp[,1:(ncol(tmp)-2)], split=tmp[,"Cluster"],
              col = viridis::viridis(21, begin=0, end=1), cluster_rows=FALSE,
              heatmap_legend_param = list(title = "weights"))
      print(h)
      } else {
        if("merged_clusters" %in% names(obj)) {
          df <- data.frame(cbind(seq(from = max(as.numeric(obj$hclust$labels)), to = 2),
                                 seq(1, ncol(obj$weights)), t(obj$weights)),
                           row.names=NULL) 
        } else {
          df <- data.frame(cbind(apply(obj$split_clusters, 2, function(x) length(unique(x))),
                                 seq(1, ncol(obj$weights)), t(obj$weights)),
                           row.names=NULL) 
        }
      rownames(df) <- NULL
      colnames(df) <- c("Clusters", "iteration", as.character(seq(1:nrow(obj$weights))))
      df <- tidyr::gather_(df, key_col="view", value_col="value", 
                           gather_cols=c(as.character(seq(1:nrow(obj$weights)))))
      
      if(!is.null(mv_names)) {
        df$view <- factor(df$view)
        if(length(mv_names) != length(levels(df$view))) stop("mv_names must be the same length as the number of views")
        levels(df$view) <- mv_names
      }
      
      g3 <- ggplot(df, aes_string("Clusters", "value")) +
        geom_area(aes_string(fill = "view", group = "view")) +
        ylab("weights") +
        viridis::scale_fill_viridis(discrete=TRUE) + theme_bw()
      if("merged_clusteers" %in% names(obj)) {
        g3 <- g3 + scale_x_reverse()
      }
      g[["weights"]] <- g3 
    }
  }
  if("weights_line" %in% type) {
    if(class(obj$weights) == "list") {
      message("-- weights line plot only available when per-cluster weights are not used.")
    } else {
      if("merged_clusters" %in% names(obj)) {
        df <- data.frame(cbind(seq(from = max(as.numeric(obj$hclust$labels)), to = 2),
                               seq(1, ncol(obj$weights)), t(obj$weights)),
                         row.names=NULL) 
      } else {
        df <- data.frame(cbind(apply(obj$split_clusters, 2, function(x) length(unique(x))),
                               seq(1, ncol(obj$weights)), t(obj$weights)),
                         row.names=NULL) 
      }
      rownames(df) <- NULL
      colnames(df) <- c("Clusters", "iteration", as.character(seq(1:nrow(obj$weights))))
      df <- tidyr::gather_(df, key_col="view", value_col="value", 
                           gather_cols=c(as.character(seq(1:nrow(obj$weights)))))
      
      if(!is.null(mv_names)) {
        df$view <- factor(df$view)
        if(length(mv_names) != length(levels(df$view))) stop("mv_names must be the same length as the number of views")
        levels(df$view) <- mv_names
      }
      
      g4 <- ggplot(df, aes_string("Clusters", "value", group = "view", colour = "view")) +
        geom_point() +
        geom_line(aes_string(lty = "view")) +
        ylab("weights") +
        viridis::scale_color_viridis(discrete=TRUE) + theme_bw()
      if("merged_clusters" %in% names(obj)) {
        g4 <- g4 + scale_x_reverse()
      }
      g[["weights_line"]] <- g4
    }
  }
  if("tree" %in% type & !"merged_clusters" %in% names(obj)) {
    if("merged_clusters" %in% names(obj)) {
      message("-- tree plot is only possible for the splitting algorithm.")
    } else {
      if(tree_type == "all") {
        aux <- obj$split_clusters
      } else if(tree_type == "final_K") {
        aux <- obj$split_clusters
        index <- which(apply(aux, 2, max) <= obj$final_K)
        aux <- obj$split_clusters[,index]
      } else {
        if(!is.numeric(tree_type)) stop("tree_type must be either 'all', 'final_K', or a number.")
        aux <- obj$split_clusters
        index <- which(apply(aux, 2, max) <= tree_type)
        aux <- obj$split_clusters[,index]
      }
      colnames(aux) <- paste0("K", apply(aux, 2, max))
      
      aux <- data.frame(aux)
      aux$color <- 0
      c <- clustree::clustree(aux, prefix = "K", 
                              edge_arrow=arg.user$edge_arrow,
                              scale_node_text=FALSE,
                              node_size_range=c(5,5),
                              edge_width = 0.5,
                              node_colour = "color",
                              node_colour_aggr="mean")  +
        guides(edge_colour = FALSE, edge_alpha = FALSE) +
        theme(legend.position = "none") +
        scale_edge_color_continuous(low = "grey10", high = "grey30") 
      
      cmod <- c
      cmod_data <- cmod$data
      cmod_data$K <- as.numeric(as.character(cmod_data$K))
      for(Kval in cmod_data$K) {
        if(Kval == min(cmod_data$K)) next;
        index  <- which(cmod_data$K == Kval)
        index_prev  <- which(cmod_data$K == Kval-1)
        drop_index <- which(cmod_data$x[index] %in% cmod_data$x[index_prev])
        cmod_data[index[drop_index],]$size <- 0
      }
      cmod_data$K <- factor(cmod_data$K)
      cmod_data$mean_color <- ifelse(cmod_data$size == 0, 0, 1)
      cmod_data$mean_color <- factor(cmod_data$mean_color)
      cmod$data <- cmod_data
      cmod <- cmod +
        scale_color_manual(values = c("white", "#21908CFF"))

      g[["tree"]] <- cmod
    }
  }
  if("probapost_boxplot" %in% type) {
    if(is.null(obj$final_probapost)) {
      message("-- probapost_boxplot is only possible for posterior probabilities from fuzzy aggregation/splitting algorithms")
    } else {
      pbox <- probapost_boxplot(obj$final_probapost)
      g[["probapost_boxplot"]] <- pbox
    }
  }
  
  if("tree_perClusterWeights" %in% type) {
    if("merged_clusters" %in% names(obj)) {
      message("-- tree plot is only possible for the splitting algorithm.")
    } else if(!is.list(obj$weights)) {
      message("-- tree plot is only possible for the splitting algorithm with per-cluster weights.")
    } else {
      ## First plot tree
      if(tree_type == "all") {
        aux <- obj$split_clusters
        index0 <- 1:ncol(aux)
      } else if(tree_type == "final_K") {
        aux <- obj$split_clusters
        index0 <- which(apply(aux, 2, max) <= obj$final_K)
        aux <- obj$split_clusters[,index0]
      } else {
        if(!is.numeric(tree_type)) stop("tree_type must be either 'all', 'final_K', or a number.")
        aux <- obj$split_clusters
        index0 <- which(apply(aux, 2, max) <= tree_type)
        aux <- obj$split_clusters[,index0]
        if(class(aux) == "numeric") stop("No splits at the provided value of tree type.")
      }
      colnames(aux) <- paste0("K", apply(aux, 2, max))
      aux <- data.frame(aux)
      aux$color <- 0
      c <- clustree::clustree(aux, prefix = "K", 
                              edge_arrow=arg.user$edge_arrow,
                              scale_node_text=FALSE,
                              node_size_range=c(5,5),
                              edge_width = 0.5,
                              node_colour = "color",
                              node_colour_aggr="mean")  +
        guides(edge_colour = FALSE, edge_alpha = FALSE) +
        theme(legend.position = "none") +
        scale_edge_color_continuous(low = "grey10", high = "grey30") 
      cmod <- c
      cmod_data <- cmod$data
      cmod_data$K <- as.numeric(as.character(cmod_data$K))
      for(Kval in cmod_data$K) {
        if(Kval == min(cmod_data$K)) next;
        index  <- which(cmod_data$K == Kval)
        index_prev  <- which(cmod_data$K == Kval-1)
        drop_index <- which(cmod_data$x[index] %in% cmod_data$x[index_prev])
        cmod_data[index[drop_index],]$size <- 0
      }
      cmod_data$K <- factor(cmod_data$K)
      cmod_data$mean_color <- ifelse(cmod_data$size == 0, 0, 1)
      cmod_data$mean_color <- factor(cmod_data$mean_color)
      cmod$data <- cmod_data
      cmod <- cmod +
        scale_color_manual(values = c("white", "#21908CFF"))
      
      ## Now plot weights of split clusters
      w <- obj$weights
      names(w) <- levels(cmod$data$K)
      w <- w[index0[-1]]
      for(i in 1:length(w)) {
        tmp <- cmod_data[which(cmod_data$K == names(w)[i]),]
        clus_choice <- tmp[which(tmp$size != 0),]$cluster
        w[[i]] <- w[[i]][clus_choice,]
        rownames(w[[i]]) <- paste0("K=",clus_choice)
      }
      w <- do.call("rbind", w)
      colnames(w) <- 1:ncol(w)
      wdf <- data.frame(w, row.names=NULL, check.names=FALSE)
      wdf$cluster <- rownames(w)
      wdf$level <- as.numeric(rep(unique(cmod$data$y)[-1], each = 2))
      wdf$ymin <- wdf$level + c(-.4, 0)
      wdf$ymax <- wdf$level + c(0, 0.4)
      wdf$ymid <- wdf$level + c(-0.2, 0.2)


      wdf_long <- gather_(wdf, key_col="View", value_col="weight", 
                          gather_cols = c(as.character(1:ncol(w))))
      wdf_long$View <- as.numeric(as.character(wdf_long$View))
      wdf_long$xmin <- wdf_long$View - 0.5
      wdf_long$xmax <- wdf_long$View + 0.5
      wdf_long$Viewname <- paste0("View ", wdf_long$View)
      if(!is.null(mv_names)) {
        if(length(mv_names) != ncol(wdf)-5) stop("mv_names must be the same length as the number of views")
        wdf_long$Viewname <- factor(wdf_long$Viewname)
        levels(wdf_long$Viewname) <- mv_names
      }
      ym <- max(wdf_long$ymax) + 0.5
      wp <- ggplot(wdf_long, aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax")) +
        geom_rect(aes_string(fill="weight")) +
        geom_text(aes_string(label="cluster", x=0, y="ymid")) +
        geom_text(aes_string(label="Viewname", x="View", y=ym)) +
        theme(axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              axis.line = element_blank()) + 
        ylim(c(min(cmod$data$y)-0.5, max(cmod$data$y))) +
        scale_fill_viridis() + 
          theme_void()
      cmod2 <- cowplot::plot_grid(plotlist=list(cmod+ 
                                         ylim(c(min(cmod$data$y)-0.5, max(cmod$data$y))), wp))
      
      g[["tree_perClusterWeights"]] <- cmod2
    }
    
  }
  
  if(length(g) > 1) {
    print(cowplot::plot_grid(plotlist=g))
  } else {
    print(g)
  }
  return(g)
}


## Not exported:
probapost_boxplot <- function(probapost) {
  aux <- data.frame(probamax = apply(probapost, 1, max),
                    lab = as.factor(apply(probapost, 1, which.max)))
  p <- ggplot(aux, aes_string(x = "lab", y = "probamax")) +
    xlab("Cluster") +
    ylab("Maximum posterior probability") +
    geom_boxplot() +
    theme_bw()
    return(p)
}

## Not exported:
probapost_threshold <- function(obj, probapost, threshold = 0.8) {
  aux <- NULL
  for (k in 2:length(obj$criterion)) {
    aux <- c(aux, sum(apply(
      cutreeNewProbapost(obj$hclust, K=k, probapost)$probapost, 1, max) > threshold) * 100 / 
        nrow(probapost))
  }
  p <- ggplot(data.frame(k = 2:length(obj$criterion), aux = aux), aes_string("k", "aux")) +
    geom_line() +
    geom_point() +
    xlab("Number of clusters") +
    ylab("Percent of max posterior prob > threshold")
  print(p)
}


##!!  I have not implemented the following functions:
# PlotClassif <- function(data, classif) {
#   don <- data.frame(classif = as.factor(classif), data = data)
#   ggpairs(don,
#           mapping = aes(color = classif),
#           columns = colnames(don)[-1])
# }
#############################################
##   exploitation graphique de la sortie splittingClustersbis
#############################################
# PlotResSplittingbis <-
#   function(Res,
#            weights_step = FALSE,
#            weights_clust = FALSE) {
#     
#     # plot des poids
#     # graphes de l'évolution des poids
#     if (weights_step != FALSE) {
#       for (h in 1:length(weights_step)) {
#         j <- weights_step[h]
#         #!!!!!!!
#         A <- melt(Res$weights[[j]])
#         A[, 1] <- as.factor(A[, 1])
#         A[, 2] <- as.factor(A[, 2])
#         colnames(A) = c("clust", "mv", "weights")
#         g <- ggplot(A, aes(fill = mv, y = weights, x = clust)) +
#           geom_bar(stat = "identity") +
#           ggtitle(paste("step ", j - 1, " - split ", Res$ksplit[j]))
#         print(g)
#       }
#     }
#     
#     if (weights_clust != FALSE) {
#       #!!!!!!
#       A <- melt(Res$weights[[length(Res$weights)]])
#       colnames(A) <- c("clust", "mv", "weights")
#       I <- NULL
#       for (h in 1:length(weights_clust))
#         I <- c(I, which(A$clust == weights_clust[h]))
#       A[, 1] <- as.factor(A[, 1])
#       A[, 2] <- as.factor(A[, 2])
#       g2 <- ggplot(A[I, ], aes(fill = mv, y = weights, x = clust)) +
#         geom_bar(stat = "identity")
#       print(g2)
#     }
#     # plot du splitting   ????
#   }
