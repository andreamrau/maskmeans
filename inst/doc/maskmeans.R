## ---- include=FALSE------------------------------------------------------
knitr::opts_chunk$set(dev='png')

## ------------------------------------------------------------------------
library(fclust)
library(mclust)

## ------------------------------------------------------------------------
set.seed(12345)
sim <- mv_simulate(type = "D6", delta=5, n=200, K=9, sigma=1.5)
X <- sim$data         ## Simulated data
mv <- c(2,2,2,1,1,2)  ## Size of each view
all.equal(sum(mv), ncol(X))

## ---- fig.width = 10, fig.height=6, dpi=36, out.width="600px", out.height="600px"----
mv_plot(mv_data=X, mv=mv, labels=sim$labels[,1])

## ------------------------------------------------------------------------
aggregate_hard_init <- kmeans(X[,1:2], centers=15)$cluster
aggregate_fuzzy_init <- FKM(X[,1:2], k=15, maxit=10, RS=10)$U
split_hard_init <- kmeans(X[,1:2], centers=2)$cluster

## ------------------------------------------------------------------------
hard_aggreg <- maskmeans(mv_data=X, mv=mv, clustering_init=aggregate_hard_init, 
                         type = "aggregation", gamma=2) 

## ---- fig.height = 8, fig.width=10, dpi=36, out.width="600px", out.height="600px"----
hard_aggreg$final_K
p <- maskmeans_plot(hard_aggreg, mv_data = X)

## ---- fig.width = 10, fig.height=6, dpi=36, out.width="600px", out.height="600px"----
# mv_plot(mv_data=X, mv=mv, labels=hard_aggreg$final_classification)
adjustedRandIndex(hard_aggreg$final_classification, sim$labels[,1])
adjustedRandIndex(aggregate_hard_init, sim$labels[,1])

## ------------------------------------------------------------------------
fuzzy_aggreg <- maskmeans(mv_data=X, mv=mv, clustering_init=aggregate_fuzzy_init, 
                         type = "aggregation", gamma=2) 

## ---- fig.height = 8, fig.width=10, dpi=36, out.width="600px", out.height="600px"----
fuzzy_aggreg$final_K
p <- maskmeans_plot(fuzzy_aggreg, mv_data = X)

## ---- fig.width = 10, fig.height=6, dpi=36, out.width="600px", out.height="600px"----
# mv_plot(mv_data=X, mv=mv, labels=fuzzy_aggreg$final_classification)
adjustedRandIndex(fuzzy_aggreg$final_classification, sim$labels[,1])
adjustedRandIndex(apply(aggregate_fuzzy_init, 1, which.max), sim$labels[,1])

## ------------------------------------------------------------------------
hard_split <- maskmeans(mv_data=X, mv=mv, clustering_init=split_hard_init, Kmax=20,
                         type = "splitting", gamma=2, perCluster_mv_weights=FALSE) 

## ---- fig.height = 10, fig.width=10, dpi=36, out.width="600px", out.height="600px"----
hard_split$final_K
p <- maskmeans_plot(hard_split, mv_data = X)

## ---- fig.width = 10, fig.height=6, dpi=36, out.width="600px", out.height="600px"----
# mv_plot(mv_data=X, mv=mv, labels=hard_split$final_classification)
adjustedRandIndex(hard_split$final_classification, sim$labels[,1])
adjustedRandIndex(split_hard_init, sim$labels[,1])

## ------------------------------------------------------------------------
hard_split_perCluster <- maskmeans(mv_data=X, mv=mv, clustering_init=split_hard_init, Kmax=20,
                         type = "splitting", gamma=2, perCluster_mv_weights=TRUE) 

## ---- fig.height = 8, fig.width=10, dpi=36, out.width="600px", out.height="600px"----
hard_split_perCluster$final_K
p <- maskmeans_plot(hard_split_perCluster, mv_data = X, type = c("criterion", "tree"))

## ----  fig.height=8, fig.width=4, dpi=36, out.width="600px", out.height="600px"----
p <- maskmeans_plot(hard_split_perCluster, mv_data = X, type = c("weights"))

## ---- fig.width = 10, fig.height=6, dpi=36, out.width="600px", out.height="600px"----
# v_plot(mv_data=X, mv=mv, labels=hard_split_perCluster$final_classification)
adjustedRandIndex(hard_split_perCluster$final_classification, sim$labels[,1])
adjustedRandIndex(split_hard_init, sim$labels[,1])

## ------------------------------------------------------------------------
Xlist <- list(X[,1:2], X[,3:4], X[,5:6], matrix(X[,7], ncol=1), 
              matrix(X[,8], ncol=1), X[,9:10])
hard_aggreg_list <- maskmeans(mv_data=Xlist, 
                              clustering_init=aggregate_hard_init, 
                              type = "aggregation", gamma=gamma) 

## ------------------------------------------------------------------------
sessionInfo()

