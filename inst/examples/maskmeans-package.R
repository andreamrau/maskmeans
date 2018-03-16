# @example /inst/examples/maskmeans-package.R

## Simulate data
set.seed(12345)
sim_1 <- mv_simulate(type = "D1")
sim_2 <- mv_simulate(type = "D2")
sim_3 <- mv_simulate(type = "D3")
sim_4 <- mv_simulate(type = "D4")
sim_5 <- mv_simulate(type = "D5")
sim_6a <- mv_simulate(type = "D6")
sim_6b <- mv_simulate(type = "D6", delta=7, n=200, K=5, sigma=0.5)

X <- sim_6a$data
mv <- c(2,2,2,1,1,2)
gamma <- 2 
Xlist <- list(X[,1:2], X[,3:4], X[,5:6], matrix(X[,7], ncol=1), matrix(X[,8], ncol=1), X[,9:10])
X_scale <- scaleview(X, mv)


#-------------------------------------------------------------------
## Double-check that all functions provide the same result as before
#-------------------------------------------------------------------


#**************************************
## Test 1: hard clustering aggregation
#**************************************
cluster_init <- sim_6a$labels[,1]
set.seed(12345)
hard_agglom <- maskmeans(mv_data=X, mv=mv, clustering_init=cluster_init, 
                         type = "aggregation", gamma=gamma) 
  
set.seed(12345)
hard_agglom_old <- hmv1(X_scale, mv=mv, gamma=gamma, cluster.init=cluster_init, 
                        weightsopt = TRUE)
  
all.equal(hard_agglom$weights, hard_agglom_old$weights,
          check.attributes = FALSE)              
all.equal(hard_agglom$criterion, hard_agglom_old$CRIT)               
all.equal(hard_agglom$merged_clusters, hard_agglom_old$merge)

#**************************************
## Test 2: fuzzy clustering aggregation
#**************************************
proba_init <- matrix(runif(nrow(X)*5), ncol=5)
proba_init <- proba_init / rowSums(proba_init)
set.seed(12345)
fuzzy_agglom <- maskmeans(mv_data=X, mv=mv, clustering_init=proba_init, 
                          type = "aggregation", gamma=gamma) 
set.seed(12345)
fuzzy_agglom_old <- hmvprobapost(X_scale, mv=mv, gamma=gamma, probapost.init=proba_init)

all.equal(fuzzy_agglom$weights, fuzzy_agglom_old$weights, check.attributes=FALSE)
all.equal(fuzzy_agglom$criterion, fuzzy_agglom_old$CRIT)            ## THESE ARE NOT EQUAL!!!
all.equal(fuzzy_agglom$merged_clusters, fuzzy_agglom_old$merge)
  
#**************************************  
## Test 3: hard clustering splitting
#**************************************
set.seed(12345)
hard_split <- maskmeans(mv_data=X, mv=mv, clustering_init=cluster_init, type = "splitting", Kmax=20,
                        perCluster_mv_weights = FALSE)  

set.seed(12345)
hard_split_old <- splittingClusters(X=X_scale, mv=mv, gamma=gamma, Kmax=20, cluster.init=cluster_init,
                             weightsopt = TRUE, testkmeans = TRUE) 

all.equal(hard_split$weights, hard_split_old$weights, check.attributes = FALSE) ## THESE ARE NOT EQUAL!!!
all.equal(hard_split$criterion, hard_split_old$CRIT)                            ## THESE ARE NOT EQUAL!!!
all.equal(hard_split$split_clusters, hard_split_old$clustersplithist)           ## THESE ARE NOT EQUAL!!!
all.equal(hard_split$ksplit, hard_split_old$ksplit)                             ## THESE ARE NOT EQUAL!!!
all.equal(hard_split$withinss, hard_split_old$withinss)                         ## THESE ARE NOT EQUAL!!!


#**************************************
## Test 4: hard clustering splitting with per-weights
#**************************************
hard_split_perCluster <- maskmeans(mv_data=X, mv=mv, clustering_init=cluster_init, type = "splitting", 
                                   Kmax=20, perCluster_mv_weights=TRUE, gamma=1)  ## ERROR

set.seed(12345)
hard_split_old_perCluster <- splittingClustersbis(X=X_scale, mv=mv, gamma=1, 
                                                  Kmax=20, cluster.init=cluster_init) 

mapply(all.equal, hard_split_perCluster$weights, hard_split_old_perCluster$weights, 
       check.attributes = FALSE)               ## THESE ARE NOT EQUAL!!!
all.equal(hard_split_perCluster$criterion, hard_split_old_perCluster$CRIT)                        ## THESE ARE NOT EQUAL!!!
all.equal(hard_split_perCluster$split_clusters, hard_split_old_perCluster$clustersplithist)       ## THESE ARE NOT EQUAL!!!
all.equal(hard_split_perCluster$ksplit, hard_split_old_perCluster$ksplit)                        
all.equal(hard_split_perCluster$withinss, hard_split_old_perCluster$withinss)                     ## THESE ARE NOT EQUAL!!!


#**************************************
## Other testing idea: using Xlist instead
#**************************************