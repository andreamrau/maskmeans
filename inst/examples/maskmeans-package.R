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

## Double-check that all functions provide the same result as before
#-------------------------------------------------------------------

## Test 1: hard clustering agglomeration
cluster_init <- sim_6a$labels[,1]
hard_agglom <- maskmeans(mv_data=X, mv=mv, clustering_init=cluster_init, type = "agglomeration") 
  
## Test 2: fuzzy clustering agglomeration
proba_init <- matrix(runif(nrow(X)*5), ncol=5)
proba_init <- proba_init / rowSums(proba_init)
fuzzy_agglom <- maskmeans(mv_data=X, mv=mv, clustering_init=proba_init, type = "agglomeration") 

## Test 3: hard clustering splitting
hard_split <- maskmeans(mv_data=X, mv=mv, clustering_init=cluster_init, type = "splitting", Kmax=20) 

## Test 4: hard clustering splitting with per-weights
hard_split_perClsuter <- maskmeans(mv_data=X, mv=mv, clustering_initcluster_init, type = "splitting", 
                                   Kmax=20, perCluster_mv_weights=TRUE) 

