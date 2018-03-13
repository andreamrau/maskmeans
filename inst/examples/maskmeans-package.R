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
cluster_init <- sim_6a$labels[,1]
gamma=2 