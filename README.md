# maskmeans: Multi-view aggregation/splitting K-means clustering algorithm

The *maskmeans* package can be installed as follows:

```
library(devtools)
install_github("andreamrau/maskmeans")
library(maskmeans)
```
The primary functions of this package are as follows:

- `maskmeans`, which itself calls one of the two following functions:
    * `mv_aggregation`
    * `mv_splitting`: Note that this algorithm allows either fixed multi-view weights across clusters (`perCluster_mv_weights = FALSE`) or per-cluster multi-view weights (`perCluster_mv_weights = TRUE`).
- `maskmeans_cutree`:  cut an aggregation tree for a specified number of clusters
- `mv_simulate` to simulate data types `"D1"`, ... `"D6"`

There are also two plotting functions:

- `mv_plot`, to provide a plotting overview of multi-view data. Univariate views are plotted as density plots, bivariate views as scatterplots, and multivariate views as scatterplots of the first two principal components. A vector of cluster labels can be added to color the points according to a unique partition (e.g., the labels of the first view).
- `maskmeans_plot`, to plot results of the `maskmeans` function. Plot types provided through this function include `type =  c("dendrogram", "heights", "weights_line", "weights", "criterion", "tree")`

See the package [vignette](https://github.com/andreamrau/maskmeans/blob/master/doc/maskmeans.html) for a full example.