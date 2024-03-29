# maskmeans: Multi-view aggregation/splitting K-means clustering algorithm



[![DOI](https://zenodo.org/badge/122631797.svg)](https://zenodo.org/badge/latestdoi/122631797)


The *maskmeans* package can be installed as follows:

```
library(devtools)
devtools::install_github("andreamrau/maskmeans")
library(maskmeans)
```
To also build the vignette, you can use the following (note that this will require the installation of some extra packages):
```
devtools::install_github("andreamrau/maskmeans", build_vignettes=TRUE)
library(maskmeans)
```

*maskmeans* incorporates algorithms for aggregating or splitting an existing hard or soft classification using multi-view data. The primary functions of this package are as follows:

- `maskmeans`, which itself calls one of the two following functions:
    * `mv_aggregation`
    * `mv_splitting`: Note that this algorithm allows either fixed multi-view weights across clusters (`perCluster_mv_weights = FALSE`) or per-cluster multi-view weights (`perCluster_mv_weights = TRUE`).
- `maskmeans_cutree`:  cut an aggregation tree for a specified number of clusters
- `mv_simulate` to simulate data types `"D1"`, ... `"D6"`

There are also two main plotting functions:

- `mv_plot`, to provide a plotting overview of multi-view data. Univariate views are plotted as density plots, bivariate views as scatterplots, and multivariate views as scatterplots of the first two principal components. A vector of cluster labels can be added to color the points according to a unique partition (e.g., the labels of the first view).
- `maskmeans_plot`, to plot results of the `maskmeans` function. Plot types provided through this function include `type =  c("dendrogram", "heights", "weights_line", "weights", "criterion", "tree")`

See the package [vignette](https://htmlpreview.github.io/?https://github.com/andreamrau/maskmeans/blob/master/doc/maskmeans.html) for a full example and description. If the package was installed with the built vignette above, it may be accessed after loading the package via  `vignette("maskmeans")`.

If you use *maskmeans* in your research, please cite our work:

- Godichon-Baggioni, A., Maugis-Rabusseau, C. and Rau, A. (2020) Multi-view cluster aggregation and splitting, with an application to multi-omic breast cancer data. Annals of Applied Statistics, 14:2, 752-767. [DOI: 10.1214/19-AOAS1317](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-14/issue-2/Multiview-cluster-aggregation-and-splitting-with-an-application-to-multiomic/10.1214/19-AOAS1317.short)