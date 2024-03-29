% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Splitting.R
\name{mv_splitting}
\alias{mv_splitting}
\title{Splitting of hard or soft clusters based on multi-view data}
\usage{
mv_splitting(
  X,
  mv,
  clustering_init,
  Kmax,
  gamma = 2,
  use_mv_weights = TRUE,
  delta = 2,
  perCluster_mv_weights = TRUE,
  verbose = TRUE,
  parallel = TRUE,
  BPPARAM = bpparam()
)
}
\arguments{
\item{X}{Matrix of multi-view data, where the first view corresponds to the 
principal data used to obtain the partition or soft clustering in \code{cluster_init}}

\item{mv}{(Optional unless \code{X} is a matrix.) If \code{X} is a matrix, vector 
corresponding to the size of each data view. 
The sum of \code{mv} should correspond to the number of columns in \code{X}.}

\item{clustering_init}{Either a vector of available cluster labels (for hard clustering) 
or a matrix
of soft classification labels (for soft clustering)}

\item{Kmax}{Maximum number of clusters for splitting}

\item{gamma}{Parameter that controls the distribution of view weights. Default value is 2.}

\item{use_mv_weights}{If \code{TRUE}, run algorithm in weighted multi-view mode; 
if FALSE, the
weight for each view is set to be equal. This option is only used for hard clustering.}

\item{delta}{Parameter that controls the weights on the soft classifications pi(i,k)}

\item{perCluster_mv_weights}{If \code{TRUE}, use cluster-specific multi-view weights.
Otherwise use classic multi-view weights.}

\item{verbose}{If \code{TRUE}, provide verbose output}

\item{parallel}{If \code{FALSE}, no parallelization. If \code{TRUE}, parallel
execution using BiocParallel (see next argument \code{BPPARAM}) for the soft splitting algorithm. A note on running
in parallel using BiocParallel: it may be advantageous to remove large, unneeded objects
from the current R environment before calling the function, as it is possible that R's
internal garbage collection will copy these files while running on worker nodes.}

\item{BPPARAM}{Optional parameter object passed internally to \code{bplapply} when
\code{parallel=TRUE}. If not specified, the parameters last registered with \code{register}
will be used.}
}
\value{
\item{split_clusters }{Matrix providing the history of each cluster splitting at 
each iteration of the algorithm}
\item{weights }{Matrix of dimension \code{v} x \code{niterations}, where \code{v} 
is the number of
views and \code{niterations} is the number of successive splits, providing the 
multi-view weights}
\item{criterion }{Value taken on by the splitting criterion at each iteration}
\item{withnss }{The within sum-of-squares for each cluster at the last iteration}
\item{ksplit }{Vector identifying which cluster was split at each iteration of the 
algorithm}
\item{all_probapost }{List of conditional probabilities for each split for the 
soft algorithm}
}
\description{
Splitting of hard or soft clusters based on multi-view data
}
