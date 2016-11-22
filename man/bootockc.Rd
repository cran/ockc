\name{bootockc}
\alias{bootockc}
\title{
  Bootstrap Order Constrained k-means Clustering
}
\description{
  Runs ockc for different numbers of clusters on bootstrap replica of the
  original data (maintaining the supplied order) and returns
  corresponding cluster assignments, centroids and Rand indices
  comparing pairs of partitions.
}
\usage{
bootockc(x, k, nboot = 100, order = NULL, correct = TRUE, seed = NULL,
         multicore = TRUE, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, k, \ldots}{
    Passed to \code{\link{ockc}}
}
  \item{nboot}{
    Number of bootstrap pairs (maintaining order).
}
  \item{order}{
    Order restriction of \code{x}. If \code{NULL} an initial run of \code{ockc}
    with \code{order=NULL} is run to calculate an order with \code{seriate} from
    package \code{seriation}
}
  \item{correct}{
    Logical, correct the index for agreement by chance?
}
  \item{seed}{
    If not \code{NULL}, a call to \code{set.seed()} is made
    before any clustering is done.
}
  \item{multicore}{
      Use parallelization, if available. For examples and additional
      documentation see \code{\link{bootFlexclust}}.
}
  \item{verbose}{
      Logical, show progress information during computations.
      Ignored if \code{multicore=TRUE}.
}
}
\value{
  Returns an object of class \code{"bootFlexclust"}.
}
\author{
  Sebastian Krey
}

\seealso{
  \code{\link{ockc}}, \code{\link{bootFlexclust}}, \code{\link{stepFlexclust}}
}
\examples{
x <- rbind(cbind(rnorm(10, mean=0), rnorm(10, mean=0,), rnorm(10, mean=0)),
           cbind(rnorm(10, mean=10), rnorm(10, mean=10), rnorm(10, mean=0)),
           cbind(rnorm(10, mean=10), rnorm(10, mean=0), rnorm(10, mean=10)),
           cbind(rnorm(10, mean=10), rnorm(10, mean=10), rnorm(10, mean=10))
           )

bockc <- bootockc(x, 2:4, nboot=4, order=c(1:10, 21:40, 11:20),
                  multicore=FALSE, verbose=FALSE)
bockc
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}% use one of  RShowDoc("KEYWORDS")
