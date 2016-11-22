\name{ockc}
\alias{ockc}
\alias{ockc-class}
\alias{show,ockc-method}
\title{
Order Contrained Solutions in k-Means Clustering
}
\description{
  Calculates an order constrained clustering solution (default k-means) on a
  data matrix.
}
\usage{
ockc(x, k, family = kccaFamily("kmeans"), order = NULL, control = NULL,
     save.data = FALSE, multicore = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A numeric matrix of data.
}
  \item{k}{
    An integer vector of number of clusters. For each element of k a clustering
    solution is computed (reusage of intermediate results makes this more
    efficient than individual calls of ockc). 
}
  \item{family}{
    Object of class \code{kccaFamily}.
}
  \item{order}{
    Order restriction of \code{x}. If \code{NULL} an order is caluclated with
    \code{\link[seriation]{seriate}} from package \code{seriation}
}
  \item{control}{
    An object of class \code{flexclustControl}.
}
  \item{save.data}{
    Save a copy of \code{x} in the return object?
}
  \item{multicore}{
      Use parallelization, if available. For examples and additional
      documentation see \code{\link{bootFlexclust}}.
}
  \item{\dots}{
    Additional options for \code{\link[seriation]{seriate}} for \code{order} calculation.
}
}
\references{
  Steinley, D. and Hubert, L. (2008).
  Order-Constrained Solutions in K-Means Clustering: Even Better Than Being
  Globally Optimal.
  Psychometrika, 73 (4), pp. 647-664.
}
\author{
  Sebastian Krey, Friedrich Leisch, Sebastian Hoffmeister
}
\seealso{
  \code{\link{kcca}}
}
\examples{
x <- rbind(cbind(rnorm(10, mean=0), rnorm(10, mean=0,), rnorm(10, mean=0)),
           cbind(rnorm(10, mean=10), rnorm(10, mean=10), rnorm(10, mean=0)),
           cbind(rnorm(10, mean=10), rnorm(10, mean=0), rnorm(10, mean=10)),
           cbind(rnorm(10, mean=10), rnorm(10, mean=10), rnorm(10, mean=10))
           )

res <- ockc(x, k=4, nboot=4, order=c(1:10, 21:40, 11:20))
res
}
\keyword{cluster}
