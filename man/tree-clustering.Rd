\name{tree.clustering}
\alias{tree.clustering}
\title{Tree Clustering}
\description{
  Perform tree clustering of a data frame.
}
\usage{
tree.clustering(
    ds,
    nsamples=1,
    burnin=100,
    maxiter=NULL,
    limits=NULL,
    fragment.size=NULL,
    max.segments=0,
    nruns=1
)
}
\arguments{
  \item{ds}{Data frame containing the data the be clustered. Only numeric
  columns are supported at the moment.}
  \item{nsamples}{Number of samples to generate using the Metropolis-Hastings
  algorithm.}
  \item{burnin}{Burn-in period, i.e. the number of initial samples to be
  ignored.}
  \item{maxiter}{Maximum number of iterations of the Metropolis-Hastings
  algorithm.  Use \code{maxiter} to ensure \code{tree.clustering} finishes
  in a reasonable time even when the accept rate is low. When 0, the number
  of iterations is unlimited. Default: \code{(nsamples + burnin)*1000}.
  }
  \item{limits}{Parameter space bounds. When NULL, defaults to the range of data
  in \code{ds}.}
  \item{fragment.size}{Vector specifying fragment size of parameters.
  When set, segment ranges are constrained to be a multiple of
  fragment size in their respective parameters. The
  vector must have the same length as the number of columns in \code{ds}.
  }
  \item{max.segments}{Maximum number of segments. When 0, the number
  of segments is unlimited.}
  \item{nruns}{Number of runs. The algorithm is run \code{nruns} times,
  and segmentations are concatenated in the result.}
}
\details{
  This function performs clustering by constructing a tree approximating
  the probability density function from which data has been drawn.
  The tree partitions the parameter space hierarchically into a number of
  rectangular segments, on which the probability density function is
  assumed to have a uniform distribution.
  Possible trees are sampled using the Metropolis-Hastings algorithm.
}
\value{
  \code{tree.clustering} returns a list of segmentations from the MCMC sampler.
  Segmentation is a list with the following elements:

  \item{likelihood}{Segmentation likelihood, i.e. the likelihood of the
  data being drawn from the segmentation.}
  \item{segments}{List of \emph{segments}.}

  \emph{Segment} is a list with the following elements:

  \item{NX}{Number of data elements in the segment.}
  \item{V}{Volume of the segment calculated by multiplying the segment lengths
  over each parameters.}
  \item{ranges}{List of segment ranges over parameters. Range
  is a vector of two values: the lower and upper bound (resp.)
    of the segment in the respective parameter.}
}
\examples{
data(trees)
tree.clustering(trees, fragment.size=c(0.1,1,0.1))
}
\keyword{file}
