\name{gradientdescent}
\alias{gradientdescent}
\title{Perform a gradient descent}
\description{
Internal function of the genetic algorithm, similar to a gradient descent algorithm, that computes the best triangular matrix associated to a permutation. This triangular matrix minimizes the penalized log-likelihood of the DAG learning problem, P being fixed.
}
\usage{
gradientdescent(P, n, XtX, L, lambda, maxite, tolobj)
}
\arguments{
  \item{P}{
Permutation from [1,p] in a matrix form.
}
\item{n}{
Number of samples of the design matrix X (>0).
}
  \item{XtX}{
Cross-product of X. Should be computed (crossprod(X)) before running the gradientdescent() function.
}
\item{L}{
Lispchitz constant (>0). Internally computed from X in the evalutation() function.
}
  \item{lambda}{
Parameter of penalization (>0).
}
  \item{maxite}{
Maximal number of iterations for the inner optimization (>0).
}
  \item{tolobj}{
Tolerance for the inner optimization, i.e. the gradient descent (>0).
}
}
\value{
A pxp triangular matrix T in a vector form.
}
\details{
  These functions are written using the Armadillo C++ classes.
}
\author{
  \packageAuthor{GADAG}
}
\seealso{
	\code{\link{chrom}}, \code{\link{evaluation}}, \code{\link{GADAG}}, \code{\link{GADAG_Run}}.
}
\examples{
########################################################
# Loading toy data
########################################################
data(toy_data)

p <- ncol(toy_data$X)
n <- nrow(toy_data$X)
P <- sample(p)
P <- chrom(P)

gradientdescent(P=P, n=n, XtX=crossprod(toy_data$X), L=1, lambda=0.1, maxite=50, tolobj=1e-6)
}
