\name{GADAG-package}
\alias{GADAG}
\docType{package}
\title{
\packageTitle{GADAG}
}
\description{
\packageDescription{GADAG}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{GADAG}
}
\author{
\packageAuthor{GADAG}

Maintainer: \packageMaintainer{GADAG}
}
\keyword{Directed Acyclic Graphs - Convex program - Optimization}
\examples{
  #############################################################
  # Loading toy data
  #############################################################
  data(toy_data)
  # toy_data is a list of two matrices corresponding to a "star"
  # DAG (node 1 activates all other nodes):
  # - toy_data$X is a 100x10 design matrix
  # - toy_data$G is the 10x10 adjacency matrix (ground trough)

  #############################################################
  # Running GADAG
  #############################################################
  # Simple run, with only the penalty term specified
  GADAG_results <- GADAG_Run(X=toy_data$X, lambda=0.1)

  # Expensive run with many evaluations if we refine the
  # termination conditions
  \dontrun{
  n.gen <- 1e10 # we allow a very large number of iterations
  tol.Shannon <- 1e-10 # the entropy of Shannon of the population
                       # has to be very small
  max.eval <- 1e10 # we allow a very large number of nested
                   # evaluation
  GADAG_results <- GADAG_Run(X=toy_data$X, lambda=0.1,
       GADAGcontrol=list(n.gen=n.gen, tol.Shannon=tol.Shannon, max.eval=max.eval))
  }

  # Expensive run if we also increase the population size
  \dontrun{
  pop.size <- 5*ncol(toy_data$G) # this is usually a good
                                 # population size
  GADAG_results <- GADAG_Run(X=toy_data$X, lambda=0.1,
      GADAGcontrol=list(pop.size=pop.size))
  }

  # You can have more information about the evolution of the
  # algorithm by turning return.level on
  \dontrun{
  return.level <- 1
  GADAG_results <- GADAG_Run(X=toy_data$X, lambda=0.1, return.level = return.level)
  }
}
