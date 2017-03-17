#################################################################
##' @title Compute the fitness of a potential solution
##' @description Internal function of the genetic algorithm that evaluates the fitness (penalized log-likelihood) of a potential solution, given as a pair of a permutation (P) and a triangular matrix (T).
##' @usage fitness(P,X,T,lambda)
##' @param P A permutation from [1,p] in a matrix form.
##' @param X Design matrix, with samples (n) in rows and variables (p) in columns.
##' @param T A pxp triangular matrix.
##' @param lambda Parameter of penalization (>0).
##' @return A numeric value corresponding to the fitness of the potential solution.
##' @author \packageAuthor{GADAG}
##' @seealso \code{\link{chrom}}
##' @examples
##'  ########################################################
##'  # Loading toy data
##'  ########################################################
##'  data(toy_data)
##'
##'  p <- ncol(toy_data$X)
##'  P <- sample(p)
##'  P <- chrom(P)
##'  T <- matrix(rnorm(p),p,p)
##'  T[upper.tri(T)] <- 0
##'  T <- T-diag(diag(T))
##'
##'  ########################################################
##'  # Computing the fitness of a potential solution
##'  ########################################################
##'  Fitness <- fitness(P=P, X=toy_data$X, T=T, lambda=0.1)

fitness <- function(P,X,T,lambda){
  # INPUTS
  # P: permutation matrix (p*p)
  # X: observation matrix (n*p)
  # T: triangular matrix (p*p)
  # lambda: penalization term (scalar)

  # OUTPUTS
  # f: fitness value

  n <- dim(X)[1]
  f <- (1/n)*norm(X - X %*% P %*% T %*% t(P),'f')^2 + lambda*(sum(abs(T)))
  return(f)
}
