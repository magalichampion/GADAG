#################################################################
#### evaluating a solution ####
##' @keywords internal
# @title Evaluate the fitness of a population
# @description Internal function of the genetic algorithm that evaluates the fitness (penalized log-likelihood) of a population of permutations. It internally computes the best triangular matrix associated to each permutation.
# @param Pop Population of permutations from [1,p] (output from create.population() function for example).
# @param X Design matrix, with samples (n) in rows and variables (p) in columns.
# @param XtX Cross-product of X. Should be computed (crossprod(X)) before running the evaluation() function.
# @param lambda Parameter of penalization (>0).
# @param Gradcontrol A list containing the parameters for controlling the inner optimization, i.e. the gradient descent
# \itemize{
# \item{tol.obj.inner}{ tolerance (>0),}
# \item{max.ite.inner}{ maximum number of iterations (>0).}
# }
# @rawNamespace export(evaluation)
# @seealso \code{\link{GADAG}}, \code{\link{GADAG_Run}}.
# @param ncores Number of cores (>1, depending on your computer).
# @return A list with the following elements:
# \itemize{
# \item{Tpop}{ Matrix with pxp columns, each column corresponding to the best triangular matrix (in a vector form) associated to each permutation of the population.}
# \item{f}{ Fitness of the population.}
# }
# @author \packageAuthor{GADAG}
# @examples
#  ########################################################
#  # Loading toy data
#  ########################################################
#  data(toy_data)
#
#  ########################################################
#  # Creating a population of permutations
#  ########################################################
#  Population <- create.population(p=ncol(toy_data$X),pop.size=20)
#
#  ########################################################
#  # Evaluating the fitness of the population
#  ########################################################
#  Evaluation <- evaluation(Pop=Population,X=toy_data$X,
#                           XtX=crossprod(toy_data$X),lambda=0.1)

evaluation <- function(Pop, X, XtX, lambda, Gradcontrol = list(tol.obj=1e-6, max.ite=50), ncores=1){
  # Pop: population of permutations (pop.size*p)
  # X: observation matrix (n*p)
  # XtX: t(X)%*%X matrix, precomputed for speed
  # lambda: penalization term
  # tol.obj: tolerance for the gradient descent (on the norm of T)
  # max.ite: maximum number of iterations for the gradient descent
  #
  # OUTPUTS
  # List of two with
  # Tpop: optimal T values for the population (pop.size*p^2 matrix, one row for each individual)
  # f: fitness of the population
  tol.obj <- Gradcontrol$tol.obj
  max.ite <- Gradcontrol$max

  n <- dim(X)[1]
  p <- dim(X)[2]
  L <- (2/n) * norm(XtX,'f') # Lispchitz constant
  if (max.ite < 0){
    stop("max.ite should be non-negative.")
  }
  if (is.vector(Pop)==TRUE){
    Pop <- t(as.matrix(Pop,ncol=length(Pop),nrow=1))
  }
  pop.size <- dim(Pop)[1]
  if (ncol(Pop)!=ncol(X)){
    stop("The number of variables of Pop does not correspond to the number of variables of X.")
  }
  my.gradientdescent <- function(i){
    gradientdescent(P=chrom(Pop[i,]), n=n, XtX=XtX, L=L, lambda=lambda, maxite=max.ite, tolobj=tol.obj)
  }

  Tpop <- matrix(unlist(mclapply(X=1:pop.size, FUN=my.gradientdescent, mc.cores = ncores, mc.preschedule=FALSE)), nrow=pop.size, byrow=TRUE)

  my.fitness <- function(i){
    fitness(P=chrom(Pop[i,]), X, matrix(Tpop[i,], p, p), lambda=lambda)
  }

  f <- unlist(mclapply(X=1:pop.size, FUN=my.fitness))

  return(list(Tpop = Tpop, f = f))
}
