#########################################
### Tuning the parameter of penalization by CV
##' @title Cross-validation for GADAG
##' @description Function to run k-fold cross-validation for GADAG to optimally tune the parameter of penalization.
##' @details The function runs \code{GADAG} \code{n.folds} times for each of the tested lambdas to compute the best solution associated to each omitted fold. The error is accumulated and the averaged error over the folds is computed. The best lambda \code{lambda.min} corresponds to the one that minimizes the error.
##' @param X Design matrix, with samples (n) in rows and variables (p) in columns.
##' @param Lambdas Optional user-supplied lambda sequence. Default is null and GADAG chooses its own sequence.
##' @param n.folds Number of folds for cross-validation (10 by default). Can be as large as the sample size (leave-one-out cross validation) but not recommended for large data sets.
##' @param threshold Thresholding value for the estimated edges.
##' @param GADAG.control A list containing parameters for controlling GADAG (termination conditions and inherent parameters of the Genetic Algortihm).
##' Some parameters (n.gen, max.eval and pop.size) are particularly critical for reducing the computational time.
##' \itemize{
##' \item{\code{n.gen}}{ maximal number of population generations (>0),}
##' \item{\code{pop.size}}{ initial population size for the genetic algorithm (>0),}
##' \item{\code{max.eval}}{ overall maximal number of calls of the evaluation function (>0, should be of the order of \code{n.gen}*\code{pop.size}),}
##' \item{\code{tol.Shannon}}{ threshold for the Shannon entropy (>0),}
##' \item{\code{p.xo}}{ crossover probability of the genetic algorithm (between 0 and 1),}
##' \item{\code{p.mut}}{ mutation probability of the genetic algorithm (between 0 and 1).}
##' }
##' @rawNamespace export(GADAG_CV)
##' @rawNamespace import(igraph)
##' @rawNamespace import(MASS)
##' @rawNamespace import(Rcpp)
##' @rawNamespace import(parallel)
##' @rawNamespace import(cvTools)
##' @rawNamespace importFrom(Rcpp, evalCpp)
##' @rawNamespace useDynLib(GADAG)
##' @param grad.control A list containing the parameters for controlling the inner optimization, i.e. the gradient descent.
##' \itemize{
##' \item{\code{tol.obj.inner}}{ tolerance (>0),}
##' \item{\code{max.ite.inner}}{ maximum number of iterations (>0).}
##' }
##' @param ncores Number of cores (>0, depending on your computer).
##' @param plot.CV If 1, plots the averaged cross validation error given sequence of lambdas.
##' @return A list with the following elements:
##' \itemize{
##' \item{\code{lambda.min}}{ Value of \code{lambda} that minimizes the averaged cross validation error \code{error.CV}.}
##' \item{\code{lambda.1se}}{ Largest value of \code{lambda} such that \code{error.CV} is within 1% of the minimum.}
##' \item{\code{nzero}}{ Number of non-zero coefficients at each \code{lambda}.}
##' \item{\code{Lambdas}}{ The values of \code{lambda} used in the fits.}
##' \item{\code{error.CV}}{ The averaged cross validation error.}
##' }
##' @seealso \code{\link{GADAG}}, \code{\link{GADAG_Run}}, \code{\link{GADAG_Analyze}}.
##' @author \packageAuthor{GADAG}
##'
##'@references
##' M. Champion, V. Picheny, M. Vignes, Inferring large graphs using l-1 penalized likelihood,
##' Statistics and Computing (2017).
##'
##' @examples
##'  #############################################################
##'  # Loading toy data
##'  #############################################################
##'  data(toy_data)
##'  # toy_data is a list of two matrices corresponding to a "star"
##'  # DAG (node 1 activates all other nodes):
##'  # - toy_data$X is a 100x10 design matrix
##'  # - toy_data$G is the 10x10 adjacency matrix (ground trough)
##'
##'  #############################################################
##'  # Tuning the parameter of penalization
##'  #############################################################
##'  # Simple run, whithout specifying GADAG parameters
##'  \dontrun{
##'  GADAG_CV_results <- GADAG_CV(X=toy_data$X)
##'  print(GADAG_CV_results$lambda.1se) # best lambda
##'  }
##'  # If desired, additional plot for the averaged cross validation
##'  # error
##'  \dontrun{
##'  GADAG_CV_results <- GADAG_CV(X=toy_data$X,plot.CV=1)
##'  }
##'  # Given the best lambda, re-run GADAG
##'  \dontrun{
##'  GADAG_results <- GADAG_Run(X=toy_data$X, lambda=GADAG_CV_results$lambda.1se)
##'  print(GADAG_results$G.best) # optimal adjacency matrix graph
##'  }

GADAG_CV <- function(X, Lambdas = NULL, n.folds = 10, threshold = 0.1,
                      GADAG.control = list(n.gen=250, tol.Shannon=1e-6, max.eval=1e7,pop.size=5*ncol(X), p.xo=.25, p.mut=.05),
                      grad.control = list(tol.obj.inner=1e-6, max.ite.inner=50),
                      ncores=1, plot.CV=0) {

  #############################################################
  # INPUTS:
  # X: design matrix (size n*p)
  # Lambdas: lambda sequence
  # n.folds:number of folds for cross-validation
  # pop.size: initial population size for GA
  # p.xo: crossover probability
  # p.mut: mutation probability
  # n.gen: maximal number of population generations
  # tol.Shannon: termination threshold for the Shannon entropy
  # max.eval: maximal number of objective function evaluations (i.e. gradient descents)
  # tol.obj.inner: tolerance for the inner optimization (gradient descent)
  # max.ite.inner: maximum number of iterations for the inner optimization (gradient descent)
  # plot.CV: 1/0, if 1 plots the averaged cross validation error
  #
  # OUTPUTS:
  # lambda.min: value of lambda that gives the minimum value of the cross validation error
  # nzero: number of non-zero coefficient at each lambda
  # Lambdas: the value of lambda used in the fits
  # error.CV: the averaged cross-validation error
  #############################################################
  if (is.null(Lambdas)){
    lambda.log <- seq(-6,0,length.out=100)
  } else {
    lambda.log <- log(Lambdas)
  }
  if (is.null(GADAG.control$n.gen)){
    n.gen <- 100
  } else {
    n.gen <- GADAG.control$n.gen
  }
  if (is.null(GADAG.control$max.eval)){
    max.eval <- 1e4
  } else {
    max.eval <- GADAG.control$max.eval
  }
  if (is.null(GADAG.control$tol.Shannon)){
    tol.Shannon <- 1e-6
  } else {
    tol.Shannon <- GADAG.control$tol.Shannon
  }
  if (is.null(GADAG.control$pop.size)){
    pop.size <- 10
  } else {
    pop.size <- GADAG.control$pop.size
  }
  if (is.null(GADAG.control$p.xo)){
    p.xo <- 0.25
  } else {
    p.xo <- GADAG.control$p.xo
  }
  if (is.null(GADAG.control$p.mut)){
    p.mut <- 0.05
  } else {
    p.mut <- GADAG.control$p.mut
  }
  if (is.null(grad.control$tol.obj.inner)){
    tol.obj.inner <- 1e-6
  } else {
    tol.obj.inner <- grad.control$tol.obj.inner
  }
  if (is.null(grad.control$max.ite.inner)){
    max.ite.inner <- 50
  } else {
    max.ite.inner <- grad.control$max.ite.inner
  }
  GADAG.control <- list(n.gen=n.gen,max.eval=max.eval,tol.Shannon,pop.size=pop.size,p.xo=p.xo,p.mut=p.mut)
  grad.control <- list(tol.obj.inner=tol.obj.inner,max.ite.inner=max.ite.inner)

  ##### Define the groups for cross validation #####
  n <- dim(X)[1]
  p <- dim(X)[2]

  if (n.folds<=nrow(X)){
    Groups <- cvFolds(n,K=n.folds)
  } else {
    cat(paste0("The number of samples being too small, the number of groups for cross validation is set to ",nrow(X)))
    n.folds <- nrow(X)
    Groups <- cvFolds(n,K=n.folds)
  }
  Groups <- Groups$which

  ##### Run GADAG for each fold #####
  Fold.CV <- function(i,l){
    Test <- which(Groups==i)
    Xtest <- X[Test,]
    Train <- c(1:n)
    Train <- Train[-Test]
    Xtrain <- X[Train,] - matrix(1,length(Train),1) %*% colMeans(X[Train,])
    Xtrain = Xtrain / sqrt( (1/length(Train)) * matrix(1,length(Train),1) %*% colSums(Xtrain^2))

    GADAG_results <- GADAG_Run(X = Xtrain, lambda = l, threshold = threshold, GADAG.control = GADAG.control, grad.control = grad.control, ncores = ncores)
    error.CV <- (1/(length(Test)*p)) * sum((Xtest-Xtest%*%GADAG_results$G.best)^2)

    nnzero <- length(which(abs(GADAG_results$G.best)>0))

    return(c(error=error.CV,nnzero=nnzero))
  }

  Lambda.CV <- function(j){
    lambda <- exp(lambda.log[j])
    results <- mclapply(X=1:n.folds, FUN=Fold.CV,l=lambda)
    error.CV <- mean(do.call(rbind, results)[,1])
    nnzero <- round(mean(do.call(rbind, results)[,2]))

    return(c(error=error.CV,nnzero=nnzero))
  }

  results <- mclapply(X=1:length(lambda.log), FUN=Lambda.CV)
  error.CV <- do.call(rbind, results)[,1]
  nnzero <- do.call(rbind, results)[,2]

  error.min <- error.CV[error.CV<min(error.CV)+0.01*min(error.CV)]
  lambda.min <- exp(lambda.log[which.min(error.CV)])
  lambda.1se <- exp(lambda.log[which(error.CV==error.min[length(error.min)])])

  if (plot.CV == 1){
    plot(lambda.log,error.CV,pch=16,xlab="log.lambda",ylab="Averaged cross validation error")
    abline(h=min(error.CV)+0.01*min(error.CV),lty=2)
    abline(h=min(error.CV)-0.01*min(error.CV),lty=2)
    abline(v=log(lambda.1se),lty=2,col="red")
  }
  return(list(lambda.min=lambda.min, lambda.1se = lambda.1se, nzero=nnzero, lambda=exp(lambda.log), error.CV=error.CV))
}
