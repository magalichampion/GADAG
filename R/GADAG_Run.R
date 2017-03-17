##' @title Run GADAG
##' @description Function to run GADAG, an algorithm that aims at inferring large sparse directed acyclic graphs by minimizing the negative log-likelihood with a convex program embedded in a genetic algorithm.
##' @param X Design matrix, with samples (n) in rows and variables (p) in columns.
##' @param lambda Parameter of penalization (>0).
##' @param threshold Thresholding value for the estimated edges.
##' @param GADAGcontrol A list containing the conditions for terminating GADAG:
##' \itemize{
##' \item{n.gen}{ maximal number of population generations (>0),}
##' \item{tol.Shannon}{ threshold for the Shannon entropy (>0),}
##' \item{max.eval}{ maximal number of iterations for the inner optimization (>0).}
##' }
##' @param GAcontrol A list containing the genetic algorithm inherent parameters:
##' \itemize{
##' \item{pop.size}{ initial population size for the genetic algorithm (>0),}
##' \item{p.xo}{ crossover probability of the genetic algorithm (between 0 and 1),}
##' \item{p.mut}{ mutation probability of the genetic algorithm (between 0 and 1).}
##' }
##' @param Gradcontrol A list containing the parameters for controlling the inner optimization, i.e. the gradient descent
##' \itemize{
##' \item{tol.obj.inner}{ tolerance (>0),}
##' \item{max.ite.inner}{ maximum number of iterations (>0).}
##' }
##' @param ncores Number of cores (>1, depending on your computer).
##' @param print.level 0 no print, 1 some info on the genetic algorithm behaviour are printed.
##' @param return.level 0 only best solution is returned, 1 evolution of the current best solution and statistics on the population fitness values are also returned.
##' @return A list with the following elements:
##' \itemize{
##' \item{f.best}{ Best fitness value.}
##' \item{P.best}{ Best node order (vector of length p).}
##' \item{T.best}{ Corresponding best edges values (vector of length p).}
##' \item{G.best}{ Best graph (matrix form).}
##' \item{f.best.evol}{ Evolution of the best fitness value across the iterations (if return.level=1).}
##' \item{P.best.evol}{ Evolution of the best node order across the iterations (if return.level=1).}
##' \item{T.best.evol}{ Evolution of the best edges values across the iterations (if return.level=1).}
##' \item{fmin.evol}{ Evolution of the minimal fitness value of the population across the iterations (if return.level=1).}
##' \item{fmean.evol}{ Evolution of the averaged fitness value of the population across the iterations (if return.level=1).}
##' \item{fp10.evol}{ Evolution of the quantiles of the fitness value across the iterations (if return.level=1).}
##' \item{fp90.evol}{ Evolution of the quantiles of the fitness value across the iterations (if return.level=1).}
##' \item{Shannon.evol}{ Evolution of the Shannon entropy of the population across the iterations (if return.level=1).}
##' }
##' @author \packageAuthor{GADAG}
##' @examples
##'  ########################################################
##'  # Loading toy data
##'  ########################################################
##'  data(toy_data)
##'
##'  ########################################################
##'  # Running GADAG
##'  ########################################################
##'  GADAG_results <- GADAG_Run(X=toy_data$X, lambda=0.1)

GADAG_Run <- function(X, lambda, threshold=0.1,
        GADAGcontrol = list(n.gen=100, tol.Shannon=1e-6, max.eval=1e4),
        GAcontrol = list(pop.size=10, p.xo=.25, p.mut=.05),
        Gradcontrol = list(tol.obj.inner=1e-6, max.ite.inner=50),
        ncores=1,print.level=0, return.level=0) {

  #############################################################
  # INPUTS:
  # X: design matrix (size n*p)
  # lambda: penalty
  # threshold: thresholding value for the edges
  # pop.size: initial population size for GA
  # p.xo: crossover probability
  # p.mut: mutation probability
  # n.gen: maximal number of population generations
  # tol.Shannon: termination threshold for the Shannon entropy
  # max.eval: maximal number of objective function evaluations (i.e. gradient descents)
  # tol.obj.inner: tolerance for the inner optimization (gradient descent)
  # max.ite.inner: maximum number of iterations for the inner optimization (gradient descent)
  # print.level: 0 no print, 1 some info on the GA behaviour
  # return.level: 0 only best solution is returned,
  #               1 evolution of the current best solution and statistics on the population fitness values are also returned
  #
  # OUTPUTS:
  # G.best: best graph
  # P.best: best individual (node order)
  # T.best: corresponding best triangular matrix (edges values)
  # f.best: best fitness value
  # More outputs can be provided (see return.level)
  #############################################################

  ############## Initialisation step ###############
  n.gen <- GADAGcontrol$n.gen
  tol.Shannon <- GADAGcontrol$tol.Shannon
  max.eval <- GADAGcontrol$max.eval
  pop.size <- GAcontrol$pop.size
  p.xo <- GAcontrol$p.xo
  p.mut <- GAcontrol$p.mut

  ##### Create and evaluate initial population #####
  if (Gradcontrol$max.ite.inner < 0){
    stop("Gradcontrol$max.ite.inner should be non-negative.")
  }

  n <- dim(X)[1]
  p <- dim(X)[2]
  XtX <- crossprod(X)

  Pop       <- create.population(p, pop.size)
  evalTandf <- evaluation(Pop, X, XtX, lambda, Gradcontrol = Gradcontrol, ncores=ncores)
  f.pop     <- evalTandf$f
  T.pop     <- evalTandf$Tpop
  f.count   <- pop.size

  ##### Initialize all indicators ####
  Shannon.val <- rep(0,p)
  for (j in (1:p)){
    frac <- colMeans(Pop==j)
    h <- frac*log(frac)
    h[is.nan(h)] <- 0
    Shannon.val[j] <- -sum(h)
  }

  pop.best        <- which.min(f.pop)
  f.best.alltimes <- min(f.pop)
  P.best.alltimes <- Pop[pop.best,]
  T.best.alltimes <- T.pop[pop.best,]

  if (return.level != 0) {
    fmin.pop.save  <- rep(0, n.gen)
    fmean.pop.save <- rep(0, n.gen)
    fp10.pop.save  <- rep(0, n.gen)
    fp90.pop.save  <- rep(0, n.gen)

    fmin.pop.save[1]  <- min(f.pop)
    fmean.pop.save[1] <- mean(f.pop)
    fp10.pop.save[1]  <- quantile(f.pop, probs = .1)
    fp90.pop.save[1]  <- quantile(f.pop, probs = .9)

    f.best.alltimes.save  <- rep(0, n.gen)
    P.best.alltimes.save  <- matrix(data = 0, nrow = n.gen, ncol = p)
    T.best.alltimes.save  <- matrix(data = 0, nrow = n.gen, ncol = p^2)

    f.best.alltimes.save[1]  <- f.best.alltimes
    P.best.alltimes.save[1,] <- P.best.alltimes
    T.best.alltimes.save[1,] <- T.best.alltimes

    Shannon.save   <- matrix(data = 0, nrow = n.gen, ncol = p)
    Shannon.save[1,]  <- Shannon.val
  }

  ##### Main Loop #####
  k <- 1
  if ((print.level != 0) ) cat("Nb gen | best f |  average (sd) f in pop | Shannon | f count \n")

  while (k < n.gen && sum(Shannon.val) > tol.Shannon && f.count<max.eval){
    if (print.level != 0) cat(k, f.best.alltimes, mean(f.pop), "(", sd(f.pop), ")", sum(Shannon.val), f.count, "\n")

    ##### Selection operator #####
    I <- selection(Pop, f.pop)

    ##### Cross-over + mutation + evaluation on selected Population ####
    Results.crossover   <- crossover(Pop=Pop[I,],p.xo=p.xo)
    Children <- Results.crossover$Children
    I.cross <- Results.crossover$I.cross

    if (length(I.cross)>1){
      Children   <- mutation(Children, p.mut=p.mut)
      evalTandfe <- evaluation(Pop=Children, X=X, XtX=XtX, lambda=lambda, Gradcontrol=Gradcontrol, ncores=ncores)
      I.cross <- I[I.cross]
      Pop[I.cross,]   <- Children
      f.pop[I.cross]  <- evalTandfe$f
      T.pop[I.cross,] <- evalTandfe$Tpop
    }

    #### Replace worst by best parent ####
    pop.worst         <- which.max(f.pop)
    f.pop[pop.worst]  <- f.best.alltimes
    Pop[pop.worst,]   <- P.best.alltimes
    T.pop[pop.worst,] <- T.best.alltimes

    #### Update counters ####
    f.count <- f.count + length(I.cross) + p-1
    k <- k + 1

    Shannon.val <- rep(0,p)
    for (j in (1:p)){
      frac <- colMeans(Pop==j)
      h <- frac*log(frac)
      h[is.nan(h)] <- 0
      Shannon.val[j] <- -sum(h)
    }

    pop.best <- which.min(f.pop)
    if (min(f.pop) < f.best.alltimes) {
      f.best.alltimes <- min(f.pop)
      P.best.alltimes <-  Pop[pop.best,]
      T.best.alltimes <- T.pop[pop.best,]
    }

    if (return.level != 0) {
      fmin.pop.save[k]  <- min(f.pop)
      fmean.pop.save[k] <- mean(f.pop)
      fp10.pop.save[k]  <- quantile(f.pop, probs = .1)
      fp90.pop.save[k]  <- quantile(f.pop, probs = .9)

      f.best.alltimes.save[k]  <- f.best.alltimes
      P.best.alltimes.save[k,] <- P.best.alltimes
      T.best.alltimes.save[k,] <- T.best.alltimes

      Shannon.save[k,]  <- Shannon.val
    }
  }
  if (print.level != 0) cat(k, f.best.alltimes, mean(f.pop), "(", sd(f.pop), ")", sum(Shannon.val), f.count, "\n")

  # END OF MAIN LOOP
  P <- chrom(P.best.alltimes)
  T <- matrix(T.best.alltimes,p,p)
  T[T<threshold] <- 0
  G.best <- P %*% T %*% t(P)
  Gbest.bin <- G.best*0
  Gbest.bin[(abs(G.best)>0)] <- 1

  ###############################################################
  ###############################################################
  #if (plot.level > 0 && return.level>0){
  #highlight.node <- c(1,2,3,4,5,6)
  #col <- rep("black", p)
  #col[highlight.node] <- c("blue","cyan","green", "magenta", "red", "grey")
  #lwd <- rep(1, p)
  #lwd[highlight.node] <- 2
  #lwd[1] <- 3
  #lty <- rep(1, p)
  #lty[1:6] <- 1:6

  # Fpop, min, mean and quantiles
  #ylim=c(min(c(fmin, GAresults$fmin.pop.save)), max(GAresults$fp90.pop.save))
  #plot(GAresults$f.best.evol, type="l", main="Cost function", xlab="generation #", ylab="-logLikelihood", lwd=2, col="red",
  #     ylim=ylim)
  #polygon(x=c(1:length(GAresults$fp10.pop.save),rev(1:length(GAresults$fp10.pop.save))), y=c(GAresults$fp10.pop.save, rev(GAresults$fp90.pop.save)), border=NA,col="lightgrey")
  #lines(GAresults$fp10.pop.save)
  #lines(GAresults$fp90.pop.save)
  #legend(x=n.gen/2, y=max(ylim), legend=c("Current best", "Population"), col=c("red","black"), lwd=c(2,1), text.width=60)

  # Shannon
  #plot(GAresults$Shannon.save[,1], type="l", main="Shannon entropy", xlab="generation #", ylab="entropy", ylim=c(0, max(GAresults$Shannon.save)),
  #     col=col[1], lwd=lwd[1], lty=lty[1])
  #for (i in 2:p) lines(GAresults$Shannon.save[,i], col=col[i], lwd=lwd[i], lty=lty[i])

  # Nodes paths - bestever
  #A <- which(GAresults$P.best.evol==1, arr.ind=T)
  #a <- sort(A[,1], index.return=TRUE)$ix
  #plot(A[a,1], A[a,2], ylim=c(0,p), type="l", main="Best permutation", xlab="generation #", ylab="Node path", col=col[1], lwd=lwd[1], lty=lty[1])

  #for (i in 2:6){
  #  A <- which(GAresults$P.best.evol==i, arr.ind=T)
  #  a <- sort(A[,1], index.return=TRUE)$ix
  #  lines(A[a,1], A[a,2], col=col[i], lwd=lwd[i], lty=lty[i])
  #}
  #legend(x=n.gen/2, y=p/2, legend=paste0("node ",1:6), lwd=lwd[1:6], lty=lty[1:6], col=col[1:6], text.width=30)

#} else {
 # cat("No convergence results to plot. \n")
#}
#}

  #############################################################
  #############################################################
  if (return.level > 0) {
    return(list(f.best=f.best.alltimes, P.best=P.best.alltimes, T.best=T.best.alltimes,G.best=Gbest.bin,
                f.best.evol=f.best.alltimes.save[1:k], P.best.evol=P.best.alltimes.save[1:k,], T.best.evol=T.best.alltimes.save[1:k,],
                fmin.evol=fmin.pop.save[1:k], fmean.evol=fmean.pop.save[1:k], fp10.evol=fp10.pop.save[1:k],
                fp90.evol=fp90.pop.save[1:k],
                Shannon.evol=Shannon.save[1:k,]))
  } else {
    return(list(f.best=f.best.alltimes, P.best=P.best.alltimes, T.best=T.best.alltimes, G.best=Gbest.bin))
  }
}
