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
########################################################
# Loading toy data
########################################################
data(toy_data)

########################################################
# Running GADAG
########################################################
GADAG_results <- GADAG_Run(X=toy_data$X,lambda=0.1)

########################################################
# Evaluating GADAG Results
########################################################    
GADAG_analysis <- GADAG_Analyse(Results=GADAG_results,
	G=toy_data$G,X=toy_data$X)
}