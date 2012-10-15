\name{KTSP.Train3TSP}
\alias{KTSP.Train3TSP}
\title{Training the 3-TSP classifier }
\description{KTSP.Train3TSP trains the 3TSP classifier. }
\usage{
classifier = KTSP.Train3TSP( data, situation);
}
\arguments{
\item{data}{ the matrix of the expression values whose columns represents samples and the rows represents the genes.}
\item{situation}{ a vector containing the training labels. Its elements should be one or zero.}
\item{classifier}{The classifier which contains the 3-TSPs, their names and their scores. It should be passed to KTSP.Classify as for classification on test samples.}

}
\author{Bahman Afsari}
\keyword{3TSP Training}
\seealso{KTSP.Classify, KTSP.Train, KTSP.CalculateSignedScore }
\examples{
### load("matETABM115.rda");
### classifier = KTSP.Train3TSP( matETABM115, ttmGroupETABM115);
### pred = KTSP.Classify(data,classifier);
### 
### mean(pred[ttmGroupETABM115==1]==1)#Type 1 error
### #[1] 0.8478261
### mean(pred[ttmGroupETABM115==0]==0)#Type 0 error
### #[1] 0.6206897
### 
}
