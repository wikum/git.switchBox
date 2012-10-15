\name{KTSP.Train}
\alias{KTSP.Train}
\title{Training the K-TSP classifier }
\description{KTSP.Train trains the K-TSP classifier. }
\usage{
classifier = KTSP.Train( data, situation, n);
}
\arguments{
\item{data}{ the matrix of the expression values whose columns represents samples and the rows represents the genes.}
\item{situation}{ a vector containing the training labels. Its elements should be one or zero.}
\item{n}{The number of disjoint TSP used for classification. If before n pairs, the score drops to zero, the TSP with zero score are ignored.}
\item{classifier}{The classifier which contains the K-TSPs, their names and their scores. It should be passed to KTSP.Classify as for classification on test samples.}
}
\author{Bahman Afsari}
\keyword{KTSP Training}
\seealso{KTSP.Classify, KTSP.Train3TSP, KTSP.CalculateSignedScore }
\examples{
### load("matETABM115.rda");
### classifier = KTSP.Train( matETABM115, ttmGroupETABM115, 3);
### pred = KTSP.Classify(data,classifier);
### 
### mean(pred[ttmGroupETABM115==1]==1)#Type 1 error
### #[1] 0.8478261
### mean(pred[ttmGroupETABM115==0]==0)#Type 0 error
### #[1] 0.6206897
### 
}
