\name{KTSP.Classifiy}
\alias{KTSP.Classify}
\title{Classify test samples with KTSP classifier }
\description{KTSP.Classify classifies new test samples using KTSP coming out of the function KTSP.Train. }
\usage{
  KTSP.Classify(test, classifier, combineFunc);
}
\arguments{
\item{test}{ the new test data. If a matrix the rows should represent genes and the columns the samples}
\item{classifier}{ The outcome of KTSP.Train or KTSP.Train3TSP}
\item{combineFunc}{ A user define function to combine the prediction of
  the individual K TSP. If missing the consensus classification among the
  majority of the TSP will be used.}
}
\author{Bahman Afsari, Luigi Marchionni}
\keyword{KTSP Classifying}
\seealso{ KTSP.Train3TSP, KTSP.Train, KTSP.CalculateSignedScore }
\examples{
### load("matETABM115.rda");
### classifier <- KTSP.Train( matETABM115, ttmGroupETABM115, 3);
### pred = KTSP.Classify(data, classifier);
### 
### mean(pred[ttmGroupETABM115==1]==1)#Type 1 error
### #[1] 0.8478261
### mean(pred[ttmGroupETABM115==0]==0)#Type 0 error
### #[1] 0.6206897

}
