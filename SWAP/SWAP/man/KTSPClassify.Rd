\name{SWAP.KTSP.Classifiy}
\alias{SWAP.KTSP.Classify}
\title{Classify test samples with KTSP classifier }
\description{SWAP.KTSP.Classify classifies new test samples using KTSP coming out of the function SWAP.KTSP.Train. }
\usage{
pred = SWAP.KTSP.Classify( test,classifier,thresh=0);
}
\arguments{
\item{test}{ the new test data. If a matrix the rows should represent genes and the columns the samples}
\item{classifier}{ The outcome of SWAP.KTSP.Train}
\item{thresh}{The threshold applied to kTSP statistic. Playing with this, the user can find the specificity and sensitivity needed. The threshold zero, corresponds to the balanced prior for both classes.}

}
\value{
\item{pred}{ The prediction of the classifier. the row numbers are equal to number of classifiers and the column is equal to the number of samples.}
}
\author{Bahman Afsari}
\keyword{KTSP Classifying}
\seealso{ KTSP.Train3TSP, KTSP.Train, KTSP.CalculateSignedScore }
\examples{
load("matETABM115.rda");
classifier = SWAP.KTSP.Train.Plain( matETABM115, ttmGroupETABM115, krange=c(3,5));
pred = SWAP.KTSP.Classify(data,classifier);

mean(pred[ttmGroupETABM115==1]==1)#Type 1 error
#[1] 0.8478261
mean(pred[ttmGroupETABM115==0]==0)#Type 0 error
#[1] 0.6206897

}
