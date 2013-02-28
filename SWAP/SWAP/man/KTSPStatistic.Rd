\name{SWAP.KTSP.Statitsics}
\alias{SWAP.KTSP.Statitsics}
\title{kTSP statistic: The number of pair votes in favor one class minus the number of pair votes against the other class. }
\description{KTSP.Train trains the K-TSP classifier. }
\usage{
statistic = SWAP.KTSP.Statitsics(data, classifiers)
}
\arguments{
\item{classifiers}{ kTSP classifiers. The outcome of the SWAP.KTSP.Statistics}
\item{data}{ a matrix of column vectors, we want to calculate the kTSP statistics.}
}
\details{
To classify, we can only compare kTSP statistic with zero. If bigger than zero then classifier[[i]]$label[1] o.w. classifier[[i]]$label[2]
}
\value{
an n by m matrix of kTSP staitistics. statistic[i,j] represents the kTSP statistics of the classifier[[i]] on the j-th column of data.
}

\author{Bahman Afsari}
\keyword{KTSP statistic}
\seealso{KTSP.Classify, KTSP.CalculateSignedScore KTSP Training}
\examples{
load("matETABM115.rda");
classifier = SWAP.KTSP.Train.Plain( matETABM115, ttmGroupETABM115, krange=c(3,5));
kTSP = SWAP.KTSP.Statitsics(matETABM115,classifier);


}
