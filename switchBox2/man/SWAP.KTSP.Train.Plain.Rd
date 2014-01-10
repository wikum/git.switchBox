\name{SWAP.KTSP.Train.Plain}

\alias{SWAP.KTSP.Train.Plain}

\title{Trains the K-TSP classifier without choosing dynamically the number of genes.  }

\description{KTSP.Train trains the K-TSP classifier. }

\usage{
SWAP.KTSP.Train.Plain( data, situation, krange= c(3,5:10), genes)
}

\arguments{
\item{data}{ the matrix of the expression values whose columns represents samples and the rows represents the genes.}
\item{situation}{ a bi-value vector containing the training labels.}
\item{krange}{The range for k, the number of disjoint TSP used for classification.}
\item{genes}{The subset of genes from data which must be investigated.}
}

\details{
The core of the algorithm is to find the pair of genes whose comparisons can distinguish between two classes. kTSP algorithm finds the most discriminative pair and then find the next discriminative disjoint pair. This process contrinues till k-pairs has been collected. This function provides kTSP for a range of k's. For more details read the references.
}

\value{
classifier is a vector of list with length equal to krange. classifier can be sent to SWAP.KTSP.Classify for classification.The list contains:
\item{$name}{The classifier name, e.g. 3-TSP which the classifier uses the top three comparison for the classification. }
\item{$TSPs}{Top Scoring Pair:an k by 2 matrix, containing the gene names of the k top pairs or their indices depending on the availability.}
\item{$score}{Top k scores: the scores of the pairs in $TSP}
\item{$label}{contains the class labels. If labels[1] corresponds to large kTSP.}
}

\author{Bahman Afsari}

\keyword{KTSP Training}

\seealso{SWAP.KTSP.Classify, SWAP.CalculateSignedScore}


\examples{
data(matETABM115);
classifier = SWAP.KTSP.Train.Plain( matETABM115, ttmGroupETABM115, krange=c(3,5),genes = 1:30);
pred = SWAP.KTSP.Classify(matETABM115,classifier);

mean(pred[ttmGroupETABM115==1]==1)#Type 1 error
#[1] 0.8478261
mean(pred[ttmGroupETABM115==0]==0)#Type 0 error
#[1] 0.6206897

}
