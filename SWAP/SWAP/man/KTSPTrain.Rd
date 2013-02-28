\name{SWAP.KTSP.Train}
\alias{SWAP.KTSP.Train}
\title{Training the K-TSP classifier }
\description{KTSP.Train trains the K-TSP classifier. }
\usage{
classifier = KTSP.Train( data, situation, krange= c(3,5:10), RestrictedPairs=c() , Filter=KTSP.Filter.Wilcoxon, FilterParam=list(UpDown=TRUE,geneNo=100),chooseKDyn=TRUE)
}
\arguments{
\item{data}{ the matrix of the expression values whose columns represents samples and the rows represents the genes.}
\item{situation}{ a bi-value vector containing the training labels.}
\item{krange}{The range for k, the number of disjoint TSP used for classification.}
\item{RestrictedPair}{The list candidate pairs from which the algorithm chooses its pairs. If it is empty, there is no restriction.  }
\item{Filter}{The function applied to filter the genes. If empty, the genes are filtered.}
\item{FilterParam}{The parameters must be passed to Filter function.}
\item{chooseKDyn}{If true, a unique classifier is chosen by dynamic algorithm. see: SWAP.ChoooseDynk }
}
\value{classifier a list or a vector of lists, which contains the classifier. The list contains of
\item{$name}{The classifier name, e.g. 3-TSP which the classifier uses the top three comparison for the classification. }
\item{$TSPs}{Top Scoring Pair:an k by 2 matrix, containing the gene names of the k top pairs or their indices depending on the availability.}
\item{$score}{Top k scores: the scores of the pairs in $TSP}
\item{$label}{contains the class labels. If labels[1] corresponds to large kTSP.}
}
\author{Bahman Afsari}
\keyword{KTSP Training}
\seealso{SWAP.KTSP.Classify,SWAP.KTSP.Train.Plain, KTSP.CalculateSignedScore }
\examples{
load("matETABM115.rda");

classifier = SWAP.KTSP.Train( matETABM115, ttmGroupETABM115,chooseKDyn=FALSE);
pred = SWAP.KTSP.Classify(data,classifier);

classifier = SWAP.KTSP.Train( matETABM115, ttmGroupETABM115,chooseKDyn=TRUE);
pred = SWAP.KTSP.Classify(data,classifier);

TSps=classifier$TSPs[1:6,];

classifier = SWAP.KTSP.Train( matETABM115, ttmGroupETABM115,chooseKDyn=FALSE,RestrictedPairs=TSps);



mean(pred[ttmGroupETABM115==1]==1)#Type 1 error
#[1] 0.8478261
mean(pred[ttmGroupETABM115==0]==0)#Type 0 error
#[1] 0.6206897

}
