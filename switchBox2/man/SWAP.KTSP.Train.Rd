\name{SWAP.KTSP.Train}

\alias{SWAP.KTSP.Train}

\title{Training the K-TSP classifier }

\description{\code{SWAP.KTSP.Train} trains the K-TSP classifier.}

\usage{
SWAP.KTSP.Train(data, situation, krange= c(3,5:10),
    RestrictedPairs=c(), Filter=SWAP.Filter.Wilcoxon,
    FilterParam=list(UpDown=TRUE, geneNo=100), chooseKDyn=TRUE)
}

\arguments{
  \item{data}{The matrix of values (e.g., gene expression)
    in which columns represents samples and rows features
    (e.g., genes).}
  \item{situation}{A binary vector containing the training labels.}
  \item{krange}{The range for \code{k}, the number of disjoint TSP
    used for classification.}
  \item{RestrictedPairs}{The list of candidate pairs from which
    the algorithm chooses its pairs.
    If it is empty, there is no restriction.}
  \item{Filter}{The function applied to filter the genes.
    If empty, the genes are not filtered.
    The default is using a Wilcoxon rank-sum test,
    using KTSP.Filter.Wilcoxon}
  \item{FilterParam}{A parameters list that must be passed
    to the \code{Filter} function.
    Default values are used by \code{SWAP.Filter.Wilcoxon},
    This has been shown a useful filter}
  \item{chooseKDyn}{Logical. If \code{TRUE}, a unique classifier
    is chosen by dynamic algorithm (see: \code{SWAP.ChoooseDynk})}
}

\value{The KTSP classifier, in the form of a list or a vector of lists,
  which contains the following components:
  \item{$name}{The classifier name.}
  \item{$TSPs}{A k by 2 matrix, containing
    the gene names of the k TSPs, or their indices
    depending on their availability.}
  \item{$score}{Top k scores for the top TSPs.}
  \item{$label}{The class labels: \dots}
}


%\note{ }

\seealso{
  \code{\link{SWAP.KTSP.Classify}},
  \code{\link{SWAP.KTSP.Train.Plain}},
  \code{\link{SWAP.CalculateSignedScore}}
}

\references{
  Marchionni, L.; Afsari, B.; Geman, D. and Leek, T.J.
  A simple and reproducible breast cancer prognostic test.
  BMC Genomics, 2013, 14:336
}

\author{Bahman Afsari \email{bahman@jhu.edu},  Luigi Marchionni \email{marchion@jhu.edu}}

\keyword{KTSP Training}

\examples{

### Load the example data
data(matETABM115)
### Train the classifier
classifier <- SWAP.KTSP.Train(matETABM115, ttmGroupETABM115, chooseKDyn=FALSE)
### Apply the classifier
pred <- SWAP.KTSP.Classify(matETABM115, classifier)

### Train the classifier selecting the best KTSP combination
classifier <- SWAP.KTSP.Train( matETABM115, ttmGroupETABM115,chooseKDyn=TRUE)
### Apply the Classifier
pred <-  SWAP.KTSP.Classify(matETABM115,classifier)

### Extract the TSP from the classifier
TSps <- classifier$TSPs[ 1:6 , ]

### Train a classifier using a restrcted set of TSP
classifier <- SWAP.KTSP.Train(matETABM115, ttmGroupETABM115,
    chooseKDyn=FALSE, RestrictedPairs=TSps);

### Compute Type 1 error
mean(pred[ttmGroupETABM115==1]==1)
#[1] 0.8478261

### Compute Type 0 error
mean(pred[ttmGroupETABM115==0]==0)
#[1] 0.6206897

### Train a classifier using the KTSP.Filter.Wilcoxon() to filter the genes
classifier <- SWAP.KTSP.Train(matETABM115, ttmGroupETABM115,
     chooseKDyn=TRUE, Filter=KTSP.Filter.Wilcoxon,
     FilterParam=list(UpDown=TRUE, geneNo=6))
### Apply the classifier
pred <- SWAP.KTSP.Classify(matETABM115,classifier)

}
