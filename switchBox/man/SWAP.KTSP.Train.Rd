\name{SWAP.KTSP.Train}

\alias{SWAP.KTSP.Train}

\title{Function for  training the K-TSP classifier.}

\description{
  \code{SWAP.KTSP.Train} trains a K-TSP classifier
  for the specific phenotype of interest.
  The classifiers resulting from using this function can be
  passed to \code{\link{SWAP.KTSP.Classify}}
  for samples classification.
}

\usage{
SWAP.KTSP.Train(inputMat, phenoGroup, krange = c(3, 5, 7:10),
  FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs, ...)
}

\arguments{
  \item{inputMat}{is a numerical matrix containing the
    measurements (\emph{e.g.}, gene expression data)
    to be used to build the K-TSP classifier.
    The columns represent samples and the
    rows represent the features (\emph{e.g.}, genes).
    The number of columns must agree
    with the length of \code{phenoGroup}.
    Note that \code{rownames(inputMat)} will be used
    as must feature names (\emph{e.g.}, gene symbols)
    in all subsequen analyses.}
  \item{phenoGroup}{ is a factor containing the training phenotypes
    with two levels.}
  \item{krange}{an integer (or a vector of integers)
    defining the number of possible disjoint 
    Top Scoring Pairs (TSPs) from which to chose
    to build the final classifier.}
  \item{FilterFunc}{is a filtering function to reduce the
    starting number of features to be used to identify the
    Top Scoring Pairs (TSP). The deault filter is based on
    the Wilcoxon rank-sum test
    and alternative filtering functions can be passed too
    (see \code{SWAP.Filter.Wilcoxon} for details).}
  \item{RestrictedPairs}{is a character matrix with two columns
    containing the feature pairs to be considered for score calculations.
    Each row should contain a pair of feature names matching the
    rownames of \code{inputMat}.
    If RestrictedPairs is missing all available feature pairs
    will be considered.}
  \item{...}{Additional argument passed to the filtering
  function \code{FilterFunc}.}
}

\value{
  The KTSP classifier, in the form of a list,
  which contains the following components:
  \item{name}{The classifier name.}
  \item{TSPs}{A \code{k} by 2 matrix, containing
    the feature names for each TSP. These names
  correspond to the \code{rownames(inputData)}.}
  \item{$score}{TSP scores for the top \code{k} TSPs.}
  \item{$label}{The class labels.}
}

\seealso{
  \code{\link{SWAP.KTSP.Classify}},
  \code{\link{SWAP.Filter.Wilcoxon}},
  \code{\link{SWAP.CalculateSignedScore}}
}


\author{
  Bahman Afsari \email{bahman.afsari@gmail.com},
  Luigi Marchionni \email{marchion@jhu.edu}
}

\source{
  \url{http://www.ncbi.nlm.nih.gov/pubmed/23682826}

  \url{http://www.ncbi.nlm.nih.gov/pubmed/16105897}

  \url{http://www.ncbi.nlm.nih.gov/pubmed/16131522}

  \url{http://www.ncbi.nlm.nih.gov/pubmed/16646797}
}

\references{
  Marchionni et al.,
  "A simple and reproducible breast cancer prognostic test.",
 \emph{BMC Genomics}, 2013, \bold{14}(1):336-342

 Tan et al.,
 "Simple decision rules for classifying human cancers from gene expression profiles.",
\emph{Bioinformatics} (2005) \bold{21}(20), 3896-3904.
 
 Xu et al., 
 "Robust prostate cancer marker genes emerge from direct integration
 of inter-study microarray data" 
\emph{Bioinformatics} (2005) \bold{21}(20), 3905-3911.

Geman et al.
"Classifying gene expression profiles from pairwise mRNA comparisons"
\emph{Statistical applications in genetics and molecular biology}
(2004) \bold{3}.1 : 1071.

}

\examples{

##################################################
### Load gene expression data for the training set
data(trainingData)


### Show group variable for the TRAINING set
table(trainingGroup)


##################################################
### Train a classifier using default filtering function based on the Wilcoxon test
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup, krange=c(3:15))

### Show the classifier
classifier


##################################################
### Train another classifier from the top 4 best features 
### according to the deafault  filtering function
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup,
			      FilterFunc=SWAP.Filter.Wilcoxon, featureNo=4)

### Show the classifier
classifier


##################################################
### To use all features "FilterFunc" must be set to NULL
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup, FilterFunc=NULL)

### Show the classifier
classifier


##################################################
### Train a classifier using and alternative filtering function.
### For instance we can use the a  "t.test"  to selec the features 
### with an absolute t-statistics larger than a specified quantile
topRttest <- function(situation, data, quant = 0.75) {
	out <- apply(data, 1, function(x, ...) t.test(x ~ situation)$statistic )
	names(out[ abs(out) > quantile(abs(out), quant) ])
}

### Show the top features selected 
topRttest(trainingGroup, matTraining, quant=0.95)

### Train a classifier using the alternative filtering function
### and also define the maximum number of TSP using "krange"
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup,
			      FilterFunc = topRttest, quant = 0.75, krange=c(15:30) )

### Show the classifier
classifier


##################################################
### Training with restricted pairs

### Define a set of specific pairs to be used for classifier development
### For this example we will a random set of features
### In a real example these pairs should  be provided by the user.
set.seed(123)
somePairs <- matrix(sample(rownames(matTraining), 6^2, replace=FALSE), ncol=2)
head(somePairs, n=3)
dim(somePairs)

### Train a classifier using the restricted feature pairs and the default filtering
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup,
			      RestrictedPairs = somePairs, krange=3:16)

### Show the classifier
classifier

}

\keyword{KTSP, classification, prediction}
