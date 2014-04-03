\name{KTSP.Classifiy}

\alias{KTSP.Classify}

\title{Function to classify samples using a KTSP classifier.}

\description{
  \code{KTSP.Classify} classifies new test samples
  using KTSP coming out of the function \code{\link{KTSP.Train}}.
  This function was used in Marchionni et al, 2013, BMC Genomics,
  and it is maintained only for backward compatibility.
  It has been replaced by \code{\link{SWAP.KTSP.Classify}}.
}

\usage{
  KTSP.Classify(data, classifier, combineFunc)
}

\arguments{
  \item{data}{the test data: a matrix in which the rows
    represent the genes and the columns the samples.}
  \item{classifier}{The output of  \code{\link{KTSP.Train}},
    a KTSP classifier.}
  \item{combineFunc}{A user defined function to combine
    the predictions of the individual K TSP.
    If missing the consensus classification among the
    majority of the TSPs will be used.}
}


\seealso{
  \code{\link{KTSP.Train}},
  \code{\link{SWAP.KTSP.Classify}},
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

### Turn into a numeric vector with values equal to 0 and 1
trainingGroupNum <- as.numeric(trainingGroup) - 1

### Show group variable for the TRAINING set
table(trainingGroupNum)


##################################################
### Train a classifier using default filtering function based on the Wilcoxon test
classifier <- KTSP.Train(matTraining, trainingGroupNum, n=8)

### Show the classifier
classifier


##################################################
### Testing on new data

### Load the example data for the TEST set
data(testingData)

### Turn into a numeric vector with values equal to 0 and 1
testingGroupNum <- as.numeric(testingGroup) - 1

### Show group variable for the TEST set
table(testingGroupNum)

### Apply the classifier to one sample of the TEST set using
### sum of votes grearter than 2
testPrediction <- KTSP.Classify(matTesting, classifier,
     combineFunc = function(x) sum(x) > 2)

### Show prediction
table(testPrediction, testingGroupNum)

}

\keyword{KTSP, classification, prediction}
