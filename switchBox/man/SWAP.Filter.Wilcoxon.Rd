\name{SWAP.Filter.Wilcoxon}

\alias{SWAP.Filter.Wilcoxon}

\title{Statistical feature filtering based on theWilcoxon test on the ranks of expressions.}

\description{\code{SWAP.Filter.Wilcoxon} filters the features to
  identify the top different ones between phenotypes
  to be used for KTSP classifier implementation.
}

\usage{
SWAP.Filter.Wilcoxon(phenoGroup, inputMat, featureNo = 100, UpDown = TRUE)
}

\arguments{
  \item{phenoGroup}{a factor with levels containing training labels
    for the phenotype of interest.}
  \item{inputMat}{a numerical matrix containing  feature
    measurements to be used to implement the classifiier
    (\emph{e.g.}, the set of gene expression values).
    The columns of this matrix correspond to samples and
    must correspond to \code{phenoGroup}.
    The rows represents the features and
    \code{rownames(inputMat)} will be used as feature names.}
  \item{featureNo}{an integer specifcying the number
    of different features to be returned.}
  \item{UpDown}{logical value specifying whether
    an equal proportion of features displaying opposite
    change across the two phenotypes should be returned
    (\emph{e.g.}an equal number of up- and down-regulated genes).}
}

\value{The names of the features that survided statistical filtering.}

 
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

### Load gene expression data for the training set
data(trainingData)

### Return the top 10 features (both up and down)
SWAP.Filter.Wilcoxon(trainingGroup, matTraining, featureNo=10)

### Return the top 10 up-regulated features
SWAP.Filter.Wilcoxon(trainingGroup, matTraining, featureNo=10)

}
