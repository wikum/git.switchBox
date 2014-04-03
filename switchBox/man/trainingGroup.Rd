\name{trainingGroup}

\alias{trainingGroup}

\docType{data}

\title{Training set phenotypes}

\description{A factor with two levels
  describing the phenotypes for the training data
  (Glas et al cohort,
  see the \code{\link[mammaPrintData]{mammaPrintData}} package).
}

\usage{
data(trainingData)
}

\format{
  The \code{trainingGroup} factor contains phenotypic information
  for the 78 samples of the training dataset.
}

\details{
This phenotype factor corresponds to the breast cancer patients' cohort
published by Glas and collleagues in BMC Genomics (2006).
The information was obtained from
the \code{mammaPrintData} package as described 
by Marchionni and colleagues in BMC Genomics (2013).
}


\seealso{
  \code{\link{testingGroup}}
}

\author{
  Bahman Afsari \email{bahman.afsari@gmail.com},
  Luigi Marchionni \email{marchion@jhu.edu}
}

\source{
  \url{http://www.ncbi.nlm.nih.gov/pubmed/?term=17074082}

  \url{http://www.ncbi.nlm.nih.gov/pubmed/?term=23682826}
}

\references{
 Marchionni et al.,
 "A simple and reproducible breast cancer prognostic test", 
 \emph{BMC Genomics}, 2013, \bold{14}(1):336-242
 
 Glas et al et al., "Converting a breast cancer microarray signature
 into a high-throughput diagnostic test",
 {BMC Genomics}, 2006, \bold{7}:278-288
}

\examples{

### Load gene expression data for the training set
data(trainingData)

### Show the class of the ``trainingGroup'' object
class(trainingGroup)

### Show group variable
table(trainingGroup)

}

\keyword{datasets}
