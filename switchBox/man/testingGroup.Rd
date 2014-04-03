\name{testingGroup}

\alias{testingGroup}

\docType{data}

\title{Testing set phenotypes}

\description{A factor with two levels
  describing the phenotypes for the testing data
  (Buyse et al cohort,
  (see the \code{\link[mammaPrintData]{mammaPrintData}} package).
}

\usage{
data(testingData)
}

\format{
  The \code{matTesting} factor contains phenotypic information
  for the 307 samples of the testing dataset.
}

\details{
This phenotype factor corresponds to the breast cancer patients' cohort
published by Buyse and collleagues in JNCI (2006).
The gene expression matrix was obtained from
the \code{mammaPrintData} package as described 
by Marchionni and colleagues in BMC Genomics (2013).
}


\seealso{
  \code{\link{trainingGroup}}
}

\author{
  Bahman Afsari \email{bahman.afsari@gmail.com},
  Luigi Marchionni \email{marchion@jhu.edu}
}

\source{
  \url{http://www.ncbi.nlm.nih.gov/pubmed/?term=16954471}
  
  \url{http://www.ncbi.nlm.nih.gov/pubmed/?term=23682826}
}

\references{
 Marchionni et al.,
 "A simple and reproducible breast cancer prognostic test", 
 \emph{BMC Genomics}, 2013, \bold{14}(1):336-342
 
 Buyse et al., "Validation and clinical utility of a 70-gene prognostic
 signature for women with node-negative breast cancer",
 \emph{JNCI}, 2006, \bold{98}(17):1183-1192
}

\examples{

### Load gene expression data for the test set
data(testingData)

### Show the class of the ``testingGroup'' object
class(testingGroup)

### Show group variable
table(testingGroup)

}

\keyword{datasets}
