\name{matTesting}

\alias{matTesting}

\docType{data}

\title{Gene expression matrix for test set data}

\description{A numerical matrix containing gene expression matrix
  for 70 genes and 307 breast cancer patients (test set data)
  from the Buyse et al cohort
  (see the \code{\link[mammaPrintData]{mammaPrintData}} package).
}

\usage{
data(matTesting)
}

\format{
  The \code{matTesting} matrix contains normalized
  expression values for the 70 gene signature (rows)
  across 307 samples (columns).
  Group information (emph{``bad'' versus ``good'' }prognosis)
  is shown in \code{colnames(matTesting)}.
}

\details{
This dataset corresponds to the breast cancer patients' cohort
published by Buyse and collleagues in JNCI (2006).
The gene expression matrix was obtained from
the \code{mammaPrintData} package as described 
by Marchionni and colleagues in BMC Genomics (2013).
}


\seealso{
  \code{\link{matTraining}}
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
data(matTesting)

### Show the class of the ``matTesting'' object
class(matTesting)

### Show the dimentions of the ``matTesting'' matrix
dim(matTesting)

### Show the first  10 sample names of the ``matTest'' matrix
head(colnames(matTesting), n=10)

### Create a phenotypic group variable for the 307 samples
testGroup <- factor(gsub(".+\\\\.", "", colnames(matTesting)))

### Show group variable
table(testGroup)

}

\keyword{datasets}
