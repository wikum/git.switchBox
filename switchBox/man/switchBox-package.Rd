\name{switchBox-package}
\alias{switchBox-package}
\alias{switchBox}
\docType{package}


\title{A package to train and validate K-Top-Scoring-Pair classifiers.}

\description{
  The switchBox package allows to train and validate
  a K-Top-Scoring-Pair (KTSP) classifier,
  as used by Marchionni et al (BMC Genomics, 2013).
  KTSP is an extension of the TSP classifier described
  by Geman and colleagues (Bioinformatics, 2005).
  The TSP algorithm is a simple binary classifier based on
  the ordering of two measurements.
}


\section{switchBox package features}{
  The switchBox package contains several utilities enabling to:
  \describe{
    {A)  Filter the features to be used to develop the classifier
    (\emph{i.e.}, differentially expressed genes);}

  {B)  Compute the scores for all available feature pairs
    to identify the top performing TSP;}
  
  {C)  Compute the scores for selected feature pairs
    to identify the top performing TSP;}

  {D)  Identify the number of $K$ TSP to be used in the final classifier;}
  
  {E)  Compute individual TSP votes for one class or the other 
    and aggregate the votes based on various methods;}
  
  {F)  Classify new samples based on the top KTSP 
    based on various methods.}
}}


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


\keyword{ package }
