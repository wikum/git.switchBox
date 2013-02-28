\name{SWAP.Filter.Wilcoxon}
\alias{SWAP.Filter.Wilcoxon}
\title{Wilcoxon Filter for DE Genes used for SWAP functions. }
\description{SWAP.Filter.Wilcoxon filters genes using Wilcoxon Test on the ranks of expressions.}
\usage{index = SWAP.Filter.Wilcoxon(situation,data,FiltParam=list(genesNo=100,UpDown=TRUE));}
\arguments{
\item{data}{ the set of expression values. The matrix of the expression values whose columns represents samples and the rows represents the genes.}
\item{situation}{ a vector containing the training labels. Its elements should be one or zero.}
\item{FiltParam$genesNo}{The number of genes to be filtered.}
\item{FiltParam$UpDown}{If true and pvalue negative, the half of 'genesNo' genes are selected from up regulated and another half are selected from down regulated. Otherwise, they are chosen from only the two sided Wilcoxon test is considered.}
}
\value{index contains the filtered genes. If the data has row names, it contains gene names otherwise it cointains the index.}
\author{Bahman Afsari}
\keyword{Wilcoxon Differential Expressed Genes}
\seealso{SWAP.Classify, SWAP.Train}
\examples{
load("matETABM115.rda");
index = SWAP.Filter.Wilcoxon(ttmGroupETABM115,matETABM115,FiltParam=list(genesNo=5,UpDown=TRUE));
}
