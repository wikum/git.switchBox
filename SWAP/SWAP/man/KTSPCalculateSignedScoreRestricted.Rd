\name{SWAP.KTSP.CalculateSignedScore}
\alias{SWAP.KTSP.CalculateSignedScore}
\title{Calculas the pair-wise score}
\description{SWAP.KTSP.CalculateSignedScore calculates the pairwise score between specific pairs of genes. }
\usage{pairwisescore = SWAP.KTSP.CalculateSignedScore.Restricted(  situation , data1, data2, pairs );}
\arguments{
\item{data1}{ the first set of expression values. The matrix of the expression values whose columns represents samples and the rows represents the genes.}
\item{data2}{ the second set of expression values.}
\item{situation}{ a vector containing the training labels. It must be a binary vector.}
\item{pairs}{a n by 2 matrix containing pair of genes. pairs can be gene names or only index.}
}
\value{
The output is a list containing the following items:
\item{labels}{The labels in the vector situation. }
\item{pairwisescore$P}{ P[i] = P(data_{pairs[1,i]}<data_{pairs[2,i]}|situation==$label[1]) }
\item{pairwisescore$Q}{ Q[i] = P(data_{pairs[1,i]}<data_{pairs[2,i]}|situation==$label[2]) }
\item{pairwisescore$score}{A matrix with pairwise scores. Its rows is equal to the rows of the data1 (genes in data1) and its columns is equal to the rows in data2 (genes in data2). The row names of $score are the gene names in data1 and $score are gene names in data2.
score[i,j] is signed because  score[i,j]=P[i,j]-Q[i,j] + C. The C term is the tie breaker and proportion to the secondary score to avoid the ties.  If R_i is rank of data_i in the whole gene profile. C = (E(R_i-R_j|1)-E(R_i-R_j|0))/K where K is big enough to make sure that the secondary score does not intervene with the primary score.}

}
\author{Bahman Afsari}
\keyword{Pairwise score}
\seealso{SWAP.KTSP.Classify, SWAP.KTSP.Train}
\examples{
load("matETABM115.rda");
Restrictedpairs = rbind(c(1,2),c(3,4),c(10,12));
pairwisescore = SWAP.KTSP.CalculateSignedScore.Restricted(  ttmGroupETABM115, matETABM115, matETABM115,Restrictedpairs);
RestrictedGenes= matrix(rownames(matETABM115)[Restrictedpairs],3,2);
pairwisescore = SWAP.KTSP.CalculateSignedScore.Restricted(  ttmGroupETABM115, matETABM115, matETABM115,RestrictedGenes);

}
