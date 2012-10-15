\name{KTSP.CalculateSignedScore}
\alias{KTSP.CalculateSignedScore}
\title{Calculating the pair-wise score}
\description{KTSP.CalculateSignedScore calculates the pairwise score between ang gene from a set of genes with any gene from another set of genes. }
\usage{pairwisescore = KTSP.CalculateSignedScore(  situation , data1, data2 );}
\arguments{
\item{data1}{ the first set of expression values. The matrix of the expression values whose columns represents samples and the rows represents the genes.}
\item{data2}{ the second set of expression values.}
\item{situation}{ a vector containing the training labels. Its elements should be one or zero.}
\item{pairwisescore$score}{a matrix with pairwise score. Its rows is equal to the rows of the data1 (genes in data1) and its columns is equal to the rows in data2 (genes in data2). It is signed because    score(i,j)=P(X_i<X_j|1)-P(X_i<X_j|0) - P(X_i>X_j|1) +P(X_i>X_j|0) + C The third and the forth terms are for avoiding the equlities. The C term is proportion to the secondary score to avoid the ties.  C = (E(X_j-X_i|1)-E(X_j-X_i|0))/K where K is big enough to make sure that the secondary score does not intervene with the primary score.}
\item{pairwisescore$a}{a[i,j]= sum(data1[i,situation==0]<data2[j,situation==0]) used for Fisher exact test. }
\item{pairwisescore$M}{ is the number of samples with label 0}
\item{pairwisescore$k}{ k[i,j] = sum(data1[i,]<data2[j,]) used for Fisher exact test. }
\item{pairwisescore$N}{ is the number of samples. }
\item{pairwisescore$P}{ P[i,j] = P(data1[i]<data[j]|situation==0) = a_ij/#(situation==0) }
\item{pairwisescore$Q}{ Q[i,j] = P(data1[i]<data[j]|situation==1) }

}
\author{Bahman Afsari}
\keyword{Pairwise score}
\seealso{KTSP.Classify, KTSP.Train, KTSP.Train3TSP}
\examples{
### load("matETABM115.rda");
### pairwisescore = KTSP.CalculateSignedScore(  ttmGroupETABM115, matETABM115, matETABM115);
}
