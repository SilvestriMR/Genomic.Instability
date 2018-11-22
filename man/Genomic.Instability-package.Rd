\name{Genomic.Instability-package}
\alias{Genomic.Instability-package}
\alias{Genomic.Instability}
\docType{package}
\title{
\packageTitle{Genomic.Instability}
}
\description{
\packageDescription{Genomic.Instability}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{Genomic.Instability}
\packageIndices{Genomic.Instability}
R package for the evaluation of Genomic Instability starting from Copy number alterations data. In particular it computes Large-scale state transitions (LST).
LST is defined as the number of chromosomal breaks between adjacent regions of at least 10 Mb. 
}
\author{
\packageAuthor{Genomic.Instability}

Maintainer: \packageMaintainer{Genomic.Instability}
}
\references{
S. B. Greene et al., “Chromosomal Instability Estimation Based on Next Generation Sequencing and Single Cell Genome Wide Copy Number Variation Analysis,” PLoS One, vol. 11, no. 11, p. e0165089, Nov. 2016. 

T. Popova et al., “Ploidy and Large-Scale Genomic Instability Consistently Identify Basal-like Breast Carcinomas with BRCA1/2 Inactivation,” Cancer Res November 1 2012 (72) (21) 5454-5462
}

\keyword{ package }
\seealso{
}
\examples{
data(ratio)
MeasureLST(data=ratio, window = 639926, ID = "AA", workflow = "SS")

}
