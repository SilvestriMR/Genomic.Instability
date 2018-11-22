\name{ratio}
\alias{ratio}
\docType{data}
\title{
Data.frame with CNA results (1 WGS sample) from control-FREEC
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("ratio")}
\format{
  A data frame with 4667 observations on the following 4 variables.
  \describe{
    \item{\code{V1}}{factor-chromosome ID}
    \item{\code{V2}}{numeric-start of genomic regions}
    \item{\code{V3}}{numeric-end of genomic regions}
    \item{\code{cn}}{numeric-copy number values for each genomic region}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
 (http://boevalab.com/FREEC/tutorial.html)
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(ratio)
str(ratio)

}
\keyword{datasets}
