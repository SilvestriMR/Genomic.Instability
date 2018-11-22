\name{MeasureLST}
\alias{MeasureLST}

\title{
Function to compute Large-scale state transitions (LST)
}
\description{
Starting from Copy Number Alteration data derived from Whole Genome Sequencing experiments, the function compute Large-scale state transitions (LST).
}
\usage{
MeasureLST(data, window, ID, workflow = c("SS", "MS"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
data.frame object. First column must contain chromosome name (e.g chr1 - factor), the second and third column the start and end position of copy number event respectively (integer), the fourth the number of copy detected in each regions (integer - For a reliable results, regions without coverage or with no significant CNA calls should be set as -1).
}
  \item{window}{
window used for the copy number analysis estimation (i.e end-start). List of data.frame in Multi-sample case (see workflow parameter).
}
  \item{ID}{
Name of the Output. It will be used as Matrix rownames. List of character in Multi-sample case (see workflow parameter).
}
  \item{workflow}{
"SS" for single sample analysis, "MS" for multiple samples (one "data" object for each samples).
}
}
\details{
Chromosome Y will be excluded from the analysis.
}
\value{
In case of "SS" workflow, the function will returns a data.frame with LST score for each chromosome and a global score (column = "score") for one sample
In case of "MS" workflow, the function will returns a data.frame with LST score for each chromosome and a global score (column = "score") for all sample considered.
}
\references{
S. B. Greene et al., “Chromosomal Instability Estimation Based on Next Generation Sequencing and Single Cell Genome Wide Copy Number Variation Analysis,” PLoS One, vol. 11, no. 11, p. e0165089, Nov. 2016. 

T. Popova et al., “Ploidy and Large-Scale Genomic Instability Consistently Identify Basal-like Breast Carcinomas with BRCA1/2 Inactivation,” Cancer Res November 1 2012 (72) (21) 5454-5462
}
\author{
Marco Silvestri
}
\note{

}


\seealso{

}
\examples{
data(ratio)
MeasureLST(data=ratio, window=ratio[1,3],ID="AA",workflow="SS")

## The function is currently defined as
function (data, window, ID, workflow = c("SS", "MS")) 
{
    if (workflow == "SS") {
        for (x in data[, 1]) {
            par <- 1e+07/window
            ch_name <- x
            idx <- data[data[, 1] \%in\% x, ]
            idx1 <- rle(idx[, 4])
            idx2 <- grep("-1", idx1$values)
            idx1$lengths[idx2] <- 0
            c <- as.vector(idx1$lengths)
            c <- append(c, 0)
            c[c < par] <- 0
            c[c > par] <- 1
            count <- as.data.frame(sum(rollsum(c, 2) > 1))
            colnames(count) <- x
            assign(x, count)
        }
        TabLST <- do.call(cbind, mget(ls(pattern = "^chr")))
        TabLST$Score <- rowSums(TabLST)
        rownames(TabLST) <- ID
        TabLST <- TabLST[, c("chr1", "chr2", "chr3", "chr4", 
            "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
            "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
            "chrX", "Score")]
        TabLST
    }
    else {
        for (y in ID) {
            data <- data[[1]]
            wind <- window[1]
            for (x in data[, 1]) {
                par <- 1e+07/wind
                ch_name <- x
                idx <- data[data[, 1] \%in\% x, ]
                idx1 <- rle(idx[, 4])
                idx2 <- grep("-1", idx1$values)
                idx1$lengths[idx2] <- 0
                c <- as.vector(idx1$lengths)
                c <- append(c, 0)
                c[c < par] <- 0
                c[c > par] <- 1
                count <- as.data.frame(sum(rollsum(c, 2) > 1))
                colnames(count) <- x
                assign(x, count)
            }
            TabLST <- do.call(cbind, mget(ls(pattern = "^chr")))
            TabLST$Score <- rowSums(TabLST)
            rownames(TabLST) <- y
            assign(paste("Bind_", y, sep = ""), TabLST)
            data <- data[-1]
            window <- window[-1]
        }
        Tab_summ <- do.call(rbind, mget(ls(pattern = "^Bind_")))
        Tab_summ <- Tab_summ[, c("chr1", "chr2", "chr3", "chr4", 
            "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
            "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
            "chrX", "Score")]
        Tab_summ
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
