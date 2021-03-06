% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deltaCaptureC.R
\name{getSizeFactorsSE}
\alias{getSizeFactorsSE}
\title{Get the size factors for SummarizedExperiment}
\usage{
getSizeFactorsSE(se)
}
\arguments{
\item{se}{A SummarizedExperiment with an assay counts}
}
\value{
The same SummarizedExperiment with an additional column in its colData
giving the size factors for counts
}
\description{
This function takes a SummarizedExperiment with an assay counts and returns
this object with a column sizeFactors added to its colData
}
\examples{
miniSEWithSizeFactors = getSizeFactorsSE(miniSE)
}
