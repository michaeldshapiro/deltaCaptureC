% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deltaCaptureC.R
\name{plotSignificantRegions}
\alias{plotSignificantRegions}
\title{This produces a plot of the region of interest showing regions of significance.}
\usage{
plotSignificantRegions(
  significantRegions,
  significanceType = "abs",
  title = "Significant Regions",
  xLabel = "viewpoint",
  legend = TRUE
)
}
\arguments{
\item{significantRegions}{a GRanges object as produced by
getSignificantRegions}

\item{significanceType}{= 'abs' a variable indicating whether to plot
significance according to min, max or abs.}

\item{title}{a title for the plot}

\item{xLabel}{= 'viewpoint' supplies an xlabel}

\item{legend}{= TRUE whether or not to show the legend}
}
\value{
a ggplot object
}
\description{
This function takes a input the GRanges object produced by
getSignificant regions and produces a ggplot of significant features
}
\examples{
plotOfSignificantRegions = plotSignificantRegions(significantRegions)
}
