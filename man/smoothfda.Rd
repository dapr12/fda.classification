% Generated by roxygen2 (4.0.1.99): do not edit by hand
\name{smoothfda}
\alias{smoothfda}
\title{Smoothfda - Computes the functional spatial median of the functional data}
\usage{
smoothfda(fdaobj, bandwidth, degree = 1)
}
\arguments{
\item{fdaobj}{An object of the FDA Class.}

\item{bandwitdh}{Bandwidth}

\item{degree}{Degree of the polynomial}
}
\value{
A component list with:
\itemize{
\item{"Data"}{}
\item{"argvals"}{}
\item{"rangevals"}{}
\item{YSmooth}{}
\item{"CVScore"}{}
\item{Mean}{}
}
}
\description{
Smoothfda - Computes the functional spatial median of the functional data
}
\examples{
\dontrun{
#Computing the Median
matplot(time,CVScores$YSmooth, type="l")
lines(time, Median$median, col="red",lwd = 4)
}
}
\references{
Ramsay, James O., and Silverman, Bernard W. (2002), Applied Functional Data Analysis, Springer,New York.
}

