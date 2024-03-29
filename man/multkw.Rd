% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multkw.R
\name{multkw}
\alias{multkw}
\title{Multivariate Kruskal-Wallis test}
\usage{
multkw(group, y, print = TRUE)
}
\arguments{
\item{group}{The factorial variable that subsets the dataset in groups. Can be a character vector, a factorial vector or an integer/numeric vector.}

\item{y}{The dataset of n numeric(or integer) variables.}

\item{print}{Whether the test should be printed (\code{TRUE}, the default) or not (e.g., to be stored in an object)}
}
\value{
Output is either a list (with \code{"simplify=FALSE"}) or a vector (with \code{"simplify=TRUE"}) containing the results of the multivariate Kruskal-Wallis test.
}
\description{
This function computes a multivariate Kruskal-Wallis test for n numeric variables relative to one factorial variable (that subsets the dataset in groups)
}
\details{
A "standard" multivariate Kruskal-Wallis test is computed, deleting all missing data.
}
\examples{
data(airquality)
datamkw<-airquality[,1:4]
multkw(y=datamkw,airquality$Month)

}
\references{
\insertRef{He.etal.2017}{ULT}
}
\seealso{
See chapter 2.2.2 and 4.2 of the \href{http://d-scholarship.pitt.edu/19411/1/Fanyin_ETD_draft_08-06-2013.pdf}{PhD manuscript of Fanyin He} and 'Methodology' of \insertCite{He.etal.2017;textual}{ULT} for more details.
}
\author{
Fanyin He (most of the statistical function)

Jacob Maugoust (packaging)
}
