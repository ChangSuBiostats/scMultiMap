% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_pair}
\alias{example_pair}
\title{A simulated independent peak-gene pair}
\format{
\subsection{\code{example_pair}}{

Two lists (one named 'gene' and the other named 'peak'), each with 2 elements:
\describe{
\item{counts}{A length 2000 vector of gene/peak counts}
\item{seq_depths}{A length 2000 vector of sequencing depths of scRNA-seq/scATAC-seq}
}
}
}
\usage{
example_pair
}
\description{
This data set is created only for illustrative purposes and is used to test
scMultiMap_IRLS() and scMultiMap_WLS().
The source code for generating this data set is in data-raw/ on Github.
}
\keyword{datasets}
