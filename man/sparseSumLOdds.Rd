% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparseSumLOdds.R
\name{sparseSumLOdds}
\alias{sparseSumLOdds}
\title{sparseSumLOdds}
\usage{
sparseSumLOdds(stab_lor_list, exprs_data, intercept)
}
\arguments{
\item{stab_lor_list}{A list of odds ratios annotated with Entrez IDs. Can simply use the output of the coef_stabilizer() function.}

\item{exprs_data}{A matrix with rows representing genes and columns representing samples.}

\item{intercept}{A boolean whether to include the intercept}
}
\description{
This function calculates the sum of the log odds per class, even if not all classifier genes are present.
}
