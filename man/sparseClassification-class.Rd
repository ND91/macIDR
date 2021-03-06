% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparseClassification-class.R
\docType{class}
\name{sparseClassification-class}
\alias{sparseClassification-class}
\title{An S4 class to represent the sparseClassification object.}
\description{
An S4 class to represent the sparseClassification object.
}
\section{Slots}{

\describe{
\item{\code{gene_lodds}}{A list of vectors, where each vector represents the log odds ratios for the classifier gene per macrophage subset.}

\item{\code{lodds}}{A matrix of summed log odds per macrophage subset.}

\item{\code{odds}}{A matrix of summed odds per macrophage subset.}

\item{\code{proportions}}{A matrix of responses [0-1] per macrophage subset.}

\item{\code{predictions}}{A data frame of predictions per sample.}

\item{\code{presence_stats}}{A data frame containing the presence and absence statistics of the classifier genes per macrophage subset.}

\item{\code{present_genes}}{A list of genes present per macrophage subset.}

\item{\code{absent_genes}}{A list of genes absent per macrophage subset.}
}}

