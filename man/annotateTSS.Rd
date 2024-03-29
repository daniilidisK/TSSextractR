% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtractKnowledge.R
\name{annotateTSS}
\alias{annotateTSS}
\title{Annotate TSS based and assign them to categories}
\usage{
annotateTSS(TSSdataframe, GeneList, plot = T)
}
\arguments{
\item{TSSdataframe}{A Dataframe containing start, strand and RRS information for each TSS}

\item{GeneList}{A Dataframe containing gene names, left & right coordinates and strand information}

\item{plot}{A logical value to show (TRUE) or not show (FALSE) plots}
}
\value{
A full Dataframe containing gene TSS categorization, gene names, etc
}
\description{
Annotate TSS based and assign them to categories
}
