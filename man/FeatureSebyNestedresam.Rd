% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PepFeatureSelect.R
\name{FeatureSebyNestedresam}
\alias{FeatureSebyNestedresam}
\title{Title}
\usage{
FeatureSebyNestedresam(
  dataset = MLtestData,
  target = "judge",
  taskid = "FeatureSebyNestedresam",
  filtermethodlrn = "classif.rpart",
  nevals = 20,
  fselector = "random_search"
)
}
\arguments{
\item{dataset}{A dataframe of Machine Learing test dataset}

\item{target}{A character of column name in dataset}

\item{taskid}{A character of taskid}

\item{filtermethodlrn}{A filter selection by which lrn}

\item{nevals}{A number of nevals}

\item{fselector}{A character of fselector}
}
\value{

}
\description{
Title
}
\examples{
rr = FeatureSebyNestedresam(MLtestData)
extract_inner_fselect_archives(rr) #https://mlr3fselect.mlr-org.com/
}
