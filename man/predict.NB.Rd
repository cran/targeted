% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.NB.R
\name{predict.NB}
\alias{predict.NB}
\title{Predictions for Naive Bayes Classifier}
\usage{
\method{predict}{NB}(object, newdata, expectation = NULL, threshold = c(0.001, 0.001), ...)
}
\arguments{
\item{object}{density object}

\item{newdata}{new data on which to make predictions}

\item{expectation}{Variable to calculate conditional expectation wrt
probabilities from NB classifier}

\item{threshold}{Threshold parameters. First element defines the threshold
on the probabilities and the second element the value to set those
truncated probabilities to.}

\item{...}{Additional arguments to lower level functions}
}
\description{
Naive Bayes Classifier predictions
}
\author{
Klaus K. Holst
}
