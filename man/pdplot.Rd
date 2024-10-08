% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pdplot.R
\name{pdplot}
\alias{pdplot}
\title{Posterior density plots}
\usage{
pdplot(
  object,
  interval = 0.9,
  display = "b1",
  order_predictor = NULL,
  order_group = NULL,
  top_label = 0.05,
  left_label = 0.1,
  title_size = 3,
  label_size = 1.6,
  line_width = 0.2,
  point_size = 0.5,
  err_bar_lwd = 0.2,
  xlab = "Parameter estimate",
  ylab = "Posterior density",
  xylab_size = 3,
  xtext_size = 8
)
}
\arguments{
\item{object}{An object generated by the meta function}

\item{interval}{What interval to display}

\item{display}{Which model parameter to display b1 or b0 (currently only b1)}

\item{order_predictor}{The order of the predictors}

\item{order_group}{The order of the groups}

\item{top_label}{The space for the top label}

\item{left_label}{The space for the left labels}

\item{title_size}{Size of the titles}

\item{label_size}{Size of the labels in the plots}

\item{line_width}{The line width}

\item{point_size}{The size of the point}

\item{err_bar_lwd}{The line width of the error bars}

\item{xlab}{The x-label text}

\item{ylab}{The y-label text}

\item{xylab_size}{x and y label size}

\item{xtext_size}{size of the x-axis text}
}
\description{
Posterior density plots
}
