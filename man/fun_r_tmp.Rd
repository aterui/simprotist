% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_r_tmp.R
\name{fun_r_tmp}
\alias{fun_r_tmp}
\title{Internal function: function for temperature-dependent population growth}
\usage{
fun_r_tmp(
  a,
  e_a,
  e_d,
  k_b = 8.6 * 1e-05,
  tmp,
  tmp_r = 298.15,
  tmp_h,
  scale_factor = 0
)
}
\arguments{
\item{a}{Normalization constant.}

\item{e_a}{Activation energy}

\item{e_d}{Deactivation energy}

\item{k_b}{Boltzmann's constant}

\item{tmp}{Temperature}

\item{tmp_r}{Reference temperature}

\item{tmp_h}{Half-r temperature}

\item{scale_factor}{Scale factor}
}
\description{
Internal function: function for temperature-dependent population growth
}
