\name{SIVCMTest}
\alias{SIVCMTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model checking test for single-index varying coefficient models (SIVCM) with functional response
}
\description{
test whether the SIVCM with functional response is adequate or not
}
\usage{
SIVCMTest(n, p, m, tm, x, ally, B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size}
  \item{p}{dimension of covariates x}
  \item{m}{number of time points}
  \item{tm}{time points}
  \item{x}{covariates (n*p matrix)}
  \item{ally}{response functions (n*m matrix)}
  \item{B}{number of bootstrap replications}
}
\value{
A list of outputs
\item{TestSTat }{test statistic}
\item{Pvalue }{p-value}
\item{Cri95 }{critical value at significance level 0.05}
}

%\author{Qing Jiang}

\examples{
##---- Step 1: generate data ----

n <- 30
m <- 15
a <- 0
p <- 3
B <- 500

data <- GenData.Sa(n,m,a)
x <- data$x  # n*p
ally <- data$ally  # n*m
beta <- data$beta  # p*m
tm <- data$tm 

##---- Step 2: model test ----

SIVCMTest(n,p,m,tm,x,ally,B)

}


