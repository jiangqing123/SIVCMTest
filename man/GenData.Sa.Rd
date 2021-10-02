\name{GenData.Sa}
\alias{GenData.Sa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for generating data
}
\description{
Generate data for Example 1(a) in "Model Checks on Single-index Varying Coefficient Models with Functional Response"
}
\usage{
GenData.Sa(n, m, a)
}
\arguments{
  \item{n}{sample size}
  \item{m}{number of time points}
  \item{a}{distance away from the null}
}

\value{
A list of outputs
\item{x }{covariates (n*p matrix)}
\item{ally }{response functions (n*m matrix)}
\item{tm }{time points (m vector)}
\item{beta }{true coefficient functions (p*m matrix)}
}

%\author{Qing Jiang}

\examples{
##---- Generate data ----
n <- 30
m <- 15
a <- 0
GenData.Sa(n,m,a)
}
