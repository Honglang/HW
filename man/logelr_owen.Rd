\name{logelr_owen}
\alias{logelr_owen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
logelr_owen function used in empirical likelihood function elm
}
\description{
computing empirical likelihood for a mean vector
}
\usage{
logelr_owen(x, mu, lam)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{lam}{
%%     ~~Describe \code{lam} here~~
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, mu, lam)
{
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    if (n <= p)
        stop("Need more observations than variables in logelr.")
    mu <- as.vector(mu)
    if (length(mu) != p)
        stop("Length of mean doesn't match number of variables in logelr.")
    z <- t(t(x) - mu)
    arg <- 1 + z \%*\% lam
    -sum(llog(arg, 1/n))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
