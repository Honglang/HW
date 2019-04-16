\name{llogpp}
\alias{llogpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
llogpp function used in empirical likelihood function elm
}
\description{
computing empirical likelihood for a mean vector
}
\usage{
llogpp(z, eps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (z, eps)
{
    ans <- z
    lo <- (z < eps)
    ans[lo] <- -1/eps^2
    ans[!lo] <- -1/z[!lo]^2
    ans
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
