\name{elm}
\alias{elm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
computing empirical likelihood for a mean vector
}
\description{
computing empirical likelihood for a mean vector
}
\usage{
elm(x, mu, lam, maxit = 25, gradtol = 1e-07, svdtol = 1e-09, itertrace = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
data matrix
}
  \item{mu}{
a mean vector for the el evaluation
}
  \item{lam}{
missing is ok
}
  \item{maxit}{
maximum iteration
}
  \item{gradtol}{
tolerance for gd
}
  \item{svdtol}{
tolerance for svd
}
  \item{itertrace}{
itertrace
}
}
\details{
el
}
\value{
\item{logelr}{log elr}
\item{lambda}{lambda}
\item{grad}{grad}
\item{hess}{hess}
\item{wts}{wts}
\item{nits}{nits}
}
\references{
http://statweb.stanford.edu/~owen/empirical/
}
\author{
Honglang Wang
}
\note{
https://math.iupui.edu/~hlwang/
}
\examples{
## The function is currently defined as
function (x, mu, lam, maxit = 25, gradtol = 1e-07, svdtol = 1e-09,
    itertrace = F)
{
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    mu <- as.vector(mu)
    if (length(mu) != p)
        stop("Mu must have same dimension as observation vectors.")
    if (n <= p)
        stop("Need more observations than variables in elm.")
    z <- t(t(x) - mu)
    TINY <- sqrt(.Machine$double.xmin)
    scale <- mean(abs(z)) + TINY
    z <- z/scale
    if (!missing(lam)) {
        lam <- as.vector(lam)
        lam <- lam * scale
        if (logelr_owen(z, rep(0, p), lam) > 0)
            lam <- rep(0, p)
    }
    if (missing(lam))
        lam <- rep(0, p)
    if (svdtol < TINY)
        svdtol <- TINY
    if (gradtol < TINY)
        gradtol <- TINY
    nwts <- c(3^-c(0:3), rep(0, 12))
    gwts <- 2^(-c(0:(length(nwts) - 1)))
    gwts <- (gwts^2 - nwts^2)^0.5
    gwts[12:16] <- gwts[12:16] * 10^-c(1:5)
    nits <- 0
    gsize <- gradtol + 1
    while (nits < maxit && gsize > gradtol) {
        arg <- 1 + z \%*\% lam
        wts1 <- as.vector(llogp(arg, 1/n))
        wts2 <- as.vector(-llogpp(arg, 1/n))^0.5
        grad <- as.matrix(-z * wts1)
        grad <- as.vector(apply(grad, 2, sum))
        gsize <- mean(abs(grad))
        hess <- z * wts2
        svdh <- svd(hess)
        if (min(svdh$d) < max(svdh$d) * svdtol)
            svdh$d <- svdh$d + max(svdh$d) * svdtol
        nstep <- svdh$v \%*\% (t(svdh$u)/svdh$d)
        nstep <- as.vector(nstep \%*\% matrix(wts1/wts2, n, 1))
        gstep <- -grad
        if (sum(nstep^2) < sum(gstep^2))
            gstep <- gstep * sum(nstep^2)^0.5/sum(gstep^2)^0.5
        ologelr <- -sum(llog(arg, 1/n))
        ninner <- 0
        for (i in 1:length(nwts)) {
            nlogelr <- logelr_owen(z, rep(0, p), lam + nwts[i] *
                nstep + gwts[i] * gstep)
            if (nlogelr < ologelr) {
                lam <- lam + nwts[i] * nstep + gwts[i] * gstep
                ninner <- i
                break
            }
        }
        nits <- nits + 1
        if (ninner == 0)
            nits <- maxit
        if (itertrace)
            print(c(lam, nlogelr, gsize, ninner))
    }
    list(logelr = nlogelr, lambda = lam * scale, grad = grad *
        scale, hess = t(hess) \%*\% hess * scale^2, wts = wts1,
        nits = nits)
  }
}

\keyword{ ~Empirical Likelihood }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
