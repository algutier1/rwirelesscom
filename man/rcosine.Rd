\name{rcosine}
\alias{rcosine}
\title{Raised Cosine}
\usage{
rcosine(x,B,Ns)
}
\description{
rcosine(x,B,Ns) = [sin(pi*x[n]/Ns)/(pi*x[n]/Ns)]*cos(pi*B*x[n]/Ns) /(1-(2*B*x[n]/Ns)^2)

Reference:  G. Proakis, Digital Communications, 3rd ed., New York: McGraw-Hill, 1995.
}
\arguments{
\item{x}{- input vector}

\item{B}{- roll-off factor}

\item{Ns}{- Ns samples per symbol}
}
\value{
Response (Double, Numeric vector) of the rcosine function applied to the input vector x with roll-off factor B and Ns samples per symbol.
}

\examples{
B=0
Ns=8
hx=seq(-3*Ns,3*Ns,by=1)
h=rcosine(hx,B,Ns)
plot(hx/Ns,h, ylim=c(-0.2,1), xlim=c(-3,3), pch=19, cex=0.1, ylab="h", xlab="T", type="l")
grid( col = "grey50", lty = "dotted")
}
\seealso{
Other rwirelesscom functions: \code{\link{f16pskdemod}};
  \code{\link{f16pskmod}}; \code{\link{f16qamdemod}};
  \code{\link{f16qammod}}; \code{\link{f64qamdemod}};
  \code{\link{f64qammod}}; \code{\link{f8pskdemod}};
  \code{\link{f8pskmod}}; \code{\link{fbpskdemod}};
  \code{\link{fqpskdemod}}; \code{\link{fqpskmod}};
  \code{\link{rcosine}}; \code{\link{sqrtrcosine}};
  \code{\link{iqdensityplot}}; \code{\link{iqscatterplot}};
  \code{\link{stemplot}}; \code{\link{eyediagram}}; \code{\link{sinc}}}

