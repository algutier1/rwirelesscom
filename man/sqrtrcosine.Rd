\name{sqrtrcosine}
\alias{sqrtrcosine}
\title{Square Root Raised Cosine}
\usage{
sqrtrcosine(x,B,Ns)
}
\description{
sqrtrcosine(x,B,Ns) = [sin((1-B)*pi*x[n]/Ns)/(pi*x[n]/Ns)+(4*B*x[n]/Ns)*cos((1+B)*pi*x[n]/Ns)]/((pi*x[n]/Ns)(1-(4*B*x[n]/Ns)^2))

Reference:

 1. S. Daumont, R. Basel, Y. Louet, "Root-Raised Cosine filter influences on PAPR distribution of single carrier signals", ISCCSP 2008, Malta, 12-14 March 2008.

 2. 3GPP TS 25.104 V6.8.0 (2004-12)
}
\arguments{
\item{x}{- input vector}

\item{B}{- roll-off factor}

\item{Ns}{- Ns samples per symbol}
}
\value{
Response (Double, Numeric vector) of the suare root raised rcosine function applied to the input vector x with roll-off factor B and Ns samples per symbol.
}

\examples{
Ns=16
B=0.5
hx=seq(-5*Ns,5*Ns,by=1)
h1=sqrtrcosine(hx,B,Ns=Ns)
plot(hx/Ns,h1, ylim=c(-0.2,1.1), xlim=c(-3,3), pch=19, cex=0.1, ylab="h", xlab="T", type="l")
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

