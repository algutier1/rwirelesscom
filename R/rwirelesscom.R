#' R Wireless Communications Package
#'
#' A basic wireless communications package for the simulation of digital communications in R. The package includes
#' modulation functions for BPSK, QPSK, 8-PSK, 16-QAM, 64-QAM, and a AWGN noise generation function. Additionally,
#' the package includes functions to plot an I (in phase) and Q (quadrature) scatter diagram, or density plot:
#' \itemize{
#' \item fNo(),
#' \item fbpskmod(), fbpskdemod(),
#' \item f8pskmod(), f8pskdemod(),
#' \item f16qammod(), f16qamdemod(),
#' \item f64qammod(), f64qamdemod,
#' \item iqscatterplot(), iqdensityplot()
#' }
#' Together these functions
#' enable the evaluation of respective bit error and symbol rates in an AWGN channel and for
#' easily viewing the respective signals and noise in a scatter plot or density plot.
#'
#' @docType package
#' @name rwirelesscom
NULL

#' fNo
#'
#' Generates a vector of normally distributed noise samples with mean of zero and noise spectral density (No/2), a.k.a. AWGN.
#' @param N: number of noise samples
#' @param No: noise spectral density
#' @param type: ="real" or "complex" defualts to real
#'
#' @export
fNo <- function(N,No,type="real") {
  if (type=="real") n = sqrt(No/2)*rnorm(N, 0, 1)
  else if(type=="complex")  n=sqrt(No/2)*complex(real=rnorm(N,0,1),imaginary=rnorm(N,0,1))
  else n=c(rep(0,N))
  return(n)
}

#' fbpskmod
#'
#' BPSK modulator receives a vector of bits, each with value 0 or 1, and outputs a
#' vector with values 1 and -1. An output value of -1 corresponds to an input value
#' of 0 and an output value of 1 corresponds to an input of 1.
#'  @name fbpskmod
#' @param bits: received vector of bits (0's and 1's).
#' @export
fbpskmod <- function(bits) {
  r <- sapply(bits, function(x) if (x==0) r=-1 else r=x)
  return(r)
}


#' fbpskdemod
#'
#' Receives a vector of real values, corresponding to a
#' BPSK modulated signal transmitted through a communications channel
#' (e.g., signal plus noise). An input value < 1 is mapped to an
#' output value of 0, otherwise to a value of 1.
#' @param r: received signal plus noise.
#' @export
fbpskdemod <- function(r) {
    r <- sapply(r,function(x) if (x>=0) r=1 else r=0)
    return(r)
}

#' fqpskmod
#'
#' Receives a vector of bits (1's and 0's). The 1's and 0's are
#' mapped to in-phase (real) and quadrature (imaginary) components.
#' Correspondingly, a bit of 1 is mapped to +1/sqrt(2), otherwise to -1/sqrt(2)
#' as follows:
#'  \tabular{cc}{
#' input \tab output \cr
#' 00 \tab  (-1 - 1j) / sqrt(2)  \cr
#' 01 \tab  (-1 + 1j) / sqrt(2) \cr
#' 10 \tab  (+1 - 1j) / sqrt(2) \cr
#' 11 \tab  (+1 + 1j) / sqrt(2)
#' }
#' @param bits: received vector of bits (0's and 1's).
#' @export
fqpskmod <- function(bits) {
  bi <- bits[seq(1,length(bits),2)]
  bq <- bits[seq(2,length(bits),2)]
  si <- (1/sqrt(2))*sapply(bi, function(x) if (x==0) r=-1 else r=x)
  sq <- (1/sqrt(2))*sapply(bq, function(x) if (x==0) r=-1 else r=x)
  s <- complex(real=si, imaginary=sq)
  return(s)
}

#' fqpskdemod
#'
#' Receives a vector of complex values, corresponding to a
#' BPSK modulated signal transmitted through a communications channel
#' (e.g., signal plus noise). The received signal, r, is mapped to its in-phase (real part) and
#' quadrature parts (imaginary part). If the in phase part is > 0 then the corresponding bit value output is 1,
#' otherwise 0. Similarly, if the quadrature part (imaginary) > 0 then the corresponding bit value
#' is 1, otherwise 0.
#' @param r: received signal plus noise.
#' @export
fqpskdemod <- function(r) {
  ri <- Re(r)
  rq <- Im(r)
  r2 <- matrix(rbind(ri,rq),length(ri)+length(rq),byrow=TRUE)
  r <- sapply(r2,function(x) if (x>=0) r=1 else r=0)
  return(r)
}

#' f8pskmod
#'
#' Receives a vector of bits (1's and 0's), which are mapped
#' to in-phase (real) and quadrature (imaginary) components according to
#' a Binary Reflective Gray Code (BRGC, see reference). Each received pair of bits
#' are are mapped to 8-PSK symbols,
#' with \eqn{E_{s}}{Es} (symbol energy) = 1.The BRGC mapping is as follows:
#' \tabular{cc}{
#' input \tab output \cr
#' 000 \tab \eqn{  0 }  \cr
#' 001 \tab \eqn{  \pi/4 } \cr
#' 011 \tab \eqn{  \pi/2} \cr
#' 010 \tab \eqn{  3 \pi/4} \cr
#' 110 \tab  \eqn{ pi} \cr
#' 111 \tab  \eqn{ -3 \pi/4} \cr
#' 101 \tab  \eqn{ - \pi/2} \cr
#' 100 \tab  \eqn{ - \pi/4} \cr
#' }
#' Reference: E. Agrell, J Lassing, E. Strom, and T. Ottosson, Gray Coding for Multilevel Constellations In Gaussian Noise, IEEE Transactions on Communications, Vol. 53, No. 1, January 2007
#' @param bits: received vector of bits (0's and 1's).
#' @export
#'
f8pskmod <- function(bits) {
  # receive symbolbits matrix (Nsym rows x Log2(M) cols )
  # transform to symbolcodes
  # transform to xreal, yimag
  M=8
  Nsymbols = length(bits) %/% log2(M)
  symbolbits<-matrix(bits,Nsymbols,log2(M),byrow=TRUE)
  m8=cbind(c(0,0,1),c(0,2,0),c(4,0,0))
  symbolbitsm8 = symbolbits %*% m8
  symbolcodes <- apply(symbolbitsm8,1,sum)
  s<-sapply(symbolcodes,ft8pskbitmap)
  return(s)
}

ft8pskbitmap <- function(x) {                                        #b2b1b0
  if (x == 0) r = complex(real=1,imaginary=0)                         # 000
  else if (x == 1) r = complex(real=1/sqrt(2),imaginary= 1/sqrt(2))   # 001
  else if (x == 3) r = complex(real=0, imaginary=1)                   # 011
  else if (x == 2) r = complex(real=-1/sqrt(2),imaginary=1/sqrt(2))   # 010
  else if (x == 6) r = complex(real=-1,imaginary=0)                   # 110
  else if (x == 7) r = complex(real=-1/sqrt(2), imaginary=-1/sqrt(2)) # 111
  else if (x == 5) r = complex(real=0, imaginary=-1)                  # 101
  else if (x == 4) r = complex(real=1/sqrt(2),imaginary=-1/sqrt(2))   # 100
  else {print(x)
    r=100}
  return(r)
}

#' f8pskdemod
#'
#' Receives a vector of complex values, r, corresponding to an
#' 8-PSK modulated signal transmitted through a communications channel
#' (e.g., signal plus noise). Three bits are output for each received symbol
#' according to the following decision rules
#' \tabular{cc}{
#' input \tab output \cr
#'  \eqn{ -\pi/8 \gee Arg(r) < \pi/8} \tab 000  \cr
#'  \eqn{  \pi/8 \ge Arg(r) < 3 \pi/8} \tab 001 \cr
#'  \eqn{ 3 \pi/8 \ge Arg(r) < 5 \pi/8} \tab 011 \cr
#'  \eqn{ 5 \pi/8 \ge Arg(r) < 7 \pi/8} \tab 010 \cr
#'  \eqn{ 7 \pi/8 \ge Arg(r) < 9 \pi/8} \tab 110 \cr
#'  \eqn{ -7 \pi/8 \ge Arg(r) < -5 \pi/8} \tab 111 \cr
#'  \eqn{ -5 \pi/8 \ge Arg(r) < -3 \pi/8} \tab 101 \cr
#'  \eqn{ -3 \pi/8 \ge Arg(r) < - \pi/8}\tab  101 \cr
#' }
#' @param r: received signal.
#' @export
f8pskdemod <- function(r) {
  r2<-sapply(r,fr8pskbitmap)
  bits2 <- as.vector(matrix(r2,1,length(r2)))
  return(bits2)
}

#' f8pskdemod
#'
#' Receives a vector of complex values, r, corresponding to an
#' 8-PSK modulated signal transmitted through a communications channel
#' (e.g., signal plus noise). Three bits are output for each received symbol
#' according to the following decision rules
#' \tabular{cc}{
#' input \tab output \cr
#'  \eqn{ -\pi/8 \gee Arg(r) < \pi/8} \tab 000  \cr
#'  \eqn{  \pi/8 \ge Arg(r) < 3 \pi/8} \tab 001 \cr
#'  \eqn{ 3 \pi/8 \ge Arg(r) < 5 \pi/8} \tab 011 \cr
#'  \eqn{ 5 \pi/8 \ge Arg(r) < 7 \pi/8} \tab 010 \cr
#'  \eqn{ 7 \pi/8 \ge Arg(r) < 9 \pi/8} \tab 110 \cr
#'  \eqn{ -7 \pi/8 \ge Arg(r) < -5 \pi/8} \tab 111 \cr
#'  \eqn{ -5 \pi/8 \ge Arg(r) < -3 \pi/8} \tab 101 \cr
#'  \eqn{ -3 \pi/8 \ge Arg(r) < - \pi/8}\tab  101 \cr
#' }
#' @param r: received signal.
#' @export
fr8pskbitmap <- function(r) {
  argr <- Arg(r)
  if ( argr >= -pi/8 && argr < pi/8 ) symbolbits = c(0,0,0)
  else if ( argr >= pi/8 && argr < 3*pi/8 ) symbolbits = c(0,0,1)
  else if ( argr >= 3*pi/8 && argr < 5*pi/8) symbolbits = c(0,1,1)
  else if ( argr >= 5*pi/8 && argr < 7*pi/8) symbolbits = c(0,1,0)
  else if ( (argr >= 7*pi/8 && argr < 9*pi/8 ) || (argr > -9*pi/8 && argr < -7*pi/8 )) symbolbits = c(1,1,0)
  else if ( argr >= -7*pi/8 && argr  < -5*pi/8) symbolbits = c(1,1,1)
  else if ( argr >= -5*pi/8 && argr  < -3*pi/8) symbolbits = c(1,0,1)
  else if ( argr >= -3*pi/8 && argr  < -pi/8) symbolbits = c(1,0,0)
  else symbolbits=c(-1,-1,-1)
  return(symbolbits)
}




#' f16qammod
#'
#' Receives a vector of bits
#' @param bits: vector of bits
#' @export
f16qammod <- function(bits) {
  # constellation quadrant x = 1+1j, 3+1j, 1+3j, 3+3j ...
  # Es = Sum((x * Conj(x)))/4 = 10
  # receive bits
  # transform to symbolbits (Nsym rows) length(bits)/4 x Log2(M) cols
  # transform to symbolcodes
  # transform to xreal(inphase), yimag(quatrature)
  M=16
  Nsymbols = length(bits) %/% log2(M)
  symbolbits<-matrix(bits,Nsymbols,log2(M),byrow=TRUE)
  m16=cbind(c(0,0,0,1),c(0,0,2,0),c(0,4,0,0),c(8,0,0,0))
  symbolbitsm16 = symbolbits %*% m16
  symbolcodes <- apply(symbolbitsm16,1,sum)
  s<-sapply(symbolcodes,ft16qambitmap)
  return(s)
}
ft16qambitmap <- function(x) {
  # b1b0  R := 00 -> -3, 01 -> -1, 11-> +1, 10 -> +3
  # b3b2  I := 00 -> -3, 01 -> -1, 11-> +1, 10 -> +3             b3b2 b1b0
  if (x == 0)      r = complex(real=-3, imaginary=-3)            #00 00
  else if (x == 1) r = complex(real=-1, imaginary=-3)            #00 01
  else if (x == 2) r = complex(real=+3, imaginary=-3)            #00 10
  else if (x == 3) r = complex(real=+1, imaginary=-3)            #00 11
  else if (x == 4) r = complex(real=-3, imaginary=-1)            #01 00
  else if (x == 5) r = complex(real=-1, imaginary=-1)            #01 01
  else if (x == 6) r = complex(real=+3, imaginary=-1)            #01 10
  else if (x == 7) r = complex(real=+1, imaginary=-1)            #01 11
  else if (x == 8) r = complex(real=-3, imaginary=+3)            #10 00
  else if (x == 9) r = complex(real=-1, imaginary=+3)            #10 01
  else if (x == 10) r = complex(real=+3,imaginary=+3)            #10 10
  else if (x == 11) r = complex(real=+1,imaginary=+3)            #10 11
  else if (x == 12) r = complex(real=-3,imaginary=+1)            #11 00
  else if (x == 13) r = complex(real=-1,imaginary=+1)            #11 01
  else if (x == 14) r = complex(real=+3,imaginary=+1)            #11 10
  else if (x == 15) r = complex(real=+1,imaginary=+1)            #11 11
  else {print(x)
    r=100}
  return(r)
}

#' f16qamdemod
#'
#' Receives a vector of bits
#' @param bits: vector of bits
#' @export



f16qamdemod <- function(r) {
  r2 <- sapply(r,fr16qambitmap)
  bitsr <- as.vector(matrix(r2,1,length(r2)))
  return(bitsr)
}
fr16qambitmap <- function(r) {
  if (Re(r) < -2 ) {
    if (Im(r) < -2) bits=c(0,0,0,0)
    else if (Im(r) < 0) bits=c(0,1,0,0)
    else if (Im(r) < 2) bits=c(1,1,0,0)
    else bits=c(1,0,0,0)
  }
  else if (Re(r) < 0 ) {
    if (Im(r) < -2) bits=c(0,0,0,1)
    else if (Im(r) < 0) bits=c(0,1,0,1)
    else if (Im(r) < 2) bits=c(1,1,0,1)
    else bits=c(1,0,0,1)
  }
  else if (Re(r) < 2) {
    if (Im(r) < -2) bits=c(0,0,1,1)
    else if (Im(r) < 0) bits=c(0,1,1,1)
    else if (Im(r) < 2) bits=c(1,1,1,1)
    else bits=c(1,0,1,1)
  }
  else if (Re(r) >=2 ) {
    if (Im(r) < -2) bits=c(0,0,1,0)
    else if (Im(r) < 0) bits=c(0,1,1,0)
    else if (Im(r) < 2) bits=c(1,1,1,0)
    else bits=c(1,0,1,0)
  }
  return(bits)
}

#' f64qammod
#'
#' Receives a vector of bits
#' @param bits: vector of bits
#' @export

f64qammod <- function(bits) {
  # x quad1 = 1+1i 1+3i 1+5i 1+7i 3+1i 3+3i 3+5i 3+7i 5+1i 5+3i 5+5i 5+7i 7+1i 7+3i 7+5i 7+7i
  # Es = sum(Re(x*Conj(x)))/16 = 42
  # receive bits
  # transform to symbolbits (Nsym rows) length(bits)/4 x Log2(M) cols
  # transform to symbolcodes
  # transform to xreal, yimag
  M=64
  Nsymbols = length(bits) %/% log2(M)
  symbolbits<-matrix(bits,Nsymbols,log2(M),byrow=TRUE)
  m64=cbind(c(0,0,0,0,0,1),c(0,0,0,0,2,0),c(0,0,0,4,0,0),c(0,0,8,0,0,0),c(0,16,0,0,0,0),c(32,0,0,0,0,0))
  symbolbits = symbolbits %*% m64
  symbolcodes <- apply(symbolbits,1,sum)
  s <- sapply(symbolcodes,ft64qambitmap)
  return(s)
}

ft64qambitmap <- function(x) {
                                                      #   imag      real
                                                      # b5b6b4b3    b2b1b0
  if (x == 0)      r = complex(real=-7, imaginary=-7)   #  000     000
  else if (x == 1) r = complex(real=-5, imaginary=-7)   #  000     001
  else if (x == 3) r = complex(real=-3, imaginary=-7)   #  000     011
  else if (x == 2) r = complex(real=-1, imaginary=-7)   #  000     010
  else if (x == 6) r = complex(real=+1, imaginary=-7)   #  000     110
  else if (x == 7) r = complex(real=+3, imaginary=-7)   #  000     111
  else if (x == 5) r = complex(real=+5, imaginary=-7)   #  000     101
  else if (x == 4) r = complex(real=+7, imaginary=-7)   #  000     100

  else if (x == 8) r = complex(real=-7, imaginary=-5)   #  001     000
  else if (x == 9) r = complex(real=-5, imaginary=-5)   #  001     001
  else if (x == 11) r = complex(real=-3, imaginary=-5)  #  001     011
  else if (x == 10) r = complex(real=-1, imaginary=-5)  #  001     010
  else if (x == 14) r = complex(real=+1, imaginary=-5)  #  001     110
  else if (x == 15) r = complex(real=+3, imaginary=-5)  #  001     111
  else if (x == 13) r = complex(real=+5, imaginary=-5)  #  001     101
  else if (x == 12) r = complex(real=+7, imaginary=-5)  #  001     100

  else if (x == 24) r = complex(real=-7, imaginary=-3)  #  011     000
  else if (x == 25) r = complex(real=-5, imaginary=-3)  #  011     001
  else if (x == 27) r = complex(real=-3, imaginary=-3)  #  011     011
  else if (x == 26) r = complex(real=-1, imaginary=-3)  #  011     010
  else if (x == 30) r = complex(real=+1, imaginary=-3)  #  011     110
  else if (x == 31) r = complex(real=+3, imaginary=-3)  #  011     111
  else if (x == 29) r = complex(real=+5, imaginary=-3)  #  011     101
  else if (x == 28) r = complex(real=+7, imaginary=-3)  #  011     100

  else if (x == 16) r = complex(real=-7, imaginary=-1)  #  010     000
  else if (x == 17) r = complex(real=-5, imaginary=-1)  #  010     001
  else if (x == 19) r = complex(real=-3, imaginary=-1)  #  010     011
  else if (x == 18) r = complex(real=-1, imaginary=-1)  #  010     010
  else if (x == 22) r = complex(real=+1, imaginary=-1)  #  010     110
  else if (x == 23) r = complex(real=+3, imaginary=-1)  #  010     111
  else if (x == 21) r = complex(real=+5, imaginary=-1)  #  010     101
  else if (x == 20) r = complex(real=+7, imaginary=-1)  #  010     100

  else if (x == 48) r = complex(real=-7, imaginary=+1)  #  110     000
  else if (x == 49) r = complex(real=-5, imaginary=+1)  #  110     001
  else if (x == 51) r = complex(real=-3, imaginary=+1)  #  110     011
  else if (x == 50) r = complex(real=-1, imaginary=+1)  #  110     010
  else if (x == 54) r = complex(real=+1, imaginary=+1)  #  110     110
  else if (x == 55) r = complex(real=+3, imaginary=+1)  #  110     111
  else if (x == 53) r = complex(real=+5, imaginary=+1)  #  110     101
  else if (x == 52) r = complex(real=+7, imaginary=+1)  #  110     100

  else if (x == 56) r = complex(real=-7, imaginary=+3)  #  111     000
  else if (x == 57) r = complex(real=-5, imaginary=+3)  #  111     001
  else if (x == 59) r = complex(real=-3, imaginary=+3)  #  111     011
  else if (x == 58) r = complex(real=-1, imaginary=+3)  #  111     010
  else if (x == 62) r = complex(real=+1, imaginary=+3)  #  111     110
  else if (x == 63) r = complex(real=+3, imaginary=+3)  #  111     111
  else if (x == 61) r = complex(real=+5, imaginary=+3)  #  111     101
  else if (x == 60) r = complex(real=+7, imaginary=+3)  #  111     100

  else if (x == 40) r = complex(real=-7, imaginary=+5)  #  101     000
  else if (x == 41) r = complex(real=-5, imaginary=+5)  #  101     001
  else if (x == 43) r = complex(real=-3, imaginary=+5)  #  101     011
  else if (x == 42) r = complex(real=-1, imaginary=+5)  #  101     010
  else if (x == 46) r = complex(real=+1, imaginary=+5)  #  101     110
  else if (x == 47) r = complex(real=+3, imaginary=+5)  #  101     111
  else if (x == 45) r = complex(real=+5, imaginary=+5)  #  101     101
  else if (x == 44) r = complex(real=+7, imaginary=+5)  #  101     100

  else if (x == 32) r = complex(real=-7, imaginary=+7)  #  100     000
  else if (x == 33) r = complex(real=-5, imaginary=+7)  #  100     001
  else if (x == 35) r = complex(real=-3, imaginary=+7)  #  100     011
  else if (x == 34) r = complex(real=-1, imaginary=+7)  #  100     010
  else if (x == 38) r = complex(real=+1, imaginary=+7)  #  100     110
  else if (x == 39) r = complex(real=+3, imaginary=+7)  #  100     111
  else if (x == 37) r = complex(real=+5, imaginary=+7)  #  100     101
  else if (x == 36) r = complex(real=+7, imaginary=+7)  #  100     100

  else {print(x)
    r=100}
  return(r)
}


#' f64qamdemod
#'
#' Receives a vector of bits
#' @param bits: vector of bits
#' @export
f64qamdemod <- function(r) {
  r2 <- sapply(r,fr64qambitmap)
  bitsr <- as.vector(matrix(r2,1,length(r2)))
  return(bitsr)
}

fr64qambitmap <- function(r) {
  if (Im(r) < -6 ) {
         if (Re(r) < -6) bits=c(0,0,0,0,0,0)
    else if (Re(r) < -4) bits=c(0,0,0,0,0,1)
    else if (Re(r) < -2) bits=c(0,0,0,0,1,1)
    else if (Re(r) < 0)  bits=c(0,0,0,0,1,0)
    else if (Re(r) < 2)  bits=c(0,0,0,1,1,0)
    else if (Re(r) < 4)  bits=c(0,0,0,1,1,1)
    else if (Re(r) < 6)  bits=c(0,0,0,1,0,1)
    else                 bits=c(0,0,0,1,0,0)
  }
  else if (Im(r) < -4 ) {
         if (Re(r) < -6) bits=c(0,0,1,0,0,0)
    else if (Re(r) < -4) bits=c(0,0,1,0,0,1)
    else if (Re(r) < -2) bits=c(0,0,1,0,1,1)
    else if (Re(r) < 0)  bits=c(0,0,1,0,1,0)
    else if (Re(r) < 2)  bits=c(0,0,1,1,1,0)
    else if (Re(r) < 4)  bits=c(0,0,1,1,1,1)
    else if (Re(r) < 6)  bits=c(0,0,1,1,0,1)
    else                 bits=c(0,0,1,1,0,0)
  }
  else if (Im(r) < -2 ) {
         if (Re(r) < -6) bits=c(0,1,1,0,0,0)
    else if (Re(r) < -4) bits=c(0,1,1,0,0,1)
    else if (Re(r) < -2) bits=c(0,1,1,0,1,1)
    else if (Re(r) < 0)  bits=c(0,1,1,0,1,0)
    else if (Re(r) < 2)  bits=c(0,1,1,1,1,0)
    else if (Re(r) < 4)  bits=c(0,1,1,1,1,1)
    else if (Re(r) < 6)  bits=c(0,1,1,1,0,1)
    else                 bits=c(0,1,1,1,0,0)
  }
  else if (Im(r) < 0 ) {
         if (Re(r) < -6) bits=c(0,1,0,0,0,0)
    else if (Re(r) < -4) bits=c(0,1,0,0,0,1)
    else if (Re(r) < -2) bits=c(0,1,0,0,1,1)
    else if (Re(r) < 0)  bits=c(0,1,0,0,1,0)
    else if (Re(r) < 2)  bits=c(0,1,0,1,1,0)
    else if (Re(r) < 4)  bits=c(0,1,0,1,1,1)
    else if (Re(r) < 6)  bits=c(0,1,0,1,0,1)
    else                 bits=c(0,1,0,1,0,0)
  }
  else if (Im(r) < 2 ) {
         if (Re(r) < -6)bits=c(1,1,0,0,0,0)
    else if (Re(r) < -4) bits=c(1,1,0,0,0,1)
    else if (Re(r) < -2) bits=c(1,1,0,0,1,1)
    else if (Re(r) < 0)  bits=c(1,1,0,0,1,0)
    else if (Re(r) < 2)  bits=c(1,1,0,1,1,0)
    else if (Re(r) < 4)  bits=c(1,1,0,1,1,1)
    else if (Re(r) < 6)  bits=c(1,1,0,1,0,1)
    else                 bits=c(1,1,0,1,0,0)
  }
  else if (Im(r) < 4 ) {
         if (Re(r) < -6) bits=c(1,1,1,0,0,0)
    else if (Re(r) < -4) bits=c(1,1,1,0,0,1)
    else if (Re(r) < -2) bits=c(1,1,1,0,1,1)
    else if (Re(r) < 0)  bits=c(1,1,1,0,1,0)
    else if (Re(r) < 2)  bits=c(1,1,1,1,1,0)
    else if (Re(r) < 4)  bits=c(1,1,1,1,1,1)
    else if (Re(r) < 6)  bits=c(1,1,1,1,0,1)
    else                 bits=c(1,1,1,1,0,0)
  }
  else if (Im(r) < 6 ) {
        if (Re(r) < -6)  bits=c(1,0,1,0,0,0)
    else if (Re(r) < -4) bits=c(1,0,1,0,0,1)
    else if (Re(r) < -2) bits=c(1,0,1,0,1,1)
    else if (Re(r) < 0)  bits=c(1,0,1,0,1,0)
    else if (Re(r) < 2)  bits=c(1,0,1,1,1,0)
    else if (Re(r) < 4)  bits=c(1,0,1,1,1,1)
    else if (Re(r) < 6)  bits=c(1,0,1,1,0,1)
    else                 bits=c(1,0,1,1,0,0)
  }
  else if (Im(r) >=6 ) {
         if (Re(r) < -6) bits=c(1,0,0,0,0,0)
    else if (Re(r) < -4) bits=c(1,0,0,0,0,1)
    else if (Re(r) < -2) bits=c(1,0,0,0,1,1)
    else if (Re(r) < 0)  bits=c(1,0,0,0,1,0)
    else if (Re(r) < 2)  bits=c(1,0,0,1,1,0)
    else if (Re(r) < 4)  bits=c(1,0,0,1,1,1)
    else if (Re(r) < 6)  bits=c(1,0,0,1,0,1)
    else                 bits=c(1,0,0,1,0,0)
  }
  return(bits)
}



iqscatterplot <- function(r) {
  ggplot(data.frame(cbind(Re(r),Im(r))), aes(x=X1, y=X2)) + geom_point(size=1)+ xlab("In-Phase") + ylab("Quadrature") + theme(axis.text=element_text(size=14),axis.title=element_text(size=12,face="bold"))
}

iqdensityplot <- function(r,iq="r") {
    if (iq!="i") { # Real Part
      ggplot(data.frame(r), aes(x=Re(r))) +  geom_histogram(binwidth=0.05, position="identity",aes(y=..density..))
    } else {    #Imaginary Part
      ggplot(data.frame(r), aes(x=Im(r))) +  geom_histogram(binwidth=0.05, position="identity",aes(y=..density..))
    }
 }



