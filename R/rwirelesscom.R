#' fNo function
#'
#' Generate a vector of normally distributed noise samples with mean of zero and noise spectral density (No/2), a.k.a. AWGN.
#' @param N number of noise samples
#' @param No: noise spectral density
#' @param Type: "real" or "complex" defualts to real

fNo <- function(N,No,Type="real") {
  if (Type=="real") n = sqrt(No/2)*rnorm(N, 0, 1)
  else if(Type=="complex")  n=sqrt(No/2)*complex(real=rnorm(N,0,1),imaginary=rnorm(N,0,1))
  else n=c(rep(0,N))
  return(n)
}

fbpskmod <- function(bits) {
  r <- sapply(bits, function(x) if (x==0) r=-1 else r=x)
  return(r)
}

fbpskdemod <- function(r) {
    r <- sapply(r,function(x) if (x>=0) r=1 else r=0)
    return(r)
}


fqpskmod <- function(bits) {
  bi <- bits[seq(1,length(bits),2)]
  bq <- bits[seq(2,length(bits),2)]
  si <- sapply(bi, function(x) if (x==0) r=-1 else r=x)
  sq <- sapply(bq, function(x) if (x==0) r=-1 else r=x)
  s <- complex(real=si, imaginary=sq)
  return(s)
}

fqpskdemod <- function(r) {
  ri <- Re(r)
  rq <- Im(r)
  r2 <- matrix(rbind(ri,rq),length(ri)+length(rq),byrow=TRUE)
  r3 <- sapply(r2,function(x) if (x>=0) r=1 else r=0)
  return(r3)
}

# 8 psk modulator w/ symbol mapping
# m8 symbol bit2symbol mapping
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

f8pskdemod <- function(r) {
  r2<-sapply(r,fr8pskbitmap)
  bits2 <- as.vector(matrix(r2,1,length(r2)))
  return(bits2)
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
f16qammod <- function(bits) {
  # receive bits
  # transform to symbolbits (Nsym rows) length(bits)/4 x Log2(M) cols
  # transform to symbolcodes
  # transform to xreal, yimag
  M=16
  Nsymbols = length(bits) %/% log2(M)
  symbolbits<-matrix(bits,Nsymbols,log2(M),byrow=TRUE)
  m16=cbind(c(0,0,0,1),c(0,0,2,0),c(0,4,0,0),c(8,0,0,0))
  symbolbitsm16 = symbolbits %*% m16
  symbolcodes <- apply(symbolbitsm16,1,sum)
  s<-sapply(symbolcodes,ft16qambitmap)
  return(s)
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


f16qamdemod <- function(r) {
  r2 <- sapply(r,fr16qambitmap)
  bitsr <- as.vector(matrix(r2,1,length(r2)))
  return(bitsr)
  }
