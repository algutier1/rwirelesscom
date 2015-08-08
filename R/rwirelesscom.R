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
  si <- (1/sqrt(2))*sapply(bi, function(x) if (x==0) r=-1 else r=x)
  sq <- (1/sqrt(2))*sapply(bq, function(x) if (x==0) r=-1 else r=x)
  s <- complex(real=si, imaginary=sq)
  return(s)
}

fqpskdemod <- function(r) {
  ri <- Re(r)
  rq <- Im(r)
  r2 <- matrix(rbind(ri,rq),length(ri)+length(rq),byrow=TRUE)
  r <- sapply(r2,function(x) if (x>=0) r=1 else r=0)
  return(r)
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

f64qamdemod <- function(r) {
  r2 <- sapply(r,fr64qambitmap)
  bitsr <- as.vector(matrix(r2,1,length(r2)))
  return(bitsr)
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
