library("rwirelesscom")
library("testthat")
library(ggplot2)

context("BPSK error rate, Eb/No = 4 dB")
test_that("Test BPSK Eb/No = 4 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  M=2
  Eb=1
  Es = log2(M)*Eb
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=T)
  s <- sqrt(Eb)*fbpskmod(bits)

  EbNodB=4
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No)
  r <- s+n
  bitsr <- fbpskdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: BPSK EbNo_dB = %d, Bits = %g, bit errors = %g, Pberr=%f",EbNodB, length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

  expect_true(Pberr < 0.03, info="BPSK EbNodb=4, Pberr should be < 0.03")
  expect_true(Pberr > 0.008, info="PBSK EbNodb=4, Pberr should be > 0.008")

} )

context("BPSK error rate, Eb/No = 8 dB")
test_that("Test BPSK EbNo_dB= 4, Modulator and Demodulator in AWGN, Bit Error Rate", {
  skip_on_cran()
  M=2
  Eb=1
  Es = log2(M)*Eb
  Nsymbols=100000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=T)
  s <- sqrt(Eb)*fbpskmod(bits)

  EbNodB=8
    No = Eb/(10^(EbNodB/10))
    n <- fNo(Nsymbols,No)
    r <- s+n
    bitsr <- fbpskdemod(r)
    biterrs<-bits[bitsr!=bits]
    b<-factor(bits)
    Pberr=length(biterrs)/length(bits)

    # str<-sprintf("Test: BPSK EbNo_dB = %d, Bits = %g, bit errors = %g, Pberr=%f",EbNodB, length(bits), length(biterrs),Pberr)
    #print("",quote=FALSE)
    #print(str,quote=FALSE)
    expect_true(Pberr < 0.00032, info="BPSK EbNodb=8, Pberr should be < 0.0032")
    expect_true(Pberr > 0.0001, info="PBSK EbNodb=8, Pberr should be > 0.001")

} )

context("QPSK error rate, Eb/No = 4 dB")
test_that("Test BPSK Eb/No = 4 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  M=4
  Es=1
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=T)
  s <- sqrt(Es)*fqpskmod(bits)

  EbNodB=4
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- fqpskdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: BPSK EbNo_dB = %d, Bits = %g, bit errors = %g, Pberr=%f",EbNodB, length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

  expect_true(Pberr < 0.02, info="BPSK EbNodb=4, Pberr should be < 0.015")
  expect_true(Pberr > 0.007, info="PBSK EbNodb=4, Pberr should be > 0.012")

} )

context("QPSK error rate, Eb/No = 8 dB")
test_that("Test BPSK EbNo_dB= 4, Modulator and Demodulator in AWGN, Bit Error Rate", {
  skip_on_cran()
  M=4
  Es=1
  Eb = Es/log2(M)
  Nsymbols=100000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=T)

  s <- sqrt(Es)*fqpskmod(bits)

  EbNodB=8
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- fqpskdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: BPSK EbNo_dB = %d, Bits = %g, bit errors = %g, Pberr=%f",EbNodB, length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)
  expect_true(Pberr < 0.00032, info="BPSK EbNodb=8, Pberr should be < 0.0032")
  expect_true(Pberr > 0.0001, info="PBSK EbNodb=8, Pberr should be > 0.001")

} )


context("8-PSK error rate, Eb/No = 7 dB")
test_that("Test 8-PSK Eb/No = 11 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  M=8
  Es=1
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=T)
  s <- f8pskmod(bits)

  EbNodB=7
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f8pskdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  # str<-sprintf("Test: %d-PSK EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  # print("",quote=FALSE)
  # print(str,quote=FALSE)

  expect_true(Pberr < 0.015, info="8-PSK EbNodb=7, Pberr should be < 0.015")
  expect_true(Pberr > 0.01, info="8-BSK EbNodb=7, Pberr should be > 0.01")

} )

context("8-PSK error rate, Eb/No = 10 dB")
test_that("Test 8-PSK Eb/No = 11 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  skip_on_cran()
  M=8
  Es=1
  Eb = Es/log2(M)
  Nsymbols=100000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=T)
  s <- f8pskmod(bits)

  EbNodB=10
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f8pskdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: %d-PSK EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

  expect_true(Pberr < 0.0011, info="8-PSK EbNodb=7, Pberr should be < 0.015")
  expect_true(Pberr > 0.0007, info="8-BSK EbNodb=7, Pberr should be > 0.01")

} )

context("16-QAM error rate, Es/No = 8 dB")
test_that("Test 16-QAM Eb/No = 12 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  M=16
  Es=10
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=T)

  s <- f16qammod(bits)
  EbNodB=8
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f16qamdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

   #str<-sprintf("Test: %d-QAM EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
   #print("",quote=FALSE)
   #print(str,quote=FALSE)

  expect_true(Pberr < 0.012, info="8-PSK EbNodb=8, Pberr should be < 0.012")
  expect_true(Pberr > 0.008, info="8-BSK EbNodb=8, Pberr should be > 0.008")

} )

context("16-QAM error rate, Es/No = 10 dB")
test_that("Test 16-QAM Eb/No = 10 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  skip_on_cran()
  M=16
  Es=10
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=T)

  s <- f16qammod(bits)
  EbNodB=10
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f16qamdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: %d-QAM EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

  expect_true(Pberr < 0.003, info="8-PSK EbNodb=10, Pberr should be < 0.005")
  expect_true(Pberr > 0.001, info="8-BSK EbNodb=10, Pberr should be > 0.002")

} )

context("64-QAM error rate, Es/No = 12 dB")
test_that("Test 64-QAM Eb/No = 12 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  M=64
  Es=42
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=T)

  #  Nbits=24
  #  Nsymbols=24/log2(M)
  #  bits=c(0,1,0,1,0,0,0,0,1,0,1,1, 1,0,0,0,1,1, 0,1,0,1,0,1)

  s <- f64qammod(bits)

  EbNodB=12
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f64qamdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  # str<-sprintf("Test: %d-PSK EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  # print("",quote=FALSE)
  # print(str,quote=FALSE)

  expect_true(Pberr < 0.012, info="8-PSK EbNodb=12, Pberr should be < 0.012")
  expect_true(Pberr > 0.008, info="8-BSK EbNodb=12, Pberr should be > 0.008")

} )

context("64-QAM error rate, Es/No = 14 dB")
test_that("Test 64-QAM Eb/No = 14 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  skip_on_cran()
  M=64
  Es=42
  Eb = Es/log2(M)
  Nsymbols=100000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=T)

  #  Nbits=24
  #  Nsymbols=24/log2(M)
  #  bits=c(0,1,0,1,0,0,0,0,1,0,1,1, 1,0,0,0,1,1, 0,1,0,1,0,1)

  s <- f64qammod(bits)

  EbNodB=14
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f64qamdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: %d-PSK EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

  expect_true(Pberr < 0.003, info="8-PSK EbNodb=14, Pberr should be < 0.003")
  expect_true(Pberr > 0.001, info="8-BSK EbNodb=14, Pberr should be > 0.001")

} )
