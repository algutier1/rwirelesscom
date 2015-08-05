library("rwirelesscom")

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

  expect_true(Pberr < 0.02, info="BPSK EbNodb=4, Pberr should be < 0.015")
  expect_true(Pberr > 0.01, info="PBSK EbNodb=4, Pberr should be > 0.012")

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
  s <- sqrt(Es)*f8pskmod(bits)

  EbNodB=7
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,Type="complex")
  r <- s+n
  bitsr <- f8pskdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: %d-PSK EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

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
  s <- sqrt(Es)*f8pskmod(bits)

  EbNodB=10
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,Type="complex")
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
