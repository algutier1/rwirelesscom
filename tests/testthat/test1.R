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
  bits <- sample(0:1,Nbits, replace=TRUE)
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
  bits <- sample(0:1,Nbits, replace=TRUE)
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
  bits <- sample(0:1,Nbits, replace=TRUE)
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
  bits <- sample(0:1,Nbits, replace=TRUE)

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
  bits <- sample(0:1,Nbits, replace=TRUE)
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
  bits <- sample(0:1,Nbits, replace=TRUE)
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

context("16-PSK error rate, Eb/No = 12 dB")
test_that("Test 16-PSK Eb/No = 12 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  skip_on_cran()
  M=16
  Es=1
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)
  s <- f16pskmod(bits)

  EbNodB=12
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f16pskdemod(r)
  biterrs<-bits[bitsr!=bits]
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: %d-PSK EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

  expect_true(Pberr < 0.01, info="16-PSK EbNodb=12, Pberr should be < 0.01")
  expect_true(Pberr > 0.006, info="16-BSK EbNodb=12, Pberr should be > 0.006")

} )

context("16-QAM error rate, Es/No = 8 dB")
test_that("Test 16-QAM Eb/No = 12 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  M=16
  Es=10
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)

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
  bits <- sample(0:1,Nbits, replace=TRUE)

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
  bits <- sample(0:1,Nbits, replace=TRUE)

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
  bits <- sample(0:1,Nbits, replace=TRUE)

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

context("Modulation w/ rcosine pulse shape")
test_that("Modulation pulse shaping with rcosine unit test", {
  skip_on_cran()

  h1=rcosine(hx,B=1,Ns=Ns)
  h0_25=rcosine(hx,B=0.25,Ns=Ns)

  M=2
  Ns=16
  hx=seq(-4*Ns,4*Ns,by=1)
  h0_5=rcosine(hx,B=0.5,Ns=Ns)
  Nsymbols=5
  Nbits=log2(M)*Nsymbols
  bpsk_bits <- sample(0:1,Nbits, replace=TRUE)
  s_bpsk <- fbpskmod(bpsk_bits,Ns,h0_5)
  x_bpsk=c(1.238326e-16,-1.082141e-03,-2.199782e-03,-3.279449e-03,-4.245149e-03,-5.024635e-03,-5.556007e-03,-5.794124e-03,-5.716294e-03,-5.326684e-03,-4.658960e-03,-3.776734e-03,-2.771533e-03,-1.758167e-03,-8.675462e-04,-2.372053e-04,0.000000e+00,8.105301e-04,1.062116e-03,6.380569e-04,-5.276463e-04,-2.435708e-03,-5.010069e-03,-8.090811e-03,-1.143259e-02,-1.470946e-02,-1.752667e-02,-1.943937e-02,-1.997811e-02,-1.868001e-02,-1.512443e-02,-8.971330e-03,0.000000e+00,1.293616e-02,2.866636e-02,4.677912e-02,6.663068e-02,8.734697e-02,1.078389e-01,1.268317e-01,1.429074e-01,1.545597e-01,1.602603e-01,1.585311e-01,1.480229e-01,1.275938e-01,9.638459e-02,5.388686e-02,-1.110223e-16,-6.600899e-02,-1.422704e-01,-2.274097e-01,-3.195888e-01,-4.165516e-01,-5.156895e-01,-6.141269e-01,-7.088205e-01,-7.966706e-01,-8.746376e-01,-9.398586e-01,-9.897597e-01,-1.022157e+00,-1.035342e+00,-1.028150e+00,-1.000000e+00,-9.520002e-01,-8.837286e-01,-7.963053e-01,-6.913630e-01,-5.709814e-01,-4.376048e-01,-2.939499e-01,-1.429074e-01,1.255843e-02,1.695056e-01,3.251034e-01,4.767094e-01,6.219323e-01,7.586777e-01,8.851772e-01,1.000000e+00,1.102048e+00,1.190537e+00,1.264967e+00,1.325081e+00,1.370825e+00,1.402298e+00,1.419713e+00,1.423357e+00,1.413562e+00,1.390678e+00,1.355069e+00,1.307105e+00,1.247170e+00,1.175683e+00,1.093111e+00,1.000000e+00,8.969968e-01,7.848741e-01,6.645488e-01,5.370937e-01,4.037421e-01,2.658814e-01,1.250378e-01,-1.714888e-02,-1.589589e-01,-2.986331e-01,-4.344185e-01,-5.646161e-01,-6.876283e-01,-8.020038e-01,-9.064770e-01,-1.000000e+00,-1.081766e+00,-1.151222e+00,-1.208070e+00,-1.252263e+00,-1.283984e+00,-1.303620e+00,-1.311736e+00,-1.309031e+00,-1.296310e+00,-1.274438e+00,-1.244311e+00,-1.206825e+00,-1.162849e+00,-1.113212e+00,-1.058689e+00,-1.000000e+00,-9.387933e-01,-8.745571e-01,-8.078457e-01,-7.392132e-01,-6.692286e-01,-5.984883e-01,-5.276234e-01,-4.573035e-01,-3.882333e-01,-3.211438e-01,-2.567783e-01,-1.958730e-01,-1.391342e-01,-8.721302e-02,-4.067990e-02,3.245267e-16,3.546979e-02,6.440063e-02,8.671748e-02,1.025236e-01,1.120989e-01,1.158890e-01,1.144874e-01,1.086096e-01,9.906186e-02,8.670679e-02,7.242628e-02,5.708509e-02,4.149633e-02,2.639103e-02,1.239293e-02,3.416071e-16,-9.445741e-03,-1.685953e-02,-2.219634e-02,-2.552118e-02,-2.699284e-02,-2.684459e-02,-2.536283e-02,-2.286518e-02,-1.967906e-02,-1.612208e-02,-1.248498e-02,-9.017944e-03,-5.920840e-03,-3.337447e-03,-1.353753e-03,-5.124106e-17,-2.372053e-04,-8.675462e-04,-1.758167e-03,-2.771533e-03,-3.776734e-03,-4.658960e-03,-5.326684e-03,-5.716294e-03,-5.794124e-03,-5.556007e-03,-5.024635e-03,-4.245149e-03,-3.279449e-03,-2.199782e-03,-1.082141e-03,0.000000e+00,-2.732857e-16,0.000000e+00,-2.732857e-16,-2.732857e-16,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.832142e-17,6.832142e-17,6.832142e-17,1.024821e-16,3.416071e-17,7.686159e-17)
  bpsk_err=max(abs(s_bpsk-x_bpsk))

  M=4
  Nsymbols=1
  Ns=4
  hx=seq(-2*Ns,2*Ns,by=1)
  h1=rcosine(hx,B=1,Ns=Ns)
  Nbits=log2(M)*Nsymbols
  qpsk_bits <- sample(0:1,Nbits, replace=TRUE)
  s_qpsk <- fqpskmod(qpsk_bits,Ns,h1)
  x_qpsk=c(-2.220446e-17+2.220446e-17i,-5.716294e-03-5.716294e-03i,-1.059389e-16-4.636454e-17i,1.714888e-02+1.714888e-02i,-6.345947e-17-1.056949e-16i,-1.200422e-01-1.200422e-01i,-3.535535e-01-3.535535e-01i,-6.002109e-01-6.002109e-01i,-7.071068e-01-7.071068e-01i,-6.002109e-01-6.002109e-01i,-3.535535e-01-3.535535e-01i,-1.200422e-01-1.200422e-01i,-7.169683e-17-1.312711e-16i,1.714888e-02+1.714888e-02i,-9.197175e-17-4.973636e-17i,-5.716294e-03-5.716294e-03i,4.440892e-17+0.000000e+00i,-8.881784e-17-1.332268e-16i,0.000000e+00+0.000000e+00i,-8.881784e-17+0.000000e+00i)
  qpsk_err=max(abs(s_qpsk-x_qpsk))

  expect_true(bpsk_err  < 1e-6, info="BPSK w/ rcosine unit test error should be < 1e-6")
  expect_true(qspsk_err < 1e-7, info="QPSK w/ rcosine unit test error should be < 1e-7")


} )
