#include <Rcpp.h>
using namespace Rcpp;

// removed from default build options
// --no-multiarch

// Learn more about how to use Rcpp at:
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//   https://support.rstudio.com/hc/en-us/articles/200486088-Using-Rcpp-with-RStudio

// Hello Example
List rcpp_hello() {
  CharacterVector x = CharacterVector::create("foo", "bar");
  NumericVector y   = NumericVector::create(0.0, 1.0);
  List z            = List::create(x, y);
  return z;
}

// Convolution Example
//          inf
//y(n) =    sum[x(k)h(n-k)]
//       k=-inf
// [[Rcpp::export]]
NumericVector conv(NumericVector h, NumericVector x) {
  int nh = h.size(), nx = x.size();
  int ny = nh + nx - 1;
  NumericVector y(ny);
   for (int n = 0; n < ny; n++)
     for (int k = 0; k <= n; k++)  y[n] += x[k]*h[n-k];
  return y;
}


/* sinc function http://www.inside-r.org/node/175318
 package phonTools */
/* error is shown below, no matching function for call to sin,
 however, this seems to compile just fine. I believe because
 it uses Rcpp sugar, not completely sure, which RStudio does
 not seem to know about*/
// [[Rcpp::export]]
NumericVector sinc(NumericVector x) {
  int nx = x.size();
  NumericVector y(nx);
  for (int n = 0; n < nx; n++)
    if (x[n]==0) y[n] = 1;
    else y[n] = sin(x[n])/x[n];
    return y;
}

/* Raised Cosine */
// [[Rcpp::export]]
NumericVector rcosine(NumericVector x, NumericVector B, NumericVector Ns) {
  /*int nx = rctype[1]; */
  int nx = x.size();
  NumericVector error(1);
  NumericVector y(nx);

  error[0] = 0;
  if (Ns.size() !=1 || B.size() != 1
          || B[0] < 0 || B[0] > 1 || Ns[0] <= 0 ) return error;
  double b=B[0];
  double y1=0, y2=0,x1=1,x2=1;
  int N = Ns[0];
  for (int n = 0; n < nx ; n++ ) {
      if (x[n]==0) {  y[n] = 1; }
      else if  ( (int)x[n] % (int)N == 0 ) { y[n]=0; }
      else if  ( x[n] ==  N/(2*b) || x[n] == -N/(2*b))  {
        x1=x[n]+0.0000000001;
        y1=(sin(M_PI*x1/N)/(M_PI*x1/N))*(cos(M_PI*b*x1/N))/(1-(2*b*x1/N)*(2*b*x1/N));
        x2=x[n]-0.0000000001;
        y2=(sin(M_PI*x2/N)/(M_PI*x2/N))*(cos(M_PI*b*x2/N))/(1-(2*b*x2/N)*(2*b*x2/N));
        y[n]= (y1+y2)/2;
         }
      else y[n]=(sin(M_PI*x[n]/N)/(M_PI*x[n]/N))*(cos(M_PI*b*x[n]/N))/(1-(2*b*x[n]/N)*(2*b*x[n]/N));
  }
  return y;
}

/* Root Raised Cosine */
// [[Rcpp::export]]
NumericVector sqrtrcosine(NumericVector x, NumericVector B, NumericVector Ns) {
  /*int nx = rctype[1]; */
  int nx = x.size();
  NumericVector error(1);
  NumericVector y(nx);

  error[0] = 0;
  if (Ns.size() !=1 || B.size() != 1
        || B[0] < 0 || B[0] > 1 || Ns[0] <= 0 ) return error;
  double b=B[0];
  int N = Ns[0];
  for (int n = 0; n < nx ; n++ ) {
    if (x[n]==0) {  y[n] = 1-b+4*b/M_PI; }
    else if  ( (x[n] == N/(4*b)) || (x[n] == -N/(4*b)) )  {  y[n]= (b/sqrt(2))*((1+2/M_PI)*sin(M_PI/(4*b))+(1-2/M_PI)*cos(M_PI/(4*b))); }
    else { y[n]=(sin(M_PI*(1-b)*x[n]/N)+(4*b*x[n]/N)*cos(M_PI*(1+b)*x[n]/N))/((M_PI*x[n]/N)*(1-(4*b*x[n]/N)*(4*b*x[n]/N)));  }
  }
  return y;
}

