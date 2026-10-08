/* Exact multi-case conditional logistic recursion; no approximated denominator. */
#include <stdlib.h>
#include <math.h>

int conditional_eval(int nrows, int p, int groups, const double* x, const int* y,
    const int* starts, const double* beta, const double* weights,
    double* loss, double* gradient, double* scores, double* stratum_ll) {
  *loss=0.0;
  for(int j=0;j<p;++j) gradient[j]=0.0;
  for(int s=0;s<groups;++s) {
    int begin=starts[s],end=starts[s+1],n=end-begin,ncase=0;
    for(int i=begin;i<end;++i) ncase+=y[i];
    if(ncase<=0||ncase>=n) return 2;
    int complement=ncase>n/2,m=complement?n-ncase:ncase;
    double sign=complement?-1.0:1.0;
    double *eta=calloc(n,sizeof(double)), *d=calloc(m+1,sizeof(double));
    double *e=calloc((m+1)*p,sizeof(double)), *observed=calloc(p,sizeof(double));
    if(!eta||!d||!e||!observed) {free(eta);free(d);free(e);free(observed);return 3;}
    double shift=-INFINITY,numerator=0.0,logscale=0.0;
    for(int i=begin;i<end;++i) {
      double value=0.0;
      for(int j=0;j<p;++j) value+=x[i*p+j]*beta[j];
      eta[i-begin]=sign*value;
      if(sign*value>shift) shift=sign*value;
      if((complement&&!y[i])||(!complement&&y[i])) {
        numerator+=sign*value;
        for(int j=0;j<p;++j) observed[j]+=sign*x[i*p+j];
      }
    }
    d[0]=1.0;
    for(int i=0;i<n;++i) {
      double w=exp(eta[i]-shift);
      int upper=m<i+1?m:i+1;
      for(int k=upper;k>=1;--k) {
        for(int j=0;j<p;++j)
          e[k*p+j]+=w*(e[(k-1)*p+j]+d[k-1]*sign*x[(begin+i)*p+j]);
        d[k]+=w*d[k-1];
      }
      double maximum=0.0;
      for(int k=0;k<=m;++k) if(d[k]>maximum) maximum=d[k];
      if(maximum>1e100) {
        for(int k=0;k<=m;++k) d[k]/=maximum;
        for(int k=0;k<(m+1)*p;++k) e[k]/=maximum;
        logscale+=log(maximum);
      }
    }
    if(!(d[m]>0.0)||!isfinite(d[m])) {free(eta);free(d);free(e);free(observed);return 1;}
    double ll=numerator-log(d[m])-logscale-m*shift;
    stratum_ll[s]=ll; *loss-=weights[s]*ll;
    for(int j=0;j<p;++j) {
      double score=observed[j]-e[m*p+j]/d[m];
      scores[s*p+j]=score;gradient[j]-=weights[s]*score;
    }
    free(eta);free(d);free(e);free(observed);
  }
  return 0;
}
