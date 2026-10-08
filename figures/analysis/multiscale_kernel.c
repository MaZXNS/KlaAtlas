#include <math.h>
#include <stdint.h>
#include <stdlib.h>

static uint64_t random64(uint64_t *s) {
    uint64_t z=(*s+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}
static int uniform_int(uint64_t *s,int n) {
    uint64_t limit=UINT64_MAX-(UINT64_MAX%(uint64_t)n),v;
    do {v=random64(s);} while(v>=limit);
    return (int)(v%(uint64_t)n);
}
static int icmp(const void*a,const void*b) {int x=*(const int*)a,y=*(const int*)b;return(x>y)-(x<y);}
static int dcmp(const void*a,const void*b) {double x=*(const double*)a,y=*(const double*)b;return(x>y)-(x<y);}
static int contains(const int*x,int n,int v) {for(int i=0;i<n;i++)if(x[i]==v)return 1;return 0;}
static void sample(const int*k,int n,int m,int*selected,int*scratch,uint64_t*s) {
    int size=m<=n/2?m:n-m,count=0;
    while(count<size) {int j=uniform_int(s,n);if(!contains(scratch,count,j))scratch[count++]=j;}
    if(m<=n/2)for(int i=0;i<m;i++)selected[i]=k[scratch[i]];
    else {int q=0;for(int i=0;i<n;i++)if(!contains(scratch,size,i))selected[q++]=k[i];}
    qsort(selected,(size_t)m,sizeof(int),icmp);
}
static double nearest(const int*x,int m) {
    double sum=x[1]-x[0]+x[m-1]-x[m-2];
    for(int i=1;i<m-1;i++){int a=x[i]-x[i-1],b=x[i+1]-x[i];sum+=a<b?a:b;}
    return sum/m;
}
static int scan(const int*x,int m,int width) {
    int right=0,best=0;
    for(int left=0;left<m;left++){
        while(right<m&&x[right]<=x[left]+width-1)right++;
        if(right-left>best)best=right-left;
    }
    return best;
}
/*
 * Layout: observed NN, null mean(log2 NN), null sd(log2 NN), NN p,
 * arithmetic null NN, geometric median null NN;
 * W observed maxima, W null-mean maxima, W upper-tail p;
 * B rows of centered null [log2 NN, maximum fraction at each width].
 * Width is an inclusive residue interval of exactly width aa.
 */
int multiscale_eval(const int*k,int n,const int*kla,int m,const int*widths,int w,
                    int bcount,uint64_t seed,double*out) {
    if(n<=m||m<2||w<1||bcount<2)return 1;
    int *chosen=malloc((size_t)m*sizeof(int)),*scratch=malloc((size_t)n*sizeof(int));
    double *logs=malloc((size_t)bcount*sizeof(double));
    if(!chosen||!scratch||!logs){free(chosen);free(scratch);free(logs);return 2;}
    int base=6+3*w,stride=w+1;uint64_t state=seed;
    out[0]=nearest(kla,m);
    for(int j=0;j<w;j++)out[6+j]=scan(kla,m,widths[j]);
    double logsum=0,arithmetic=0;int nn_extreme=0;
    for(int b=0;b<bcount;b++){
        sample(k,n,m,chosen,scratch,&state);
        double nn=nearest(chosen,m);logs[b]=log2(nn);
        out[base+b*stride]=logs[b];logsum+=logs[b];arithmetic+=nn;
        if(nn<=out[0])nn_extreme++;
        for(int j=0;j<w;j++)out[base+b*stride+1+j]=(double)scan(chosen,m,widths[j])/m;
    }
    out[1]=logsum/bcount;out[4]=arithmetic/bcount;double ss=0;
    for(int b=0;b<bcount;b++){double d=logs[b]-out[1];ss+=d*d;out[base+b*stride]=d;}
    out[2]=sqrt(ss/(bcount-1));out[3]=(1.0+nn_extreme)/(bcount+1);
    qsort(logs,(size_t)bcount,sizeof(double),dcmp);
    out[5]=exp2(bcount%2?logs[bcount/2]:.5*(logs[bcount/2-1]+logs[bcount/2]));
    for(int j=0;j<w;j++){
        double mean=0;int extreme=0;
        for(int b=0;b<bcount;b++){double value=out[base+b*stride+1+j];mean+=value;if(value>=out[6+j]/m)extreme++;}
        mean/=bcount;out[6+w+j]=mean*m;out[6+2*w+j]=(1.0+extreme)/(bcount+1);
        for(int b=0;b<bcount;b++)out[base+b*stride+1+j]-=mean;
    }
    free(chosen);free(scratch);free(logs);return 0;
}
