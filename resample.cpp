#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
// Independent per-scenario xorshift64* stream; unbiased bounded integers.
struct RNG {
 uint64_t s;
 uint32_t next(){s^=s>>12;s^=s<<25;s^=s>>27;return (s*2685821657736338717ULL)>>32;}
 uint32_t bounded(uint32_t b){uint64_t m=uint64_t(next())*b;uint32_t l=m;if(l<b){uint32_t t=-b%b;while(l<t){m=uint64_t(next())*b;l=m;}}return m>>32;}
};
extern "C" void resample(const double* x,const double* y,int R,int nx,int ny,int B,uint64_t seed,
 double* pd,double* pt,double* bd,double* bt){
 RNG rng{seed?seed:1};int N=nx+ny;std::vector<double> z(N);std::vector<int> ix(N);
 for(int r=0;r<R;r++){
  auto a=x+r*nx;auto b=y+r*ny;double sx=0,sy=0,qx=0,qy=0;
  for(int i=0;i<nx;i++){sx+=a[i];qx+=a[i]*a[i];z[i]=a[i];}
  for(int i=0;i<ny;i++){sy+=b[i];qy+=b[i]*b[i];z[nx+i]=b[i];}
  double d=sx/nx-sy/ny, se=std::sqrt((qx-sx*sx/nx)/(nx-1)/nx+(qy-sy*sy/ny)/(ny-1)/ny);
  int cd=0,ct=0;
  for(int k=0;k<B;k++){
   // Independent groupwise nonparametric bootstrap.
   double ax=0,ay=0,vx=0,vy=0;
   for(int i=0;i<nx;i++){double v=a[rng.bounded(nx)];ax+=v;vx+=v*v;}
   for(int i=0;i<ny;i++){double v=b[rng.bounded(ny)];ay+=v;vy+=v*v;}
   double ds=ax/nx-ay/ny;
   double ses=std::sqrt(std::max(0.,(vx-ax*ax/nx)/(nx-1)/nx+(vy-ay*ay/ny)/(ny-1)/ny));
   bd[r*B+k]=ds;bt[r*B+k]=ses>0?(ds-d)/ses:std::numeric_limits<double>::quiet_NaN();
   // Draw an equiprobable subset without replacement by partial Fisher-Yates.
   for(int i=0;i<N;i++)ix[i]=i;
   ax=0;vx=0;
   for(int i=0;i<nx;i++){int j=i+rng.bounded(N-i);std::swap(ix[i],ix[j]);double v=z[ix[i]];ax+=v;vx+=v*v;}
   ay=sx+sy-ax;vy=qx+qy-vx;ds=ax/nx-ay/ny;
   ses=std::sqrt(std::max(0.,(vx-ax*ax/nx)/(nx-1)/nx+(vy-ay*ay/ny)/(ny-1)/ny));
   cd+=std::abs(ds)>=std::abs(d)-1e-12;
   ct+=std::abs(ds/ses)>=std::abs(d/se)-1e-12;
  }
  pd[r]=(1.+cd)/(B+1.);pt[r]=(1.+ct)/(B+1.);
 }
}
