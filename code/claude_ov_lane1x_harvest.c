// claude_ov_lane1x_harvest.c -- Lane 1 (overnight 2026-07-25): wide harvester for
// rational points of the contact-7 THREE-ROOT surface R(s,t,u)=0.
// Derived from code/claude_c7_harvest.c, with the tier-1 "has a root mod p"
// bitmap extended to ALL primes 7..113 (27 primes), which cuts the statistical
// false-positive rate from (2/3)^12 to (2/3)^27 = 5.1e-5 BEFORE the expensive
// tier-2 loop, and adds 40 tier-2 primes.  Sweeps pairs (v1,v2) of height <= H;
// reports pairs whose residual cubic has a rational root (=> a third rational
// Weierstrass root => quintic type [1,1,1,2] => 2-rank 3 => [2,2,14]) and,
// separately, SPLITALL candidates (all five roots rational => order 112).
//   usage: ./harvest H NTHREADS THREAD
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define NSM 27
static int SM[NSM]={7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113};
static uint8_t *HASROOT[NSM];
static int32_t *SINV[NSM]; static int32_t SI2[NSM];
#define NBG 40
static int32_t BG[NBG]; static int32_t *BINV[NBG]; static uint8_t *BQR[NBG]; static int32_t BI2[NBG];

static int32_t pmod(int64_t a,int32_t p){int64_t r=a%p;if(r<0)r+=p;return (int32_t)r;}
static int32_t mmul(int32_t a,int32_t b,int32_t p){return (int32_t)(((int64_t)a*b)%p);}
static int32_t mpw(int32_t a,int32_t e,int32_t p){int64_t r=1,x=a;while(e){if(e&1)r=r*x%p;x=x*x%p;e>>=1;}return (int32_t)r;}
static int isp(int n){if(n<2)return 0;for(int d=2;(int64_t)d*d<=n;d++)if(n%d==0)return 0;return 1;}

static inline int resid(int32_t n1,int32_t d1,int32_t n2,int32_t d2,int32_t p,int32_t*inv,int32_t i2,
                        int32_t*A,int32_t*B,int32_t*C){
  int32_t a1=pmod(d1,p), a2=pmod(d2,p); if(!a1||!a2) return 0;
  int32_t u=mmul(pmod(n1,p),inv[a1],p), w=mmul(pmod(n2,p),inv[a2],p);
  int32_t up=pmod(u+1,p), wp=pmod(w+1,p); if(!up||!wp) return 0;
  int32_t u2=mmul(u,u,p), w2=mmul(w,w,p); if(u2==w2) return 0;
  int32_t u3=mmul(u2,u,p),u5=mmul(u3,u2,p),w3=mmul(w2,w,p),w5=mmul(w3,w2,p);
  int32_t G1=mmul(pmod(-(int64_t)u5+u3+mmul(u2,i2,p),p), inv[mmul(up,up,p)], p);
  int32_t G2=mmul(pmod(-(int64_t)w5+w3+mmul(w2,i2,p),p), inv[mmul(wp,wp,p)], p);
  int32_t c4=mmul(pmod(G1-G2,p), inv[pmod(u2-w2,p)], p);
  int32_t c0=pmod(G1-(int64_t)mmul(c4,u2,p),p);
  int32_t q4=c4,q3=pmod(2LL*c4-1,p),q2=pmod((int64_t)c4+c0-i2,p);
  int32_t e1=pmod((int64_t)u+w,p), e2=mmul(u,w,p);
  *A=pmod((int64_t)q4+e1,p);
  *B=pmod((int64_t)q3-e2+(int64_t)mmul(e1,*A,p),p);
  *C=pmod((int64_t)q2-(int64_t)mmul(e2,*A,p)+(int64_t)mmul(e1,*B,p),p);
  return 1;
}
static inline int32_t cdisc(int32_t A,int32_t B,int32_t C,int32_t p){
  int64_t A2=(int64_t)A*A%p,A3=A2*A%p,B2=(int64_t)B*B%p,B3=B2*B%p,C2=(int64_t)C*C%p;
  int64_t d=(18LL*A%p)*B%p*C%p; d-=(4*A3%p)*C%p; d+=A2*B2%p; d-=4*B3%p; d-=27*C2%p; return pmod(d,p);
}
int main(int argc,char**argv){
  int H=argc>1?atoi(argv[1]):100, nthr=argc>2?atoi(argv[2]):1, thr=argc>3?atoi(argv[3]):0;
  for(int k=0;k<NSM;k++){int32_t p=SM[k];
    SINV[k]=malloc(4*p); SINV[k][0]=0; for(int32_t t=1;t<p;t++)SINV[k][t]=mpw(t,p-2,p); SI2[k]=SINV[k][2%p];
    HASROOT[k]=calloc((size_t)p*p*p,1);
    // mark by construction: for each (v,A,B) the unique C = -(v^3+Av^2+Bv)
    for(int32_t A=0;A<p;A++)for(int32_t B=0;B<p;B++)for(int32_t v=0;v<p;v++){
      int64_t val=((int64_t)v*v%p*v + (int64_t)A*v%p*v + (int64_t)B*v)%p;
      int32_t C=pmod(-val,p);
      HASROOT[k][((size_t)A*p+B)*p+C]=1; } }
  { int n=0,q=127; while(n<NBG){ if(isp(q)) BG[n++]=q; q+=2; } }
  for(int k=0;k<NBG;k++){int32_t p=BG[k]; BINV[k]=malloc(4*p);BINV[k][0]=0;for(int32_t t=1;t<p;t++)BINV[k][t]=mpw(t,p-2,p);
    BI2[k]=BINV[k][2]; BQR[k]=calloc(p,1); for(int32_t t=0;t<p;t++)BQR[k][mmul(t,t,p)]=1; }
  int cap=4*H*H+16; int32_t*NU=malloc(4*cap),*DE=malloc(4*cap); int nv=0;
  for(int den=1;den<=H;den++)for(int num=-H;num<=H;num++){int a=num<0?-num:num,b=den;while(b){int t=a%b;a=b;b=t;}
    if(a!=1)continue; if(num==0||num==den||num==-den)continue; NU[nv]=num;DE[nv]=den;nv++;}
  if(thr==0) fprintf(stderr,"H=%d #v=%d pairs=%.4g\n",H,nv,(double)nv*(nv-1)/2);
  long long pairs=0,t1=0,r3=0,r5=0;
  for(int i=thr;i<nv;i+=nthr) for(int j=i+1;j<nv;j++){
    if(DE[j]==DE[i]&&(NU[j]==NU[i]||NU[j]==-NU[i])) continue;
    pairs++;
    int32_t A,B,C; int ok=1;
    for(int k=0;k<NSM;k++){ int32_t p=SM[k];
      if(!resid(NU[i],DE[i],NU[j],DE[j],p,SINV[k],SI2[k],&A,&B,&C)) continue;
      if(!HASROOT[k][((size_t)A*p+B)*p+C]){ ok=0; break; } }
    if(!ok) continue;
    t1++;
    int hasroot=1, split=1;
    for(int k=0;k<NBG;k++){ int32_t p=BG[k];
      if(!resid(NU[i],DE[i],NU[j],DE[j],p,BINV[k],BI2[k],&A,&B,&C)) continue;
      int found=0; for(int32_t v=0;v<p;v++){ int64_t t=((int64_t)v*v%p*v+(int64_t)A*v%p*v+(int64_t)B*v+C)%p; if(t==0){found=1;break;} }
      if(!found){ hasroot=0; split=0; break; }
      int32_t D=cdisc(A,B,C,p); if(D!=0&&!BQR[k][D]) split=0; }
    if(hasroot){ r3++; printf("R3 %d/%d %d/%d%s\n",NU[i],DE[i],NU[j],DE[j],split?"  SPLITALL_CAND":""); fflush(stdout);
      if(split) r5++; }
  }
  printf("HARVEST_DONE thr=%d pairs=%lld tier1=%lld r3=%lld r5=%lld\n",thr,pairs,t1,r3,r5);
  return 0;
}
