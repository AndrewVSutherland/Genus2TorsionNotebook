// Lane 2 (claude_ov_c7z): SPLIT-ALL sieve for the contact-7 order-112 target [2,2,2,14].
//
// Reformulation (derived + verified 2026-07-25, see notes/claude_ov_c7z_2026_07_25.md):
//   the chart quintic Q(v) satisfies, identically in (c4,c0),
//        Q(-1) = -1/2,   Q'(-1) = 3,   Q'(0) = 2 Q(0),
//   so with z_i = 1/(1+v_i) the five roots obey the SYMMETRIC system
//        sum z_i = SIG(=6),   prod z_i = PRD(=2),   sum 1/(1-z_i) = TAU(=3).
//   k rational z_i <-> quintic factor type with k rational roots <-> 2-rank k-1 (k<5)/4 (k=5).
//   k=5 (all z rational)  ==>  torsion [2,2,2,14], order 112, a new record.
//
// Given a pair (z1,z2), the residual cubic t^3 - A t^2 + B t - C is FORCED:
//    A = SIG - z1 - z2,  C = PRD/(z1 z2),  T = TAU - 1/(1-z1) - 1/(1-z2),
//    u = 1 - T,  Bn = T*(1-A-C) - 3 + 2A,  B = Bn/u.
// SPLIT-ALL <=> that cubic splits completely over Q
//           <=> disc > 0 over R  AND  disc is a rational square  AND  it has a rational root.
// disc = 18ABC - 4A^3C + A^2B^2 - 4B^3 - 27C^2 ;  with W = u^3*disc (division free)
//    W = 18*A*Bn*C*u^2 - 4*A^3*C*u^3 + A^2*Bn^2*u - 4*Bn^3 - 27*C^2*u^3,
//    sign(disc) = sign(W)*sign(u),  legendre(disc,p) = legendre(W,p)*legendre(u,p).
//
// Funnel:  (0) real test sign(disc)>0   [division-free doubles, cache-tiled]
//          (1) Legendre chain on disc over up to NPR primes (each halves)
//          (2) survivors printed as SPLITCAND for exact rational re-verification.
// Two independent necessary conditions -> false-positive rate ~ 2^-NPR per pair.
//
// Usage: ./splitsieve H nthreads [SIGn SIGd PRDn PRDd TAUn TAUd]   (defaults 6/1 2/1 3/1)
//        ./splitsieve --check n1 d1 n2 d2 [consts]   -> print A,B,C,disc data for one pair
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define NPR 60
static int32_t PR[NPR];
static int32_t *INVT[NPR];      // full inverse table mod p
static uint8_t *QRT[NPR];       // QR indicator mod p
#define NTAB 3                  // primes kept in hot streaming arrays
#define NCOLD (NPR-NTAB)         // remaining primes, interleaved cold block per value

static int64_t SIGn=6,SIGd=1,PRDn=2,PRDd=1,TAUn=3,TAUd=1;

static int isp(int32_t n){ if(n<2) return 0; for(int32_t d=2; d*d<=n; d++) if(n%d==0) return 0; return 1; }
static int32_t pmod(int64_t a,int32_t p){ int64_t r=a%p; if(r<0) r+=p; return (int32_t)r; }
static inline int32_t mm(int32_t a,int32_t b,int32_t p){ return (int32_t)(((int64_t)a*b)%p); }
static int32_t mpw(int32_t a,int32_t e,int32_t p){ int64_t r=1,x=a; while(e){ if(e&1) r=r*x%p; x=x*x%p; e>>=1;} return (int32_t)r; }

static void setup_primes(void){
  int n=0; int32_t q=101;
  while(n<NPR){ if(isp(q)) PR[n++]=q; q+=2; }
  for(int k=0;k<NPR;k++){
    int32_t p=PR[k];
    INVT[k]=malloc(4*(size_t)p); QRT[k]=calloc(p,1);
    INVT[k][0]=0; for(int32_t t=1;t<p;t++) INVT[k][t]=mpw(t,p-2,p);
    for(int32_t t=1;t<p;t++) QRT[k][mm(t,t,p)]=1;
  }
}

// value list
static int32_t *VN,*VD; static double *VZ,*VIZ,*VIT;
static int32_t *TZ[NTAB],*TIZ[NTAB],*TIT[NTAB];
static int32_t *COLD;   // COLD[(size_t)i*3*NCOLD + 3*(k-NTAB) + {0,1,2}]
static int nv=0;

// modular residual-cubic data for the pair; returns 0 if prime unusable
static inline int modpair(int k,int32_t n1,int32_t d1,int32_t n2,int32_t d2,
                          int32_t z1,int32_t iz1,int32_t it1,
                          int32_t z2,int32_t iz2,int32_t it2,
                          int32_t *pW,int32_t *pu){
  int32_t p=PR[k]; int32_t *iv=INVT[k];
  int32_t SG=mm(pmod(SIGn,p),iv[pmod(SIGd,p)],p);
  int32_t PD=mm(pmod(PRDn,p),iv[pmod(PRDd,p)],p);
  int32_t TA=mm(pmod(TAUn,p),iv[pmod(TAUd,p)],p);
  int32_t A=pmod((int64_t)SG-z1-z2,p);
  int32_t C=mm(PD,mm(iz1,iz2,p),p);
  int32_t T=pmod((int64_t)TA-it1-it2,p);
  int32_t u=pmod(1LL-T,p); if(!u) return 0;
  int32_t Bn=pmod((int64_t)mm(T,pmod(1LL-A-C,p),p)-3+2LL*A,p);
  int64_t A2=(int64_t)A*A%p, A3=A2*A%p, u2=(int64_t)u*u%p, u3=u2*u%p;
  int64_t Bn2=(int64_t)Bn*Bn%p, Bn3=Bn2*Bn%p, C2=(int64_t)C*C%p;
  int64_t W = 18LL*A%p*Bn%p*C%p*u2%p;
  W -= 4LL*A3%p*C%p*u3%p;
  W += A2*Bn2%p*u%p;
  W -= 4LL*Bn3%p;
  W -= 27LL*C2%p*u3%p;
  *pW=pmod(W,p); *pu=u; return 1;
}
static inline void getmod(int k,int32_t n,int32_t d,int32_t*z,int32_t*iz,int32_t*it){
  int32_t p=PR[k]; int32_t *iv=INVT[k];
  int32_t nn=pmod(n,p), dd=pmod(d,p);
  int32_t zz = dd? mm(nn,iv[dd],p) : -1;
  if(zz<0){ *z=-1; return; }
  int32_t o=pmod(1LL-zz,p);
  *z=zz; *iz = nn? iv[zz] : -1; *it = o? iv[o] : -1;
  if(!nn) *iz=-1;
  if(!o)  *it=-1;
}

static long long g_pairs=0,g_real=0,g_cand=0;

int main(int argc,char**argv){
  if(argc>1 && !strcmp(argv[1],"--check")){
    setup_primes();
    int32_t n1=atoi(argv[2]),d1=atoi(argv[3]),n2=atoi(argv[4]),d2=atoi(argv[5]);
    if(argc>=12){ SIGn=atoll(argv[6]);SIGd=atoll(argv[7]);PRDn=atoll(argv[8]);PRDd=atoll(argv[9]);TAUn=atoll(argv[10]);TAUd=atoll(argv[11]); }
    double S=(double)SIGn/SIGd,P=(double)PRDn/PRDd,TT=(double)TAUn/TAUd;
    double z1=(double)n1/d1,z2=(double)n2/d2;
    double A=S-z1-z2, C=P/(z1*z2), T=TT-1.0/(1-z1)-1.0/(1-z2), u=1-T;
    double Bn=T*(1-A-C)-3+2*A, B=Bn/u;
    double disc=18*A*B*C-4*A*A*A*C+A*A*B*B-4*B*B*B-27*C*C;
    printf("A=%.15g B=%.15g C=%.15g u=%.15g disc=%.15g\n",A,B,C,u,disc);
    double W=18*A*Bn*C*u*u-4*A*A*A*C*u*u*u+A*A*Bn*Bn*u-4*Bn*Bn*Bn-27*C*C*u*u*u;
    printf("W=%.15g  sign(W)*sign(u)=%d\n",W,((W>0)==(u>0))?1:-1);
    int nz=0,ok=1;
    for(int k=0;k<NPR;k++){
      int32_t z1m,iz1,it1,z2m,iz2,it2,Wm,um;
      getmod(k,n1,d1,&z1m,&iz1,&it1); getmod(k,n2,d2,&z2m,&iz2,&it2);
      if(z1m<0||iz1<0||it1<0||z2m<0||iz2<0||it2<0) continue;
      if(!modpair(k,n1,d1,n2,d2,z1m,iz1,it1,z2m,iz2,it2,&Wm,&um)) continue;
      if(!Wm) continue;
      nz++;
      int leg = (QRT[k][Wm]?1:-1)*(QRT[k][um]?1:-1);
      if(leg<0){ ok=0; printf("  disc is a NON-residue mod %d -> not a square\n",PR[k]); break; }
    }
    printf("legendre-chain: nonzero primes tested=%d passed=%d\n",nz,ok);
    return 0;
  }
  int H=argc>1?atoi(argv[1]):100;
  int nthr=argc>2?atoi(argv[2]):1;
  if(argc>=9){ SIGn=atoll(argv[3]);SIGd=atoll(argv[4]);PRDn=atoll(argv[5]);PRDd=atoll(argv[6]);TAUn=atoll(argv[7]);TAUd=atoll(argv[8]); }
  omp_set_num_threads(nthr);
  setup_primes();
  double S=(double)SIGn/SIGd,P=(double)PRDn/PRDd,TT=(double)TAUn/TAUd;
  printf("SPLITSIEVE H=%d threads=%d  SIG=%lld/%lld PRD=%lld/%lld TAU=%lld/%lld  primes=%d\n",
         H,nthr,SIGn,SIGd,PRDn,PRDd,TAUn,TAUd,NPR); fflush(stdout);

  // enumerate z = n/d, gcd(n,d)=1, 1<=d<=H, |n|<=H, n!=0, z!=1, sorted by height=max(|n|,d)
  long cap=(long)4*H*H+64;
  VN=malloc(4*cap);VD=malloc(4*cap);
  int *lvstart=malloc(4*(H+2));
  for(int h=1;h<=H;h++){
    lvstart[h]=nv;
    for(int d=1;d<=h;d++) for(int n=-h;n<=h;n++){
      if(n==0) continue;
      if(d!=h && (n!=h && n!=-h)) continue;      // height exactly h
      int a=n<0?-n:n,b=d; while(b){int t=a%b;a=b;b=t;} if(a!=1) continue;
      if(n==d) continue;                          // z=1 degenerate
      VN[nv]=n;VD[nv]=d;nv++;
    }
  }
  lvstart[H+1]=nv;
  printf("#values=%d\n",nv); fflush(stdout);
  VZ=malloc(8*(size_t)nv);VIZ=malloc(8*(size_t)nv);VIT=malloc(8*(size_t)nv);
  for(int k=0;k<NTAB;k++){ TZ[k]=malloc(4*(size_t)nv);TIZ[k]=malloc(4*(size_t)nv);TIT[k]=malloc(4*(size_t)nv); }
  COLD=malloc(4*(size_t)nv*3*NCOLD);
  if(!COLD){ printf("COLD alloc failed (%zu bytes)\n",4*(size_t)nv*3*NCOLD); return 1; }
  printf("cold table %.2f GB\n", 4.0*nv*3*NCOLD/1e9); fflush(stdout);
  #pragma omp parallel for schedule(static)
  for(int i=0;i<nv;i++){
    double z=(double)VN[i]/VD[i];
    VZ[i]=z;VIZ[i]=1.0/z;VIT[i]=1.0/(1.0-z);
    for(int k=0;k<NTAB;k++){ int32_t a,b,c; getmod(k,VN[i],VD[i],&a,&b,&c); TZ[k][i]=a;TIZ[k][i]=b;TIT[k][i]=c; }
    for(int k=NTAB;k<NPR;k++){ int32_t a,b,c; getmod(k,VN[i],VD[i],&a,&b,&c);
      COLD[(size_t)i*3*NCOLD+3*(k-NTAB)+0]=a; COLD[(size_t)i*3*NCOLD+3*(k-NTAB)+1]=b; COLD[(size_t)i*3*NCOLD+3*(k-NTAB)+2]=c; }
  }
  printf("tables built\n"); fflush(stdout);

  int JT=8192;
  if(argc>9) JT=atoi(argv[9]);
  double t0=omp_get_wtime();
  for(int h=1;h<=H;h++){
    int i0=lvstart[h], i1=lvstart[h+1];
    if(i1<=i0) continue;
    long long lp=0,lr=0,lc=0;
    int nnew=i1-i0, ntile=(i1+JT-1)/JT;
    long long ntask=(long long)nnew*ntile;
    #pragma omp parallel reduction(+:lp,lr,lc)
    {
      #pragma omp for schedule(dynamic,8)
      for(long long tk=0;tk<ntask;tk++){
        int jt=(int)(tk/nnew), i=i0+(int)(tk%nnew);
        int j0=jt*JT, j1=j0+JT; if(j1>i1) j1=i1;
        {
          double z1=VZ[i],iz1=VIZ[i],it1=VIT[i];
          double a1=S-z1, c1=P*iz1, t1=TT-it1;
          int jend = (j1<i)?j1:i;
          for(int j=j0;j<jend;j++){
            double A=a1-VZ[j];
            double C=c1*VIZ[j];
            double T=t1-VIT[j];
            double u=1.0-T;
            double Bn=T*(1.0-A-C)-3.0+2.0*A;
            double A2=A*A,Bn2=Bn*Bn,u2=u*u,u3=u2*u;
            double w1=18.0*A*Bn*C*u2, w2=-4.0*A2*A*C*u3, w3=A2*Bn2*u, w4=-4.0*Bn2*Bn, w5=-27.0*C*C*u3;
            double W=w1+w2+w3+w4+w5;
            lp++;
            if((W>0.0)!=(u>0.0)){
              double mag=fabs(w1)+fabs(w2)+fabs(w3)+fabs(w4)+fabs(w5);
              if(fabs(W) > 1e-9*mag) continue;   // real test (with cancellation guard)
            }
            lr++;
            // Legendre chain
            int nz=0,ok=1;
            for(int k=0;k<NPR;k++){
              int32_t za,ia,ta,zb,ib,tb;
              if(k<NTAB){ za=TZ[k][i];ia=TIZ[k][i];ta=TIT[k][i];zb=TZ[k][j];ib=TIZ[k][j];tb=TIT[k][j]; }
              else { size_t oi=(size_t)i*3*NCOLD+3*(k-NTAB), oj=(size_t)j*3*NCOLD+3*(k-NTAB);
                     za=COLD[oi];ia=COLD[oi+1];ta=COLD[oi+2]; zb=COLD[oj];ib=COLD[oj+1];tb=COLD[oj+2]; }
              if(za<0||ia<0||ta<0||zb<0||ib<0||tb<0) continue;
              int32_t Wm,um;
              if(!modpair(k,VN[i],VD[i],VN[j],VD[j],za,ia,ta,zb,ib,tb,&Wm,&um)) continue;
              if(!Wm) continue;
              nz++;
              int leg=(QRT[k][Wm]?1:-1)*(QRT[k][um]?1:-1);
              if(leg<0){ ok=0; break; }
            }
            if(ok && nz>=40){
              lc++;
              #pragma omp critical
              { printf("SPLITCAND %d/%d %d/%d  nz=%d\n",VN[i],VD[i],VN[j],VD[j],nz); fflush(stdout); }
            }
          }
        }
      }
    }
    g_pairs+=lp; g_real+=lr; g_cand+=lc;
    if(h%25==0 || h==H){
      printf("PROGRESS H<=%d  pairs=%lld realpass=%lld (%.4f%%) cands=%lld  %.1fs\n",
             h,g_pairs,g_real,100.0*g_real/(g_pairs?g_pairs:1),g_cand,omp_get_wtime()-t0);
      fflush(stdout);
    }
  }
  printf("SEARCH_DONE H=%d pairs=%lld realpass=%lld cands=%lld time=%.1fs\n",
         H,g_pairs,g_real,g_cand,omp_get_wtime()-t0);
  return 0;
}
