// scan230.c — diagnostic for target (2,30) on Elkies' A_1(5) threefold
// Universal curve: y^2 + (Q'-xQ) y = Q^2 Q',  Q = q2 x^2 + q1 x + q0, Q' = Q - x.
// Completed square: y^2 = F(x), F = (Q'-xQ)^2 + 4 Q^2 Q'  (deg <= 6).
// Built-in: rational 5-torsion class {Q=0,y=0} - K.
// (2,30) = Z/2 x Z/2 x Z/3 x Z/5 needs at every good odd prime p:
//   2-rank(J(F_p)) >= 2  and  3 | #J(F_p)   (5 | #J automatic; self-tested).
// Modes:
//   scan  : full F_p^3 enumeration for p in {7,11,13,17,19,23}, count allowed residues
//   base  : random-sextic baseline at same primes
//   sieve : rational triples q_i = n/d, |n|<=H, d<=D, sieved at primes 3..47
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int64_t P;                 // current prime
static int32_t chitab[64];        // quadratic character mod P: -1,0,1
static int64_t nonres;            // a quadratic nonresidue mod P

static inline int64_t addm(int64_t a, int64_t b){ int64_t r=a+b; return r>=P?r-P:r; }
static inline int64_t subm(int64_t a, int64_t b){ int64_t r=a-b; return r<0?r+P:r; }
static inline int64_t mulm(int64_t a, int64_t b){ return (a*b)%P; }

static int64_t powm(int64_t a, int64_t e){
  int64_t r=1; a%=P; if(a<0)a+=P;
  while(e){ if(e&1) r=mulm(r,a); a=mulm(a,a); e>>=1; }
  return r;
}
static int64_t invm(int64_t a){ return powm(a, P-2); }

static void set_prime(int64_t p){
  P=p;
  for(int64_t i=0;i<p;i++) chitab[i]=-1;
  chitab[0]=0;
  for(int64_t i=1;i<p;i++) chitab[mulm(i,i)]=1;
  nonres=2; while(chitab[nonres]!=-1) nonres++;
}

// polynomial degree with trimming; coefficients c[0..6]
static int pdeg(int64_t *c, int n){
  for(int i=n;i>=0;i--) if(c[i]) return i;
  return -1;
}

// squarefree test of poly c (deg d>=1) via gcd(F,F') mod P
static int is_squarefree(int64_t *c, int d){
  int64_t a[8], b[8]; int da=d, db;
  for(int i=0;i<=d;i++) a[i]=c[i];
  for(int i=1;i<=d;i++) b[i-1]=mulm(c[i], i%P);
  db=pdeg(b,d-1);
  if(db<0) return 0; // derivative zero (char divides all exps) -> not squarefree here
  while(db>=0){
    // a mod b
    int64_t lcinv=invm(b[db]);
    while(da>=db){
      int64_t t=mulm(a[da],lcinv);
      if(t){ for(int i=0;i<=db;i++) a[da-db+i]=subm(a[da-db+i], mulm(t,b[i])); }
      da--; while(da>=0 && !a[da]) da--;
    }
    if(da<0) return db==0; // b | a : squarefree iff gcd is constant
    // swap
    int64_t tmp[8]; int dt=da;
    for(int i=0;i<=da;i++) tmp[i]=a[i];
    da=db; for(int i=0;i<=db;i++) a[i]=b[i];
    db=dt; for(int i=0;i<=dt;i++) b[i]=tmp[i];
  }
  return 1;
}

// Build F coefficients from q0,q1,q2 (already reduced mod P). Returns degree (trimmed), F into c[0..6].
static int buildF(int64_t q0, int64_t q1, int64_t q2, int64_t *c){
  // h = -q2 x^3 + (q2-q1) x^2 + (q1-1-q0) x + q0
  int64_t h[4];
  h[3]=subm(0,q2); h[2]=subm(q2,q1); h[1]=subm(subm(q1,1),q0); h[0]=q0;
  // Q = [q0,q1,q2], Qp = Q - x
  int64_t Q[3]={q0,q1,q2}, Qp[3]={q0,subm(q1,1),q2};
  int64_t h2[7]={0,0,0,0,0,0,0}, Q2[5]={0,0,0,0,0}, f[7]={0,0,0,0,0,0,0};
  for(int i=0;i<4;i++) for(int j=0;j<4;j++) h2[i+j]=addm(h2[i+j], mulm(h[i],h[j]));
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) Q2[i+j]=addm(Q2[i+j], mulm(Q[i],Q[j]));
  for(int i=0;i<5;i++) for(int j=0;j<3;j++) f[i+j]=addm(f[i+j], mulm(Q2[i],Qp[j]));
  for(int i=0;i<7;i++) c[i]=addm(h2[i], mulm(4%P, f[i]));
  return pdeg(c,6);
}

// evaluate F at x in F_p
static inline int64_t evalF(int64_t *c, int d, int64_t x){
  int64_t r=c[d];
  for(int i=d-1;i>=0;i--) r=addm(mulm(r,x), c[i]);
  return r;
}

// result of local analysis
typedef struct { int good; int64_t nJ; int rank2; int rank; } locres;

// analyze curve y^2=F over F_P: good?, #J, 2-rank
static locres analyze(int64_t *c, int d){
  locres R={0,0,0,0};
  if(d<5) return R;            // boundary/degenerate for this chart
  if(!is_squarefree(c,d)) return R;
  R.good=1;
  int64_t n1=0, n2=0, r1=0, r2=0;
  for(int64_t x=0;x<P;x++){
    int64_t v=evalF(c,d,x);
    n1 += 1 + chitab[v];
    n2 += (v==0)?1:2;          // every element of F_p* is a square in F_{p^2}
    if(v==0) r1++;
  }
  r2=r1;
  // x = u + v*sqrt(nonres), v in 1..(p-1)/2 (conjugate pairs, double contribution)
  int64_t s=nonres;
  for(int64_t v=1;v<=(P-1)/2;v++){
    int64_t v2s=mulm(mulm(v,v),s);
    for(int64_t u=0;u<P;u++){
      // evaluate F(u+v*w), w^2=s: keep (a + b*w)
      int64_t a=c[d], b=0;
      for(int i=d-1;i>=0;i--){
        int64_t na=addm(addm(mulm(a,u), mulm(mulm(b,v),s)), c[i]);
        int64_t nb=addm(mulm(a,v), mulm(b,u));
        a=na; b=nb;
      }
      if(a==0 && b==0){ r2+=2; n2+=2; }
      else {
        int64_t nrm=subm(mulm(a,a), mulm(mulm(b,b),s)); (void)v2s;
        n2 += 2*(1 + chitab[nrm]);
      }
    }
  }
  if(d==6){
    n1 += (chitab[c[6]]==1)?2:0;
    n2 += 2;
  } else { // d==5
    n1 += 1; n2 += 1;
  }
  R.nJ = (n1*n1 + n2)/2 - P;
  // 2-rank from factorization type via (r1, r2, d)
  int64_t a1=r1, a2=(r2-r1)/2;
  int64_t m = d - a1 - 2*a2;
  int rank;
  if(d==6){
    int64_t k=a1+a2;
    if(m==0)      rank = (a1==0)? (int)k-1 : (int)k-2;
    else if(m==3) rank = (int)k-1;           // k+1 factors, odd part -> (k+1)-2
    else if(m==4) rank = (a1==0)? (int)k : (int)k-1;  // (k+1)-1 or (k+1)-2
    else if(m==5) rank = (int)k-1;
    else          rank = 0;                  // m==6: (6) or (3,3), both rank 0
  } else { // d==5: extra rational Weierstrass point at infinity (odd orbit)
    int64_t k=a1+a2 + ((m>0)?1:0);
    rank = (int)k-1;
  }
  if(rank<0) rank=0;
  R.rank=rank;
  R.rank2 = (rank>=2);
  return R;
}

static uint64_t rng_state=88172645463325252ULL;
static inline uint64_t rnd(void){ rng_state^=rng_state<<13; rng_state^=rng_state>>7; rng_state^=rng_state<<17; return rng_state; }

int main(int argc, char **argv){
  if(argc<2){ fprintf(stderr,"usage: scan230 scan|base|sieve [H D]\n"); return 1; }
  int64_t primes[]={7,11,13,17,19,23};
  if(!strcmp(argv[1],"scan")){
    printf("p, total, good, n5div(selftest), n3div, nrank2, nboth, frac_both, allowed(=bad+both), allowed_frac\n");
    for(int pi=0;pi<6;pi++){
      set_prime(primes[pi]);
      int64_t tot=0,good=0,n5=0,n3=0,nr2=0,nboth=0;
      int64_t c[8];
      for(int64_t q0=0;q0<P;q0++)for(int64_t q1=0;q1<P;q1++)for(int64_t q2=0;q2<P;q2++){
        tot++;
        int d=buildF(q0,q1,q2,c);
        locres R=analyze(c,d);
        if(!R.good) continue;
        good++;
        if(R.nJ%5==0) n5++;
        int c3=(R.nJ%3==0);
        if(c3) n3++;
        if(R.rank2) nr2++;
        if(c3&&R.rank2) nboth++;
      }
      int64_t allowed=(tot-good)+nboth;
      printf("%lld, %lld, %lld, %lld, %lld, %lld, %lld, %.4f, %lld, %.4f\n",
        (long long)P,(long long)tot,(long long)good,(long long)n5,(long long)n3,
        (long long)nr2,(long long)nboth,(double)nboth/good,(long long)allowed,(double)allowed/tot);
    }
  } else if(!strcmp(argv[1],"base")){
    printf("baseline: random sextics, 20000 per prime\n");
    printf("p, good, b5div, b3div, brank2, bboth\n");
    for(int pi=0;pi<6;pi++){
      set_prime(primes[pi]);
      int64_t good=0,b5=0,b3=0,br2=0,bb=0, c[8];
      for(int it=0;it<20000;it++){
        for(int i=0;i<7;i++) c[i]=rnd()%P;
        int d=pdeg(c,6);
        locres R=analyze(c,d);
        if(!R.good) continue;
        good++;
        if(R.nJ%5==0) b5++;
        int c3=(R.nJ%3==0);
        if(c3) b3++;
        if(R.rank2) br2++;
        if(c3&&R.rank2) bb++;
      }
      printf("%lld, %lld, %lld, %lld, %lld, %lld\n",(long long)P,(long long)good,
        (long long)b5,(long long)b3,(long long)br2,(long long)bb);
    }
  } else if(!strcmp(argv[1],"one")){
    set_prime(atoll(argv[2]));
    int64_t c[8];
    int d=buildF(atoll(argv[3]),atoll(argv[4]),atoll(argv[5]),c);
    locres R=analyze(c,d);
    printf("p=%lld deg=%d good=%d #J=%lld rank=%d\n",(long long)P,d,R.good,(long long)R.nJ,R.rank);
  } else if(!strcmp(argv[1],"sieve5")){
    // 2-parameter sieve on the boundary surface q2 = -1/4 (deg-5 chart)
    int64_t H = (argc>2)? atoll(argv[2]) : 60;
    int64_t Dm = (argc>3)? atoll(argv[3]) : 8;
    int64_t sprimes[]={3,5,7,11,13,17,19,23,29,31,37,41,43,47};
    int ns=14;
    typedef struct { int64_t n,d; } ratv;
    ratv *vals=malloc(sizeof(ratv)*8000); int nv=0;
    for(int64_t d=1;d<=Dm;d++)for(int64_t n=-H;n<=H;n++){
      int64_t aa=n<0?-n:n; if(aa==0){ if(d>1) continue; }
      int64_t x1=aa,y1=d,t; while(y1){t=x1%y1;x1=y1;y1=t;}
      if(aa>0 && x1>1) continue;
      vals[nv].n=n; vals[nv].d=d; nv++;
    }
    fprintf(stderr,"sieve5: %d values/coord, %lld pairs\n",nv,(long long)nv*nv);
    int64_t *red=malloc(sizeof(int64_t)*nv*ns);
    for(int vi=0;vi<nv;vi++)for(int si=0;si<ns;si++){
      int64_t p=sprimes[si];
      if(vals[vi].d%p==0){ red[vi*ns+si]=-1; continue; }
      P=p;
      int64_t n=vals[vi].n%p; if(n<0)n+=p;
      red[vi*ns+si]=mulm(n, invm(vals[vi].d%p));
    }
    int64_t survA=0, survB=0;
    for(int i0=0;i0<nv;i0++)for(int i1=0;i1<nv;i1++){
      int okA=1, okB=1, anygood=0;
      for(int si=0;si<ns && (okA||okB);si++){
        int64_t p=sprimes[si];
        int64_t r0=red[i0*ns+si], r1=red[i1*ns+si];
        if(r0<0||r1<0) continue;
        set_prime(p);
        int64_t q2m = mulm(P-1, invm(4%P));   // -1/4 mod p
        int64_t c[8];
        int d=buildF(r0,r1,q2m,c);
        locres R=analyze(c,d);
        if(!R.good) continue;
        anygood=1;
        if(!R.rank2){ okA=0; okB=0; break; }
        if(p!=3 && R.nJ%3!=0) okB=0;
      }
      if(okA && anygood){
        survA++;
        printf("%s q0=%lld/%lld q1=%lld/%lld q2=-1/4\n", okB?"SURV_B":"SURV_A",
          (long long)vals[i0].n,(long long)vals[i0].d,(long long)vals[i1].n,(long long)vals[i1].d);
        if(okB) survB++;
      }
    }
    fprintf(stderr,"DONE sieve5 survA=%lld survB=%lld\n",(long long)survA,(long long)survB);
    free(vals); free(red);
  } else if(!strcmp(argv[1],"sieve")){
    int64_t H = (argc>2)? atoll(argv[2]) : 15;
    int64_t Dm = (argc>3)? atoll(argv[3]) : 5;
    int64_t sprimes[]={3,5,7,11,13,17,19,23,29,31,37,41,43,47};
    int ns=14;
    // build list of rational values n/d as pairs
    typedef struct { int64_t n,d; } ratv;
    ratv *vals=malloc(sizeof(ratv)*4000); int nv=0;
    for(int64_t d=1;d<=Dm;d++)for(int64_t n=-H;n<=H;n++){
      int64_t g, aa=n<0?-n:n, bb=d; if(aa==0){ if(d>1) continue; }
      g=aa; int64_t t; int64_t x1=g,y1=bb; while(y1){t=x1%y1;x1=y1;y1=t;} g=x1?x1:1;
      if(aa>0 && g>1) continue;
      vals[nv].n=n; vals[nv].d=d; nv++;
    }
    fprintf(stderr,"sieve: %d values/coord, %lld triples, %d primes\n",nv,(long long)nv*nv*nv,ns);
    int64_t survA=0, survB=0, tried=0;
    // precompute reductions of each value mod each prime (-1 = denominator divisible by p)
    int64_t *red=malloc(sizeof(int64_t)*nv*ns);
    for(int vi=0;vi<nv;vi++)for(int si=0;si<ns;si++){
      int64_t p=sprimes[si];
      if(vals[vi].d%p==0){ red[vi*ns+si]=-1; continue; }
      P=p;
      int64_t n=vals[vi].n%p; if(n<0)n+=p;
      red[vi*ns+si]=mulm(n, invm(vals[vi].d%p));
    }
    for(int i0=0;i0<nv;i0++)for(int i1=0;i1<nv;i1++)for(int i2=0;i2<nv;i2++){
      tried++;
      int okA=1, okB=1, anygood=0;
      for(int si=0;si<ns && (okA||okB);si++){
        int64_t p=sprimes[si];
        int64_t r0=red[i0*ns+si], r1=red[i1*ns+si], r2v=red[i2*ns+si];
        if(r0<0||r1<0||r2v<0) continue; // bad denominator -> allowed
        set_prime(p);
        int64_t c[8];
        int d=buildF(r0,r1,r2v,c);
        locres R=analyze(c,d);
        if(!R.good) continue;          // bad reduction -> allowed
        anygood=1;
        if(!R.rank2){ okA=0; okB=0; break; }
        if(p!=3 && R.nJ%3!=0) okB=0;
      }
      if(okA && anygood){
        survA++;
        if(okB){
          survB++;
          printf("SURV_B q=(%lld/%lld, %lld/%lld, %lld/%lld)\n",
            (long long)vals[i0].n,(long long)vals[i0].d,(long long)vals[i1].n,
            (long long)vals[i1].d,(long long)vals[i2].n,(long long)vals[i2].d);
        } else {
          printf("SURV_A q=(%lld/%lld, %lld/%lld, %lld/%lld)\n",
            (long long)vals[i0].n,(long long)vals[i0].d,(long long)vals[i1].n,
            (long long)vals[i1].d,(long long)vals[i2].n,(long long)vals[i2].d);
        }
      }
      if(tried%20000000==0) fprintf(stderr,"progress %lld tried, A=%lld B=%lld\n",
        (long long)tried,(long long)survA,(long long)survB);
    }
    fprintf(stderr,"DONE tried=%lld survA=%lld survB=%lld\n",(long long)tried,(long long)survA,(long long)survB);
    free(vals); free(red);
  }
  return 0;
}
