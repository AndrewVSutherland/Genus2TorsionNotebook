/* Lane 8 (overnight 2026-07-25): Z/63 sieve on the contact-7 and contact-9 charts.
 *
 *  contact-9 (1 parameter a):   f = x^5 + (a^2-9)x^4 + (-105/8 a + 36)x^3
 *                                   + (63/4 a - 10479/256)x^2 + (-9a + 1449/64)x
 *                                   + (2a - 315/64)
 *  contact-7 (2 parameters a,b):f = x^5 + (b^2-7)x^4 + (2ab+21)x^3 + (a^2-7b-35)x^2
 *                                   + (-7a+2b+35)x + (2a - 35/4)
 *
 * Both carry a marked rational class of order 9 resp. 7 for free.  A rational
 * point of order 63 needs 63 | #J(F_p) at EVERY good prime; the chart landscape
 * (claude_ov_lane8_c79.gp) measures that density at ~1/7 per prime with no
 * obstruction, so a k-prime sieve cuts by ~7^-k.
 *
 * build: gcc -O3 -march=native -o c79sieve claude_ov_lane8_c79sieve.c
 * usage: ./c79sieve 9 <Hnum> <Hden>          (contact-9, a = n/d)
 *        ./c79sieve 7 <Hnum> <Hden>          (contact-7, a = n1/d1, b = n2/d2)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define NP 12
static int32_t primes[NP] = {11,13,17,19,23,29,31,37,41,43,47,53};

static int32_t P;
static int8_t  chi[64];
static int32_t NU;

static inline int32_t md(int64_t a){ a %= P; if(a<0) a+=P; return (int32_t)a; }
static inline int32_t mul(int32_t a,int32_t b){ return (int32_t)(((int64_t)a*b)%P); }
static int32_t inv(int32_t a){ int32_t r=1,e=P-2,b=a; while(e){ if(e&1) r=mul(r,b); b=mul(b,b); e>>=1;} return r; }

/* #J(F_p) for the monic quintic with coefficients c[0..4] (c5 = 1); 0 if singular */
static int64_t nJ(int32_t *c)
{
    int32_t i,j,N1,dg;
    /* squarefree check via gcd(f,f') on degree-5 monic f */
    { int32_t a[8],b[8],da=5,db,k;
      a[5]=1; for(i=0;i<5;i++) a[i]=c[i];
      for(i=0;i<5;i++) b[i]=md((int64_t)(i+1)*a[i+1]);
      db=4; while(db>=0 && b[db]==0) db--;
      if(db<0) return 0;
      while(db>=0){
        int32_t iv=inv(b[db]);
        while(da>=db){
          int32_t co=mul(a[da],iv);
          if(co) for(k=0;k<=db;k++) a[da-db+k]=md(a[da-db+k]-mul(co,b[k]));
          a[da]=0; da--; while(da>=0&&a[da]==0) da--;
          if(da<db) break;
        }
        for(k=0;k<8;k++){ int32_t t=a[k]; a[k]=b[k]; b[k]=t; }
        { int32_t t=da; da=db; db=t; }
        while(db>=0&&b[db]==0) db--;
      }
      if(da>0) return 0;
    }
    /* N1 over F_p : 1 point at infinity (odd degree) */
    N1 = 1;
    for(i=0;i<P;i++){
        int64_t v = 1;
        for(j=4;j>=0;j--) v = ( v*i + c[j] ) % P;
        N1 += 1 + chi[v];
    }
    /* N2 over F_{p^2} = F_p[w]/(w^2 - NU) : 1 point at infinity */
    { int64_t N2 = 1;
      for(i=0;i<P;i++) for(j=0;j<P;j++){
        int32_t A=1,B=0,k;
        for(k=4;k>=0;k--){
            int32_t nA = md((int64_t)A*i + (int64_t)NU*((int64_t)B*j%P));
            int32_t nB = md((int64_t)A*j + (int64_t)B*i);
            A = md(nA + c[k]); B = nB;
        }
        { int32_t nn = md((int64_t)A*A - (int64_t)NU*((int64_t)B*B%P));
          N2 += 1 + chi[nn]; }
      }
      { int64_t s1 = P+1-N1;
        int64_t s2n = N2 - (int64_t)P*P - 1 + s1*s1;
        if(s2n & 1) return -1;
        int64_t s2 = s2n/2;
        return (int64_t)P*P + 1 - ((int64_t)P+1)*s1 + s2; }
    }
}

static void setupP(int32_t p)
{
    int32_t i;
    P = p;
    for(i=0;i<P;i++) chi[i]=0;
    for(i=1;i<P;i++) chi[(int32_t)((int64_t)i*i%P)] = 1;
    for(i=1;i<P;i++) if(!chi[i]) chi[i] = -1;
    chi[0]=0;
    for(NU=2;NU<P;NU++) if(chi[NU]==-1) break;
}

static void c9coef(int32_t a, int32_t *c)
{
    int32_t i8=inv(8), i4=inv(4), i256=inv(256), i64=inv(64);
    c[4] = md((int64_t)a*a - 9);
    c[3] = md(-(int64_t)mul(md(105),mul(i8,a)) + 36);
    c[2] = md((int64_t)mul(md(63),mul(i4,a)) - mul(md(10479%P),i256));
    c[1] = md(-9LL*a + mul(md(1449%P),i64));
    c[0] = md(2LL*a - mul(md(315%P),i64));
}

static void c7coef(int32_t a, int32_t b, int32_t *c)
{
    int32_t i4=inv(4);
    c[4] = md((int64_t)b*b - 7);
    c[3] = md(2LL*mul(a,b) + 21);
    c[2] = md((int64_t)a*a - 7LL*b - 35);
    c[1] = md(-7LL*a + 2LL*b + 35);
    c[0] = md(2LL*a - mul(md(35),i4));
}

int main(int argc, char **argv)
{
    if(argc < 4){ fprintf(stderr,"usage: %s {7|9} Hnum Hden\n", argv[0]); return 1; }
    int32_t which = atoi(argv[1]);
    int64_t HN = atoll(argv[2]), HD = atoll(argv[3]);
    int32_t k, a, b, c[5];

    /* build the admissible-residue bitmaps */
    uint8_t *bm[NP];
    for(k=0;k<NP;k++){
        setupP(primes[k]);
        if(which==9){
            bm[k] = calloc(P,1);
            int32_t ng=0;
            for(a=0;a<P;a++){ c9coef(a,c); int64_t n=nJ(c); if(n>0 && n%63==0){ bm[k][a]=1; ng++; } }
            fprintf(stderr,"p=%d admissible a: %d/%d\n",P,ng,P);
        } else {
            bm[k] = calloc((size_t)P*P,1);
            int32_t ng=0;
            for(a=0;a<P;a++) for(b=0;b<P;b++){ c7coef(a,b,c); int64_t n=nJ(c);
                if(n>0 && n%63==0){ bm[k][(size_t)a*P+b]=1; ng++; } }
            fprintf(stderr,"p=%d admissible (a,b): %d/%d\n",P,ng,P*P);
        }
    }
    /* modular inverse tables for denominators */
    static int32_t invtab[NP][64];
    for(k=0;k<NP;k++){ setupP(primes[k]); for(a=1;a<P;a++) invtab[k][a]=inv(a); }

    int64_t checked=0, surv=0;
    if(which==9){
        for(int64_t d=1; d<=HD; d++){
            for(int64_t n=-HN; n<=HN; n++){
                int64_t g,x1=n<0?-n:n,y1=d; while(y1){int64_t t=x1%y1;x1=y1;y1=t;} g=x1?x1:d;
                if(g!=1) continue;
                checked++;
                int ok=1;
                for(k=0;k<NP;k++){
                    int32_t p=primes[k];
                    int32_t dm=(int32_t)(d%p); if(dm==0) continue;
                    int32_t nm=(int32_t)(((n%p)+p)%p);
                    int32_t r=(int32_t)(((int64_t)nm*invtab[k][dm])%p);
                    if(!bm[k][r]){ ok=0; break; }
                }
                if(ok){ surv++; printf("SURVIVOR c9 a = %lld/%lld\n",(long long)n,(long long)d); fflush(stdout); }
            }
            if(d%200==0){ printf("PROGRESS c9 d=%lld checked=%lld surv=%lld\n",(long long)d,(long long)checked,(long long)surv); fflush(stdout); }
        }
        printf("SEARCH_DONE c9 HN=%lld HD=%lld checked=%lld survivors=%lld\n",
               (long long)HN,(long long)HD,(long long)checked,(long long)surv);
    } else {
        for(int64_t d1=1; d1<=HD; d1++) for(int64_t n1=-HN; n1<=HN; n1++){
            { int64_t x1=n1<0?-n1:n1,y1=d1; while(y1){int64_t t=x1%y1;x1=y1;y1=t;} if((x1?x1:d1)!=1) continue; }
            for(int64_t d2=1; d2<=HD; d2++) for(int64_t n2=-HN; n2<=HN; n2++){
                { int64_t x1=n2<0?-n2:n2,y1=d2; while(y1){int64_t t=x1%y1;x1=y1;y1=t;} if((x1?x1:d2)!=1) continue; }
                checked++;
                int ok=1;
                for(k=0;k<NP;k++){
                    int32_t p=primes[k];
                    int32_t dm1=(int32_t)(d1%p), dm2=(int32_t)(d2%p);
                    if(dm1==0||dm2==0) continue;
                    int32_t am=(int32_t)((((n1%p)+p)%p*(int64_t)invtab[k][dm1])%p);
                    int32_t bm2=(int32_t)((((n2%p)+p)%p*(int64_t)invtab[k][dm2])%p);
                    if(!bm[k][(size_t)am*p+bm2]){ ok=0; break; }
                }
                if(ok){ surv++; printf("SURVIVOR c7 a = %lld/%lld  b = %lld/%lld\n",
                                       (long long)n1,(long long)d1,(long long)n2,(long long)d2); fflush(stdout); }
            }
            if(n1==HN && d1%5==0){ printf("PROGRESS c7 d1=%lld checked=%lld surv=%lld\n",(long long)d1,(long long)checked,(long long)surv); fflush(stdout); }
        }
        printf("SEARCH_DONE c7 HN=%lld HD=%lld checked=%lld survivors=%lld\n",
               (long long)HN,(long long)HD,(long long)checked,(long long)surv);
    }
    return 0;
}
