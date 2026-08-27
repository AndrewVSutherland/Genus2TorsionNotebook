/* Lane 8 (2026-07-25, third session): a MUCH bigger Z/63 sieve on the contact
 * charts.  Rewrite of claude_ov_lane8_c79sieve.c with
 *   - up to 16 primes instead of 12,
 *   - the coprime (n,d) pairs and their residues precomputed once,
 *   - a staged prime-by-prime filter (pass 1 over all j, then only survivors),
 *   - pthreads over the outer parameter.
 * The previous run reached |n|,d <= 100 on contact-7 (1.48e8 points, 0
 * survivors); with the same machine this version reaches |n|,d <= 1000
 * (1.4e12 points), a factor of ~10^4.
 *
 *   contact-9 (1 parameter a):   f = x^5 + (a^2-9)x^4 + (-105/8 a + 36)x^3
 *                                    + (63/4 a - 10479/256)x^2 + (-9a + 1449/64)x
 *                                    + (2a - 315/64)
 *   contact-7 (2 parameters a,b):f = x^5 + (b^2-7)x^4 + (2ab+21)x^3 + (a^2-7b-35)x^2
 *                                    + (-7a+2b+35)x + (2a - 35/4)
 * Both carry a marked rational class of order 9 resp. 7 for free; the sieve
 * imposes 63 | #J(F_p) at every prime.
 *
 * build: gcc -O3 -march=native -pthread -o /tmp/c63sieve code/claude_ov_lane8_c63sieve.c
 * usage: ./c63sieve {7|9} <Hnum> <Hden> <nthreads> [nprimes]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>

#define MAXNP 16
static int32_t allprimes[MAXNP] = {11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71};
static int32_t NP = 12;

static int32_t P;
static int8_t  chi[128];
static int32_t NU;

static inline int32_t md(int64_t a){ a %= P; if(a<0) a+=P; return (int32_t)a; }
static inline int32_t mul(int32_t a,int32_t b){ return (int32_t)(((int64_t)a*b)%P); }
static int32_t inv(int32_t a){ int32_t r=1,e=P-2,b=a; while(e){ if(e&1) r=mul(r,b); b=mul(b,b); e>>=1;} return r; }

/* #J(F_p) for the monic quintic with coefficients c[0..4] (c5 = 1); 0 if singular */
static int64_t nJ(int32_t *c)
{
    int32_t i,j,N1;
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
    N1 = 1;
    for(i=0;i<P;i++){
        int64_t v = 1;
        for(j=4;j>=0;j--) v = ( v*i + c[j] ) % P;
        N1 += 1 + chi[v];
    }
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

/* ---------------- pair table ---------------- */
static int64_t  NPAIR;
static int32_t *pn, *pd;          /* numerator, denominator */
static uint8_t *res[MAXNP];       /* residue of n/d mod p, or 255 when p | d */
static uint8_t *bm[MAXNP];        /* admissible bitmaps */

static int64_t HNg, HDg;
static int32_t NTH;
static int64_t g_checked = 0, g_surv = 0;
static pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;

typedef struct { int32_t id; } targ;

static void *worker7(void *v)
{
    targ *T = (targ*)v;
    int64_t i, j, k;
    int64_t checked = 0, surv = 0;
    int64_t *cand = malloc(sizeof(int64_t)*NPAIR);
    int64_t *cand2 = malloc(sizeof(int64_t)*NPAIR);
    for(i = T->id; i < NPAIR; i += NTH){
        int64_t nc = 0;
        /* prime 0 pass over all j */
        { uint8_t ai = res[0][i];
          if(ai == 255){ for(j=0;j<NPAIR;j++) cand[nc++] = j; }
          else { uint8_t *row = bm[0] + (size_t)ai*allprimes[0];
                 uint8_t *rj = res[0];
                 for(j=0;j<NPAIR;j++){ uint8_t b = rj[j]; if(b==255 || row[b]) cand[nc++] = j; } }
        }
        for(k = 1; k < NP && nc; k++){
            uint8_t ai = res[k][i];
            if(ai == 255) continue;
            uint8_t *row = bm[k] + (size_t)ai*allprimes[k];
            uint8_t *rj = res[k];
            int64_t m = 0, t;
            for(t = 0; t < nc; t++){ int64_t jj = cand[t]; uint8_t b = rj[jj];
                if(b==255 || row[b]) cand2[m++] = jj; }
            memcpy(cand, cand2, sizeof(int64_t)*m);
            nc = m;
        }
        checked += NPAIR;
        for(j = 0; j < nc; j++){
            int64_t jj = cand[j];
            surv++;
            pthread_mutex_lock(&mtx);
            printf("SURVIVOR c7 a = %d/%d  b = %d/%d\n", pn[i], pd[i], pn[jj], pd[jj]);
            fflush(stdout);
            pthread_mutex_unlock(&mtx);
        }
        if((i / NTH) % 4096 == 0 && T->id == 0){
            printf("PROGRESS c7 i=%lld/%lld checked=%lld surv=%lld\n",
                   (long long)i,(long long)NPAIR,(long long)checked,(long long)surv);
            fflush(stdout);
        }
    }
    free(cand); free(cand2);
    pthread_mutex_lock(&mtx); g_checked += checked; g_surv += surv; pthread_mutex_unlock(&mtx);
    return 0;
}

static void *worker9(void *v)
{
    targ *T = (targ*)v;
    int64_t i, k;
    int64_t checked = 0, surv = 0;
    for(i = T->id; i < NPAIR; i += NTH){
        int ok = 1;
        for(k = 0; k < NP; k++){
            uint8_t a = res[k][i];
            if(a == 255) continue;
            if(!bm[k][a]){ ok = 0; break; }
        }
        checked++;
        if(ok){ surv++;
            pthread_mutex_lock(&mtx);
            printf("SURVIVOR c9 a = %d/%d\n", pn[i], pd[i]); fflush(stdout);
            pthread_mutex_unlock(&mtx);
        }
    }
    pthread_mutex_lock(&mtx); g_checked += checked; g_surv += surv; pthread_mutex_unlock(&mtx);
    return 0;
}

static int64_t gcd64(int64_t a, int64_t b){ if(a<0) a=-a; while(b){ int64_t t=a%b; a=b; b=t; } return a; }

int main(int argc, char **argv)
{
    if(argc < 5){ fprintf(stderr,"usage: %s {7|9} Hnum Hden nthreads [nprimes]\n", argv[0]); return 1; }
    int32_t which = atoi(argv[1]);
    HNg = atoll(argv[2]); HDg = atoll(argv[3]); NTH = atoi(argv[4]);
    if(argc > 5) NP = atoi(argv[5]);
    if(NP > MAXNP) NP = MAXNP;
    int32_t k, a, b, c[5];

    for(k=0;k<NP;k++){
        setupP(allprimes[k]);
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

    /* build the coprime pair table */
    { int64_t cap = 4*(HNg+1)*HDg/2 + 16, cnt = 0;
      pn = malloc(sizeof(int32_t)*cap); pd = malloc(sizeof(int32_t)*cap);
      for(int64_t d=1; d<=HDg; d++) for(int64_t n=-HNg; n<=HNg; n++){
          if(gcd64(n,d)!=1) continue;
          if(cnt>=cap){ fprintf(stderr,"pair table overflow\n"); return 1; }
          pn[cnt]=(int32_t)n; pd[cnt]=(int32_t)d; cnt++;
      }
      NPAIR = cnt;
      fprintf(stderr,"pairs: %lld\n",(long long)NPAIR);
    }
    for(k=0;k<NP;k++){
        int32_t p = allprimes[k];
        setupP(p);
        int32_t it[128]; for(a=1;a<p;a++) it[a]=inv(a);
        res[k] = malloc(NPAIR);
        for(int64_t i=0;i<NPAIR;i++){
            int32_t dm = pd[i]%p;
            if(dm==0){ res[k][i]=255; continue; }
            int32_t nm = ((pn[i]%p)+p)%p;
            res[k][i] = (uint8_t)(((int64_t)nm*it[dm])%p);
        }
    }

    pthread_t th[256]; targ ta[256];
    if(NTH>256) NTH=256;
    for(k=0;k<NTH;k++){ ta[k].id=k; pthread_create(&th[k],0, which==9?worker9:worker7, &ta[k]); }
    for(k=0;k<NTH;k++) pthread_join(th[k],0);

    printf("SEARCH_DONE c%d HN=%lld HD=%lld nprimes=%d pairs=%lld checked=%lld survivors=%lld\n",
           which,(long long)HNg,(long long)HDg,NP,(long long)NPAIR,(long long)g_checked,(long long)g_surv);
    return 0;
}
