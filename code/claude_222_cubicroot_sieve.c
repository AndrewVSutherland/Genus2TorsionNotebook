// claude_222_cubicroot_sieve.c — bitmap sieve for a rational root of the
// residual cubic over the parametrized [22] families (the one remaining
// extra-2 condition => [2,22] candidate).  Families:
//   0 = famA- : branch-point family, cubic coeffs A0..A3 in s
//   1 = famA+ : same with s -> -s
//   2 = famB  : Flynn43-parametrized family, cubic coeffs B0..B3 in tau
// A parameter n/d survives iff for every sieve prime p (p ∤ d) the cubic at
// r = n*d^{-1} mod p has a root mod p.  Survivors go to stdout for exact
// checking in PARI.  Usage: ./sieve FAM H NP   (NP = 30 validation / 60 prod)
// Build: gcc -O3 -march=native -fopenmp -o sieve claude_222_cubicroot_sieve.c
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

static int64_t A0[10] = {0,0,-1,2,-3,2,-1,0,0,0};
static int64_t A1[10] = {1,-2,2,-2,3,-2,1,0,0,0};
static int64_t A2[10] = {1,-4,5,-4,2,0,0,0,0,0};
static int64_t A3[10] = {1,-2,1,0,0,0,0,0,0,0};
static int64_t B0[10] = {-121265525824,54070679104,-10154594160,1022522080,-57942460,1742244,-21609,0,0,0};
static int64_t B1[10] = {127538165824,-59392919104,12012194160,-1363322080,92542460,-3582244,61609,0,0,0};
static int64_t B2[10] = {61432588800,-41020972800,11217152000,-1615856000,130052000,-5562800,98800,0,0,0};
static int64_t B3[10] = {6272640000,-5322240000,1857600000,-340800000,34600000,-1840000,40000,0,0,0};
static int32_t PRIMES60[60] = {211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,
 307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,
 449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577};

static inline int64_t md(int64_t a, int64_t p){ a %= p; if (a < 0) a += p; return a; }

static int64_t inv_mod(int64_t a, int64_t p){
  int64_t t = 0, nt = 1, r = p, nr = md(a,p);
  while (nr != 0){ int64_t q = r/nr, tmp;
    tmp = t - q*nt; t = nt; nt = tmp;
    tmp = r - q*nr; r = nr; nr = tmp; }
  return md(t, p);
}

static int64_t horner(int64_t *c, int64_t r, int64_t p){
  int64_t v = 0;
  for (int j = 9; j >= 0; j--) v = md(v*r + md(c[j],p), p);
  return v;
}

int main(int argc, char **argv){
  if (argc < 4){ fprintf(stderr, "usage: sieve FAM H NP\n"); return 1; }
  int fam = atoi(argv[1]);
  int64_t H = atoll(argv[2]);
  int np = atoi(argv[3]); if (np > 60) np = 60;
  int64_t *c[4];
  static int64_t Am0[10], Am1[10], Am2[10], Am3[10];
  if (fam == 2){ c[0]=B0; c[1]=B1; c[2]=B2; c[3]=B3; }
  else if (fam == 0){ c[0]=A0; c[1]=A1; c[2]=A2; c[3]=A3; }
  else { // famA+ : s -> -s  (negate odd-degree coefficients)
    for (int j=0;j<10;j++){ int sg = (j&1)? -1 : 1;
      Am0[j]=sg*A0[j]; Am1[j]=sg*A1[j]; Am2[j]=sg*A2[j]; Am3[j]=sg*A3[j]; }
    c[0]=Am0; c[1]=Am1; c[2]=Am2; c[3]=Am3;
  }
  // build tables: tab[i][r] = 1 iff cubic at parameter r has a root mod p (or degenerates)
  uint8_t *tab[60];
  for (int i = 0; i < np; i++){
    int64_t p = PRIMES60[i];
    tab[i] = calloc(p, 1);
    for (int64_t r = 0; r < p; r++){
      int64_t c3 = horner(c[3], r, p), c2 = horner(c[2], r, p),
              c1 = horner(c[1], r, p), c0 = horner(c[0], r, p);
      if (c3 == 0){ tab[i][r] = 1; continue; }          // degenerate: pass
      int ok = 0;
      for (int64_t x = 0; x < p; x++){
        int64_t v = md(((c3*x + c2) % p) * x + c1, p);
        v = md(v*x + c0, p);
        if (v == 0){ ok = 1; break; }
      }
      tab[i][r] = (uint8_t) ok;
    }
  }
  int64_t nsurv_total = 0;
  #pragma omp parallel reduction(+:nsurv_total)
  {
    int32_t *buf = malloc(sizeof(int32_t)*(2*H+2));
    int32_t *buf2 = malloc(sizeof(int32_t)*(2*H+2));
    #pragma omp for schedule(dynamic, 16)
    for (int64_t d = 1; d <= H; d++){
      // round 0: prime 0 incremental
      int64_t nb = 0;
      { int64_t p = PRIMES60[0];
        if (d % p == 0){ for (int64_t n = -H; n <= H; n++) buf[nb++] = (int32_t)n; }
        else { int64_t di = inv_mod(d, p);
          int64_t r = md(-H*di, p);
          for (int64_t n = -H; n <= H; n++){
            if (tab[0][r]) buf[nb++] = (int32_t)n;
            r += di; if (r >= p) r -= p; } } }
      // later rounds
      for (int i = 1; i < np && nb > 0; i++){
        int64_t p = PRIMES60[i];
        if (d % p == 0) continue;
        int64_t di = inv_mod(d, p);
        int64_t nb2 = 0;
        for (int64_t k = 0; k < nb; k++){
          int64_t n = buf[k];
          int64_t r = md(n % p, p); r = (r*di) % p;
          if (tab[i][r]) buf2[nb2++] = (int32_t)n;
        }
        int32_t *tmp = buf; buf = buf2; buf2 = tmp; nb = nb2;
      }
      if (nb > 0){
        #pragma omp critical
        { for (int64_t k = 0; k < nb; k++) printf("SURV %d %lld\n", buf[k], (long long)d); }
        nsurv_total += nb;
      }
    }
    free(buf); free(buf2);
  }
  fprintf(stderr, "fam=%d H=%lld np=%d survivors(unreduced)=%lld\n", fam, (long long)H, np, (long long)nsurv_total);
  return 0;
}
