// claude_sym_sieve.c — orbit-level search for (2,2,2,12) curves via the multigrade model:
//   X1^2+X2^2+X3^2 = B1^2+B2^2,  X1^4+X2^4+X3^4 = B1^4+B2^4   (X1<X2<X3, B1>B2>0)
// Facts: dd = e1^2-4e2 = (x3-x1-x2)^2 - 4x1x2 = w^2 > 0 requires X3 > X1+X2 (largest dominates).
// Enumerate (X1<X2<=H2); factor 4x1x2 = r*r' (r<=r', same parity); x3 = x1+x2+(r+r')/2 must be
// a square X3^2 (X3 <= X3MAX); then u=(e1+w)/2, v=(e1-w)/2 must be squares B1^2, B2^2.
// Compile: gcc -O3 -march=native -fopenmp -o sym_sieve claude_sym_sieve.c -lm
// Usage: ./sym_sieve H2 X3MAX [threads]
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

typedef unsigned __int128 u128;

static int32_t *spf;   // smallest prime factor table up to H2

static inline uint64_t isqrt64(uint64_t n) {
    uint64_t r = (uint64_t)sqrtl((long double)n);
    while (r*r > n) r--;
    while ((r+1)*(r+1) <= n) r++;
    return r;
}
static inline int is_square64(uint64_t n, uint64_t *root) {
    static const uint64_t good = 0x2030213022121200ULL; // not used; fall through
    (void)good;
    uint64_t m = n & 63;
    if (!((m==0)|(m==1)|(m==4)|(m==9)|(m==16)|(m==17)|(m==25)|(m==33)|(m==36)|(m==41)|(m==49)|(m==57))) return 0;
    uint64_t r = isqrt64(n);
    if (r*r == n) { *root = r; return 1; }
    return 0;
}

// merge factorization of X (via spf table) into (pr, ex) with multiplicity mult
static int add_factors(int64_t X, int mult, uint64_t *pr, int *ex, int np) {
    while (X > 1) {
        int64_t p = spf[X];
        int e = 0;
        while (X % p == 0) { X /= p; e++; }
        int found = -1;
        for (int i = 0; i < np; i++) if (pr[i] == (uint64_t)p) { found = i; break; }
        if (found >= 0) ex[found] += e*mult;
        else { pr[np] = p; ex[np] = e*mult; np++; }
    }
    return np;
}
// factorization of N = 4*X1^2*X2^2 from spf factorizations of X1, X2
static int factor_N(int64_t X1, int64_t X2, uint64_t *pr, int *ex) {
    pr[0] = 2; ex[0] = 2;   // the factor 4
    int np = 1;
    np = add_factors(X1, 2, pr, ex, np);
    np = add_factors(X2, 2, pr, ex, np);
    if (ex[0] == 0) { /* cannot happen */ }
    return np;
}

int64_t H2 = 30000, X3MAX = 1000000;
long long nhits = 0;

static __thread uint64_t divs[1<<17];

static void handle_pair(int64_t X1, int64_t X2) {
    uint64_t x1 = (uint64_t)X1*X1, x2 = (uint64_t)X2*X2;
    u128 N = (u128)4*x1*x2;            // = (2*X1*X2)^2: exceeds u64 once X1*X2 > 2^31, so keep in u128
    uint64_t rcap = 2*(uint64_t)X1*(uint64_t)X2;   // = sqrt(N) exactly (<= 5e9 at H2=5e4)
    uint64_t pr[24]; int ex[24];
    int np = factor_N(X1, X2, pr, ex);
    // generate only divisors r <= rcap = sqrt(N): each pair (r, N/r) is visited from its smaller
    // member, and a product > rcap only has multiples > rcap, so capping mid-generation is lossless.
    int nd = 0;
    divs[nd++] = 1;
    int overflow = 0;
    for (int i = 0; i < np; i++) {
        int base = nd;
        u128 pe = 1;
        for (int e = 1; e <= ex[i]; e++) {
            pe *= pr[i];
            if (pe > rcap) break;
            for (int j = 0; j < base; j++) {
                u128 dv = divs[j]*pe;
                if (dv > rcap) continue;
                if (nd < (1<<17)) divs[nd++] = (uint64_t)dv; else overflow = 1;
            }
        }
    }
    if (overflow) {
        #pragma omp critical
        { fprintf(stderr, "WARN divisor overflow at X1=%lld X2=%lld\n", (long long)X1, (long long)X2); }
    }
    uint64_t x12 = x1 + x2;
    uint64_t x3cap = (uint64_t)X3MAX*(uint64_t)X3MAX;
    for (int i = 0; i < nd; i++) {
        uint64_t r = divs[i];
        u128 rq128 = N / r;                    // cofactor can exceed u64 (up to N when r=1)
        if (rq128 > (u128)2*x3cap) continue;   // D1 >= rq/2 would already blow past x3cap
        uint64_t rq = (uint64_t)rq128;         // now rq <= 2*x3cap = 8e12: safe in u64
        if (((r ^ rq) & 1)) continue;          // same parity
        uint64_t D1 = (r + rq) / 2;            // = x3 - x1 - x2 > 0
        uint64_t w  = (rq - r) / 2;            // dd = w^2
        uint64_t x3 = x12 + D1;
        if (x3 > x3cap) continue;
        uint64_t X3r;
        if (!is_square64(x3, &X3r)) continue;
        if (X3r <= (uint64_t)X2) continue;      // ordering X3 largest (auto, but guard)
        uint64_t e1 = x12 + x3;
        if ((e1 ^ w) & 1) continue;             // u,v integral
        uint64_t u = (e1 + w)/2, v = (e1 - w)/2;
        if (w == 0) continue;   // B1=B2 degenerate family (the t=+-1 conics)
        uint64_t B1r, B2r;
        if (!is_square64(u, &B1r)) continue;
        if (!is_square64(v, &B2r)) continue;
        if (v == 0) continue;
        // gcd(X1,X2,X3) primitive check
        int64_t a = X1, b = X2 % X1; while (b) { int64_t t = a % b; a = b; b = t; }
        int64_t g = a; b = (int64_t)(X3r % (uint64_t)g); while (b) { int64_t t = g % b; g = b; b = t; }
        #pragma omp critical
        {
            nhits++;
            printf("ORBIT X={%lld,%lld,%llu} B={%llu,%llu} gcd=%lld %s\n",
                (long long)X1, (long long)X2, (unsigned long long)X3r,
                (unsigned long long)B1r, (unsigned long long)B2r, (long long)g,
                w == 0 ? "DEGEN(B1=B2)" : "");
            fflush(stdout);
        }
    }
}

int main(int argc, char **argv) {
    if (argc > 1) H2 = atoll(argv[1]);
    if (argc > 2) X3MAX = atoll(argv[2]);
    int nt = argc > 3 ? atoi(argv[3]) : omp_get_max_threads();
    omp_set_num_threads(nt);
    spf = calloc(H2+1, sizeof(int32_t));
    for (int64_t i = 2; i <= H2; i++) {
        if (spf[i] == 0) for (int64_t j = i; j <= H2; j += i) if (spf[j] == 0) spf[j] = (int32_t)i;
    }
    fprintf(stderr, "sym_sieve H2=%lld X3MAX=%lld threads=%d\n",
        (long long)H2, (long long)X3MAX, nt);
    double t0 = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 4)
    for (int64_t X2 = 2; X2 <= H2; X2++) {
        for (int64_t X1 = 1; X1 < X2; X1++) handle_pair(X1, X2);
        if (X2 % 2000 == 0) {
            #pragma omp critical
            { fprintf(stderr, "PROGRESS X2=%lld t=%.1fs hits=%lld\n",
                (long long)X2, omp_get_wtime()-t0, nhits); }
        }
    }
    fprintf(stderr, "SEARCH_DONE H2=%lld X3MAX=%lld hits=%lld t=%.1fs\n",
        (long long)H2, (long long)X3MAX, nhits, omp_get_wtime()-t0);
    return 0;
}
