// claude_abc_sieve.c — rational points on Ari's (2,2,2,12) surface
//   y^2 = (B^2-A^2)(B^2-C^2)(B^2-A^2-C^2)
//   z^2 = (B^2-A^2-C^2)(B^2(A^2+C^2) - (A^4+A^2C^2+C^4))
// Enumerate 1 <= C < A <= H, B in (C,A) or (sqrt(A^2+C^2), H], both F,G > 0 there.
// All-even in each variable => positive orthant suffices; A<->C symmetry => A > C.
// Compile: gcc -O3 -march=native -funroll-loops -fopenmp -o sieve_abc claude_abc_sieve.c -lm
// Usage: ./sieve_abc H [nthreads]
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

typedef __int128 i128;
typedef unsigned __int128 u128;

static uint64_t qm64[64], qm63[63], qm65[65], qm11[11];

static void build_tables(void) {
    for (int i = 0; i < 64; i++) qm64[(uint64_t)(i*i) & 63] = 1;
    for (int i = 0; i < 63; i++) qm63[(uint64_t)(i*i) % 63] = 1;
    for (int i = 0; i < 65; i++) qm65[(uint64_t)(i*i) % 65] = 1;
    for (int i = 0; i < 11; i++) qm11[(uint64_t)(i*i) % 11] = 1;
}

static inline uint64_t u128_mod(u128 n, uint64_t m) { return (uint64_t)(n % m); }

// exact square test for nonnegative 128-bit n; returns root+1 (nonzero) if square, else 0
static inline uint64_t issq128(u128 n, uint64_t *root) {
    if (!qm64[(uint64_t)n & 63]) return 0;
    if (!qm63[u128_mod(n, 63)]) return 0;
    if (!qm65[u128_mod(n, 65)]) return 0;
    if (!qm11[u128_mod(n, 11)]) return 0;
    long double d = sqrtl((long double)n);
    uint64_t r = (uint64_t)d;
    for (int64_t off = -2; off <= 2; off++) {
        uint64_t rr = r + off;
        if ((u128)rr * rr == n) { *root = rr; return 1; }
    }
    return 0;
}

static inline int64_t gcd64(int64_t a, int64_t b) {
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a < 0 ? -a : a;
}

static int64_t H;
static long long hits = 0;

static void scan_B(int64_t A_, int64_t C_, int64_t B0, int64_t B1) {
    int64_t a2 = A_*A_, c2 = C_*C_, s2 = a2 + c2;
    int64_t q4 = a2*a2 + a2*c2 + c2*c2;          // A^4+A^2C^2+C^4  (<= ~4.9e17 at H=2.2e4)
    for (int64_t B_ = B0; B_ <= B1; B_++) {
        int64_t b2 = B_*B_;
        int64_t f1 = b2 - a2, f2 = b2 - c2, f3 = b2 - s2;
        int64_t f12 = f1 * f2;                    // |f1|,|f2| < 5e8 => fits
        // F = f12 * f3 > 0 in both ranges
        uint64_t m64 = ((uint64_t)f12 * (uint64_t)f3) & 63;   // two's complement: low bits correct
        if (!qm64[m64]) continue;
        u128 F = (u128)((i128)f12 * f3);
        uint64_t ry, rz;
        if (!issq128(F, &ry)) continue;
        int64_t g2 = b2 % 2 || 1 ? (int64_t)0 : 0; (void)g2;
        int64_t gg = (int64_t)((i128)b2 * s2 - q4) ; // b2*s2 <= ~4e17, q4 <= ~5e17: fits int64
        u128 G = (u128)((i128)f3 * gg);
        if ((i128)f3 * gg < 0) continue;          // shouldn't happen in-range, guard anyway
        if (!issq128(G, &rz)) continue;
        int64_t g = gcd64(gcd64(A_, B_), C_);
        int64_t Ap = A_/g, Bp = B_/g, Cp = C_/g;
        // discriminant factors (nonzero <=> nondegenerate); zeros of linear ones excluded by ranges
        int64_t Ap2 = Ap*Ap, Bp2 = Bp*Bp, Cp2 = Cp*Cp;
        int64_t d1 = Bp2 - Ap2 - 2*Cp2;
        int64_t d2 = Bp2 - Ap2 - Cp2;
        int64_t d3 = Bp2 - Ap2 - Ap*Cp - Cp2;
        int64_t d4 = Bp2 - Ap2 + Ap*Cp - Cp2;
        i128 d5 = (i128)Bp2*Bp2 - (i128)Ap2*Bp2 + 2*(i128)Ap2*Cp2 - 2*(i128)Bp2*Cp2 + (i128)Cp2*Cp2;
        i128 d6 = (i128)Ap2*Ap2 - (i128)Ap2*Bp2 + (i128)Ap2*Cp2 - (i128)Bp2*Cp2 + (i128)Cp2*Cp2;
        i128 d7 = (i128)Ap2*Ap2 - 2*(i128)Ap2*Bp2 + (i128)Bp2*Bp2 + 2*(i128)Ap2*Cp2 - (i128)Bp2*Cp2;
        int degen = (d1==0)||(d2==0)||(d3==0)||(d4==0)||(d5==0)||(d6==0)||(d7==0);
        #pragma omp critical
        {
            hits++;
            printf("HIT %s A=%lld B=%lld C=%lld  (primitive %lld %lld %lld)  y=%llu z=%llu\n",
                   degen ? "DEGEN" : "GOOD",
                   (long long)A_, (long long)B_, (long long)C_,
                   (long long)Ap, (long long)Bp, (long long)Cp,
                   (unsigned long long)ry, (unsigned long long)rz);
            fflush(stdout);
        }
    }
}

int main(int argc, char **argv) {
    H = argc > 1 ? atoll(argv[1]) : 10000;
    int nt = argc > 2 ? atoi(argv[2]) : omp_get_max_threads();
    omp_set_num_threads(nt);
    build_tables();
    fprintf(stderr, "sieve_abc H=%lld threads=%d\n", (long long)H, nt);
    double t0 = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 4)
    for (int64_t A_ = 2; A_ <= H; A_++) {
        for (int64_t C_ = 1; C_ < A_; C_++) {
            // range 1: C < B < A
            scan_B(A_, C_, C_ + 1, A_ - 1);
            // range 2: B > sqrt(A^2+C^2)
            int64_t s2 = A_*A_ + C_*C_;
            int64_t B0 = (int64_t)sqrtl((long double)s2);
            while (B0*B0 <= s2) B0++;
            if (B0 <= H) scan_B(A_, C_, B0, H);
        }
        if (A_ % 2000 == 0) {
            #pragma omp critical
            { fprintf(stderr, "PROGRESS A=%lld t=%.1fs hits=%lld\n", (long long)A_, omp_get_wtime()-t0, hits); }
        }
    }
    fprintf(stderr, "SEARCH_DONE H=%lld hits=%lld t=%.1fs\n", (long long)H, hits, omp_get_wtime()-t0);
    return 0;
}
