// lawsweep.c — LAW-TARGETED T5 sweep for (2,2,2,12) hits on M(2,2,2,6).
// PROVEN LAW (lane C): V'1,V'2 both squares => rho'(rho'-1) is a rational square
// => rn = -a^2, rd = c^2-a^2 (coprime a<c) up to the involution
// sigma: (u,rho') -> ((4u-3)/(4u-4), 1-rho') which swaps member rho' <-> 1-rho'
// and maps hits to hits.  So EVERY hit member is sigma-equivalent to one of the
// form rho' = -a^2/(c^2-a^2).  This sweep enumerates exactly those members.
// Known hit: (a,c) = (7,17) i.e. rho' = -49/240, at u = -97/48.
// Usage: lawsweep CMAX N NTHREADS   (members a<c<=CMAX coprime; u=p/q, |p|,q<=N)
//        lawsweep --test p q rn rd
// Compile: gcc -O3 -march=native -fopenmp -o lawsweep lawsweep.c -lm
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>

typedef __int128 i128;
static uint8_t sq256[256];
static uint8_t sq45045[45045];
static uint8_t sq97[97], sq101[101];
static void init_tabs(void) {
    for (int i = 0; i < 256; i++) sq256[(i * i) & 255] = 1;
    for (int64_t i = 0; i < 45045; i++) sq45045[(i * i) % 45045] = 1;
    for (int i = 0; i < 97; i++) sq97[(i * i) % 97] = 1;
    for (int i = 0; i < 101; i++) sq101[(i * i) % 101] = 1;
}
static int is_sq128(i128 x) {
    if (x <= 0) return 0;
    if (!sq256[(uint64_t)x & 255]) return 0;
    if (!sq45045[(uint64_t)(x % 45045)]) return 0;
    if (!sq97[(uint64_t)(x % 97)]) return 0;
    if (!sq101[(uint64_t)(x % 101)]) return 0;
    long double xd = (long double)x;
    int64_t r = (int64_t)sqrtl(xd);
    for (int64_t t = r - 2; t <= r + 2; t++)
        if (t >= 0 && (i128)t * t == x) return 1;
    return 0;
}
static int64_t gcd64(int64_t x, int64_t y) {
    if (x < 0) x = -x; if (y < 0) y = -y;
    while (y) { int64_t t = x % y; x = y; y = t; }
    return x;
}
static int test_one(int64_t p, int64_t q, int64_t rn, int64_t rd, int verbose) {
    int64_t QQ = 4*p*p - 6*p*q + 3*q*q;
    int64_t C3 = 16*p*p*p*p - 40*p*p*p*q + 40*p*p*q*q - 18*p*q*q*q + 3*q*q*q*q;
    int64_t C2 = -16*p*p*p*p + 32*p*p*p*q - 28*p*p*q*q + 10*p*q*q*q - q*q*q*q;
    int64_t C1 = (8*p*p*p - 12*p*p*q + 10*p*q*q - 3*q*q*q)*q;
    int64_t C0 = (-2*p + q)*q*q*q;
    i128 X1 = (i128)rn * ((i128)QQ*rn - (i128)(2*p-q)*q*rd);
    i128 X2 = (i128)QQ*rn*rn - (i128)(4*p*p-4*p*q+2*q*q)*rn*rd + (i128)(2*p-q)*q*rd*rd;
    i128 inner = ((i128)C3*rn*rn*rn) + ((i128)C2*rn*rn)*rd + ((i128)C1*rn)*rd*rd
                 + (i128)C0*rd*rd*rd;
    i128 X3 = (i128)rn * inner;
    i128 X4 = (i128)C3*rn*rn + (i128)(-16*p*p*p+28*p*p*q-18*p*q*q+4*q*q*q)*q*rn*rd
              + (i128)(2*p-q)*(2*p-q)*q*q*rd*rd;
    int m = (is_sq128(X1)?1:0) | (is_sq128(X2)?2:0) | (is_sq128(X3)?4:0) | (is_sq128(X4)?8:0);
    if (verbose)
        printf("p=%lld q=%lld rn=%lld rd=%lld: X1sq=%d X2sq=%d X3sq=%d X4sq=%d\n",
            (long long)p,(long long)q,(long long)rn,(long long)rd, m&1,(m>>1)&1,(m>>2)&1,(m>>3)&1);
    return m;
}
int main(int argc, char *argv[]) {
    init_tabs();
    if (argc >= 6 && !strcmp(argv[1], "--test")) {
        test_one(atoll(argv[2]), atoll(argv[3]), atoll(argv[4]), atoll(argv[5]), 1);
        return 0;
    }
    if (argc < 3) { puts("usage: lawsweep CMAX N [nthreads]"); return 1; }
    int64_t CMAX = atoll(argv[1]), N = atoll(argv[2]);
    int nth = argc > 3 ? atoi(argv[3]) : 3;
    omp_set_num_threads(nth);
    // build member list rn = -a^2, rd = c^2-a^2, coprime a<c<=CMAX
    int64_t M = 0, cap = CMAX*CMAX/2 + 16;
    int64_t *mrn = malloc(cap*sizeof(int64_t)), *mrd = malloc(cap*sizeof(int64_t));
    for (int64_t a = 1; a < CMAX; a++)
        for (int64_t c = a+1; c <= CMAX; c++) {
            if (gcd64(a,c) != 1) continue;
            mrn[M] = -a*a; mrd[M] = c*c - a*a; M++;
        }
    fprintf(stderr, "lawsweep: %lld members (CMAX=%lld), N=%lld, %d threads\n",
            (long long)M, (long long)CMAX, (long long)N, nth);
    int64_t hits = 0, pass1 = 0;
    double t0 = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 2) reduction(+:hits,pass1)
    for (int64_t q = 1; q <= N; q++) {
        for (int64_t p = -N; p <= N; p++) {
            if (p == 0 || p == q || 2*p == q || gcd64(p, q) != 1) continue;
            int64_t QQ = 4*p*p - 6*p*q + 3*q*q;
            int64_t A1 = (2*p-q)*q;
            uint8_t QQ8 = (uint8_t)(QQ & 255), A18 = (uint8_t)(((A1 % 256)+256)&255);
            for (int64_t mi = 0; mi < M; mi++) {
                int64_t rn = mrn[mi], rd = mrd[mi];
                uint8_t rn8 = (uint8_t)(((rn % 256)+256)&255);
                uint8_t rd8 = (uint8_t)(rd & 255);
                uint8_t x18 = (uint8_t)((rn8 * (uint8_t)((QQ8*rn8 - A18*rd8)&255))&255);
                if (!sq256[x18]) continue;
                i128 X1 = (i128)rn * ((i128)QQ*rn - (i128)A1*rd);
                if (X1 <= 0 || !is_sq128(X1)) continue;
                pass1++;
                int m = test_one(p, q, rn, rd, 0);
                if ((m & 7) == 7) {
                    #pragma omp critical(out)
                    { printf("%s p=%lld q=%lld rn=%lld rd=%lld (X4 %s)\n",
                             (m==15)?"HIT":"NEAR", (long long)p,(long long)q,
                             (long long)rn,(long long)rd, (m&8)?"sq":"nonsq");
                      fflush(stdout); }
                    if (m == 15) hits++;
                }
            }
        }
        if (q % 100 == 0)
            fprintf(stderr, "  q=%lld done (%.1fs)\n", (long long)q, omp_get_wtime()-t0);
    }
    fprintf(stderr, "LAWSWEEP DONE CMAX=%lld N=%lld members=%lld X1-passes=%lld hits=%lld time=%.1fs\n",
            (long long)CMAX, (long long)N, (long long)M, (long long)pass1,
            (long long)hits, omp_get_wtime()-t0);
    return 0;
}
