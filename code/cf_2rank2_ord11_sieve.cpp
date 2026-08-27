// [2,22] hunt: C sieve for D_inf order on the 2-rank-2 family
//     f = (x^2-1)(x^2+a x+b)(x^2+c x+d),   a in [0..H], b,c,d in [-H..H],
// symmetry quotient (a,b) <= (c,d) lexicographic (swap of the two free
// quadratics; x -> -x handled by a >= 0).
//
// Stage 1 (this program): polynomial-Pell CF order of sqrt(f) mod the three
// primes 101, 103, 107 (exact order of D_inf in J(F_p)).  A global class of
// order N reduces to exact order N at every good prime, so no true hit is
// lost.  Survivors (all three primes agree on an INTERESTING order:
// odd 7..13, or a multiple of 11 up to 44) go to stdout/TSV for the exact
// Magma stage (code/cf_2rank2_ord11_exact.m).
// Order 5 is deliberately NOT collected (dense already at H=12; see
// results/cf_2rank2_ord11_H12.log).
//
// Self-tests at startup (abort on failure):
//   f14 = (x^2+1)(x^4+5x^2+4x+4)      -> 14 mod 101 and 103
//   f18 = (x^2-x+1)(x^4-x^3+9x^2+8x-8)-> 18 mod 101
//   f28-trap                           -> 7 mod 101
//   family point (a,b,c,d)=(8,9,10,-9) -> 5 at all three primes
//
// Build:  g++ -O3 -march=native -std=c++17 -fopenmp \
//             code/cf_2rank2_ord11_sieve.cpp -o /tmp/cf22sieve
// Run:    /tmp/cf22sieve H THREADS > data/cf_2rank2_ord11_sieve_H<H>.tsv

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

static const int NP = 3;
static const int PRIMES[NP] = {101, 103, 107};
static const int MAXORD = 44;
static const int MAXSTEPS = 60;

struct Poly {
    int deg;              // -1 for zero
    long long c[16];      // coefficients mod p, c[i] ~ x^i
};

static inline void pnorm(Poly &A, int p) {
    while (A.deg >= 0 && A.c[A.deg] % p == 0) A.deg--;
    for (int i = 0; i <= A.deg; i++) { A.c[i] %= p; if (A.c[i] < 0) A.c[i] += p; }
}

static inline long long ppow(long long b, long long e, long long m) {
    long long r = 1; b %= m;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    return r;
}
static inline long long pinv(long long a, long long p) { return ppow(a, p - 2, p); }

// A -= B * (k * x^sh)
static inline void submulshift(Poly &A, const Poly &B, long long k, int sh, int p) {
    for (int i = 0; i <= B.deg; i++) {
        int j = i + sh;
        long long v = (A.deg >= j ? A.c[j] : 0);
        v = (v - B.c[i] * k) % p; if (v < 0) v += p;
        if (j > A.deg) { for (int t = A.deg + 1; t <= j; t++) A.c[t] = 0; A.deg = j; }
        A.c[j] = v;
    }
    pnorm(A, p);
}

// quotient and remainder of A by B (B != 0)
static inline void pdivmod(const Poly &A, const Poly &B, Poly &Qo, Poly &R, int p) {
    R = A; Qo.deg = -1;
    for (int i = 0; i < 16; i++) Qo.c[i] = 0;
    if (B.deg < 0) { return; }
    long long binv = pinv(B.c[B.deg], p);
    while (R.deg >= B.deg && R.deg >= 0) {
        int sh = R.deg - B.deg;
        long long k = R.c[R.deg] * binv % p;
        if (sh > Qo.deg) { for (int t = Qo.deg + 1; t <= sh; t++) Qo.c[t] = 0; Qo.deg = sh; }
        Qo.c[sh] = (Qo.c[sh] + k) % p;
        submulshift(R, B, k, sh, p);
    }
    pnorm(Qo, p);
}

static inline void pmul(const Poly &A, const Poly &B, Poly &C, int p) {
    C.deg = (A.deg < 0 || B.deg < 0) ? -1 : A.deg + B.deg;
    for (int i = 0; i < 16; i++) C.c[i] = 0;
    if (C.deg < 0) return;
    for (int i = 0; i <= A.deg; i++)
        for (int j = 0; j <= B.deg; j++)
            C.c[i + j] = (C.c[i + j] + A.c[i] * B.c[j]) % p;
    pnorm(C, p);
}

static inline void padd(const Poly &A, const Poly &B, Poly &C, int p) {
    C.deg = A.deg > B.deg ? A.deg : B.deg;
    for (int i = 0; i <= C.deg; i++) {
        long long v = (i <= A.deg ? A.c[i] : 0) + (i <= B.deg ? B.c[i] : 0);
        C.c[i] = v % p;
    }
    pnorm(C, p);
}

static inline void psub(const Poly &A, const Poly &B, Poly &C, int p) {
    C.deg = A.deg > B.deg ? A.deg : B.deg;
    for (int i = 0; i <= C.deg; i++) {
        long long v = (i <= A.deg ? A.c[i] : 0) - (i <= B.deg ? B.c[i] : 0);
        v %= p; if (v < 0) v += p;
        C.c[i] = v;
    }
    pnorm(C, p);
}

// squarefree test: gcd(f, f') constant
static bool squarefree(const Poly &F, int p) {
    Poly A = F, B, T, Qo, R;
    B.deg = F.deg - 1;
    for (int i = 0; i < 16; i++) B.c[i] = 0;
    for (int i = 1; i <= F.deg; i++) B.c[i - 1] = (long long)i % p * F.c[i] % p;
    pnorm(B, p);
    if (B.deg < 0) return false;
    while (B.deg > 0) { pdivmod(A, B, Qo, R, p); A = B; B = R; if (B.deg < 0) return false; }
    return B.deg == 0;   // nonzero constant gcd tail => coprime
}

// CF order of sqrt(F), monic sextic F mod p; 0 = none within budget
static int cforder(const Poly &F, int p) {
    // polynomial part s of sqrt(F): s = x^3+s2x^2+s1x+s0
    long long i2 = pinv(2, p);
    long long s2 = F.c[5] * i2 % p;
    long long s1 = (F.c[4] - s2 * s2 % p + (long long)p * p) % p * i2 % p;
    long long s0 = (F.c[3] - 2 * s2 * s1 % p + (long long)p * p) % p * i2 % p;
    Poly S; S.deg = 3;
    for (int i = 0; i < 16; i++) S.c[i] = 0;
    S.c[3] = 1; S.c[2] = s2; S.c[1] = s1; S.c[0] = s0;
    Poly Pi, Qi;
    Pi.deg = -1; for (int i = 0; i < 16; i++) Pi.c[i] = 0;
    Qi.deg = 0;  for (int i = 0; i < 16; i++) Qi.c[i] = 0;
    Qi.c[0] = 1;
    int total = 0;
    Poly T, Ai, R, Pn, Qn, T2;
    for (int it = 0; it <= MAXSTEPS; it++) {
        if (Qi.deg < 0) return 0;
        padd(Pi, S, T, p);
        pdivmod(T, Qi, Ai, R, p);
        if (Ai.deg < 0) return 0;
        total += Ai.deg;
        if (total > MAXORD) return 0;
        pmul(Ai, Qi, T2, p);
        psub(T2, Pi, Pn, p);
        pmul(Pn, Pn, T2, p);
        psub(F, T2, T, p);          // F - Pn^2
        pdivmod(T, Qi, Qn, R, p);
        if (R.deg >= 0) return 0;   // non-exact division: abort
        Pi = Pn; Qi = Qn;
        if (it >= 1 && Qi.deg == 0 && Qi.c[0] != 0) return total;
    }
    return 0;
}

static void family_poly(long long a, long long b, long long c, long long d,
                        int p, Poly &F) {
    // (x^2-1)(x^2+ax+b)(x^2+cx+d)
    Poly A, B, C, T;
    A.deg = 2; for (int i = 0; i < 16; i++) A.c[i] = 0;
    A.c[2] = 1; A.c[0] = p - 1;
    B.deg = 2; for (int i = 0; i < 16; i++) B.c[i] = 0;
    B.c[2] = 1; B.c[1] = ((a % p) + p) % p; B.c[0] = ((b % p) + p) % p;
    C.deg = 2; for (int i = 0; i < 16; i++) C.c[i] = 0;
    C.c[2] = 1; C.c[1] = ((c % p) + p) % p; C.c[0] = ((d % p) + p) % p;
    pmul(A, B, T, p);
    pmul(T, C, F, p);
}

static inline bool interesting(int o) {
    if (o >= 7 && o <= 13 && (o & 1)) return true;
    return (o > 0 && o % 11 == 0);
}

int main(int argc, char **argv) {
    int H = argc > 1 ? atoi(argv[1]) : 40;
    int threads = argc > 2 ? atoi(argv[2]) : 8;
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    // ---- self-tests ----
    {
        // f14 = (x^2+1)(x^4+5x^2+4x+4) = x^6+6x^4+4x^3+9x^2+4x+4
        // f18 = (x^2-x+1)(x^4-x^3+9x^2+8x-8)
        //     = x^6-2x^5+11x^4-2x^3-7x^2+16x-8
        // f28t = x^6+2x^5-5x^4-14x^3-3x^2+24x+28
        long long f14c[7] = {4, 4, 9, 4, 6, 0, 1};
        long long f18c[7] = {-8, 16, -7, -2, 11, -2, 1};
        long long f28c[7] = {28, 24, -3, -14, -5, 2, 1};
        for (int tp = 0; tp < 2; tp++) {
            int p = PRIMES[tp];
            Poly F; F.deg = 6;
            for (int i = 0; i < 16; i++) F.c[i] = 0;
            for (int i = 0; i <= 6; i++) F.c[i] = ((f14c[i] % p) + p) % p;
            if (cforder(F, p) != 14) { fprintf(stderr, "SELFTEST FAIL f14 p=%d\n", p); return 1; }
        }
        Poly F; F.deg = 6;
        for (int i = 0; i < 16; i++) F.c[i] = 0;
        for (int i = 0; i <= 6; i++) F.c[i] = ((f18c[i] % 101) + 101) % 101;
        if (cforder(F, 101) != 18) { fprintf(stderr, "SELFTEST FAIL f18\n"); return 1; }
        for (int i = 0; i < 16; i++) F.c[i] = 0;
        F.deg = 6;
        for (int i = 0; i <= 6; i++) F.c[i] = ((f28c[i] % 101) + 101) % 101;
        if (cforder(F, 101) != 7) { fprintf(stderr, "SELFTEST FAIL f28trap\n"); return 1; }
        for (int tp = 0; tp < NP; tp++) {
            int p = PRIMES[tp];
            Poly G; family_poly(8, 9, 10, -9, p, G);
            if (!squarefree(G, p) || cforder(G, p) != 5) {
                fprintf(stderr, "SELFTEST FAIL family(8,9,10,-9) p=%d\n", p); return 1;
            }
        }
        fprintf(stderr, "SELFTEST_PASS\n");
    }
    long long tested = 0, out = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:tested, out)
    for (int a = 0; a <= H; a++) {
        char buf[1 << 16]; int blen = 0;
        for (int b = -H; b <= H; b++) {
            for (int c = -H; c <= H; c++) {
                for (int d = -H; d <= H; d++) {
                    if (a > c || (a == c && b > d)) continue;   // swap quotient
                    tested++;
                    int o0 = -2, ngood = 0;
                    bool ok = true;
                    for (int tp = 0; tp < NP && ok; tp++) {
                        int p = PRIMES[tp];
                        Poly F; family_poly(a, b, c, d, p, F);
                        if (!squarefree(F, p)) continue;   // bad reduction: try next prime
                        int o = cforder(F, p);
                        if (ngood == 0) {
                            if (!interesting(o)) { ok = false; break; }
                            o0 = o; ngood = 1;
                        } else if (o != o0) { ok = false; break; }
                        else ngood++;
                    }
                    if (!ok || ngood < 2) continue;
                    out++;
                    blen += snprintf(buf + blen, sizeof(buf) - blen,
                                     "%d %d %d %d %d\n", a, b, c, d, o0);
                    if (blen > (1 << 15)) {
#pragma omp critical
                        { fputs(buf, stdout); }
                        blen = 0;
                    }
                }
            }
        }
        if (blen) {
#pragma omp critical
            { fputs(buf, stdout); }
        }
    }
    fprintf(stderr, "SIEVE_DONE H=%d tested=%lld survivors=%lld\n", H, tested, out);
    return 0;
}
