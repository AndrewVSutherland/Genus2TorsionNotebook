// [2,22] hunt, route 4: 11-divisibility sieve on the 2-rank-2 family
//     f = (x^2-1)(x^2+a x+b)(x^2+c x+d),  a in [0..H], b,c,d in [-H..H],
// swap-symmetry quotient (a,b) <= (c,d).
//
// A rational 11-torsion class (ANY class, not just D_inf) forces
// 11 | #J(F_p) at every good prime.  Chain over p in {3,5,7,13,17,19,23,29}:
// skip bad-reduction primes, reject on the first good prime with
// 11 ndiv #J(F_p), require >= 5 good primes passed.  Survivors go to the
// exact Magma stage (code/div11_2rank2_exact.m).
// #J(F_q) = (#C(F_q)^2 + #C(F_{q^2}))/2 - q,   #C(F_q) = q + 2 + sum chi(f(x))
// (monic sextic => two rational points at infinity).  chi over F_{p^2} via
// chi_p(Norm).
//
// Self-tests: Flynn's order-11 family members t = 1, 2, 5 (11 | #J at every
// good p in the chain) and the CF vector f14 (14 | #J at p=3 must hold...
// checked as: #J(F_3) computed equals Magma's value 14 for f14).
//
// Build: g++ -O3 -march=native -std=c++17 -fopenmp \
//            code/div11_2rank2_sieve.cpp -o /tmp/div11sieve
// Run:   /tmp/div11sieve H THREADS > data/div11_2rank2_sieve_H<H>.tsv

#include <cstdio>
#include <cstdlib>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

static const int NPR = 8;
static const int PRS[NPR] = {3, 5, 7, 13, 17, 19, 23, 29};
static const int MINGOOD = 5;

struct PrimeCtx {
    int p;
    int ns;              // a non-square mod p
    signed char chi[32]; // chi[x] for x in [0,p)
};

static PrimeCtx CTX[NPR];

static void init_ctx() {
    for (int i = 0; i < NPR; i++) {
        int p = PRS[i];
        CTX[i].p = p;
        for (int x = 0; x < p; x++) CTX[i].chi[x] = -1;
        CTX[i].chi[0] = 0;
        for (int x = 1; x < p; x++) CTX[i].chi[(x * x) % p] = 1;
        for (int x = 1; x < p; x++) if (CTX[i].chi[x] == -1) { CTX[i].ns = x; break; }
    }
}

// f coefficients mod p (monic sextic), c[0..6]
// squarefree via gcd(f, f') over F_p
static bool sqfree_modp(const int *c, int p) {
    int A[8], B[8], degA = 6, degB;
    for (int i = 0; i <= 6; i++) A[i] = c[i];
    for (int i = 1; i <= 6; i++) B[i - 1] = (int)((long long)i * c[i] % p);
    degB = 5; while (degB >= 0 && B[degB] == 0) degB--;
    if (degB < 0) return false;
    while (degB > 0) {
        // A mod B
        while (degA >= degB) {
            // inverse of B[degB]
            long long binv = 1, base = B[degB], e = p - 2;
            while (e) { if (e & 1) binv = binv * base % p; base = base * base % p; e >>= 1; }
            long long k = (long long)A[degA] * binv % p;
            int sh = degA - degB;
            for (int i = 0; i <= degB; i++) {
                A[i + sh] = (int)(((long long)A[i + sh] - k * B[i]) % p);
                if (A[i + sh] < 0) A[i + sh] += p;
            }
            while (degA >= 0 && A[degA] == 0) degA--;
            if (degA < 0) break;
        }
        // swap
        int T[8], degT = degA;
        for (int i = 0; i <= (degA < 0 ? -1 : degA); i++) T[i] = A[i];
        for (int i = 0; i <= degB; i++) A[i] = B[i];
        degA = degB;
        if (degT < 0) return false;
        for (int i = 0; i <= degT; i++) B[i] = T[i];
        degB = degT;
    }
    return B[0] != 0;   // constant nonzero gcd
}

// #J(F_p) mod 11 test; assumes f squarefree mod p.  Returns J mod 11... we
// only need divisibility, but compute #J exactly (fits in long long).
static long long jac_count(const int *c, const PrimeCtx &ctx) {
    int p = ctx.p;
    // #C(F_p)
    long long s1 = 0;
    for (int x = 0; x < p; x++) {
        long long v = 0;
        for (int i = 6; i >= 0; i--) v = (v * x + c[i]) % p;
        s1 += ctx.chi[v];
    }
    long long N1 = p + 2 + s1;
    // #C(F_{p^2}): elements u + v*w, w^2 = ns; chi2(y) = chi_p(Norm(y))
    // evaluate f at t = u + v w with F_p^2 arithmetic
    long long s2 = 0;
    int ns = ctx.ns;
    for (int u = 0; u < p; u++) {
        for (int v = 0; v < p; v++) {
            // Horner in F_p^2: acc = (ar, ai)
            long long ar = 0, ai = 0;
            for (int i = 6; i >= 0; i--) {
                // acc *= t ; t = (u, v)
                long long nr = (ar * u + ai * v % p * ns) % p;
                long long ni = (ar * v + ai * u) % p;
                ar = (nr + c[i]) % p;
                ai = ni;
            }
            // Norm = ar^2 - ns*ai^2
            long long nrm = ((ar * ar - ns * ai % p * ai) % p + p) % p;
            s2 += ctx.chi[nrm];
        }
    }
    long long N2 = (long long)p * p + 2 + s2;
    return (N1 * N1 + N2) / 2 - p;
}

static void reduce_family(long long a, long long b, long long c, long long d,
                          int p, int *co) {
    // (x^2-1)(x^2+ax+b)(x^2+cx+d) expanded over Z then reduced:
    // q = x^4+(a+c)x^3+(b+d+ac)x^2+(ad+bc)x+bd; f = (x^2-1)*q
    long long e3 = a + c, e2 = b + d + a * c, e1 = a * d + b * c, e0 = b * d;
    long long f6 = 1, f5 = e3, f4 = e2 - 1, f3 = e1 - e3, f2 = e0 - e2,
              f1 = -e1, f0 = -e0;
    long long fs[7] = {f0, f1, f2, f3, f4, f5, f6};
    for (int i = 0; i <= 6; i++) co[i] = (int)(((fs[i] % p) + p) % p);
}

int main(int argc, char **argv) {
    int H = argc > 1 ? atoi(argv[1]) : 100;
    int threads = argc > 2 ? atoi(argv[2]) : 8;
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    init_ctx();
    // ---- self-tests ----
    {
        // f14 = x^6+6x^4+4x^3+9x^2+4x+4 : Jacobian torsion contains Z/14,
        // and #J(F_3): compute and check 14 | #J(F_p) at good p in chain.
        long long f14[7] = {4, 4, 9, 4, 6, 0, 1};
        int ok14 = 0;
        for (int i = 0; i < NPR; i++) {
            int co[7];
            for (int k = 0; k <= 6; k++) co[k] = (int)(((f14[k] % PRS[i]) + PRS[i]) % PRS[i]);
            if (!sqfree_modp(co, PRS[i])) continue;
            long long J = jac_count(co, CTX[i]);
            if (J % 14 != 0) { fprintf(stderr, "SELFTEST FAIL f14 p=%d J=%lld\n", PRS[i], J); return 1; }
            ok14++;
        }
        if (ok14 < 5) { fprintf(stderr, "SELFTEST FAIL f14 goodprimes=%d\n", ok14); return 1; }
        // Flynn order-11 members t = 1, 2, 5:
        // F_t = x^6+2x^5+(2t+3)x^4+2x^3+(t^2+1)x^2+2t(1-t)x+t^2
        long long ts[3] = {1, 2, 5};
        for (int j = 0; j < 3; j++) {
            long long t = ts[j];
            long long F[7] = {t * t, 2 * t * (1 - t), t * t + 1, 2, 2 * t + 3, 2, 1};
            int okf = 0;
            for (int i = 0; i < NPR; i++) {
                int co[7];
                for (int k = 0; k <= 6; k++) co[k] = (int)(((F[k] % PRS[i]) + PRS[i]) % PRS[i]);
                if (!sqfree_modp(co, PRS[i])) continue;
                long long J = jac_count(co, CTX[i]);
                if (J % 11 != 0) { fprintf(stderr, "SELFTEST FAIL flynn t=%lld p=%d J=%lld\n", t, PRS[i], J); return 1; }
                okf++;
            }
            if (okf < 5) { fprintf(stderr, "SELFTEST FAIL flynn t=%lld goodprimes=%d\n", t, okf); return 1; }
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
                    if (a > c || (a == c && b > d)) continue;
                    tested++;
                    int ngood = 0;
                    bool ok = true;
                    for (int i = 0; i < NPR; i++) {
                        int co[7];
                        reduce_family(a, b, c, d, PRS[i], co);
                        if (!sqfree_modp(co, PRS[i])) continue;
                        long long J = jac_count(co, CTX[i]);
                        if (J % 11 != 0) { ok = false; break; }
                        ngood++;
                    }
                    if (!ok || ngood < MINGOOD) continue;
                    out++;
                    blen += snprintf(buf + blen, sizeof(buf) - blen,
                                     "%d %d %d %d\n", a, b, c, d);
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
