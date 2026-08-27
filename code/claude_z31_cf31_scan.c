/* Task B3 : chart I CF-31 fiber sieve for Z/31 torsion via ord(D_inf) = 31.
 *
 * Chart I (5 parameters):  y^2 = f(x),
 *   f = x^6 + c4 x^4 + c3 x^3 + c2 x^2 + c1 x + c0        (monic, c5 = 0)
 * marked class D_inf = [inf+ - inf-]; its EXACT order is computed by the
 * polynomial continued fraction of sqrt(f) (Platonov-Petrunin / polynomial
 * Pell): quasi-period total degree, INCLUDING deg a_0 = 3, equals the order
 * (see .claude/skills/pell-cf-order/SKILL.md).  ord(D_inf) = 31 is TWO
 * conditions on the 5-dim chart, so Sigma31 = {ord = 31} is a 3-fold:
 * |Sigma31(F_P)| ~ c*P^3 and the fiber over a fixed slice (c4,c3,c2) is
 * generically FINITE (~c points, c = O(1)).  The search fixes a rational
 * slice, enumerates the (c1,c0) fiber over F_P for 3-4 primes, matches
 * fiber points across primes by CRT, rationally reconstructs (c1,c0), and
 * confirms at an independent check prime; exact verification over Q is a
 * separate (Magma/PARI) step on the printed CAND lines.
 *
 * CF core notes:
 *   - exact port of CFOrderFp in code/claude_z31_magma_ref.m (which built
 *     data/claude_z31_vectors_chartI.tsv): same recursion, same guards
 *     (Qi = 0 abort; degenerate deg a_i < 1 at i >= 1 abort; total > cap
 *     early abort; exact-division check; quasi-period stop on CONSTANT
 *     NONZERO Q_i with i >= 1), same iteration bound cap + 10.
 *   - the core takes a GENERAL monic sextic (coeffs f0..f5, c5 free): the
 *     depression x -> x - c5/6 is only a chart normalization; selftest
 *     checks that raw (c5 != 0) and depressed forms give the same order.
 *     surface/sample/fibstat/fiber/sweep fix c5 = 0 (the chart).
 *   - squarefreeness (gcd(f, f') constant) is checked only AFTER a CF hit
 *     (hits are rare, the gcd is not free); non-squarefree hits are
 *     discarded (bad reduction, the CF/order theory does not apply there).
 *
 * Modes
 *   selftest [tsv]        (1) f14/f18/f28guess raw + depressed mod 5 primes
 *                         (2) ALL rows of data/claude_z31_vectors_chartI.tsv:
 *                             CF(cap 200) must equal N_inf exactly (0 = 0)
 *   bench P niter [seed]  single-thread throughput probe at the production
 *                         target (candidates/sec/core) + density estimate
 *   surface P [dumpfile]  exhaustive |Sigma31(F_P)| over F_P^5 (codim check)
 *   sample P n seed       sampled |Sigma31(F_P)| estimate (for large P)
 *   fibstat P n seed      fiber-size histogram over n random slices
 *   fiber P n4 d4 n3 d3 n2 d2
 *                         the (c1,c0) fiber over slice (n4/d4,n3/d3,n2/d2)
 *   sweep H npr P1..Pk lo hi recbound chkP int|rat
 *                         main search; slices (c4,c3,c2) of height <= H
 *                         (int: integers in [-H,H]; rat: rationals), shard
 *                         [lo,hi) of the c4-index space; fibers at the k
 *                         primes (give them ASCENDING: empty-fiber early
 *                         exit tests the cheapest prime first), CRT +
 *                         rational reconstruction with |num|,|den| <=
 *                         recbound, confirmation at check prime chkP.
 *
 * Env: Z31_TARGET (default 31) overrides the target order in bench/surface/
 * sample/fibstat/fiber/sweep -- used to validate the pipeline end-to-end on
 * the known order-14/18 rational points.  selftest ignores it (cap 200).
 *
 * Build: gcc -O3 -march=native -fopenmp -o claude_z31_cf31_scan claude_z31_cf31_scan.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t  i64;
typedef __uint128_t u128;

static int TARGET = 31;                       /* Z31_TARGET override        */

/* ------------------------------------------------------------- modular ops */
typedef struct { u32 p; u64 minv; u32 inv2; u32 *invtab; } Mod;

static inline u32 mmul(Mod *M, u32 a, u32 b)
{
    u64 z = (u64)a * (u64)b;
    u64 q = (u64)(((u128)z * (u128)M->minv) >> 64);
    u64 r = z - q * (u64)M->p;
    if (r >= M->p) r -= M->p;
    return (u32)r;
}
static inline u32 madd(Mod *M, u32 a, u32 b) { u32 s = a + b; return s >= M->p ? s - M->p : s; }
static inline u32 msub(Mod *M, u32 a, u32 b) { return a >= b ? a - b : a + M->p - b; }

static void mod_init(Mod *M, u32 p)
{
    M->p = p;
    M->minv = (u64)(((u128)1 << 64) / p);
    M->invtab = (u32 *)malloc(sizeof(u32) * p);
    M->invtab[0] = 0;
    if (p > 1) M->invtab[1] = 1;
    for (u32 i = 2; i < p; i++)
        M->invtab[i] = (u32)((u64)(p - p / i) * (u64)M->invtab[p % i] % p);
    M->inv2 = M->invtab[2 % p];
}

/* reduce n/d mod p; sets *bad = 1 if p | d */
static inline u32 ratmod(Mod *M, i64 n, i64 d, int *bad)
{
    i64 dd = d % (i64)M->p; if (dd < 0) dd += M->p;
    if (dd == 0) { *bad = 1; return 0; }
    i64 nn = n % (i64)M->p; if (nn < 0) nn += M->p;
    return mmul(M, (u32)nn, M->invtab[dd]);
}

static inline int pdeg(u32 *a, int hi)
{
    for (int i = hi; i >= 0; i--) if (a[i]) return i;
    return -1;
}

/* --------------------------------------------------------- CF order over Fp
 * f monic degree 6 given as f6[0..5] (f6[5] = coeff of x^5, may be nonzero;
 * leading coeff 1 implicit).  Exact port of CFOrderFp (Magma reference):
 * returns the exact order of D_inf if a quasi-period closes with total
 * degree <= cap, else 0 (budget exceeded / degenerate).  Total INCLUDES
 * deg a_0 = 3.
 */
static int cford_fp(Mod *M, u32 *f6, int cap)
{
    u32 f[7], s[4], P[7], Q[7], A[7], Pn[7], R[7], T[9];
    int i, j, k, dQ, dA, tot;

    for (i = 0; i < 6; i++) f[i] = f6[i];
    f[6] = 1;

    /* polynomial part of sqrt(f): s = x^3 + s2 x^2 + s1 x + s0 */
    s[3] = 1;
    s[2] = mmul(M, f[5], M->inv2);
    s[1] = mmul(M, msub(M, f[4], mmul(M, s[2], s[2])), M->inv2);
    s[0] = mmul(M, msub(M, f[3], mmul(M, madd(M, s[2], s[2]), s[1])), M->inv2);

    memset(P, 0, sizeof(P));
    memset(Q, 0, sizeof(Q));
    Q[0] = 1;
    tot = 0;

    for (i = 0; i <= cap + 10; i++) {
        dQ = pdeg(Q, 6);
        if (dQ < 0) return 0;                    /* Q == 0 */

        /* A = (P + s) div Q ; P + s has degree <= 3 */
        for (j = 0; j < 7; j++) T[j] = P[j];
        for (j = 0; j < 4; j++) T[j] = madd(M, T[j], s[j]);
        int dT = pdeg(T, 6);
        memset(A, 0, sizeof(A));
        if (dT >= dQ) {
            u32 ilc = M->invtab[Q[dQ]];
            for (k = dT - dQ; k >= 0; k--) {
                u32 c = mmul(M, T[k + dQ], ilc);
                A[k] = c;
                if (c) for (j = 0; j <= dQ; j++)
                    T[k + j] = msub(M, T[k + j], mmul(M, c, Q[j]));
            }
        }
        dA = pdeg(A, 6);
        if (dA < 0) return 0;                    /* zero partial quotient */
        if (i >= 1 && dA < 1) return 0;          /* degenerate orbit (Magma parity) */
        tot += dA;
        if (tot > cap) return 0;                 /* early abort */

        /* Pn = A*Q - P */
        memset(Pn, 0, sizeof(Pn));
        for (j = 0; j <= dA; j++) {
            if (!A[j]) continue;
            for (k = 0; k <= dQ; k++) {
                if (j + k > 6) { if (mmul(M, A[j], Q[k])) return 0; continue; }
                Pn[j + k] = madd(M, Pn[j + k], mmul(M, A[j], Q[k]));
            }
        }
        for (j = 0; j < 7; j++) Pn[j] = msub(M, Pn[j], P[j]);
        if (pdeg(Pn, 6) > 3) return 0;

        /* R = f - Pn^2 (degree <= 6) */
        memset(T, 0, sizeof(T));
        for (j = 0; j <= 3; j++) {
            if (!Pn[j]) continue;
            for (k = 0; k <= 3; k++) T[j + k] = madd(M, T[j + k], mmul(M, Pn[j], Pn[k]));
        }
        for (j = 0; j < 7; j++) R[j] = msub(M, f[j], T[j]);

        /* Qn = R / Q, remainder must vanish (always does in valid orbits) */
        int dR = pdeg(R, 6);
        u32 Qn[7];
        memset(Qn, 0, sizeof(Qn));
        if (dR >= dQ) {
            u32 ilc = M->invtab[Q[dQ]];
            for (k = dR - dQ; k >= 0; k--) {
                u32 c = mmul(M, R[k + dQ], ilc);
                Qn[k] = c;
                if (c) for (j = 0; j <= dQ; j++)
                    R[k + j] = msub(M, R[k + j], mmul(M, c, Q[j]));
            }
        }
        if (pdeg(R, 6) >= 0) return 0;           /* nonexact division */

        for (j = 0; j < 7; j++) { P[j] = Pn[j]; Q[j] = Qn[j]; }
        int dQn = pdeg(Q, 6);
        if (i >= 1 && dQn == 0) return tot;      /* constant nonzero: period */
        if (dQn < 0) return 0;
    }
    return 0;
}

/* squarefree test: gcd(f, f') constant, f monic sextic f6[0..5] + lc 1 */
static int sqfree6(Mod *M, u32 *f6)
{
    u32 a[7], b[7], t;
    int j, k;
    for (j = 0; j < 6; j++) a[j] = f6[j];
    a[6] = 1;
    for (j = 0; j < 5; j++) b[j] = mmul(M, (u32)((j + 1) % M->p), f6[j + 1]);
    b[5] = 6 % M->p;
    b[6] = 0;
    while (1) {
        int db = pdeg(b, 6);
        if (db < 0) return 0;                    /* gcd = a, deg >= 1 */
        if (db == 0) return 1;                   /* gcd constant       */
        int da = pdeg(a, 6);
        if (da >= db) {
            u32 ilc = M->invtab[b[db]];
            for (k = da - db; k >= 0; k--) {
                u32 c = mmul(M, a[k + db], ilc);
                if (c) for (j = 0; j <= db; j++)
                    a[k + j] = msub(M, a[k + j], mmul(M, c, b[j]));
            }
        }
        for (j = 0; j < 7; j++) { t = a[j]; a[j] = b[j]; b[j] = t; }
    }
}

/* chart-I candidate test: f = x^6 + a4 x^4 + a3 x^3 + a2 x^2 + c1 x + c0 */
static inline int hit31(Mod *M, u32 a4, u32 a3, u32 a2, u32 c1, u32 c0)
{
    u32 f6[6];
    f6[0] = c0; f6[1] = c1; f6[2] = a2; f6[3] = a3; f6[4] = a4; f6[5] = 0;
    if (cford_fp(M, f6, TARGET) != TARGET) return 0;
    return sqfree6(M, f6);
}

/* ------------------------------------------------------------- rational rec */
static int ratrecon(i64 x, i64 m, i64 bound, i64 *num, i64 *den)
{
    i64 r0 = m, r1 = x % m; if (r1 < 0) r1 += m;
    i64 t0 = 0, t1 = 1;
    while (r1 > bound) {
        if (r1 == 0) return 0;
        i64 q = r0 / r1;
        i64 r2 = r0 - q * r1, t2 = t0 - q * t1;
        r0 = r1; r1 = r2; t0 = t1; t1 = t2;
    }
    if (t1 == 0) return 0;
    i64 a = r1, b = t1;
    if (b < 0) { b = -b; a = -a; }
    if (b > bound) return 0;
    i64 g0 = a < 0 ? -a : a, g1 = b;
    while (g1) { i64 t = g0 % g1; g0 = g1; g1 = t; }
    if (g0 != 1) return 0;
    *num = a; *den = b;
    return 1;
}

static i64 crt2(i64 a1, i64 m1, i64 a2, i64 m2, i64 *M12)
{
    i64 m = m1 * m2;
    i64 t = 0, nt = 1, r = m2, nr = m1 % m2;
    while (nr) { i64 q = r / nr; i64 tmp = t - q * nt; t = nt; nt = tmp; tmp = r - q * nr; r = nr; nr = tmp; }
    if (t < 0) t += m2;
    i64 x = a1 + m1 * (((a2 - a1) % m2 + m2) % m2 * t % m2);
    x %= m; if (x < 0) x += m;
    *M12 = m;
    return x;
}

/* ------------------------------------------------------------------ rng     */
static inline u64 xs64(u64 *st)
{
    u64 x = *st;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    *st = x;
    return x;
}
static inline u64 smix(u64 z)
{
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/* ------------------------------------------------------------------ selftest */
/* Taylor shift mod p: g(x) = f(x + t), f monic sextic (f6[0..5], lc 1);
   writes g7[0..6] (g7[6] must come out 1) */
static void shift6(Mod *M, u32 *f6, u32 t, u32 *g7)
{
    int C[7][7];
    u32 f[7], tp[7];
    int j, k;
    for (j = 0; j <= 6; j++) { C[j][0] = 1; for (k = 1; k <= j; k++)
        C[j][k] = (k == j) ? 1 : C[j - 1][k - 1] + C[j - 1][k]; }
    for (j = 0; j < 6; j++) f[j] = f6[j];
    f[6] = 1;
    tp[0] = 1 % M->p;
    for (j = 1; j <= 6; j++) tp[j] = mmul(M, tp[j - 1], t);
    for (k = 0; k <= 6; k++) {
        u32 acc = 0;
        for (j = k; j <= 6; j++)
            acc = madd(M, acc, mmul(M, mmul(M, (u32)(C[j][k] % M->p), f[j]), tp[j - k]));
        g7[k] = acc;
    }
}

static void do_selftest(char *tsvpath)
{
    /* (a) skill vectors, raw (general c5) AND depressed x -> x - c5/6 */
    struct { char *name; i64 c[6]; int expect; } kv[3] = {
        { "f14",      {  4,  4,  9,   4,  6,  0 }, 14 },
        { "f18",      { -8, 16, -7,  -2, 11, -2 }, 18 },
        { "f28guess", { 28, 24, -3, -14, -5,  2 },  7 },
    };
    u32 primes[5] = { 101, 211, 499, 1009, 2003 };
    int failA = 0, passA = 0;
    for (int t = 0; t < 5; t++) {
        Mod M; mod_init(&M, primes[t]);
        for (int v = 0; v < 3; v++) {
            u32 f6[6], g7[7];
            for (int j = 0; j < 6; j++) {
                i64 c = kv[v].c[j] % (i64)M.p; if (c < 0) c += M.p;
                f6[j] = (u32)c;
            }
            int oraw = cford_fp(&M, f6, 200);
            u32 tsh = mmul(&M, msub(&M, 0, f6[5]), M.invtab[6 % M.p]);  /* -c5/6 */
            shift6(&M, f6, tsh, g7);
            int okdep = (g7[6] == 1) && (g7[5] == 0);
            int odep = cford_fp(&M, g7, 200);
            int ok = okdep && oraw == kv[v].expect && odep == kv[v].expect;
            if (ok) passA++; else failA++;
            printf("SELFTEST_VEC p=%u %-9s raw=%d depressed=%d expect=%d %s\n",
                   M.p, kv[v].name, oraw, odep, kv[v].expect, ok ? "PASS" : "FAIL");
        }
        free(M.invtab);
    }

    /* (b) every row of the chart-I reference TSV: CF(cap 200) == N_inf */
    FILE *fp = fopen(tsvpath, "r");
    if (!fp) { printf("SELFTEST_TSV cannot open %s\nSELFTEST_RESULT FAIL\n", tsvpath); return; }
    Mod cache[8]; u32 cachep[8]; int ncache = 0;
    long long P, c4, c3, c2, c1, c0, N;
    int rows = 0, passB = 0, failB = 0;
    while (fscanf(fp, "%lld %lld %lld %lld %lld %lld %lld",
                  &P, &c4, &c3, &c2, &c1, &c0, &N) == 7) {
        rows++;
        Mod *M = NULL;
        for (int i = 0; i < ncache; i++) if (cachep[i] == (u32)P) M = &cache[i];
        if (!M) { mod_init(&cache[ncache], (u32)P); cachep[ncache] = (u32)P; M = &cache[ncache]; ncache++; }
        u32 f6[6];
        f6[0] = (u32)c0; f6[1] = (u32)c1; f6[2] = (u32)c2;
        f6[3] = (u32)c3; f6[4] = (u32)c4; f6[5] = 0;
        int cf = cford_fp(M, f6, 200);
        if (cf == (int)N) passB++;
        else {
            failB++;
            printf("SELFTEST_TSV_FAIL P=%lld c=[%lld,%lld,%lld,%lld,%lld] N_inf=%lld cf=%d\n",
                   P, c4, c3, c2, c1, c0, N, cf);
        }
    }
    fclose(fp);
    for (int i = 0; i < ncache; i++) free(cache[i].invtab);
    printf("SELFTEST_TSV rows=%d pass=%d fail=%d\n", rows, passB, failB);
    printf("SELFTEST_VECTORS pass=%d fail=%d (15 = 3 vectors x 5 primes)\n", passA, failA);
    printf("SELFTEST_RESULT %s\n", (failA == 0 && failB == 0 && rows == 626) ? "PASS" : "FAIL");
}

/* ------------------------------------------------------------------ bench   */
static void do_bench(u32 p, u64 niter, u64 seed)
{
    Mod M; mod_init(&M, p);
    u64 st = smix(seed ? seed : 31);
    u64 hits = 0, sf = 0;
    double t0 = omp_get_wtime();
    for (u64 it = 0; it < niter; it++) {
        u32 f6[6];
        f6[5] = 0;
        f6[4] = (u32)(xs64(&st) % p); f6[3] = (u32)(xs64(&st) % p);
        f6[2] = (u32)(xs64(&st) % p); f6[1] = (u32)(xs64(&st) % p);
        f6[0] = (u32)(xs64(&st) % p);
        if (cford_fp(&M, f6, TARGET) == TARGET) {
            hits++;
            if (sqfree6(&M, f6)) sf++;
        }
    }
    double dt = omp_get_wtime() - t0;
    printf("BENCH p=%u target=%d niter=%llu time=%.3fs  %.1f ns/cand  %.0f cand/s/core"
           "  hits=%llu sqfree=%llu density=%.3e c_impl=%.3f\n",
           p, TARGET, (unsigned long long)niter, dt, dt / (double)niter * 1e9,
           (double)niter / dt, (unsigned long long)hits, (unsigned long long)sf,
           (double)sf / (double)niter, (double)sf / (double)niter * (double)p * (double)p);
    free(M.invtab);
}

/* ------------------------------------------------------------------ surface */
static void do_surface(u32 p, char *dumpfile)
{
    Mod M; mod_init(&M, p);
    u64 total = 0, sing = 0;
    FILE *fp = dumpfile ? fopen(dumpfile, "w") : NULL;
    double t0 = omp_get_wtime();
#pragma omp parallel reduction(+:total,sing)
    {
        char *buf = fp ? (char *)malloc(1 << 20) : NULL; size_t bl = 0;
#pragma omp for collapse(2) schedule(dynamic, 1)
        for (int a4 = 0; a4 < (int)p; a4++)
            for (int a3 = 0; a3 < (int)p; a3++) {
                u32 f6[6];
                f6[5] = 0; f6[4] = (u32)a4; f6[3] = (u32)a3;
                for (u32 a2 = 0; a2 < p; a2++) {
                    f6[2] = a2;
                    for (u32 c1 = 0; c1 < p; c1++) {
                        f6[1] = c1;
                        for (u32 c0 = 0; c0 < p; c0++) {
                            f6[0] = c0;
                            if (cford_fp(&M, f6, TARGET) != TARGET) continue;
                            if (!sqfree6(&M, f6)) { sing++; continue; }
                            total++;
                            if (fp) {
                                bl += (size_t)sprintf(buf + bl, "%d %d %u %u %u\n", a4, a3, a2, c1, c0);
                                if (bl > (1 << 20) - 64) {
#pragma omp critical
                                    fwrite(buf, 1, bl, fp);
                                    bl = 0;
                                }
                            }
                        }
                    }
                }
            }
        if (fp && bl) {
#pragma omp critical
            fwrite(buf, 1, bl, fp);
        }
        if (buf) free(buf);
    }
    if (fp) fclose(fp);
    double dt = omp_get_wtime() - t0;
    double p3 = (double)p * p * p, p5 = p3 * (double)p * p;
    printf("SURFACE p=%u target=%d |Sigma(F_p)|=%llu sing31=%llu  P^3=%.0f  c=|Sigma|/P^3=%.4f"
           "  density=%.3e  (%.1f s, %.2f Mcand/s total)\n",
           p, TARGET, (unsigned long long)total, (unsigned long long)sing, p3,
           (double)total / p3, (double)total / p5, dt, p5 / dt / 1e6);
    fflush(stdout);
    free(M.invtab);
}

static void do_sample(u32 p, u64 nsamp, u64 seed)
{
    Mod M; mod_init(&M, p);
    u64 total = 0, sing = 0;
    double t0 = omp_get_wtime();
#pragma omp parallel for reduction(+:total,sing) schedule(static)
    for (i64 it = 0; it < (i64)nsamp; it++) {
        u64 st = smix(seed ^ (u64)it);
        u32 f6[6];
        f6[5] = 0;
        f6[4] = (u32)(xs64(&st) % p); f6[3] = (u32)(xs64(&st) % p);
        f6[2] = (u32)(xs64(&st) % p); f6[1] = (u32)(xs64(&st) % p);
        f6[0] = (u32)(xs64(&st) % p);
        if (cford_fp(&M, f6, TARGET) != TARGET) continue;
        if (!sqfree6(&M, f6)) { sing++; continue; }
        total++;
    }
    double dt = omp_get_wtime() - t0;
    double p3 = (double)p * p * p, p5 = p3 * (double)p * p;
    double dens = (double)total / (double)nsamp;
    double est = dens * p5, rel = total ? 1.0 / __builtin_sqrt((double)total) : 0;
    printf("SAMPLE p=%u target=%d n=%llu hits=%llu sing31=%llu  est|Sigma|=%.4e (+-%.1f%%)"
           "  c=est/P^3=%.4f  density=%.3e  (%.1f s, %.2f Mcand/s)\n",
           p, TARGET, (unsigned long long)nsamp, (unsigned long long)total,
           (unsigned long long)sing, est, 100 * rel, est / p3, dens, dt,
           (double)nsamp / dt / 1e6);
    free(M.invtab);
}

/* ------------------------------------------------------------------ fibers  */
#define FIBCAP 512
#define MAXPR  6

static int fiber_of(Mod *M, u32 a4, u32 a3, u32 a2, u32 *out, int maxout)
{
    int n = 0;
    for (u32 c1 = 0; c1 < M->p; c1++)
        for (u32 c0 = 0; c0 < M->p; c0++)
            if (hit31(M, a4, a3, a2, c1, c0)) {
                if (n < maxout) { out[2 * n] = c1; out[2 * n + 1] = c0; }
                n++;
                if (n >= maxout) return n;
            }
    return n;
}

static void do_fibstat(u32 p, u64 nslices, u64 seed)
{
    Mod M; mod_init(&M, p);
    u64 histo[17]; memset(histo, 0, sizeof(histo));
    u64 fsum = 0;
    double t0 = omp_get_wtime();
#pragma omp parallel
    {
        u64 lh[17]; memset(lh, 0, sizeof(lh));
        u64 ls = 0;
        u32 out[2 * FIBCAP];
#pragma omp for schedule(dynamic, 1)
        for (i64 it = 0; it < (i64)nslices; it++) {
            u64 st = smix(seed ^ (u64)(it + 1));
            u32 a4 = (u32)(xs64(&st) % p), a3 = (u32)(xs64(&st) % p), a2 = (u32)(xs64(&st) % p);
            int n = fiber_of(&M, a4, a3, a2, out, FIBCAP);
            lh[n > 16 ? 16 : n]++;
            ls += (u64)n;
        }
#pragma omp critical
        { for (int i = 0; i < 17; i++) histo[i] += lh[i]; fsum += ls; }
    }
    double dt = omp_get_wtime() - t0;
    printf("FIBSTAT p=%u target=%d slices=%llu avg_fiber=%.3f empty=%llu (%.1f%%)  histo:",
           p, TARGET, (unsigned long long)nslices, (double)fsum / (double)nslices,
           (unsigned long long)histo[0], 100.0 * (double)histo[0] / (double)nslices);
    for (int i = 0; i < 17; i++) if (histo[i]) printf(" %d:%llu", i, (unsigned long long)histo[i]);
    printf("  (%.1f s wall)\n", dt);
    free(M.invtab);
}

static void do_fiber(u32 p, i64 n4, i64 d4, i64 n3, i64 d3, i64 n2, i64 d2)
{
    Mod M; mod_init(&M, p);
    int bad = 0;
    u32 a4 = ratmod(&M, n4, d4, &bad), a3 = ratmod(&M, n3, d3, &bad), a2 = ratmod(&M, n2, d2, &bad);
    if (bad) { printf("FIBER p=%u bad reduction of slice\n", p); free(M.invtab); return; }
    u32 out[2 * FIBCAP];
    double t0 = omp_get_wtime();
    int n = fiber_of(&M, a4, a3, a2, out, FIBCAP);
    printf("FIBER p=%u target=%d slice=(%lld/%lld,%lld/%lld,%lld/%lld) -> (%u,%u,%u)  #fiber=%d  (%.2f s)\n",
           p, TARGET, (long long)n4, (long long)d4, (long long)n3, (long long)d3,
           (long long)n2, (long long)d2, a4, a3, a2, n, omp_get_wtime() - t0);
    for (int i = 0; i < n && i < FIBCAP; i++)
        printf("   (c1,c0) = (%u, %u)\n", out[2 * i], out[2 * i + 1]);
    free(M.invtab);
}

/* ------------------------------------------------------------------ sweep   */
typedef struct { i64 n, d; } Rat;

static int build_rats(int H, int intonly, Rat *out, int cap)
{
    int n = 0;
    for (i64 den = 1; den <= (intonly ? 1 : (i64)H); den++)
        for (i64 num = -H; num <= H; num++) {
            i64 a = num < 0 ? -num : num, b = den;
            while (b) { i64 t = a % b; a = b; b = t; }
            if (num == 0 && den != 1) continue;
            if (!(num == 0 && den == 1) && a != 1) continue;
            if (n < cap) { out[n].n = num; out[n].d = den; }
            n++;
        }
    return n;
}

/* confirmation of a reconstructed candidate at the check prime */
static int confirm_chk(Mod *MC, Rat s4, Rat s3, Rat s2, i64 n1, i64 d1, i64 n0, i64 d0)
{
    int bad = 0;
    u32 a4 = ratmod(MC, s4.n, s4.d, &bad), a3 = ratmod(MC, s3.n, s3.d, &bad);
    u32 a2 = ratmod(MC, s2.n, s2.d, &bad);
    u32 c1 = ratmod(MC, n1, d1, &bad), c0 = ratmod(MC, n0, d0, &bad);
    if (bad) return -1;
    return hit31(MC, a4, a3, a2, c1, c0);
}

static void do_sweep(int H, int npr, u32 *prs, int lo, int hi, i64 recbound,
                     u32 chkp, int intonly)
{
    int rcap = 4 * (H + 1) * (H + 1);
    Rat *rats = (Rat *)malloc(sizeof(Rat) * (size_t)rcap);
    int nr = build_rats(H, intonly, rats, rcap);
    if (hi > nr || hi <= 0) hi = nr;
    if (lo < 0) lo = 0;
    int nsh = hi - lo;

    u128 mprod = 1;
    for (int t = 0; t < npr; t++) mprod *= prs[t];
    if (mprod > (u128)9000000000000000000ULL) {
        printf("SWEEP_ABORT prime product overflows int64\n"); free(rats); return;
    }
    for (int t = 0; t < npr; t++)
        if (prs[t] == chkp) printf("SWEEP_WARN check prime %u is also a sieve prime\n", chkp);
    printf("SWEEP target=%d H=%d %s #slice-coords=%d shard=[%d,%d) slices=%lld nprimes=%d",
           TARGET, H, intonly ? "int" : "rat", nr, lo, hi,
           (long long)nsh * nr * nr, npr);
    for (int t = 0; t < npr; t++) printf(" %u", prs[t]);
    printf(" M=%.3e recbound=%lld chkP=%u threads=%d\n",
           (double)(u64)mprod, (long long)recbound, chkp, omp_get_max_threads());
    fflush(stdout);

    Mod M[MAXPR], MC;
    for (int t = 0; t < npr; t++) mod_init(&M[t], prs[t]);
    mod_init(&MC, chkp);
    u64 done = 0, skipped = 0, cand = 0, nconf = 0, tuples = 0, overflow = 0;
    u64 empty[MAXPR], fibsum[MAXPR], fibcnt[MAXPR];
    memset(empty, 0, sizeof(empty));
    memset(fibsum, 0, sizeof(fibsum));
    memset(fibcnt, 0, sizeof(fibcnt));
    u64 prog = 0;
    double t0 = omp_get_wtime();

#pragma omp parallel reduction(+:done,skipped,cand,nconf,tuples,overflow) \
                     reduction(+:empty[:MAXPR],fibsum[:MAXPR],fibcnt[:MAXPR])
    {
        u32 *fib = (u32 *)malloc(sizeof(u32) * 2 * FIBCAP * MAXPR);
        int nf[MAXPR];
#pragma omp for collapse(3) schedule(dynamic, 2)
        for (int ii = 0; ii < nsh; ii++)
            for (int j = 0; j < nr; j++)
                for (int k = 0; k < nr; k++) {
                    int i = lo + ii;
                    Rat s4 = rats[i], s3 = rats[j], s2 = rats[k];
                    int ok = 1, bad = 0;
                    for (int t = 0; t < npr && ok; t++) {
                        u32 a4 = ratmod(&M[t], s4.n, s4.d, &bad);
                        u32 a3 = ratmod(&M[t], s3.n, s3.d, &bad);
                        u32 a2 = ratmod(&M[t], s2.n, s2.d, &bad);
                        if (bad) { ok = 0; skipped++; break; }
                        nf[t] = fiber_of(&M[t], a4, a3, a2, fib + 2 * FIBCAP * t, FIBCAP);
                        fibcnt[t]++;
                        fibsum[t] += (u64)nf[t];
                        if (nf[t] == 0) { empty[t]++; ok = 0; }
                        else if (nf[t] >= FIBCAP) { overflow++; ok = 0; }
                    }
                    if (!bad) done++;
                    if (ok) {
                        int idx[MAXPR];
                        for (int t = 0; t < npr; t++) idx[t] = 0;
                        while (1) {
                            i64 mm = prs[0];
                            i64 x = fib[2 * idx[0]], y = fib[2 * idx[0] + 1];
                            for (int t = 1; t < npr; t++) {
                                i64 nm;
                                x = crt2(x, mm, fib[2 * FIBCAP * t + 2 * idx[t]], prs[t], &nm);
                                y = crt2(y, mm, fib[2 * FIBCAP * t + 2 * idx[t] + 1], prs[t], &nm);
                                mm = nm;
                            }
                            tuples++;
                            i64 n1, d1, n0, d0;
                            if (ratrecon(x, mm, recbound, &n1, &d1)
                                && ratrecon(y, mm, recbound, &n0, &d0)) {
                                int cf = confirm_chk(&MC, s4, s3, s2, n1, d1, n0, d0);
                                cand++;
                                if (cf == 1) nconf++;
#pragma omp critical
                                {
                                    printf("CAND c4=%lld/%lld c3=%lld/%lld c2=%lld/%lld "
                                           "c1=%lld/%lld c0=%lld/%lld conf=%d\n",
                                           (long long)s4.n, (long long)s4.d,
                                           (long long)s3.n, (long long)s3.d,
                                           (long long)s2.n, (long long)s2.d,
                                           (long long)n1, (long long)d1,
                                           (long long)n0, (long long)d0, cf);
                                    fflush(stdout);
                                }
                            }
                            int t = 0;
                            while (t < npr && ++idx[t] >= nf[t]) { idx[t] = 0; t++; }
                            if (t == npr) break;
                        }
                    }
                    u64 mydone;
#pragma omp atomic capture
                    mydone = ++prog;
                    if (mydone % 2000 == 0) {
#pragma omp critical
                        {
                            double el = omp_get_wtime() - t0;
                            printf("PROGRESS slices=%llu/%lld cand=%llu  %.1f s  (%.2f slices/s)\n",
                                   (unsigned long long)mydone, (long long)nsh * nr * nr,
                                   (unsigned long long)cand, el, (double)mydone / el);
                            fflush(stdout);
                        }
                    }
                }
        free(fib);
    }
    double el = omp_get_wtime() - t0;
    printf("SWEEP_DONE H=%d shard=[%d,%d) slices=%llu skipped=%llu tuples=%llu "
           "cand=%llu conf=%llu overflow=%llu  %.1f s  (%.3f s/slice avg)\n",
           H, lo, hi, (unsigned long long)done, (unsigned long long)skipped,
           (unsigned long long)tuples, (unsigned long long)cand,
           (unsigned long long)nconf, (unsigned long long)overflow, el,
           el * omp_get_max_threads() / (done ? (double)done : 1.0));
    for (int t = 0; t < npr; t++)
        printf("SWEEP_FIBSTAT p=%u fibers=%llu empty=%llu (%.1f%%) avg_size=%.3f\n",
               prs[t], (unsigned long long)fibcnt[t], (unsigned long long)empty[t],
               fibcnt[t] ? 100.0 * (double)empty[t] / (double)fibcnt[t] : 0.0,
               fibcnt[t] ? (double)fibsum[t] / (double)fibcnt[t] : 0.0);
    fflush(stdout);
    for (int t = 0; t < npr; t++) free(M[t].invtab);
    free(MC.invtab);
    free(rats);
}

int main(int argc, char **argv)
{
    char *tenv = getenv("Z31_TARGET");
    if (tenv) TARGET = atoi(tenv);
    if (argc < 2) { fprintf(stderr, "usage: %s selftest|bench|surface|sample|fibstat|fiber|sweep ...\n", argv[0]); return 1; }
    if (!strcmp(argv[1], "selftest"))
        do_selftest(argc > 2 ? argv[2] : "data/claude_z31_vectors_chartI.tsv");
    else if (!strcmp(argv[1], "bench") && argc >= 4)
        do_bench((u32)atoi(argv[2]), (u64)atoll(argv[3]), argc > 4 ? (u64)atoll(argv[4]) : 31);
    else if (!strcmp(argv[1], "surface") && argc >= 3)
        do_surface((u32)atoi(argv[2]), argc > 3 ? argv[3] : NULL);
    else if (!strcmp(argv[1], "sample") && argc >= 5)
        do_sample((u32)atoi(argv[2]), (u64)atoll(argv[3]), (u64)atoll(argv[4]));
    else if (!strcmp(argv[1], "fibstat") && argc >= 5)
        do_fibstat((u32)atoi(argv[2]), (u64)atoll(argv[3]), (u64)atoll(argv[4]));
    else if (!strcmp(argv[1], "fiber") && argc >= 9)
        do_fiber((u32)atoi(argv[2]), atoll(argv[3]), atoll(argv[4]), atoll(argv[5]),
                 atoll(argv[6]), atoll(argv[7]), atoll(argv[8]));
    else if (!strcmp(argv[1], "sweep") && argc >= 8) {
        int H = atoi(argv[2]);
        int npr = atoi(argv[3]);
        if (npr < 1 || npr > MAXPR || argc < 8 + npr) { fprintf(stderr, "bad sweep args\n"); return 1; }
        u32 prs[MAXPR];
        for (int t = 0; t < npr; t++) prs[t] = (u32)atoi(argv[4 + t]);
        int lo = atoi(argv[4 + npr]), hi = atoi(argv[5 + npr]);
        i64 rb = atoll(argv[6 + npr]);
        u32 chkp = (u32)atoi(argv[7 + npr]);
        int intonly = (argc > 8 + npr) && !strcmp(argv[8 + npr], "int");
        do_sweep(H, npr, prs, lo, hi, rb, chkp, intonly);
    }
    else { fprintf(stderr, "unknown mode or missing args\n"); return 1; }
    return 0;
}
