/* Lane 3 / route B2' : the [1,1,2,2] + ord(D_inf)=11 locus.
 *
 * MONIC chart:  y^2 = f(x) = x (x-1) (x^2 + p1 x + p2) (x^2 + p3 x + p4)
 *   factor type [1,1,2,2]  =>  2-rank 2  =>  torsion contains (Z/2)^2
 *   ord(D_inf) = 11        =>  torsion contains Z/2 x Z/22 = [2,22]
 *
 * Both known [2,22] witnesses live here:
 *   19044.h.2  (p1,p2,p3,p4) = (-3, 8, 4, 27)
 *   BLP C4corr (p1,p2,p3,p4) = (-34/15, 32/45, 31/15, -2/9)
 *
 * ord(D_inf) = 11 is TWO conditions (dim J = 2), so Sigma'' = {ord=11} is a
 * SURFACE in A^4.  A blind 4-parameter box sweep therefore hits it with
 * probability 0.  The correct search fixes (p1,p2) -- 2 free parameters --
 * and solves the resulting 0-dimensional fiber for rational (p3,p4), which we
 * do by exhaustive mod-P fiber enumeration at three primes + CRT + rational
 * reconstruction, then exact verification over Q.
 *
 * Modes
 *   selftest              CF order over F_p on the validated vectors + anchors
 *   surface P [dumpfile]  exhaustive |Sigma''(F_P)| over F_P^4 (codim check)
 *   fiber P n1 d1 n2 d2   fiber of Sigma'' over (p1,p2)=(n1/d1,n2/d2) mod P
 *   sweep H P Q R lo hi   main search; (p1,p2) rationals of height <= H,
 *                         shard [lo,hi) of the p1-index space
 *
 * Build: gcc -O3 -march=native -fopenmp -o claude_ov_b2p_scan claude_ov_b2p_scan.c
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

/* --------------------------------------------------------- CF order over Fp
 * f is monic of degree 6, given as f[0..5] (f[6] = 1 implicit).
 * Returns the exact order of D_inf if a quasi-period is found with total
 * degree <= cap, else 0.  Early-aborts as soon as the running total exceeds
 * cap -- for cap = 11 that is at most 9 partial quotients.
 */
#define PDEG 8
static inline int pdeg(u32 *a, int hi)
{
    for (int i = hi; i >= 0; i--) if (a[i]) return i;
    return -1;
}

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

    for (i = 0; i <= 12; i++) {
        dQ = pdeg(Q, 6);
        if (dQ < 0) return 0;                    /* Q == 0 */

        /* A = (P + s) div Q ; P+s has degree <= 3 */
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
        if (dA < 0) return 0;                    /* degenerate */
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

        /* R must be divisible by Q */
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
        if (pdeg(R, 6) >= 0) return 0;           /* remainder nonzero */

        for (j = 0; j < 7; j++) { P[j] = Pn[j]; Q[j] = Qn[j]; }
        int dQn = pdeg(Q, 6);
        if (i >= 1 && dQn == 0) return tot;
        if (dQn < 0) return 0;
    }
    return 0;
}

/* build f = x(x-1)(x^2+p1 x+p2)(x^2+p3 x+p4) mod p, return 0 if singular-ish */
static inline void build_quartic(Mod *M, u32 p1, u32 p2, u32 *q4)
{
    /* x(x-1) = x^2 - x  ;  times (x^2 + p1 x + p2) */
    u32 a[3]; a[0] = 0; a[1] = M->p - 1; a[2] = 1;   /* x^2 - x */
    if (a[1] == M->p) a[1] = 0;
    u32 b[3]; b[0] = p2; b[1] = p1; b[2] = 1;
    for (int i = 0; i < 5; i++) q4[i] = 0;
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
        q4[i + j] = madd(M, q4[i + j], mmul(M, a[i], b[j]));
}

static inline void build_f(Mod *M, u32 *q4, u32 p3, u32 p4, u32 *f6)
{
    /* f = q4 * (x^2 + p3 x + p4), degree 6, leading 1 */
    u32 c[3]; c[0] = p4; c[1] = p3; c[2] = 1;
    u32 t[7]; for (int i = 0; i < 7; i++) t[i] = 0;
    for (int i = 0; i < 5; i++) { if (!q4[i]) continue;
        for (int j = 0; j < 3; j++) t[i + j] = madd(M, t[i + j], mmul(M, q4[i], c[j])); }
    for (int i = 0; i < 6; i++) f6[i] = t[i];
}

/* squarefree-ness of x(x-1)q1 q2 : the four factors pairwise coprime & quadratics separable */
static inline int nondegenerate(Mod *M, u32 p1, u32 p2, u32 p3, u32 p4)
{
    u32 d1 = msub(M, mmul(M, p1, p1), mmul(M, 4 % M->p, p2));
    u32 d2 = msub(M, mmul(M, p3, p3), mmul(M, 4 % M->p, p4));
    if (!d1 || !d2) return 0;
    if (!p2 || !p4) return 0;                                     /* q_i(0) != 0 */
    u32 o1 = madd(M, madd(M, 1, p1), p2), o2 = madd(M, madd(M, 1, p3), p4);
    if (!o1 || !o2) return 0;                                     /* q_i(1) != 0 */
    /* Res(q1,q2) = (p4-p2)^2 - (p4-p2)(p3-p1)p1 + (p3-p1)^2 p2 */
    u32 u = msub(M, p4, p2), v = msub(M, p3, p1);
    u32 res = madd(M, msub(M, mmul(M, u, u), mmul(M, mmul(M, u, v), p1)), mmul(M, mmul(M, v, v), p2));
    if (!res) return 0;
    return 1;
}

/* ------------------------------------------------------------- rational rec */
/* reconstruct a/b = x mod m with |a|,|b| <= bound ; returns 1 on success */
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
    /* gcd(a,b) must be 1 */
    i64 g0 = a < 0 ? -a : a, g1 = b;
    while (g1) { i64 t = g0 % g1; g0 = g1; g1 = t; }
    if (g0 != 1) return 0;
    *num = a; *den = b;
    return 1;
}

static i64 crt2(i64 a1, i64 m1, i64 a2, i64 m2, i64 *M12)
{
    /* returns x mod m1*m2 */
    i64 m = m1 * m2;
    /* inverse of m1 mod m2 */
    i64 t = 0, nt = 1, r = m2, nr = m1 % m2;
    while (nr) { i64 q = r / nr; i64 tmp = t - q * nt; t = nt; nt = tmp; tmp = r - q * nr; r = nr; nr = tmp; }
    if (t < 0) t += m2;
    i64 x = a1 + m1 * (((a2 - a1) % m2 + m2) % m2 % m2 * t % m2);
    x %= m; if (x < 0) x += m;
    *M12 = m;
    return x;
}

/* ------------------------------------------------------------------- modes */
static void do_selftest(void)
{
    u32 primes[] = { 101, 251, 1009, 10007 };
    for (int t = 0; t < 4; t++) {
        Mod M; mod_init(&M, primes[t]);
        u32 f6[6];
        /* 19044.h.2 : x(x-1)(x^2-3x+8)(x^2+4x+27) */
        u32 q4[5];
        u32 p1 = M.p - 3, p2 = 8 % M.p, p3 = 4 % M.p, p4 = 27 % M.p;
        build_quartic(&M, p1, p2, q4);
        build_f(&M, q4, p3, p4, f6);
        int o1 = cford_fp(&M, f6, 11);
        /* C4corr : (p1,p2,p3,p4) = (-34/15, 32/45, 31/15, -2/9) */
        u32 i15 = M.invtab[15 % M.p], i45 = M.invtab[45 % M.p], i9 = M.invtab[9 % M.p];
        u32 c1 = mmul(&M, M.p - (34 % M.p), i15), c2 = mmul(&M, 32 % M.p, i45);
        u32 c3 = mmul(&M, 31 % M.p, i15), c4 = mmul(&M, M.p - (2 % M.p), i9);
        build_quartic(&M, c1, c2, q4);
        build_f(&M, q4, c3, c4, f6);
        int o2 = cford_fp(&M, f6, 11);
        printf("SELFTEST p=%u  19044.h.2 -> %d   C4corr -> %d   (both expect 11)\n",
               M.p, o1, o2);
        free(M.invtab);
    }
    /* random-density probe: how often does a random monic sextic in the chart
       have CF order exactly 11 over F_p ?  (should be ~ 1/p^2 * O(1)) */
    for (int t = 0; t < 2; t++) {
        u32 p = (t == 0) ? 251 : 1009;
        Mod M; mod_init(&M, p);
        u64 tested = 0, hits = 0;
        u64 st = 123456789;
        for (u64 it = 0; it < 4000000ULL; it++) {
            st ^= st << 13; st ^= st >> 7; st ^= st << 17;
            u32 a = st % p; st ^= st << 13; st ^= st >> 7; st ^= st << 17;
            u32 b = st % p; st ^= st << 13; st ^= st >> 7; st ^= st << 17;
            u32 c = st % p; st ^= st << 13; st ^= st >> 7; st ^= st << 17;
            u32 d = st % p;
            if (!nondegenerate(&M, a, b, c, d)) continue;
            u32 q4[5], f6[6];
            build_quartic(&M, a, b, q4);
            build_f(&M, q4, c, d, f6);
            tested++;
            if (cford_fp(&M, f6, 11) == 11) hits++;
        }
        printf("DENSITY p=%u  tested=%llu  ord11=%llu  rate=%.3e   1/p^2=%.3e\n",
               p, (unsigned long long)tested, (unsigned long long)hits,
               (double)hits / (double)tested, 1.0 / ((double)p * p));
        free(M.invtab);
    }
}

static void do_surface(u32 p, char *dumpfile)
{
    Mod M; mod_init(&M, p);
    u64 total = 0;
    FILE *fp = dumpfile ? fopen(dumpfile, "w") : NULL;
    double t0 = omp_get_wtime();
#pragma omp parallel reduction(+:total)
    {
        char *buf = (char *)malloc(1 << 20); size_t bl = 0;
#pragma omp for schedule(dynamic, 1)
        for (int a = 0; a < (int)p; a++) {
            u32 q4[5], f6[6];
            for (u32 b = 0; b < p; b++) {
                build_quartic(&M, (u32)a, b, q4);
                for (u32 c = 0; c < p; c++)
                    for (u32 d = 0; d < p; d++) {
                        if (!nondegenerate(&M, (u32)a, b, c, d)) continue;
                        build_f(&M, q4, c, d, f6);
                        if (cford_fp(&M, f6, 11) == 11) {
                            total++;
                            if (fp) {
                                bl += sprintf(buf + bl, "%d %u %u %u\n", a, b, c, d);
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
        free(buf);
    }
    if (fp) fclose(fp);
    printf("SURFACE p=%u  |Sigma''(F_p)| = %llu   p^2 = %llu   ratio = %.4f   (%.1f s)\n",
           p, (unsigned long long)total, (unsigned long long)p * p,
           (double)total / ((double)p * p), omp_get_wtime() - t0);
    fflush(stdout);
    free(M.invtab);
}

static int fiber_of(Mod *M, u32 p1, u32 p2, u32 *out, int maxout)
{
    u32 q4[5], f6[6];
    int n = 0;
    build_quartic(M, p1, p2, q4);
    for (u32 c = 0; c < M->p; c++)
        for (u32 d = 0; d < M->p; d++) {
            if (!nondegenerate(M, p1, p2, c, d)) continue;
            build_f(M, q4, c, d, f6);
            if (cford_fp(M, f6, 11) == 11) {
                if (n < maxout) { out[2 * n] = c; out[2 * n + 1] = d; }
                n++;
                if (n >= maxout) return n;
            }
        }
    return n;
}

static void do_fiber(u32 p, i64 n1, i64 d1, i64 n2, i64 d2)
{
    Mod M; mod_init(&M, p);
    i64 a = ((n1 % p) + p) % p, b = ((d1 % p) + p) % p;
    i64 c = ((n2 % p) + p) % p, d = ((d2 % p) + p) % p;
    if (b == 0 || d == 0) { printf("FIBER p=%u  bad reduction of (p1,p2)\n", p); return; }
    u32 P1 = mmul(&M, (u32)a, M.invtab[b]), P2 = mmul(&M, (u32)c, M.invtab[d]);
    u32 out[4096];
    int n = fiber_of(&M, P1, P2, out, 2048);
    printf("FIBER p=%u  (p1,p2)=(%lld/%lld,%lld/%lld)  #fiber = %d\n",
           p, (long long)n1, (long long)d1, (long long)n2, (long long)d2, n);
    for (int i = 0; i < n && i < 2048; i++)
        printf("   (p3,p4) = (%u, %u)\n", out[2 * i], out[2 * i + 1]);
    free(M.invtab);
}

/* --------------------------------------------------------------- main sweep */
typedef struct { i64 n, d; } Rat;

static int build_rats(int H, Rat *out, int cap)
{
    int n = 0;
    for (i64 den = 1; den <= H; den++)
        for (i64 num = -H; num <= H; num++) {
            i64 a = num < 0 ? -num : num, b = den;
            while (b) { i64 t = a % b; a = b; b = t; }
            if (a != 1 && !(num == 0 && den == 1)) continue;
            if (num == 0 && den != 1) continue;
            if (n < cap) { out[n].n = num; out[n].d = den; }
            n++;
        }
    return n;
}

#define MAXPR 6
#define FIBCAP 512

static void do_sweep(int H, int npr, u32 *prs, int lo, int hi, i64 recbound)
{
    Rat *rats = (Rat *)malloc(sizeof(Rat) * 4 * (H + 1) * (H + 1));
    int nr = build_rats(H, rats, 4 * (H + 1) * (H + 1));
    if (hi > nr || hi <= 0) hi = nr;
    printf("SWEEP H=%d  #rationals=%d  nprimes=%d  shard [%d,%d)  recbound=%lld\n",
           H, nr, npr, lo, hi, (long long)recbound);
    fflush(stdout);

    Mod M[MAXPR];
    for (int t = 0; t < npr; t++) mod_init(&M[t], prs[t]);
    u64 done = 0, cand = 0, skipped = 0;
    double t0 = omp_get_wtime();
    int reported = 0;

#pragma omp parallel reduction(+:done,cand,skipped)
    {
        u32 *fib = (u32 *)malloc(sizeof(u32) * 2 * FIBCAP * MAXPR);
        int nf[MAXPR];
#pragma omp for schedule(dynamic, 1)
        for (int i = lo; i < hi; i++) {
            for (int j = 0; j < nr; j++) {
                i64 n1 = rats[i].n, d1 = rats[i].d, n2 = rats[j].n, d2 = rats[j].d;
                if (n2 == 0) continue;                      /* p2 != 0 */
                int bad = 0;
                for (int t = 0; t < npr; t++)
                    if (d1 % prs[t] == 0 || d2 % prs[t] == 0) bad = 1;
                if (bad) { skipped++; continue; }
                done++;
                int ok = 1;
                for (int t = 0; t < npr && ok; t++) {
                    u32 p = prs[t];
                    u32 P1 = mmul(&M[t], (u32)(((n1 % p) + p) % p), M[t].invtab[d1 % p]);
                    u32 P2 = mmul(&M[t], (u32)(((n2 % p) + p) % p), M[t].invtab[d2 % p]);
                    nf[t] = fiber_of(&M[t], P1, P2, fib + 2 * FIBCAP * t, FIBCAP);
                    if (nf[t] == 0) ok = 0;
                    if (nf[t] >= FIBCAP) { skipped++; ok = 0; }
                }
                if (!ok) continue;

                /* iterate over all npr-tuples of fiber points */
                int idx[MAXPR]; for (int t = 0; t < npr; t++) idx[t] = 0;
                while (1) {
                    i64 mm = prs[0], x = fib[2 * idx[0]], y = fib[2 * idx[0] + 1];
                    for (int t = 1; t < npr; t++) {
                        i64 nm;
                        x = crt2(x, mm, fib[2 * FIBCAP * t + 2 * idx[t]], prs[t], &nm);
                        y = crt2(y, mm, fib[2 * FIBCAP * t + 2 * idx[t] + 1], prs[t], &nm);
                        mm = nm;
                    }
                    i64 n3, d3, n4, d4;
                    if (ratrecon(x, mm, recbound, &n3, &d3)
                        && ratrecon(y, mm, recbound, &n4, &d4)) {
#pragma omp critical
                        {
                            printf("CAND p1=%lld/%lld p2=%lld/%lld p3=%lld/%lld p4=%lld/%lld\n",
                                   (long long)n1, (long long)d1, (long long)n2, (long long)d2,
                                   (long long)n3, (long long)d3, (long long)n4, (long long)d4);
                            fflush(stdout);
                        }
                        cand++;
                    }
                    int t = 0;
                    while (t < npr && ++idx[t] >= nf[t]) { idx[t] = 0; t++; }
                    if (t == npr) break;
                }
            }
#pragma omp critical
            {
                reported++;
                if (reported % 25 == 0) {
                    printf("PROGRESS i=%d/%d cand=%llu  %.1f s\n",
                           reported, hi - lo, (unsigned long long)cand, omp_get_wtime() - t0);
                    fflush(stdout);
                }
            }
        }
        free(fib);
    }
    printf("SWEEP_DONE H=%d shard [%d,%d) pairs=%llu skipped=%llu candidates=%llu  %.1f s\n",
           H, lo, hi, (unsigned long long)done, (unsigned long long)skipped,
           (unsigned long long)cand, omp_get_wtime() - t0);
    for (int t = 0; t < npr; t++) free(M[t].invtab);
    free(rats);
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s selftest|surface|fiber|sweep ...\n", argv[0]); return 1; }
    if (!strcmp(argv[1], "selftest")) do_selftest();
    else if (!strcmp(argv[1], "surface")) do_surface((u32)atoi(argv[2]), argc > 3 ? argv[3] : NULL);
    else if (!strcmp(argv[1], "fiber"))
        do_fiber((u32)atoi(argv[2]), atoll(argv[3]), atoll(argv[4]), atoll(argv[5]), atoll(argv[6]));
    else if (!strcmp(argv[1], "sweep")) {
        /* sweep H npr p1 p2 ... lo hi recbound */
        int H = atoi(argv[2]);
        int npr = atoi(argv[3]);
        u32 prs[MAXPR];
        for (int t = 0; t < npr; t++) prs[t] = (u32)atoi(argv[4 + t]);
        int lo = atoi(argv[4 + npr]), hi = atoi(argv[5 + npr]);
        i64 rb = atoll(argv[6 + npr]);
        do_sweep(H, npr, prs, lo, hi, rb);
    }
    else { fprintf(stderr, "unknown mode\n"); return 1; }
    return 0;
}
