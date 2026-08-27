/* claude_z31_kummer_sieve.c — chart U' universal Kummer-ladder sieve for
 * rational 31-torsion on genus-2 Jacobians (task B2 of the Z31 hunt).
 *
 * Chart U': parameters (u1,u0,v1,v0,w4,w3,w2,w1,w0) with
 *   u = x^2 + u1 x + u0 (monic),  v = v1 x + v0,
 *   w = w4 x^4 + w3 x^3 + w2 x^2 + w1 x + w0,
 *   f = v^2 + u*w   (deg 6, leading coeff w4, any sign),
 * marked class D = [div(u,v) - D_inf].  Every (genus-2 curve /Q, rational
 * divisor class) pair appears in this chart.  Sieve condition: 31*D = 0.
 *
 * Kummer embedding (Flynn / Cassels-Flynn convention, pinned by B1):
 *   kappa(D) = (1 : -u1 : u0 : (F0 - 2 y1y2)/(u1^2 - 4 u0)),  worked
 *   projectively as (den : -u1 den : u0 den : F0 - 2 y1y2), den = u1^2-4u0.
 *   den == 0 mod P  =>  PASS-UNKNOWN (never resolved in C, forwarded).
 * kappa(0) = (0:0:0:1).  31*D = 0  <=>  kappa(31D) = (0:0:0:1) projectively
 * (kappa identifies +-, and D != 0 is automatic since k1 = den != 0).
 *
 * Ladder: Montgomery ladder on the Kummer surface; doubling via Flynn's
 * quartic forms delta_i, pseudo-addition via the biquadratic forms B_ij:
 * with base w = kappa(D) (w_1 != 0), and (x,y) = (kappa(mD), kappa((m+1)D)),
 *   kappa((2m+1)D)_i  ∝  2 w_1 B_i1(x,y) - B_11(x,y) w_i .
 * The delta_i / B_i1 / defining-equation tables in claude_z31_kummer_tables.h
 * are MECHANICALLY generated (code/claude_z31_gen_tables.py) from Flynn's
 * files https://people.maths.ox.ac.uk/flynn/genus2/kummer/{duplication,
 * defining.equations,biquadratic.forms}; validated against Magma Jacobian
 * arithmetic (results/claude_z31_kummer_probe.log) and B1's reference
 * vectors (data/claude_z31_vectors_chartU.tsv, .._kummer.tsv).
 *
 * Modes
 *   selftest chartU.tsv kummer.tsv   full validation vs B1 vectors
 *   polytest vecfile                 raw polynomial agreement vs Python ref
 *   fiber P u1n u1d u0n u0d v1n v1d v0n v0d w4n w4d w3n w3d w2n w2d
 *                                    enumerate (w1,w0) in F_P^2, print hits
 *   sweep H npr P1..Pnpr lo hi recbound  <14 slice ints as in fiber>  vA vB
 *                                    slice sweep + CRT + rational recon;
 *                                    vA,vB in 0..6 vary that slice param over
 *                                    rationals of height<=H (vA=-1: single
 *                                    slice, the end-to-end path)
 *
 * Build: gcc -O3 -march=native -fopenmp -o claude_z31_kummer_sieve \
 *            claude_z31_kummer_sieve.c
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
static inline u32 mred(Mod *M, i64 x) { i64 r = x % (i64)M->p; if (r < 0) r += M->p; return (u32)r; }

/* ------------------------------- generated Flynn tables + monomial evals -- */
#include "claude_z31_kummer_tables.h"

/* every prime used must be odd and > MAXCOEFF (=96): term coefficients are
 * used as raw residues, and inv2 must exist */
static void check_prime(u32 p)
{
    if (p <= MAXCOEFF || (p & 1) == 0) {
        fprintf(stderr, "prime %u invalid: need odd p > %d\n", p, MAXCOEFF);
        exit(1);
    }
}

static Term3 *D3TAB[5] = { DELTA1_T, DELTA2_T, DELTA3_T, DELTA4_T, KUMMEQN_T };
static int    D3N[5]   = { N_DELTA1, N_DELTA2, N_DELTA3, N_DELTA4, N_KUMMEQN };
static Term4 *B4TAB[4] = { BB11_T, BB21_T, BB31_T, BB41_T };
static int    B4N[4]   = { N_BB11, N_BB21, N_BB31, N_BB41 };

static inline u32 eval_t3(Mod *M, Term3 *T, int n, u32 *fm, u32 *km)
{
    u32 acc = 0;
    for (int j = 0; j < n; j++) {
        u32 t = mmul(M, fm[T[j].fi], km[T[j].ki]);
        int16_t c = T[j].c;
        if (c >= 0) acc = madd(M, acc, mmul(M, (u32)c, t));
        else        acc = msub(M, acc, mmul(M, (u32)(-c), t));
    }
    return acc;
}
static inline u32 eval_t4(Mod *M, Term4 *T, int n, u32 *fm, u32 *kx, u32 *ky)
{
    u32 acc = 0;
    for (int j = 0; j < n; j++) {
        u32 t = mmul(M, fm[T[j].fi], mmul(M, kx[T[j].ki], ky[T[j].li]));
        int16_t c = T[j].c;
        if (c >= 0) acc = madd(M, acc, mmul(M, (u32)c, t));
        else        acc = msub(M, acc, mmul(M, (u32)(-c), t));
    }
    return acc;
}

/* doubling: out = delta(x); safe for out == x */
static inline void kdbl(Mod *M, u32 *fm, u32 *x, u32 *out)
{
    u32 km[NKMON4];
    eval_kmon4(M, x, km);
    for (int i = 0; i < 4; i++) out[i] = eval_t3(M, D3TAB[i], D3N[i], fm, km);
}
/* pseudo-add of x,y with difference w (w[0] != 0 required); out != x,y,w */
static inline void kpadd1(Mod *M, u32 *fm, u32 *x, u32 *y, u32 *w, u32 *out)
{
    u32 kx[NKMON2], ky[NKMON2], B[4];
    eval_kmon2(M, x, kx);
    eval_kmon2(M, y, ky);
    for (int i = 0; i < 4; i++) B[i] = eval_t4(M, B4TAB[i], B4N[i], fm, kx, ky);
    u32 tw = madd(M, w[0], w[0]);
    for (int i = 0; i < 4; i++)
        out[i] = msub(M, mmul(M, tw, B[i]), mmul(M, B[0], w[i]));
}
static inline u32 keqn(Mod *M, u32 *fm, u32 *x)
{
    u32 km[NKMON4];
    eval_kmon4(M, x, km);
    return eval_t3(M, D3TAB[4], D3N[4], fm, km);
}

#define KZERO(k)  ((k)[0] == 0 && (k)[1] == 0 && (k)[2] == 0 && (k)[3] != 0)
#define KALLZ(k)  ((k)[0] == 0 && (k)[1] == 0 && (k)[2] == 0 && (k)[3] == 0)

/* ladder: out = kappa(n*D) from base = kappa(D); base[0] != 0 required.
 * returns 1 ok, 0 degenerate (all-zero encountered -> caller marks UNKNOWN) */
static int laddern(Mod *M, u32 *fm, u32 *base, u64 n, u32 *out)
{
    u32 R0[4], R1[4], T0[4];
    if (n == 1) { memcpy(out, base, 4 * sizeof(u32)); return 1; }
    memcpy(R0, base, sizeof(R0));
    kdbl(M, fm, base, R1);
    int hb = 63;
    while (!((n >> hb) & 1)) hb--;
    for (int b = hb - 1; b >= 0; b--) {
        kpadd1(M, fm, R0, R1, base, T0);
        if ((n >> b) & 1) {
            kdbl(M, fm, R1, R1);
            memcpy(R0, T0, sizeof(T0));
        } else {
            kdbl(M, fm, R0, R0);
            memcpy(R1, T0, sizeof(T0));
        }
        if (KALLZ(R0) || KALLZ(R1)) return 0;
    }
    memcpy(out, R0, 4 * sizeof(u32));
    return 1;
}

/* ------------------------------------------------- chart U' mod-p plumbing */
typedef struct { u32 u1, u0, v1, v0, w4, w3, w2; } Slice;

static inline void chart_f(Mod *M, Slice *S, u32 w1, u32 w0, u32 *f)
{
    f[6] = S->w4;
    f[5] = madd(M, S->w3, mmul(M, S->u1, S->w4));
    f[4] = madd(M, S->w2, madd(M, mmul(M, S->u1, S->w3), mmul(M, S->u0, S->w4)));
    f[3] = madd(M, w1, madd(M, mmul(M, S->u1, S->w2), mmul(M, S->u0, S->w3)));
    f[2] = madd(M, w0, madd(M, mmul(M, S->u1, w1),
                            madd(M, mmul(M, S->u0, S->w2), mmul(M, S->v1, S->v1))));
    u32 vv = mmul(M, S->v1, S->v0);
    f[1] = madd(M, mmul(M, S->u1, w0), madd(M, mmul(M, S->u0, w1), madd(M, vv, vv)));
    f[0] = madd(M, mmul(M, S->u0, w0), mmul(M, S->v0, S->v0));
}

/* kappa(D) projective rep; returns den = u1^2 - 4 u0 (0 => degenerate) */
static inline u32 chart_kappa(Mod *M, Slice *S, u32 *f, u32 *k)
{
    u32 s1 = msub(M, 0, S->u1), s2 = S->u0;
    u32 s2q = mmul(M, s2, s2), s2c = mmul(M, s2q, s2);
    u32 F0 = madd(M, f[0], f[0]);
    F0 = madd(M, F0, mmul(M, f[1], s1));
    F0 = madd(M, F0, madd(M, mmul(M, f[2], s2), mmul(M, f[2], s2)));
    F0 = madd(M, F0, mmul(M, f[3], mmul(M, s1, s2)));
    F0 = madd(M, F0, madd(M, mmul(M, f[4], s2q), mmul(M, f[4], s2q)));
    F0 = madd(M, F0, mmul(M, f[5], mmul(M, s1, s2q)));
    F0 = madd(M, F0, madd(M, mmul(M, f[6], s2c), mmul(M, f[6], s2c)));
    u32 y1y2 = madd(M, msub(M, mmul(M, mmul(M, S->v1, S->v1), S->u0),
                            mmul(M, mmul(M, S->v1, S->v0), S->u1)),
                    mmul(M, S->v0, S->v0));
    u32 den = msub(M, mmul(M, S->u1, S->u1), madd(M, madd(M, S->u0, S->u0),
                                                  madd(M, S->u0, S->u0)));
    k[0] = den;
    k[1] = mmul(M, msub(M, 0, S->u1), den);
    k[2] = mmul(M, S->u0, den);
    k[3] = msub(M, F0, madd(M, y1y2, y1y2));
    return den;
}

/* candidate test: 1 = 31D == 0 (HIT), 0 = no, 2 = PASS-UNKNOWN (degenerate) */
static inline int test31(Mod *M, Slice *S, u32 w1, u32 w0)
{
    u32 f[7], fm[NFMON], base[4], out[4];
    chart_f(M, S, w1, w0, f);
    u32 den = chart_kappa(M, S, f, base);
    if (den == 0) return 2;
    eval_fmons(M, f, fm);
    if (!laddern(M, fm, base, 31, out)) return 2;
    if (KZERO(out)) return 1;
    if (KALLZ(out)) return 2;
    return 0;
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

/* --------------------------------------------------------------- selftest */
static void do_polytest(char *path)
{
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    long nrows = 0, fails = 0;
    i64 v[25];
    while (1) {
        int got = 0;
        for (int i = 0; i < 25; i++)
            if (fscanf(fp, "%lld", (long long *)&v[i]) == 1) got++;
        if (got == 0) break;
        if (got != 25) { fprintf(stderr, "short row\n"); fails++; break; }
        nrows++;
        check_prime((u32)v[0]);
        Mod M; mod_init(&M, (u32)v[0]);
        u32 f[7], x[4], y[4], fm[NFMON], km[NKMON4], kx[NKMON2], ky[NKMON2];
        for (int i = 0; i < 7; i++) f[i] = mred(&M, v[1 + i]);
        for (int i = 0; i < 4; i++) x[i] = mred(&M, v[8 + i]);
        for (int i = 0; i < 4; i++) y[i] = mred(&M, v[12 + i]);
        eval_fmons(&M, f, fm);
        eval_kmon4(&M, x, km);
        eval_kmon2(&M, x, kx);
        eval_kmon2(&M, y, ky);
        for (int i = 0; i < 4; i++) {
            u32 d = eval_t3(&M, D3TAB[i], D3N[i], fm, km);
            if (d != (u32)v[16 + i]) { fails++; printf("POLYFAIL delta%d row %ld\n", i + 1, nrows); }
            u32 b = eval_t4(&M, B4TAB[i], B4N[i], fm, kx, ky);
            if (b != (u32)v[20 + i]) { fails++; printf("POLYFAIL B%d1 row %ld\n", i + 1, nrows); }
        }
        u32 e = eval_t3(&M, D3TAB[4], D3N[4], fm, km);
        if (e != (u32)v[24]) { fails++; printf("POLYFAIL kummeqn row %ld\n", nrows); }
        free(M.invtab);
    }
    fclose(fp);
    printf("POLYTEST rows=%ld fails=%ld %s\n", nrows, fails, fails ? "FAIL" : "PASS");
}

static void normalize4(Mod *M, u32 *k)
{
    for (int i = 0; i < 4; i++)
        if (k[i]) {
            u32 inv = M->invtab[k[i]];
            for (int j = 0; j < 4; j++) k[j] = mmul(M, k[j], inv);
            return;
        }
}

static void do_selftest(char *chartU_path, char *kummer_path)
{
    /* ---- part 1: kummer.tsv — kappa(D) and kappa(31D) vs Magma ---- */
    FILE *fp = fopen(kummer_path, "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", kummer_path); exit(1); }
    int krows = 0, kden0 = 0, kfail = 0, k31fail = 0;
    i64 v[20];
    while (1) {
        int got = 0;
        for (int i = 0; i < 20; i++)
            if (fscanf(fp, "%lld", (long long *)&v[i]) == 1) got++;
        if (got == 0) break;
        if (got != 20) { fprintf(stderr, "kummer short row\n"); kfail++; break; }
        krows++;
        check_prime((u32)v[0]);
        Mod M; mod_init(&M, (u32)v[0]);
        Slice S = { mred(&M, v[1]), mred(&M, v[2]), mred(&M, v[3]), mred(&M, v[4]),
                    mred(&M, v[5]), mred(&M, v[6]), mred(&M, v[7]) };
        u32 w1 = mred(&M, v[8]), w0 = mred(&M, v[9]);
        u32 f[7], fm[NFMON], base[4], out[4], ref[4];
        chart_f(&M, &S, w1, w0, f);
        u32 den = chart_kappa(&M, &S, f, base);
        if (den == 0) { kden0++; free(M.invtab); continue; }   /* PASS-UNKNOWN */
        eval_fmons(&M, f, fm);
        memcpy(out, base, sizeof(out));
        normalize4(&M, out);
        for (int i = 0; i < 4; i++) ref[i] = mred(&M, v[12 + i]);
        if (memcmp(out, ref, sizeof(ref))) {
            kfail++;
            printf("KFAIL row %d p=%lld got (%u:%u:%u:%u) want (%u:%u:%u:%u)\n",
                   krows, (long long)v[0], out[0], out[1], out[2], out[3],
                   ref[0], ref[1], ref[2], ref[3]);
        }
        if (!laddern(&M, fm, base, 31, out)) {
            k31fail++; printf("K31DEGEN row %d\n", krows);
        } else {
            normalize4(&M, out);
            for (int i = 0; i < 4; i++) ref[i] = mred(&M, v[16 + i]);
            if (memcmp(out, ref, sizeof(ref))) {
                k31fail++;
                printf("K31FAIL row %d p=%lld got (%u:%u:%u:%u) want (%u:%u:%u:%u)\n",
                       krows, (long long)v[0], out[0], out[1], out[2], out[3],
                       ref[0], ref[1], ref[2], ref[3]);
            }
        }
        free(M.invtab);
    }
    fclose(fp);
    printf("SELFTEST kummer: rows=%d den0_unknown=%d kappaD_fail=%d kappa31D_fail=%d\n",
           krows, kden0, kfail, k31fail);

    /* ---- part 2: chartU.tsv — exact-order ladder + div31 discrimination ---- */
    fp = fopen(chartU_path, "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", chartU_path); exit(1); }
    int rows = 0, unk = 0, eqfail = 0, ordfail = 0, d31fail = 0, s31pass = 0, agree = 0;
    while (1) {
        int got = 0;
        for (int i = 0; i < 12; i++)
            if (fscanf(fp, "%lld", (long long *)&v[i]) == 1) got++;
        if (got == 0) break;
        if (got != 12) { fprintf(stderr, "chartU short row\n"); ordfail++; break; }
        rows++;
        check_prime((u32)v[0]);
        Mod M; mod_init(&M, (u32)v[0]);
        Slice S = { mred(&M, v[1]), mred(&M, v[2]), mred(&M, v[3]), mred(&M, v[4]),
                    mred(&M, v[5]), mred(&M, v[6]), mred(&M, v[7]) };
        u32 w1 = mred(&M, v[8]), w0 = mred(&M, v[9]);
        u64 order = (u64)v[10];
        int flag = (int)v[11];
        u32 f[7], fm[NFMON], base[4], out[4];
        chart_f(&M, &S, w1, w0, f);
        u32 den = chart_kappa(&M, &S, f, base);
        if (den == 0) {
            unk++;
            printf("  PASS-UNKNOWN (den=0) row %d p=%lld\n", rows, (long long)v[0]);
            free(M.invtab); continue;
        }
        eval_fmons(&M, f, fm);
        if (keqn(&M, fm, base) != 0) { eqfail++; printf("  EQFAIL row %d\n", rows); }
        int rowok = 1;
        if (!laddern(&M, fm, base, order, out) || !KZERO(out)) {
            ordfail++; rowok = 0;
            printf("  ORDFAIL row %d p=%lld order=%llu\n", rows, (long long)v[0],
                   (unsigned long long)order);
        }
        if (flag) {   /* ladder must certify 31 | ord: (order/31)*D != 0 */
            if (!laddern(&M, fm, base, order / 31, out) || KZERO(out)) {
                d31fail++; rowok = 0;
                printf("  D31FAIL row %d p=%lld order=%llu\n", rows, (long long)v[0],
                       (unsigned long long)order);
            }
        }
        if (rowok) agree++;
        /* the production sieve predicate: expect PASS only if order in {1,31} */
        if (laddern(&M, fm, base, 31, out) && KZERO(out)) {
            s31pass++;
            if (order != 31 && order != 1)
                printf("  SIEVE31-SPURIOUS row %d order=%llu\n", rows,
                       (unsigned long long)order);
        }
        free(M.invtab);
    }
    fclose(fp);
    printf("SELFTEST chartU: rows=%d pass_unknown=%d eqfail=%d ordfail=%d "
           "d31fail=%d div31_agree=%d/%d sieve31_pass=%d\n",
           rows, unk, eqfail, ordfail, d31fail, agree, rows - unk, s31pass);
    int fails = kfail + k31fail + eqfail + ordfail + d31fail;
    printf("SELFTEST_DONE FAILURES=%d %s\n", fails, fails ? "FAIL" : "PASS");
}

/* ------------------------------------------------------------------ fiber */
#define FIBCAP 8192

typedef struct { u32 w1, w0, stat; } FibPt;

/* enumerate (w1,w0) in F_P^2; returns hit count (capped), hits in out[] */
static i64 fiber_collect(Mod *M, Slice *S, FibPt *out, int cap, i64 *n_unk)
{
    i64 n = 0, nu = 0;
#pragma omp parallel reduction(+:nu)
    {
        FibPt loc[256];
        int nl = 0;
#pragma omp for schedule(dynamic, 4)
        for (int w1 = 0; w1 < (int)M->p; w1++) {
            for (u32 w0 = 0; w0 < M->p; w0++) {
                int r = test31(M, S, (u32)w1, w0);
                if (r == 0) continue;
                if (r == 2) nu++;
                if (nl < 256) { loc[nl].w1 = (u32)w1; loc[nl].w0 = w0; loc[nl].stat = r; nl++; }
                if (nl == 256) {
#pragma omp critical
                    {
                        for (int i = 0; i < nl && n < cap; i++) out[n++] = loc[i];
                    }
                    nl = 0;
                }
            }
        }
#pragma omp critical
        {
            for (int i = 0; i < nl && i + n < cap; i++) out[n++] = loc[i];
        }
    }
    *n_unk = nu;
    return n;
}

static int reduce_slice(Mod *M, i64 *r14, Slice *S)
{
    u32 vals[7];
    for (int i = 0; i < 7; i++) {
        i64 num = r14[2 * i], den = r14[2 * i + 1];
        if (den % M->p == 0) return 0;
        vals[i] = mmul(M, mred(M, num), M->invtab[mred(M, den)]);
    }
    S->u1 = vals[0]; S->u0 = vals[1]; S->v1 = vals[2]; S->v0 = vals[3];
    S->w4 = vals[4]; S->w3 = vals[5]; S->w2 = vals[6];
    return 1;
}

static void do_fiber(u32 p, i64 *r14)
{
    check_prime(p);
    Mod M; mod_init(&M, p);
    Slice S;
    if (!reduce_slice(&M, r14, &S)) { printf("FIBER_BAD p=%u slice denominator divisible by p\n", p); return; }
    u32 den = msub(&M, mmul(&M, S.u1, S.u1),
                   madd(&M, madd(&M, S.u0, S.u0), madd(&M, S.u0, S.u0)));
    if (den == 0) { printf("FIBER_DEGENERATE p=%u u1^2-4u0=0 (whole fiber PASS-UNKNOWN)\n", p); return; }
    FibPt *out = (FibPt *)malloc(sizeof(FibPt) * FIBCAP);
    i64 nu = 0;
    double t0 = omp_get_wtime();
    i64 n = fiber_collect(&M, &S, out, FIBCAP, &nu);
    double dt = omp_get_wtime() - t0;
    i64 cand = (i64)p * p;
    int nth = omp_get_max_threads();
    for (i64 i = 0; i < n; i++)
        printf("HITF %u %u %u %s\n", p, out[i].w1, out[i].w0,
               out[i].stat == 2 ? "U" : "P");
    printf("FIBER_DONE p=%u hits=%lld unknown=%lld cand=%lld secs=%.2f "
           "rate=%.0f/s threads=%d rate_per_core=%.0f/s\n",
           p, (long long)n, (long long)nu, (long long)cand, dt,
           (double)cand / dt, nth, (double)cand / dt / nth);
    free(out); free(M.invtab);
}

/* ------------------------------------------------------------------ sweep */
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
#define TUPCAP 1000000LL

static void do_sweep(int H, int npr, u32 *prs, int lo, int hi, i64 recbound,
                     i64 *base14, int vA, int vB)
{
    Rat *rats = (Rat *)malloc(sizeof(Rat) * 4 * (H + 1) * (H + 1) + sizeof(Rat));
    int nr = build_rats(H, rats, 4 * (H + 1) * (H + 1));
    int nA = (vA >= 0) ? nr : 1, nB = (vB >= 0) ? nr : 1;
    if (hi > nA || hi <= 0) hi = nA;
    printf("SWEEP H=%d #rats=%d nprimes=%d shard[%d,%d)x%d recbound=%lld vary=(%d,%d)\n",
           H, nr, npr, lo, hi, nB, (long long)recbound, vA, vB);
    fflush(stdout);

    /* sort primes ascending so empty small fibers abort early */
    for (int i = 0; i < npr; i++)
        for (int j = i + 1; j < npr; j++)
            if (prs[j] < prs[i]) { u32 t = prs[i]; prs[i] = prs[j]; prs[j] = t; }
    Mod M[MAXPR];
    for (int t = 0; t < npr; t++) { check_prime(prs[t]); mod_init(&M[t], prs[t]); }
    FibPt *fib = (FibPt *)malloc(sizeof(FibPt) * FIBCAP * MAXPR);
    i64 nf[MAXPR];
    u64 slices = 0, skipped = 0, cands = 0;
    double t0 = omp_get_wtime();

    for (int ia = lo; ia < hi; ia++) {
        for (int ib = 0; ib < nB; ib++) {
            i64 r14[14];
            memcpy(r14, base14, sizeof(r14));
            if (vA >= 0) { r14[2 * vA] = rats[ia].n; r14[2 * vA + 1] = rats[ia].d; }
            if (vB >= 0) { r14[2 * vB] = rats[ib].n; r14[2 * vB + 1] = rats[ib].d; }
            slices++;
            int ok = 1;
            i64 ntup = 1;
            for (int t = 0; t < npr && ok; t++) {
                Slice S;
                if (!reduce_slice(&M[t], r14, &S)) { skipped++; ok = 0; break; }
                u32 den = msub(&M[t], mmul(&M[t], S.u1, S.u1),
                               madd(&M[t], madd(&M[t], S.u0, S.u0), madd(&M[t], S.u0, S.u0)));
                if (den == 0) { skipped++; ok = 0; break; }
                i64 nu = 0;
                nf[t] = fiber_collect(&M[t], &S, fib + FIBCAP * t, FIBCAP, &nu);
                if (nf[t] == 0) ok = 0;
                if (nf[t] >= FIBCAP) { skipped++; ok = 0; }
                ntup *= nf[t] ? nf[t] : 1;
                if (ntup > TUPCAP) { skipped++; ok = 0; }
            }
            if (!ok) continue;
            /* CRT over all tuples of fiber points, then rational reconstruction */
            int idx[MAXPR];
            for (int t = 0; t < npr; t++) idx[t] = 0;
            while (1) {
                i64 mm = prs[0];
                i64 x = fib[idx[0]].w1, y = fib[idx[0]].w0;
                for (int t = 1; t < npr; t++) {
                    i64 nm;
                    x = crt2(x, mm, fib[FIBCAP * t + idx[t]].w1, prs[t], &nm);
                    y = crt2(y, mm, fib[FIBCAP * t + idx[t]].w0, prs[t], &nm);
                    mm = nm;
                }
                i64 n3, d3, n4, d4;
                if (ratrecon(x, mm, recbound, &n3, &d3)
                    && ratrecon(y, mm, recbound, &n4, &d4)) {
                    printf("CAND u1=%lld/%lld u0=%lld/%lld v1=%lld/%lld v0=%lld/%lld "
                           "w4=%lld/%lld w3=%lld/%lld w2=%lld/%lld "
                           "w1=%lld/%lld w0=%lld/%lld\n",
                           (long long)r14[0], (long long)r14[1],
                           (long long)r14[2], (long long)r14[3],
                           (long long)r14[4], (long long)r14[5],
                           (long long)r14[6], (long long)r14[7],
                           (long long)r14[8], (long long)r14[9],
                           (long long)r14[10], (long long)r14[11],
                           (long long)r14[12], (long long)r14[13],
                           (long long)n3, (long long)d3, (long long)n4, (long long)d4);
                    fflush(stdout);
                    cands++;
                }
                int t = 0;
                while (t < npr && ++idx[t] >= nf[t]) { idx[t] = 0; t++; }
                if (t == npr) break;
            }
        }
        printf("PROGRESS ia=%d/%d slices=%llu cands=%llu %.1f s\n",
               ia - lo + 1, hi - lo, (unsigned long long)slices,
               (unsigned long long)cands, omp_get_wtime() - t0);
        fflush(stdout);
    }
    printf("SWEEP_DONE shard[%d,%d) slices=%llu skipped=%llu candidates=%llu %.1f s\n",
           lo, hi, (unsigned long long)slices, (unsigned long long)skipped,
           (unsigned long long)cands, omp_get_wtime() - t0);
    for (int t = 0; t < npr; t++) free(M[t].invtab);
    free(fib); free(rats);
}

/* ------------------------------------------------------------------- main */
int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s selftest|polytest|fiber|sweep ...\n", argv[0]);
        return 1;
    }
    if (!strcmp(argv[1], "selftest")) {
        char *cu = argc > 2 ? argv[2] : "data/claude_z31_vectors_chartU.tsv";
        char *ku = argc > 3 ? argv[3] : "data/claude_z31_vectors_kummer.tsv";
        do_selftest(cu, ku);
    } else if (!strcmp(argv[1], "polytest")) {
        do_polytest(argv[2]);
    } else if (!strcmp(argv[1], "fiber")) {
        if (argc < 17) { fprintf(stderr, "fiber P + 14 slice ints\n"); return 1; }
        i64 r14[14];
        for (int i = 0; i < 14; i++) r14[i] = atoll(argv[3 + i]);
        do_fiber((u32)atoi(argv[2]), r14);
    } else if (!strcmp(argv[1], "sweep")) {
        int H = atoi(argv[2]);
        int npr = atoi(argv[3]);
        if (npr > MAXPR) { fprintf(stderr, "npr > %d\n", MAXPR); return 1; }
        u32 prs[MAXPR];
        for (int t = 0; t < npr; t++) prs[t] = (u32)atoi(argv[4 + t]);
        int a = 4 + npr;
        int lo = atoi(argv[a]), hi = atoi(argv[a + 1]);
        i64 rb = atoll(argv[a + 2]);
        if (argc < a + 3 + 14 + 2) { fprintf(stderr, "sweep needs 14 slice ints + vA vB\n"); return 1; }
        i64 base14[14];
        for (int i = 0; i < 14; i++) base14[i] = atoll(argv[a + 3 + i]);
        int vA = atoi(argv[a + 17]), vB = atoi(argv[a + 18]);
        do_sweep(H, npr, prs, lo, hi, rb, base14, vA, vB);
    } else {
        fprintf(stderr, "unknown mode %s\n", argv[1]);
        return 1;
    }
    return 0;
}
