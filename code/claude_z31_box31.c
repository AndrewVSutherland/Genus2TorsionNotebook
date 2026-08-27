/* Task B5: broad 31-divisibility box scan (zero-assumption backstop).
 *
 * Enumerate integral sextics  y^2 = f(x),  f = c6 x^6 + ... + c0  (c6 != 0),
 * in coefficient boxes and keep the curves with  31 | #J(F_p)  for ALL good
 * primes p among {37,41,43,47,53,59,61,67}  (necessary for a rational point
 * of order 31: J(Q)_tors injects into J(F_p) for every odd prime p of good
 * reduction).  Survivors go to an exact Magma check.
 *
 * Point counts (p odd, p !| c6, f squarefree mod p):
 *   #C(F_p)   = sum_{x in F_p} (1 + chi_p(f(x)))  +  (1 + chi_p(c6))
 *   #C(F_p^2) : an element t of F_p^2^* is a square  <=>  Norm(t) is a square
 *               in F_p (kernel of Norm has odd order p+1).  With
 *               F_p^2 = F_p(w), w^2 = n (n a fixed non-residue),
 *               Norm(u+vw) = u^2 - n v^2.  Every element of F_p^* is a square
 *               in F_p^2 (its norm is a square), so x in F_p contributes
 *               1 point if f(x)=0 else 2; conjugate points x = a+bw, a-bw
 *               (b != 0) contribute equally, so only b = 1..(p-1)/2 is
 *               enumerated and doubled; infinity contributes 2 (c6 in F_p^*).
 *   #J(F_p)   = (#C(F_p)^2 + #C(F_p^2))/2 - p        [ = chi(1), chi the
 *               Frobenius charpoly; standard Newton-identity consequence ]
 *
 * Fast path (box modes): per "row" (c6..c1 fixed, c0 the innermost loop) and
 * per prime we precompute
 *   gv[x]  = f(x) - c0 mod p                       for x in F_p
 *   pk[i]  = u | (p - n*v^2) << 16                 for x = a+bw, b>0, where
 *            f(x) - c0 = u + v*w
 *   h(T)   = Res_x(f'(x), f(x) - c0 + T)  in Newton form (degree 5 in T, via
 *            6 Euclid resultants + divided differences); h(c0) = 0  <=>  p is
 *            a prime of bad reduction for this curve (p !| c6 is checked
 *            separately) -- this catches repeated roots over the closure,
 *            not just F_p-rational ones.
 * Then each c0 costs ~p(p+1)/2 add+mul+byte-lookups per prime, and the
 * quadratic character / mod-p reduction is one table byte:
 *   PT1[idx] (idx < 2p)    : points contributed by an F_p value idx mod p
 *   PT2[idx] (idx < 4p^2)  : points contributed by a b>0 point whose norm is
 *                            idx mod p  (idx = t*t + (p - n v^2), t = u+c0)
 * Primes are tested cheapest-first (37 first) with early abort on the first
 * good prime with 31 !| #J, so ~30/31 of curves cost only the p=37 stage.
 *
 * Symmetries deduped in `box` mode (documented; disable with -nodedup):
 *   y -> -y                    : trivial (same f), nothing to do.
 *   x -> -x                    : (c5,c3,c1) -> (-c5,-c3,-c1); canonical form
 *                                first nonzero of (c5,c3,c1) positive.
 *   scaling x -> x/t, y -> y/t^3  (c_i -> t^{6-i} c_i), t = 2,3:
 *                                skip if t^{6-i} | c_i for i=5..0.
 *   scaling x -> tx,  y -> t^3 y  (c_i -> t^i c_i),     t = 2,3:
 *                                skip if t^i | c_i for i=1..6.
 *   content divisible by a square s^2 (4,9,25,49): y^2 = s^2 g  ~  y^2 = g.
 *   x-translation: NOT deduped by default (normalizing |c5| <= 3|c6| can
 *                                lose in-box-only models; enable with -t5).
 *   x -> 1/x reversal (c_i <-> c_{6-i}): NOT deduped (accepted duplicates).
 * Non-squarefree f over Q have disc = 0, hence bad reduction at every test
 * prime; they are dropped by the nt >= minnt rule.
 *
 * Survivor rule: a curve survives if it passes 31 | #J(F_p) at EVERY good
 * test prime AND at least minnt primes were good (default npr-2).  Without
 * the minnt floor, survivors are dominated (~20x) by curves whose disc is
 * divisible by test primes (each bad prime is a free pass); with 8 primes
 * and minnt=6 the survivor rate is ~31^-6 x C(8,2)/p^2-ish, a few per 1e11.
 * Coverage cost: curves with >= 3 bad primes among the 8 are lost --
 * accepted and documented (the RM witness has disc divisible by 61, so it
 * passes 7 of 8: nt=7).
 *
 * Default primes: {37,41,43,47,53,59,61,67} (-np N uses the first N).
 *
 * Modes
 *   selftest [nrand]        witness + gp cross-check data + internal checks
 *   check c6 c5 c4 c3 c2 c1 c0     detailed per-prime counts for one curve
 *   box  H shard nsh [-np N] [-t5] [-nodedup]     symmetric box |c_i| <= H
 *   boxr lo6 hi6 lo5 hi5 lo4 hi4 lo3 hi3 lo2 hi2 lo1 hi1 lo0 hi0 shard nsh
 *        [-np N]            explicit ranges, NO dedupe (validation / pins)
 *
 * Sharding: the (c6,c5,c4) task list is dealt round-robin, task_index mod
 * nsh == shard; shards are independent processes (one per box / restart
 * unit), OpenMP parallelizes inside a shard.
 *
 * Build: gcc -O3 -march=native -fopenmp -o claude_z31_box31 claude_z31_box31.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

typedef uint8_t  u8;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int32_t  i32;
typedef int64_t  i64;
typedef __uint128_t u128;

#define NPMAX 8
#define PMAX  67
#define ELL   31u

static u32 PR[NPMAX] = { 37, 41, 43, 47, 53, 59, 61, 67 };

/* ------------------------------------------------------------- modular ops */
typedef struct { u32 p; u64 minv; u32 invtab[PMAX]; } Mod;

static Mod MD[NPMAX];
static u8  SQF[NPMAX][PMAX];            /* nonzero-square flag              */
static u32 NQR[NPMAX];                  /* smallest quadratic non-residue   */
static u8  PT1[NPMAX][2 * PMAX];        /* F_p contribution of idx mod p    */
static u8  ISZ[NPMAX][2 * PMAX];        /* idx mod p == 0                   */
static u8  PT2[NPMAX][4 * PMAX * PMAX]; /* F_p^2 contribution, norm = idx mod p */
static u8  MRED[NPMAX][PMAX * PMAX];    /* idx mod p for idx < p^2 (build path) */

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
static inline u32 modi(i64 v, u32 p) { i64 r = v % (i64)p; if (r < 0) r += p; return (u32)r; }

static void global_init(void)
{
    for (int t = 0; t < NPMAX; t++) {
        u32 p = PR[t];
        MD[t].p = p;
        MD[t].minv = (u64)(((u128)1 << 64) / p);
        MD[t].invtab[0] = 0; MD[t].invtab[1] = 1;
        for (u32 i = 2; i < p; i++)
            MD[t].invtab[i] = (u32)((u64)(p - p / i) * (u64)MD[t].invtab[p % i] % p);
        memset(SQF[t], 0, sizeof(SQF[t]));
        for (u32 x = 1; x < p; x++) SQF[t][(u32)((u64)x * x % p)] = 1;
        u32 n = 2; while (SQF[t][n]) n++;
        NQR[t] = n;
        for (u32 idx = 0; idx < 2 * p; idx++) {
            u32 m = idx >= p ? idx - p : idx;
            PT1[t][idx] = m == 0 ? 1 : (SQF[t][m] ? 2 : 0);
            ISZ[t][idx] = (m == 0);
        }
        for (u32 idx = 0; idx < 4 * p * p; idx++) {
            u32 m = idx % p;
            PT2[t][idx] = m == 0 ? 1 : (SQF[t][m] ? 2 : 0);
        }
        for (u32 idx = 0; idx < PMAX * PMAX; idx++) MRED[t][idx] = (u8)(idx % p);
    }
}

/* --------------------------------------------- degree-6 poly helpers mod p */
/* remainder-degree of gcd(A,B); A,B given by coeff arrays + exact degrees */
static int polgcd_deg(Mod *M, u32 *A0, int dA, u32 *B0, int dB)
{
    u32 A[8], B[8];
    memcpy(A, A0, sizeof(u32) * (dA + 1));
    memcpy(B, B0, sizeof(u32) * (dB + 1));
    while (dB >= 0) {
        /* A mod B (in place on A), then swap */
        u32 il = M->invtab[B[dB]];
        for (int k = dA - dB; k >= 0; k--) {
            u32 q = mmul(M, A[k + dB], il);
            if (q) for (int j = 0; j <= dB; j++) A[k + j] = msub(M, A[k + j], mmul(M, q, B[j]));
        }
        int dR = -1;
        for (int i = dB - 1; i >= 0; i--) if (A[i]) { dR = i; break; }
        u32 T[8];
        memcpy(T, B, sizeof(u32) * (dB + 1));
        if (dR >= 0) memcpy(B, A, sizeof(u32) * (dR + 1));
        memcpy(A, T, sizeof(u32) * (dB + 1));
        int tmp = dB; dB = dR; dA = tmp;
    }
    return dA;
}

/* resultant of A (deg dA0) and B (deg dB0), both with exact degrees */
static u32 polres(Mod *M, u32 *A0, int dA0, u32 *B0, int dB0)
{
    u32 p = M->p;
    u32 A[8], B[8];
    int dA = dA0, dB = dB0, neg = 0;
    memcpy(A, A0, sizeof(u32) * (dA + 1));
    memcpy(B, B0, sizeof(u32) * (dB + 1));
    if (dA < dB) {
        u32 T[8];
        memcpy(T, A, sizeof(u32) * (dA + 1));
        memcpy(A, B, sizeof(u32) * (dB + 1));
        memcpy(B, T, sizeof(u32) * (dA + 1));
        int tmp = dA; dA = dB; dB = tmp;
        if ((dA & 1) && (dB & 1)) neg ^= 1;
    }
    u32 res = 1;
    while (1) {
        if (dB == 0) {                       /* Res(A, const c) = c^dA */
            u32 e = 1;
            for (int i = 0; i < dA; i++) e = mmul(M, e, B[0]);
            res = mmul(M, res, e);
            break;
        }
        /* R = A mod B */
        u32 il = M->invtab[B[dB]];
        for (int k = dA - dB; k >= 0; k--) {
            u32 q = mmul(M, A[k + dB], il);
            if (q) for (int j = 0; j <= dB; j++) A[k + j] = msub(M, A[k + j], mmul(M, q, B[j]));
        }
        int dR = -1;
        for (int i = dB - 1; i >= 0; i--) if (A[i]) { dR = i; break; }
        if (dR < 0) return 0;
        /* Res(A,B) = (-1)^(dA dB) lc(B)^(dA-dR) Res(B,R) */
        u32 e = 1;
        for (int i = 0; i < dA - dR; i++) e = mmul(M, e, B[dB]);
        res = mmul(M, res, e);
        if ((dA & 1) && (dB & 1)) neg ^= 1;
        u32 T[8];
        memcpy(T, B, sizeof(u32) * (dB + 1));
        memcpy(B, A, sizeof(u32) * (dR + 1));
        memcpy(A, T, sizeof(u32) * (dB + 1));
        dA = dB; dB = dR;
    }
    if (neg && res) res = p - res;
    return res;
}

/* --------------------------------------------------------------- fast rows */
typedef struct {
    u32 p, nq, c6p, inf1;
    u32 gv[PMAX];
    u32 npk;
    u32 pk[PMAX * (PMAX - 1) / 2];
    u32 hc[6];                 /* Newton coeffs of h(T) at nodes 0..5 */
    int cbad;                  /* p | c6 : prime unusable for this row */
} Row;

/* cd[i] = coefficient of x^i (degree order); c0 = cd[0] is NOT used here */
static void row_build(Row *R, int t, i64 *cd)
{
    Mod *M = &MD[t];
    u32 p = M->p;
    R->p = p; R->nq = NQR[t];
    u32 c[7];
    for (int i = 1; i <= 6; i++) c[i] = modi(cd[i], p);
    R->c6p = c[6];
    R->cbad = (c[6] == 0);
    if (R->cbad) return;
    R->inf1 = SQF[t][c[6]] ? 2 : 0;
    for (u32 x = 0; x < p; x++) {
        u32 acc = c[6];
        for (int i = 5; i >= 1; i--) acc = madd(M, mmul(M, acc, x), c[i]);
        R->gv[x] = mmul(M, acc, x);
    }
    /* 4 independent (u,v) Horner chains per iteration hide latency; all
       reductions go through the byte table MRED (operands < p, products
       < p^2, three-term sums < 3p < p^2), avoiding 128-bit Barrett here.
       Overshooting lanes (a0+l >= p) wrap a and compute duplicate values
       that are simply not stored.  n*b is hoisted per b. */
    u32 n = R->nq, k = 0;
    u8 *MR = MRED[t];
    for (u32 b = 1; b <= (p - 1) / 2; b++) {
        u32 nb = MR[n * b];
        for (u32 a0 = 0; a0 < p; a0 += 4) {
            u32 u[4], v[4], aa[4];
            for (int l = 0; l < 4; l++) {
                u32 a = a0 + (u32)l;
                aa[l] = a >= p ? a - p : a;
                u[l] = c[6]; v[l] = 0;
            }
            for (int i = 5; i >= 0; i--)
                for (int l = 0; l < 4; l++) {
                    u32 nu = (u32)MR[u[l] * aa[l]] + (u32)MR[v[l] * nb] + (i >= 1 ? c[i] : 0);
                    u32 nv = (u32)MR[u[l] * b] + (u32)MR[v[l] * aa[l]];
                    u[l] = MR[nu];
                    v[l] = MR[nv];
                }
            for (int l = 0; l < 4 && a0 + (u32)l < p; l++) {
                u32 s = MR[(u32)MR[v[l] * v[l]] * n];
                R->pk[k++] = u[l] | ((p - s) << 16);
            }
        }
    }
    R->npk = k;
    /* h(T) = Res_x(f'(x), g(x)+T), g = f - c0 ; degree 5 in T, nodes T=0..5 */
    u32 A[6];
    for (int i = 0; i <= 5; i++) A[i] = (u32)((u64)(i + 1) * c[i + 1] % p);
    u32 d[6];
    for (u32 j = 0; j < 6; j++) {
        u32 Bp[7];
        Bp[0] = j % p;
        for (int i = 1; i <= 6; i++) Bp[i] = c[i];
        d[j] = polres(M, A, 5, Bp, 6);
    }
    for (u32 kk = 1; kk <= 5; kk++)
        for (u32 j = 5; j >= kk; j--)
            d[j] = mmul(M, msub(M, d[j], d[j - 1]), M->invtab[kk]);
    for (int i = 0; i < 6; i++) R->hc[i] = d[i];
}

/* 0 = bad reduction at p ; 1 = good, counts returned */
static inline int row_test_full(Row *R, int t, u32 c0p, u32 *pn1, u32 *pn2, u32 *pJ)
{
    Mod *M = &MD[t];
    u32 p = R->p;
    u32 hv = R->hc[5];
    for (int k = 4; k >= 0; k--) hv = madd(M, mmul(M, hv, msub(M, c0p, (u32)k)), R->hc[k]);
    if (hv == 0) return 0;
    u8 *P1 = PT1[t], *IZ = ISZ[t], *P2 = PT2[t];
    u32 acc = 0, z = 0;
    for (u32 x = 0; x < p; x++) { u32 v = R->gv[x] + c0p; acc += P1[v]; z += IZ[v]; }
    u32 n1 = acc + R->inf1;
    u32 acc2 = 0, acc3 = 0;
    u32 *pk = R->pk;
    u32 npk = R->npk;
    for (u32 i = 0; i + 1 < npk; i += 2) {
        u32 pv = pk[i], pw = pk[i + 1];
        u32 tt = (pv & 0xffffu) + c0p, tu = (pw & 0xffffu) + c0p;
        acc2 += P2[tt * tt + (pv >> 16)];
        acc3 += P2[tu * tu + (pw >> 16)];
    }
    if (npk & 1) {
        u32 pv = pk[npk - 1];
        u32 tt = (pv & 0xffffu) + c0p;
        acc2 += P2[tt * tt + (pv >> 16)];
    }
    acc2 += acc3;
    u32 n2 = 2 * p - z + 2 * acc2 + 2;
    *pn1 = n1; *pn2 = n2;
    *pJ = (n1 * n1 + n2) / 2 - p;
    return 1;
}

/* 2 = bad reduction ; 1 = 31 | #J ; 0 = 31 !| #J */
static inline int row_test(Row *R, int t, u32 c0p)
{
    u32 n1, n2, J;
    if (!row_test_full(R, t, c0p, &n1, &n2, &J)) return 2;
    return (J % ELL) == 0;
}

/* ------------------------------------------------- direct (slow) reference */
/* independent brute-force counter; 1 = good reduction (counts filled), 0 = bad */
static int direct_prime(int t, i64 *cd, u32 *pn1, u32 *pn2, u32 *pJ)
{
    Mod *M = &MD[t];
    u32 p = M->p;
    u32 f[7], fd[6];
    for (int i = 0; i <= 6; i++) f[i] = modi(cd[i], p);
    if (!f[6]) return 0;
    for (int i = 0; i <= 5; i++) fd[i] = (u32)((u64)(i + 1) * f[i + 1] % p);
    if (polgcd_deg(M, f, 6, fd, 5) > 0) return 0;
    u32 n1 = SQF[t][f[6]] ? 2 : 0;
    for (u32 x = 0; x < p; x++) {
        u32 v = f[6];
        for (int i = 5; i >= 0; i--) v = madd(M, mmul(M, v, x), f[i]);
        n1 += v == 0 ? 1 : (SQF[t][v] ? 2 : 0);
    }
    u32 n = NQR[t], n2 = 2;
    for (u32 a = 0; a < p; a++)
        for (u32 b = 0; b < p; b++) {
            u32 u = f[6], v = 0;
            for (int i = 5; i >= 0; i--) {
                u32 nu = madd(M, madd(M, mmul(M, u, a), mmul(M, n, mmul(M, v, b))), f[i]);
                u32 nv = madd(M, mmul(M, u, b), mmul(M, v, a));
                u = nu; v = nv;
            }
            if (b == 0) { n2 += u == 0 ? 1 : 2; continue; }
            u32 nm = msub(M, mmul(M, u, u), mmul(M, n, mmul(M, v, v)));
            n2 += nm == 0 ? 1 : (SQF[t][nm] ? 2 : 0);
        }
    *pn1 = n1; *pn2 = n2;
    *pJ = (n1 * n1 + n2) / 2 - p;
    return 1;
}

/* ------------------------------------------------------------------- check */
static int check_curve(i64 *cf /* c6..c0 input order */, int npr, int verbose)
{
    i64 cd[7];
    for (int i = 0; i < 7; i++) cd[i] = cf[6 - i];
    int ngood = 0, npass = 0;
    for (int t = 0; t < npr; t++) {
        u32 n1, n2, J;
        int g = direct_prime(t, cd, &n1, &n2, &J);
        if (!g) {
            if (verbose) printf("CHECKP p=%u BAD_REDUCTION\n", PR[t]);
            continue;
        }
        ngood++;
        if (J % ELL == 0) npass++;
        if (verbose)
            printf("CHECKP p=%u N1=%u N2=%u J=%u Jmod31=%u %s\n",
                   PR[t], n1, n2, J, J % ELL, (J % ELL == 0) ? "PASS" : "fail");
    }
    if (verbose)
        printf("CHECK %lld %lld %lld %lld %lld %lld %lld ngood=%d npass=%d %s\n",
               (long long)cf[0], (long long)cf[1], (long long)cf[2], (long long)cf[3],
               (long long)cf[4], (long long)cf[5], (long long)cf[6], ngood, npass,
               (ngood >= 1 && npass == ngood) ? "SURVIVES" : "REJECTED");
    return ngood >= 1 && npass == ngood;
}

/* ---------------------------------------------------------------- selftest */
static u64 rst = 88172645463325252ULL;
static inline u64 xs64(void) { rst ^= rst << 13; rst ^= rst >> 7; rst ^= rst << 17; return rst; }
static inline i64 rndi(i64 lo, i64 hi) { return lo + (i64)(xs64() % (u64)(hi - lo + 1)); }

static void do_selftest(int nrand)
{
    int npr = 8;
    /* (a) the RM witness: y^2 = -3356x^6+11364x^5-18347x^4+17202x^3-9863x^2+3264x-504 */
    i64 W[7] = { -3356, 11364, -18347, 17202, -9863, 3264, -504 };
    printf("WITNESS check (expect 31 | #J at every good prime):\n");
    int wpass = check_curve(W, npr, 1);
    printf("WITNESS %s\n", wpass ? "PASS" : "FAIL");

    /* (b) random curves, good at all 6 primes, for the PARI cross-check */
    FILE *fp = fopen("results/claude_z31_box31_chkdata.gp", "w");
    if (!fp) { fprintf(stderr, "cannot open chkdata file (run from repo root)\n"); exit(1); }
    fprintf(fp, "chkdata = List();\n");
    int made = 0;
    while (made < nrand) {
        i64 cd[7];
        for (int i = 0; i < 7; i++) cd[i] = rndi(-50, 50);
        if (cd[6] == 0) continue;
        u32 n1[NPMAX], n2[NPMAX], J[NPMAX];
        int allgood = 1;
        for (int t = 0; t < npr && allgood; t++)
            if (!direct_prime(t, cd, &n1[t], &n2[t], &J[t])) allgood = 0;
        if (!allgood) continue;
        /* consistency: fast row path must agree with direct path */
        for (int t = 0; t < npr; t++) {
            Row R;
            row_build(&R, t, cd);
            u32 m1, m2, mJ;
            int g = R.cbad ? 0 : row_test_full(&R, t, modi(cd[0], R.p), &m1, &m2, &mJ);
            if (!g || m1 != n1[t] || m2 != n2[t] || mJ != J[t]) {
                printf("SELFTEST_FASTMISMATCH curve %d prime %u\n", made, PR[t]);
                exit(1);
            }
        }
        for (int t = 0; t < npr; t++) {
            fprintf(fp, "listput(chkdata, [%lld,%lld,%lld,%lld,%lld,%lld,%lld,%u,%u,%u,%u]);\n",
                    (long long)cd[6], (long long)cd[5], (long long)cd[4], (long long)cd[3],
                    (long long)cd[2], (long long)cd[1], (long long)cd[0],
                    PR[t], n1[t], n2[t], J[t]);
            printf("CHK %lld %lld %lld %lld %lld %lld %lld p=%u N1=%u N2=%u J=%u\n",
                   (long long)cd[6], (long long)cd[5], (long long)cd[4], (long long)cd[3],
                   (long long)cd[2], (long long)cd[1], (long long)cd[0],
                   PR[t], n1[t], n2[t], J[t]);
        }
        made++;
    }
    fclose(fp);
    printf("CHKDATA %d curves x %d primes written to results/claude_z31_box31_chkdata.gp\n",
           nrand, npr);

    /* (c) h(T)-resultant bad-reduction detector vs direct gcd, exhaustively in c0 */
    u64 htested = 0, hbad = 0, hmis = 0;
    for (int it = 0; it < 3000; it++) {
        int t = it % npr;
        u32 p = PR[t];
        i64 cd[7];
        cd[0] = 0;
        do { for (int i = 1; i <= 6; i++) cd[i] = rndi(-999, 999); } while (modi(cd[6], p) == 0);
        Row R;
        row_build(&R, t, cd);
        for (u32 c0 = 0; c0 < p; c0++) {
            Mod *M = &MD[t];
            u32 hv = R.hc[5];
            for (int k = 4; k >= 0; k--) hv = madd(M, mmul(M, hv, msub(M, c0, (u32)k)), R.hc[k]);
            u32 f[7], fd[6];
            for (int i = 1; i <= 6; i++) f[i] = modi(cd[i], p);
            f[0] = c0;
            for (int i = 0; i <= 5; i++) fd[i] = (u32)((u64)(i + 1) * f[i + 1] % p);
            int bad = polgcd_deg(M, f, 6, fd, 5) > 0;
            htested++;
            if (bad) hbad++;
            if (bad != (hv == 0)) hmis++;
        }
    }
    printf("HGCD tested=%llu bad=%llu mismatch=%llu %s\n",
           (unsigned long long)htested, (unsigned long long)hbad,
           (unsigned long long)hmis, hmis == 0 ? "PASS" : "FAIL");

    /* (d) fast path vs direct path on random curves incl. bad-reduction cases */
    u64 ftested = 0, fmis = 0, fgood = 0;
    for (int it = 0; it < 2000; it++) {
        i64 cd[7];
        do { for (int i = 0; i < 7; i++) cd[i] = rndi(-200, 200); } while (cd[6] == 0);
        for (int t = 0; t < npr; t++) {
            u32 n1, n2, J, m1, m2, mJ;
            int g1 = direct_prime(t, cd, &n1, &n2, &J);
            Row R;
            row_build(&R, t, cd);
            int g2 = R.cbad ? 0 : row_test_full(&R, t, modi(cd[0], R.p), &m1, &m2, &mJ);
            ftested++;
            if (g1 != g2) { fmis++; continue; }
            if (g1) { fgood++; if (n1 != m1 || n2 != m2 || J != mJ) fmis++; }
        }
    }
    printf("FASTDIRECT tested=%llu good=%llu mismatch=%llu %s\n",
           (unsigned long long)ftested, (unsigned long long)fgood,
           (unsigned long long)fmis, fmis == 0 ? "PASS" : "FAIL");
    printf("SELFTEST_DONE\n");
}

/* ------------------------------------------------------------------ ubench */
static void do_ubench(void)
{
    i64 cd[7] = { -504, 3264, -9863, 17202, -18347, 11364, -3356 }; /* witness, degree order */
    Row R;
    for (int t = 0; t < 2; t++) {
        u32 p = PR[t];
        double t0 = omp_get_wtime();
        u32 sink = 0;
        int nb = 20000;
        for (int i = 0; i < nb; i++) { cd[1] = 3264 + (i & 7); row_build(&R, t, cd); sink += R.pk[5]; }
        double t1 = omp_get_wtime();
        u64 nt2 = 2000000;
        u32 n1 = 0, n2 = 0, J = 0, good = 0;
        for (u64 i = 0; i < nt2; i++) {
            u32 c0p = (u32)(i % p);
            good += row_test_full(&R, t, c0p, &n1, &n2, &J);
            sink += J;
        }
        double t2 = omp_get_wtime();
        printf("UBENCH p=%u build=%.2fus test=%.1fns (sink %u good %u)\n",
               p, (t1 - t0) / nb * 1e6, (t2 - t1) / nt2 * 1e9, sink, good);
    }
}

/* -------------------------------------------------------------- box engine */
typedef struct { i64 lo[7], hi[7]; int dedupe, t5, npr, minnt, shard, nsh; } BoxSpec;
typedef struct { i32 c6, c5, c4; } Task;

static u64 g_enum, g_skip, g_tested, g_surv, g_nt0, g_tasks_done;
static u64 g_testP[NPMAX], g_passP[NPMAX];

static inline i64 iabs64(i64 v) { return v < 0 ? -v : v; }
static i64 igcd(i64 a, i64 b) { a = iabs64(a); b = iabs64(b); while (b) { i64 t = a % b; a = b; b = t; } return a; }

static void do_box(BoxSpec *B)
{
    int npr = B->npr;
    /* ---- build the shard's task list over (c6, c5, c4) ---- */
    int cap = 1 << 16, ntask = 0;
    Task *tasks = (Task *)malloc(sizeof(Task) * cap);
    u64 tidx = 0;
    for (i64 c6 = B->lo[6]; c6 <= B->hi[6]; c6++) {
        if (c6 == 0) continue;
        i64 lo5 = (B->dedupe && B->lo[5] < 0) ? 0 : B->lo[5];
        for (i64 c5 = lo5; c5 <= B->hi[5]; c5++) {
            if (B->t5 && iabs64(c5) > 3 * iabs64(c6)) continue;
            for (i64 c4 = B->lo[4]; c4 <= B->hi[4]; c4++) {
                if ((i64)(tidx++ % B->nsh) != B->shard) continue;
                if (ntask == cap) { cap *= 2; tasks = (Task *)realloc(tasks, sizeof(Task) * cap); }
                tasks[ntask].c6 = (i32)c6; tasks[ntask].c5 = (i32)c5; tasks[ntask].c4 = (i32)c4;
                ntask++;
            }
        }
    }
    printf("BOX shard=%d/%d npr=%d minnt=%d dedupe=%d t5=%d tasks=%d threads=%d ranges",
           B->shard, B->nsh, npr, B->minnt, B->dedupe, B->t5, ntask, omp_get_max_threads());
    for (int i = 6; i >= 0; i--) printf(" [%lld,%lld]", (long long)B->lo[i], (long long)B->hi[i]);
    printf("\n");
    fflush(stdout);
    double t0 = omp_get_wtime();

#pragma omp parallel for schedule(dynamic, 1)
    for (int ti = 0; ti < ntask; ti++) {
        i64 cd[7];
        cd[6] = tasks[ti].c6; cd[5] = tasks[ti].c5; cd[4] = tasks[ti].c4;
        u64 l_enum = 0, l_skip = 0, l_tested = 0, l_surv = 0, l_nt0 = 0;
        u64 l_testP[NPMAX] = { 0 }, l_passP[NPMAX] = { 0 };
        Row rows[NPMAX];
        i64 lo3 = (B->dedupe && cd[5] == 0 && B->lo[3] < 0) ? 0 : B->lo[3];
        for (i64 c3 = lo3; c3 <= B->hi[3]; c3++) {
            cd[3] = c3;
            for (i64 c2 = B->lo[2]; c2 <= B->hi[2]; c2++) {
                cd[2] = c2;
                i64 lo1 = (B->dedupe && cd[5] == 0 && c3 == 0 && B->lo[1] < 0) ? 0 : B->lo[1];
                for (i64 c1 = lo1; c1 <= B->hi[1]; c1++) {
                    cd[1] = c1;
                    i64 nrow = B->hi[0] - B->lo[0] + 1;
                    int A2 = 0, A3 = 0, sq4 = 0, sq9 = 0, sq25 = 0, sq49 = 0;
                    if (B->dedupe) {
                        /* whole row is an x->tx scaling of a smaller tuple */
                        if ((c1 % 2 == 0) && (c2 % 4 == 0) && (c3 % 8 == 0) &&
                            (cd[4] % 16 == 0) && (cd[5] % 32 == 0) && (cd[6] % 64 == 0))
                            { l_enum += nrow; l_skip += nrow; continue; }
                        if ((c1 % 3 == 0) && (c2 % 9 == 0) && (c3 % 27 == 0) &&
                            (cd[4] % 81 == 0) && (cd[5] % 243 == 0) && (cd[6] % 729 == 0))
                            { l_enum += nrow; l_skip += nrow; continue; }
                        /* x -> x/t scaling: c0-condition checked inside the loop */
                        A2 = (cd[5] % 2 == 0) && (cd[4] % 4 == 0) && (c3 % 8 == 0) &&
                             (c2 % 16 == 0) && (c1 % 32 == 0);
                        A3 = (cd[5] % 3 == 0) && (cd[4] % 9 == 0) && (c3 % 27 == 0) &&
                             (c2 % 81 == 0) && (c1 % 243 == 0);
                        i64 g = igcd(igcd(igcd(cd[6], cd[5]), igcd(cd[4], c3)), igcd(c2, c1));
                        sq4 = (g % 4 == 0); sq9 = (g % 9 == 0);
                        sq25 = (g % 25 == 0); sq49 = (g % 49 == 0);
                    }
                    int built[NPMAX] = { 0 };
                    for (i64 c0 = B->lo[0]; c0 <= B->hi[0]; c0++) {
                        l_enum++;
                        if (B->dedupe &&
                            ((A2 && c0 % 64 == 0) || (A3 && c0 % 729 == 0) ||
                             (sq4 && c0 % 4 == 0) || (sq9 && c0 % 9 == 0) ||
                             (sq25 && c0 % 25 == 0) || (sq49 && c0 % 49 == 0)))
                            { l_skip++; continue; }
                        l_tested++;
                        int nt = 0, pass = 1;
                        for (int t = 0; t < npr; t++) {
                            if (!built[t]) { row_build(&rows[t], t, cd); built[t] = 1; }
                            if (rows[t].cbad) continue;
                            int r = row_test(&rows[t], t, modi(c0, PR[t]));
                            if (r == 2) continue;
                            l_testP[t]++;
                            if (r == 0) { pass = 0; break; }
                            l_passP[t]++; nt++;
                        }
                        if (!pass) continue;
                        if (nt < B->minnt) { l_nt0++; continue; }
                        l_surv++;
#pragma omp critical
                        {
                            printf("SURV %lld %lld %lld %lld %lld %lld %lld nt=%d\n",
                                   (long long)cd[6], (long long)cd[5], (long long)cd[4],
                                   (long long)cd[3], (long long)cd[2], (long long)cd[1],
                                   (long long)c0, nt);
                            fflush(stdout);
                        }
                    }
                }
            }
        }
        __atomic_fetch_add(&g_enum, l_enum, __ATOMIC_RELAXED);
        __atomic_fetch_add(&g_skip, l_skip, __ATOMIC_RELAXED);
        __atomic_fetch_add(&g_tested, l_tested, __ATOMIC_RELAXED);
        __atomic_fetch_add(&g_surv, l_surv, __ATOMIC_RELAXED);
        __atomic_fetch_add(&g_nt0, l_nt0, __ATOMIC_RELAXED);
        for (int t = 0; t < npr; t++) {
            __atomic_fetch_add(&g_testP[t], l_testP[t], __ATOMIC_RELAXED);
            __atomic_fetch_add(&g_passP[t], l_passP[t], __ATOMIC_RELAXED);
        }
        u64 dn = __atomic_add_fetch(&g_tasks_done, 1, __ATOMIC_RELAXED);
        if (dn % 32 == 0) {
#pragma omp critical
            {
                double el = omp_get_wtime() - t0;
                u64 ts = __atomic_load_n(&g_tested, __ATOMIC_RELAXED);
                printf("PROGRESS tasks=%llu/%d tested=%llu surv=%llu %.1fs eta=%.0fs\n",
                       (unsigned long long)dn, ntask, (unsigned long long)ts,
                       (unsigned long long)__atomic_load_n(&g_surv, __ATOMIC_RELAXED),
                       el, el * ((double)ntask / (double)dn - 1.0));
                fflush(stdout);
            }
        }
    }
    double el = omp_get_wtime() - t0;
    int nthreads = omp_get_max_threads();
    printf("BOX_DONE shard=%d/%d enum=%llu skip=%llu tested=%llu lowntdrop=%llu surv=%llu wall=%.2fs\n",
           B->shard, B->nsh, (unsigned long long)g_enum, (unsigned long long)g_skip,
           (unsigned long long)g_tested, (unsigned long long)g_nt0,
           (unsigned long long)g_surv, el);
    for (int t = 0; t < npr; t++)
        printf("PRIME p=%u tested=%llu passed=%llu rate=%.5f (1/31=%.5f)\n",
               PR[t], (unsigned long long)g_testP[t], (unsigned long long)g_passP[t],
               g_testP[t] ? (double)g_passP[t] / (double)g_testP[t] : 0.0, 1.0 / 31.0);
    printf("RATE tested_per_sec=%.3e per_core=%.3e (threads=%d)\n",
           (double)g_tested / el, (double)g_tested / el / nthreads, nthreads);
    fflush(stdout);
    free(tasks);
}

/* -------------------------------------------------------------------- main */
int main(int argc, char **argv)
{
    global_init();
    if (argc < 2) {
        fprintf(stderr, "usage: %s selftest [nrand] | check c6..c0 | box H shard nsh [-np N] [-t5] [-nodedup] | boxr lo6 hi6 .. lo0 hi0 shard nsh [-np N]\n", argv[0]);
        return 1;
    }
    if (!strcmp(argv[1], "selftest")) {
        do_selftest(argc > 2 ? atoi(argv[2]) : 20);
        return 0;
    }
    if (!strcmp(argv[1], "ubench")) { do_ubench(); return 0; }
    if (!strcmp(argv[1], "check")) {
        if (argc < 9) { fprintf(stderr, "check needs 7 coefficients c6..c0\n"); return 1; }
        i64 cf[7];
        for (int i = 0; i < 7; i++) cf[i] = atoll(argv[2 + i]);
        check_curve(cf, 8, 1);
        return 0;
    }
    BoxSpec B;
    memset(&B, 0, sizeof(B));
    B.npr = 8; B.minnt = -1; B.nsh = 1;
    int argp = 0;
    if (!strcmp(argv[1], "box")) {
        if (argc < 5) { fprintf(stderr, "box H shard nsh\n"); return 1; }
        i64 H = atoll(argv[2]);
        for (int i = 0; i < 7; i++) { B.lo[i] = -H; B.hi[i] = H; }
        B.shard = atoi(argv[3]); B.nsh = atoi(argv[4]);
        B.dedupe = 1;
        argp = 5;
    } else if (!strcmp(argv[1], "boxr")) {
        if (argc < 18) { fprintf(stderr, "boxr lo6 hi6 lo5 hi5 lo4 hi4 lo3 hi3 lo2 hi2 lo1 hi1 lo0 hi0 shard nsh\n"); return 1; }
        for (int i = 0; i < 7; i++) { B.lo[6 - i] = atoll(argv[2 + 2 * i]); B.hi[6 - i] = atoll(argv[3 + 2 * i]); }
        B.shard = atoi(argv[16]); B.nsh = atoi(argv[17]);
        B.dedupe = 0;
        argp = 18;
    } else {
        fprintf(stderr, "unknown mode %s\n", argv[1]);
        return 1;
    }
    for (int i = argp; i < argc; i++) {
        if (!strcmp(argv[i], "-np") && i + 1 < argc) B.npr = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-minnt") && i + 1 < argc) B.minnt = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-t5")) B.t5 = 1;
        else if (!strcmp(argv[i], "-nodedup")) B.dedupe = 0;
        else { fprintf(stderr, "unknown flag %s\n", argv[i]); return 1; }
    }
    if (B.npr < 1 || B.npr > NPMAX) { fprintf(stderr, "bad -np\n"); return 1; }
    if (B.minnt < 0) { B.minnt = B.npr - 2; if (B.minnt < 1) B.minnt = 1; }
    if (B.shard < 0 || B.nsh < 1 || B.shard >= B.nsh) { fprintf(stderr, "bad shard\n"); return 1; }
    do_box(&B);
    return 0;
}
