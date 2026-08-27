/* claude_z31_ek8_sieve.c — Task B4: Elkies-Kumar disc-8 (RM by sqrt 2) fiber
 * sieve for rational 31-torsion on genus-2 Jacobians.
 *
 * Chart: the EK rational (r,s)-plane for Y_-(8)  (arXiv:1209.3527 ancillary
 * data, directory 8/igusa8.txt; recorded in data/claude_z31_ek8_source/):
 *   I2  = -4(3s+8r-2)
 *   I4  =  4(9rs+4r^2+4r+1)
 *   I6  = -(16/3)(9rs+4r^2+4r+1)(3s+8r-2) + (4/3)(54r^2 s+81rs-16r^3-24r^2-12r-2)
 *   I10 = -8 r^3 s^2
 * Every (r,s) with I10 != 0 off the special-automorphism locus is the moduli
 * point of a genus-2 curve /Q whose Jacobian has RM by Z[sqrt2].
 *
 * Condition sieved: 31 | #J(F_P) for SOME quadratic twist of the fiber curve,
 * i.e. 31 | chi_P(1) or 31 | chi_P(-1) — twist-class invariant (the multiset
 * {chi(1), chi(-1)} does not depend on which twist Mestre reconstruction
 * happens to produce).  chi is obtained EXACTLY from N1 = #C(F_P),
 * N2 = #C(F_{P^2}):  s1 = P+1-N1, s2 = (N2-(P^2+1)+s1^2)/2,
 * chi(±1) = P^2+1 ± (P+1) s1 + s2.
 *
 * Per-point pipeline mod P (P ≡ 3 mod 4, P > 5, P != 31):
 *   (r,s) -> IC -> Clebsch (A,B,C,D) -> Mestre conic A_ij + cubic a_ijk
 *   (formulas as in Magma package Geometry/CrvG2/mestre.m; the van Wamelen
 *   U/DP rescalings only normalize the twist over Q and are omitted)
 *   -> deterministic conic point -> parametrize -> sextic f(m)
 *   -> N1, N2 by exhaustive char-sum counts (F_{P^2} = F_P[i], i^2 = -1).
 * DEGEN (I10 = 0, conic det = 0, no conic point, deg f < 5, f not
 * squarefree) is forwarded as wildcard-accept (mask 7), never dropped.
 *
 * Twist bookkeeping: local Mestre models carry NO canonical twist, so the
 * per-prime mask (bit0: 31|chi(1), bit1: 31|chi(-1)) is recorded but
 * cross-prime sign consistency is NOT enforceable in C (any sign pattern on
 * <= 6 primes is realizable by some quadratic character anyway); the global
 * twist d is recovered in the exact stage (code/claude_z31_ek8_exact.m).
 *
 * Modes
 *   selftest                        anchors vs Magma (RM witness fibers)
 *   ref FILE.tsv                    validate vs claude_z31_ek8_ref.m vectors
 *   fiber P rnum rden               accepted s-fiber over slice r mod P
 *   bench P [nrows]                 throughput measurement (single thread)
 *   sweep H npr P1..Pk nver V1..Vm lo hi recbound
 *                                   r-slices of height <= H, shard [lo,hi),
 *                                   CRT over the npr primes, ratrecon bound
 *                                   recbound, membership check at V1..Vm
 *
 * Build: gcc -O3 -march=native -fopenmp -o claude_z31_ek8_sieve claude_z31_ek8_sieve.c
 * The RM witness (1830.2.a.q, J(Q)_tors = Z/31) sits at
 *   r = 6000/96721, s = 3557520/96721        (96721 = 311^2).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

typedef uint32_t u32;
typedef uint64_t u64;
typedef int32_t  i32;
typedef int64_t  i64;
typedef int8_t   i8;
typedef __uint128_t u128;

/* ------------------------------------------------------------- modular ops */
typedef struct {
    u32 p; u64 minv; u32 *invtab; i8 *chitab;
    u32 inv2, inv3;
    u32 i120, i6750, i202500, i4556250;          /* inverses of the Clebsch dens */
    u32 i9, i27, i54, i162, i324;                /* inverses of the cubic dens   */
    u32 c720, c8640, c108000, c62208, c972000, c1620000, c3037500, c6075000;
} EK;

static inline u32 mmul(EK *M, u32 a, u32 b)
{
    u64 z = (u64)a * (u64)b;
    u64 q = (u64)(((u128)z * (u128)M->minv) >> 64);
    u64 r = z - q * (u64)M->p;
    if (r >= M->p) r -= M->p;
    return (u32)r;
}
static inline u32 madd(EK *M, u32 a, u32 b) { u32 s = a + b; return s >= M->p ? s - M->p : s; }
static inline u32 msub(EK *M, u32 a, u32 b) { return a >= b ? a - b : a + M->p - b; }
static inline u32 mneg(EK *M, u32 a) { return a ? M->p - a : 0; }

static u32 mpow(EK *M, u32 a, u32 e)
{
    u32 r = 1;
    while (e) { if (e & 1) r = mmul(M, r, a); a = mmul(M, a, a); e >>= 1; }
    return r;
}

/* sqrt mod p for p = 3 mod 4; returns 1 and *rt on success, 0 if non-residue */
static int msqrt(EK *M, u32 a, u32 *rt)
{
    if (a == 0) { *rt = 0; return 1; }
    u32 s = mpow(M, a, (M->p + 1) / 4);
    if (mmul(M, s, s) != a) return 0;
    *rt = s;
    return 1;
}

static void ek_init(EK *M, u32 p)
{
    if (p % 4 != 3 || p <= 5 || p == 31) {
        fprintf(stderr, "prime %u not allowed (need p = 3 mod 4, p > 5, p != 31)\n", p);
        exit(1);
    }
    M->p = p;
    M->minv = (u64)(((u128)1 << 64) / p);
    M->invtab = (u32 *)malloc(sizeof(u32) * p);
    M->invtab[0] = 0;
    if (p > 1) M->invtab[1] = 1;
    for (u32 i = 2; i < p; i++)
        M->invtab[i] = (u32)((u64)(p - p / i) * (u64)M->invtab[p % i] % p);
    M->chitab = (i8 *)malloc(p);
    memset(M->chitab, -1, p);
    M->chitab[0] = 0;
    for (u32 x = 1; x <= (p - 1) / 2; x++) M->chitab[(u32)((u64)x * x % p)] = 1;
    M->inv2 = M->invtab[2 % p];      M->inv3 = M->invtab[3 % p];
    M->i120 = M->invtab[120 % p];    M->i6750 = M->invtab[6750 % p];
    M->i202500 = M->invtab[202500 % p];
    M->i4556250 = M->invtab[4556250 % p];
    M->i9 = M->invtab[9 % p];        M->i27 = M->invtab[27 % p];
    M->i54 = M->invtab[54 % p];      M->i162 = M->invtab[162 % p];
    M->i324 = M->invtab[324 % p];
    M->c720 = 720 % p;       M->c8640 = 8640 % p;    M->c108000 = 108000 % p;
    M->c62208 = 62208 % p;   M->c972000 = 972000 % p;
    M->c1620000 = 1620000 % p; M->c3037500 = 3037500 % p; M->c6075000 = 6075000 % p;
}

static void ek_free(EK *M) { free(M->invtab); free(M->chitab); }

/* ------------------------------------------------- small poly helpers mod p */
static inline int pdeg(u32 *a, int hi)
{
    for (int i = hi; i >= 0; i--) if (a[i]) return i;
    return -1;
}

/* c[0..da+db] += a * b */
static inline void pmulacc(EK *M, u32 *a, int da, u32 *b, int db, u32 *c)
{
    for (int i = 0; i <= da; i++) {
        if (!a[i]) continue;
        for (int j = 0; j <= db; j++)
            c[i + j] = madd(M, c[i + j], mmul(M, a[i], b[j]));
    }
}

/* degree of gcd(f, f') for f of degree df (squarefree test); result 0 = squarefree */
static int gcd_ff_deg(EK *M, u32 *fin, int df)
{
    u32 a[8], b[8];
    int da = df, db;
    for (int i = 0; i <= df; i++) a[i] = fin[i];
    for (int i = 1; i <= df; i++) b[i - 1] = mmul(M, fin[i], i % M->p);
    db = pdeg(b, df - 1);
    if (db < 0) return da;               /* f' = 0: cannot happen, p > 6 */
    while (1) {
        /* a mod b */
        u32 ilc = M->invtab[b[db]];
        for (int k = da - db; k >= 0; k--) {
            u32 c = mmul(M, a[k + db], ilc);
            if (c) for (int j = 0; j <= db; j++)
                a[k + j] = msub(M, a[k + j], mmul(M, c, b[j]));
        }
        da = pdeg(a, db - 1);
        if (da < 0) return db;
        /* swap */
        u32 t[8]; int dt;
        memcpy(t, a, sizeof(u32) * (da + 1)); dt = da;
        memcpy(a, b, sizeof(u32) * (db + 1)); da = db;
        memcpy(b, t, sizeof(u32) * (dt + 1)); db = dt;
    }
}

/* --------------------------------------------------- the per-point pipeline */
#define EK_OK    0
#define EK_DEGEN 1

/* (r,s) in F_p  ->  sextic f[0..6] (deg 5 or 6) of a Mestre model.
 * Returns EK_OK/EK_DEGEN. */
static int ek8_sextic(EK *M, u32 r, u32 s, u32 *f, int *degout)
{
    u32 p = M->p;
    if (r == 0 || s == 0) return EK_DEGEN;                    /* I10 = 0 */

    /* Igusa-Clebsch */
    u32 t3s = madd(M, madd(M, s, madd(M, s, s)), msub(M, mmul(M, 8 % p, r), 2 % p));
                                                    /* 3s + 8r - 2 */
    u32 q1 = madd(M, mmul(M, 9 % p, mmul(M, r, s)),
             madd(M, mmul(M, 4 % p, mmul(M, r, r)),
             madd(M, mmul(M, 4 % p, r), 1)));       /* 9rs + 4r^2 + 4r + 1 */
    u32 q2 = msub(M, mmul(M, madd(M, mmul(M, 54 % p, mmul(M, r, r)), mmul(M, 81 % p, r)), s),
             madd(M, mmul(M, 16 % p, mmul(M, r, mmul(M, r, r))),
             madd(M, mmul(M, 24 % p, mmul(M, r, r)),
             madd(M, mmul(M, 12 % p, r), 2 % p))));
                       /* 54r^2 s + 81rs - 16r^3 - 24r^2 - 12r - 2 */
    u32 I2  = mneg(M, mmul(M, 4 % p, t3s));
    u32 I4  = mmul(M, 4 % p, q1);
    u32 I6  = mmul(M, M->inv3,
                   msub(M, mmul(M, 4 % p, q2), mmul(M, 16 % p, mmul(M, q1, t3s))));
    u32 rs2 = mmul(M, r, mmul(M, s, s));
    u32 I10 = mneg(M, mmul(M, 8 % p, mmul(M, mmul(M, r, r), rs2)));
    if (I10 == 0) return EK_DEGEN;

    /* Clebsch A,B,C,D */
    u32 A = mmul(M, mneg(M, I2), M->i120);
    u32 A2 = mmul(M, A, A);
    u32 B = mmul(M, madd(M, I4, mmul(M, M->c720, A2)), M->i6750);
    u32 C = mmul(M, madd(M, msub(M, I6, mmul(M, M->c8640, mmul(M, A2, A))),
                            mmul(M, M->c108000, mmul(M, A, B))), M->i202500);
    u32 Dn = madd(M, I10,
             msub(M, mmul(M, M->c62208, mmul(M, A2, mmul(M, A2, A))),
             madd(M, mmul(M, M->c972000, mmul(M, mmul(M, A2, A), B)),
             msub(M, mmul(M, M->c1620000, mmul(M, A2, C)),
             madd(M, mmul(M, M->c3037500, mmul(M, A, mmul(M, B, B))),
                     mmul(M, M->c6075000, mmul(M, B, C)))))));
    u32 D = mmul(M, mneg(M, Dn), M->i4556250);

    /* Mestre conic  A11 x1^2 + A22 x2^2 + A33 x3^2 + 2A12 x1x2 + 2A31 x1x3 + 2A23 x2x3 */
    u32 BB_AC = madd(M, mmul(M, B, B), mmul(M, A, C));
    u32 A11 = madd(M, madd(M, C, C), mmul(M, mmul(M, A, B), M->inv3));
    u32 A22 = D;
    u32 A33 = madd(M, mmul(M, mmul(M, B, D), M->inv2),
                   mmul(M, mmul(M, madd(M, C, C), BB_AC), M->i9));
    u32 A23 = mmul(M, madd(M, mmul(M, B, BB_AC), mmul(M, C, A11)), M->inv3);
    u32 A31 = D;
    u32 A12 = mmul(M, madd(M, BB_AC, BB_AC), M->inv3);

    /* conic determinant (degenerate <=> extra automorphisms) */
    u32 det = madd(M, msub(M, mmul(M, A11, msub(M, mmul(M, A22, A33), mmul(M, A23, A23))),
                              mmul(M, A12, msub(M, mmul(M, A12, A33), mmul(M, A23, A31)))),
                      mmul(M, A31, msub(M, mmul(M, A12, A23), mmul(M, A22, A31))));
    if (det == 0) return EK_DEGEN;

    /* Mestre cubic coefficients */
    u32 B2c = mmul(M, B, B), C2c = mmul(M, C, C), D2c = mmul(M, D, D);
    u32 B3 = mmul(M, B2c, B), C3 = mmul(M, C2c, C), B4 = mmul(M, B2c, B2c);
    u32 B5 = mmul(M, B4, B);
    u32 AB = mmul(M, A, B), AC = mmul(M, A, C), AD = mmul(M, A, D);
    u32 BC = mmul(M, B, C), BD = mmul(M, B, D), CD = mmul(M, C, D);
    u32 ABC = mmul(M, AB, C), ABD = mmul(M, AB, D), ACD = mmul(M, AC, D);
    u32 BCD = mmul(M, BC, D);

    u32 a111 = mmul(M, mmul(M, 2 % p,
        madd(M, msub(M, mmul(M, A2, C), mmul(M, 6 % p, BC)), mmul(M, 9 % p, D))), M->i9);
    u32 a112 = mmul(M, madd(M, madd(M, mmul(M, 2 % p, B3), mmul(M, 4 % p, ABC)),
        madd(M, mmul(M, 12 % p, C2c), mmul(M, 3 % p, AD))), M->i9);
    u32 a113 = mmul(M, madd(M, madd(M, mmul(M, 3 % p, mmul(M, A, B3)),
        mmul(M, 4 % p, mmul(M, A2, BC))),
        madd(M, mmul(M, 12 % p, mmul(M, B2c, C)),
        madd(M, mmul(M, 18 % p, mmul(M, A, C2c)), mmul(M, 9 % p, BD)))), M->i27);
    u32 a122 = a113;
    u32 a123 = mmul(M, madd(M, madd(M, mmul(M, 6 % p, B4), mmul(M, 12 % p, mmul(M, AB, BC))),
        madd(M, madd(M, mmul(M, 4 % p, mmul(M, A2, C2c)), mmul(M, 12 % p, mmul(M, B, C2c))),
        madd(M, mmul(M, 9 % p, ABD), mmul(M, 36 % p, CD)))), M->i54);
    u32 a133 = mmul(M, madd(M, madd(M, mmul(M, 3 % p, mmul(M, A, B4)),
        mmul(M, 4 % p, mmul(M, A2, mmul(M, B2c, C)))),
        madd(M, madd(M, mmul(M, 16 % p, mmul(M, B3, C)), mmul(M, 26 % p, mmul(M, ABC, C))),
        madd(M, mmul(M, 24 % p, C3),
        madd(M, mmul(M, 9 % p, mmul(M, B2c, D)), mmul(M, 6 % p, ACD))))), M->i54);
    u32 a222 = mmul(M, madd(M, madd(M, mmul(M, 9 % p, B4), mmul(M, 18 % p, mmul(M, AB, BC))),
        madd(M, madd(M, mmul(M, 8 % p, mmul(M, A2, C2c)), mmul(M, 6 % p, mmul(M, B, C2c))),
        mneg(M, mmul(M, 9 % p, CD)))), M->i27);
    u32 a223 = mmul(M, madd(M, msub(M, mneg(M, mmul(M, 2 % p, mmul(M, B3, C))),
        madd(M, mmul(M, 4 % p, mmul(M, ABC, C)), mmul(M, 12 % p, C3))),
        madd(M, mmul(M, 27 % p, mmul(M, B2c, D)), mmul(M, 24 % p, ACD))), M->i54);
    u32 a233 = mmul(M, madd(M, madd(M, mmul(M, 9 % p, B5), mmul(M, 18 % p, mmul(M, A, mmul(M, B3, C)))),
        madd(M, madd(M, mmul(M, 8 % p, mmul(M, A2, mmul(M, B, C2c))),
                        mmul(M, 6 % p, mmul(M, B2c, C2c))),
        madd(M, mneg(M, mmul(M, 9 % p, BCD)), mmul(M, 81 % p, D2c)))), M->i162);
    u32 a333 = mmul(M, madd(M, msub(M, mneg(M, mmul(M, 18 % p, mmul(M, B4, C))),
        madd(M, mmul(M, 36 % p, mmul(M, AB, mmul(M, BC, C))),
        madd(M, mmul(M, 16 % p, mmul(M, A2, C3)), mmul(M, 12 % p, mmul(M, B, C3))))),
        madd(M, mmul(M, 81 % p, mmul(M, B3, D)),
        madd(M, mmul(M, 108 % p, mmul(M, ABC, D)), mmul(M, 180 % p, mmul(M, C2c, D))))), M->i324);

    /* deterministic conic point (x1, x2, 1) */
    u32 xi = 0, eta = 0; int found = 0;
    for (u32 x1 = 0; x1 < p && !found; x1++) {
        u32 qq2 = A22;
        u32 qq1 = madd(M, madd(M, mmul(M, A12, x1), A23), madd(M, mmul(M, A12, x1), A23));
        u32 qq0 = madd(M, mmul(M, A11, mmul(M, x1, x1)),
                  madd(M, mmul(M, madd(M, A31, A31), x1), A33));
        if (qq2) {
            u32 disc = msub(M, mmul(M, qq1, qq1), mmul(M, 4 % p, mmul(M, qq2, qq0)));
            u32 rt;
            if (msqrt(M, disc, &rt)) {
                xi = x1;
                eta = mmul(M, msub(M, rt, qq1), M->invtab[mmul(M, 2 % p, qq2)]);
                found = 1;
            }
        } else if (qq1) {
            xi = x1;
            eta = mmul(M, mneg(M, qq0), M->invtab[qq1]);
            found = 1;
        }
    }
    if (!found) return EK_DEGEN;

    /* parametrize: a(m) = A11 m^2 + 2A12 m + A22,
     * b(m) = 2[(A11 xi + A12 eta + A31) m + (A12 xi + A22 eta + A23)]
     * X1 = xi a - m b, X2 = eta a - b, X3 = a   (each quadratic in m) */
    u32 am[3], bm[2], X1[3], X2[3], X3[3];
    am[0] = A22; am[1] = madd(M, A12, A12); am[2] = A11;
    {
        u32 b1 = madd(M, madd(M, mmul(M, A11, xi), mmul(M, A12, eta)), A31);
        u32 b0 = madd(M, madd(M, mmul(M, A12, xi), mmul(M, A22, eta)), A23);
        bm[0] = madd(M, b0, b0); bm[1] = madd(M, b1, b1);
    }
    X1[0] = mmul(M, xi, am[0]);
    X1[1] = msub(M, mmul(M, xi, am[1]), bm[0]);
    X1[2] = msub(M, mmul(M, xi, am[2]), bm[1]);
    X2[0] = msub(M, mmul(M, eta, am[0]), bm[0]);
    X2[1] = msub(M, mmul(M, eta, am[1]), bm[1]);
    X2[2] = mmul(M, eta, am[2]);
    X3[0] = am[0]; X3[1] = am[1]; X3[2] = am[2];

    /* f(m) = M(X1,X2,X3), degree <= 6 */
    u32 P11[5], P22[5], P33[5], P12[5], P13[5], P23[5];
    memset(P11, 0, sizeof P11); pmulacc(M, X1, 2, X1, 2, P11);
    memset(P22, 0, sizeof P22); pmulacc(M, X2, 2, X2, 2, P22);
    memset(P33, 0, sizeof P33); pmulacc(M, X3, 2, X3, 2, P33);
    memset(P12, 0, sizeof P12); pmulacc(M, X1, 2, X2, 2, P12);
    memset(P13, 0, sizeof P13); pmulacc(M, X1, 2, X3, 2, P13);
    memset(P23, 0, sizeof P23); pmulacc(M, X2, 2, X3, 2, P23);

    u32 acc[7], term[7];
    memset(acc, 0, sizeof acc);
    struct { u32 coef; u32 *pa; int da; u32 *pb; int db; } terms[10] = {
        { a111, P11, 4, X1, 2 }, { a222, P22, 4, X2, 2 }, { a333, P33, 4, X3, 2 },
        { mmul(M, 3 % p, a112), P11, 4, X2, 2 },
        { mmul(M, 3 % p, a113), P11, 4, X3, 2 },
        { mmul(M, 3 % p, a122), P22, 4, X1, 2 },
        { mmul(M, 3 % p, a133), P33, 4, X1, 2 },
        { mmul(M, 3 % p, a233), P33, 4, X2, 2 },
        { mmul(M, 3 % p, a223), P22, 4, X3, 2 },
        { mmul(M, 6 % p, a123), P12, 4, X3, 2 } };
    for (int t = 0; t < 10; t++) {
        if (!terms[t].coef) continue;
        memset(term, 0, sizeof term);
        pmulacc(M, terms[t].pa, terms[t].da, terms[t].pb, terms[t].db, term);
        for (int i = 0; i < 7; i++)
            acc[i] = madd(M, acc[i], mmul(M, terms[t].coef, term[i]));
    }
    int deg = pdeg(acc, 6);
    if (deg < 5) return EK_DEGEN;
    if (gcd_ff_deg(M, acc, deg) > 0) return EK_DEGEN;     /* not squarefree */
    for (int i = 0; i < 7; i++) f[i] = acc[i];
    *degout = deg;
    return EK_OK;
}

/* exact chi(1), chi(-1) of y^2 = f over F_p by counting C(F_p), C(F_{p^2}) */
static void ek8_count(EK *M, u32 *f, int deg, i64 *chi1, i64 *chim1)
{
    u32 p = M->p;
    /* N1 and number of F_p-roots z */
    i64 N1 = (deg == 6) ? (1 + M->chitab[f[6]]) : 1;
    i64 z = 0;
    for (u32 x = 0; x < p; x++) {
        u32 v = f[deg];
        for (int i = deg - 1; i >= 0; i--) v = madd(M, mmul(M, v, x), f[i]);
        i8 c = M->chitab[v];
        N1 += 1 + c;
        if (c == 0) z++;
    }
    /* N2 : F_{p^2} = F_p[i], i^2 = -1 ; x = a + b i, b = 1..(p-1)/2 doubled */
    i64 N2 = ((deg == 6) ? 2 : 1) + (2 * (i64)p - z);
    for (u32 b = 1; b <= (p - 1) / 2; b++) {
        /* values f(k + b i), k = 0..6 (Horner in F_{p^2}) */
        u32 vr[7], vi[7];
        for (u32 k = 0; k <= 6; k++) {
            u32 ur = f[deg], ui = 0;
            for (int i = deg - 1; i >= 0; i--) {
                u32 nr = msub(M, mmul(M, ur, k), mmul(M, ui, b));
                u32 ni = madd(M, mmul(M, ur, b), mmul(M, ui, k));
                ur = madd(M, nr, f[i]); ui = ni;
            }
            vr[k] = ur; vi[k] = ui;
        }
        /* forward differences of order 6 */
        u32 dr[7], di[7];
        for (int j = 0; j <= 6; j++) { dr[j] = vr[j]; di[j] = vi[j]; }
        for (int o = 1; o <= 6; o++)
            for (int j = 6; j >= o; j--) {
                dr[j] = msub(M, dr[j], dr[j - 1]);
                di[j] = msub(M, di[j], di[j - 1]);
            }
        /* re-pack: D[j] = j-th forward difference at a = 0 */
        i32 rowsum = 0;
        u32 D0r = dr[0], D1r = dr[1], D2r = dr[2], D3r = dr[3], D4r = dr[4], D5r = dr[5], D6r = dr[6];
        u32 D0i = di[0], D1i = di[1], D2i = di[2], D3i = di[3], D4i = di[4], D5i = di[5], D6i = di[6];
        for (u32 a = 0; a < p; a++) {
            u32 nrm = madd(M, mmul(M, D0r, D0r), mmul(M, D0i, D0i));
            rowsum += M->chitab[nrm];
            D0r = madd(M, D0r, D1r); D1r = madd(M, D1r, D2r); D2r = madd(M, D2r, D3r);
            D3r = madd(M, D3r, D4r); D4r = madd(M, D4r, D5r); D5r = madd(M, D5r, D6r);
            D0i = madd(M, D0i, D1i); D1i = madd(M, D1i, D2i); D2i = madd(M, D2i, D3i);
            D3i = madd(M, D3i, D4i); D4i = madd(M, D4i, D5i); D5i = madd(M, D5i, D6i);
        }
        N2 += 2 * ((i64)p + rowsum);
    }
    i64 s1 = (i64)p + 1 - N1;
    i64 s2 = (N2 - ((i64)p * p + 1) + s1 * s1) / 2;
    *chi1  = (i64)p * p + 1 - ((i64)p + 1) * s1 + s2;
    *chim1 = (i64)p * p + 1 + ((i64)p + 1) * s1 + s2;
}

/* full test: mask bit0 = 31|chi(1), bit1 = 31|chi(-1); DEGEN -> 7 */
#define MASK_DEGEN 7
static inline u32 ek8_mask(EK *M, u32 r, u32 s, i64 *c1out, i64 *cm1out)
{
    u32 f[7]; int deg;
    if (ek8_sextic(M, r, s, f, &deg) == EK_DEGEN) return MASK_DEGEN;
    i64 c1, cm1;
    ek8_count(M, f, deg, &c1, &cm1);
    if (c1out) *c1out = c1;
    if (cm1out) *cm1out = cm1;
    return (u32)((c1 % 31 == 0) | ((cm1 % 31 == 0) << 1));
}

/* ------------------------------------------------------------- rational rec */
static int ratrecon(i64 x, i64 m, i64 bound, i64 *num, i64 *den)
{
    i64 r0 = m, r1 = x % m; if (r1 < 0) r1 += m;
    i64 t0 = 0, t1 = 1;
    while (r1 > bound) {
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

/* ----------------------------------------------------------------- selftest */
/* RM witness r = 6000/96721, s = 3557520/96721; Magma anchor multisets
 * {chi(1), chi(-1)} from claude_z31_ek8_witness.m smoke (results/claude_z31_ek8_witness.log) */
static void do_selftest(void)
{
    struct { u32 p; i64 c1, c2; } anchors[] = {
        { 103,  10199,   11447 },
        { 211,  52898,   37634 },
        { 1019, 984064,  1098304 },
        { 1031, 1009783, 1121239 },
        { 1039, 1019900, 1144700 },
        { 1051, 1088968, 1122632 },
        { 1063, 1223204, 1044452 },
        { 1087, 1276898, 1094114 } };
    int pass = 0, fail = 0;
    for (int t = 0; t < 8; t++) {
        EK M; ek_init(&M, anchors[t].p);
        u32 p = M.p;
        u32 iden = M.invtab[96721 % p];
        u32 r0 = mmul(&M, 6000 % p, iden), s0 = mmul(&M, 3557520 % p, iden);
        i64 c1 = 0, cm1 = 0;
        u32 mask = ek8_mask(&M, r0, s0, &c1, &cm1);
        int ok = mask != MASK_DEGEN && mask != 0 &&
                 ((c1 == anchors[t].c1 && cm1 == anchors[t].c2) ||
                  (c1 == anchors[t].c2 && cm1 == anchors[t].c1));
        printf("SELFTEST p=%u r0=%u s0=%u mask=%u chi={%lld,%lld} expect={%lld,%lld} %s\n",
               p, r0, s0, mask, (long long)c1, (long long)cm1,
               (long long)anchors[t].c1, (long long)anchors[t].c2, ok ? "PASS" : "FAIL");
        if (ok) pass++; else fail++;
        ek_free(&M);
    }
    printf("SELFTEST_SUMMARY pass=%d fail=%d\n", pass, fail);
    if (fail) exit(1);
}

/* ---------------------------------------------------------------- ref check */
static void do_ref(char *fname)
{
    FILE *fp = fopen(fname, "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", fname); exit(1); }
    char line[512], flag[16];
    u32 curp = 0;
    EK M;
    int pass = 0, fail = 0, skip = 0;
    while (fgets(line, sizeof line, fp)) {
        u32 p, r0, s0; long long L1, Lm1, s1, s2;
        if (sscanf(line, "%u\t%u\t%u\t%15s\t%lld\t%lld\t%lld\t%lld",
                   &p, &r0, &s0, flag, &L1, &Lm1, &s1, &s2) != 8) continue;
        if (p != curp) { if (curp) ek_free(&M); ek_init(&M, p); curp = p; }
        if (!strcmp(flag, "ERR") || !strcmp(flag, "MISM")) { skip++; continue; }
        i64 c1 = 0, cm1 = 0;
        u32 mask = ek8_mask(&M, r0, s0, &c1, &cm1);
        int ok;
        if (!strcmp(flag, "DEGEN")) ok = (mask == MASK_DEGEN);
        else ok = (mask != MASK_DEGEN) &&
                  ((c1 == L1 && cm1 == Lm1) || (c1 == Lm1 && cm1 == L1));
        if (ok) pass++;
        else {
            fail++;
            printf("REF_FAIL p=%u r0=%u s0=%u flag=%s got mask=%u chi={%lld,%lld} want={%lld,%lld}\n",
                   p, r0, s0, flag, mask, (long long)c1, (long long)cm1, L1, Lm1);
        }
    }
    if (curp) ek_free(&M);
    fclose(fp);
    printf("REF_SUMMARY file=%s pass=%d fail=%d skip=%d\n", fname, pass, fail, skip);
    if (fail) exit(1);
}

/* -------------------------------------------------------------------- fiber */
static void do_fiber(u32 p, i64 rnum, i64 rden)
{
    EK M; ek_init(&M, p);
    i64 nn = ((rnum % p) + p) % p, dd = ((rden % p) + p) % p;
    if (dd == 0 || nn == 0) { printf("FIBER p=%u r=%lld/%lld BAD_REDUCTION\n",
                                     p, (long long)rnum, (long long)rden); return; }
    u32 r0 = mmul(&M, (u32)nn, M.invtab[dd]);
    int nacc = 0, ndeg = 0;
    for (u32 s0 = 0; s0 < p; s0++) {
        i64 c1, cm1;
        u32 mask = ek8_mask(&M, r0, s0, &c1, &cm1);
        if (mask == MASK_DEGEN) { ndeg++; printf("  s0=%u DEGEN\n", s0); }
        else if (mask) { nacc++; printf("  s0=%u mask=%u chi1=%lld chim1=%lld\n",
                                        s0, mask, (long long)c1, (long long)cm1); }
    }
    printf("FIBER p=%u r=%lld/%lld r0=%u accepted=%d degen=%d of %u\n",
           p, (long long)rnum, (long long)rden, r0, nacc, ndeg, p);
    ek_free(&M);
}

/* -------------------------------------------------------------------- bench */
static void do_bench(u32 p, int nrows)
{
    EK M; ek_init(&M, p);
    double t0 = omp_get_wtime();
    u64 cnt = 0, acc = 0, deg = 0;
    for (int row = 0; row < nrows; row++) {
        u32 r0 = 1 + (u32)((u64)(row * 2654435761u) % (p - 1));
        for (u32 s0 = 1; s0 < p; s0++) {
            u32 mask = ek8_mask(&M, r0, s0, NULL, NULL);
            cnt++;
            if (mask == MASK_DEGEN) deg++;
            else if (mask) acc++;
        }
    }
    double dt = omp_get_wtime() - t0;
    printf("BENCH p=%u rows=%d candidates=%llu accepted=%llu degen=%llu  %.2f s  "
           "%.1f cand/s/core  accept_rate=%.4f (expect ~%.4f)\n",
           p, nrows, (unsigned long long)cnt, (unsigned long long)acc,
           (unsigned long long)deg, dt, (double)cnt / dt,
           (double)acc / (double)cnt, 2.0 / 31.0 - 1.0 / (31.0 * 31.0));
    ek_free(&M);
}

/* -------------------------------------------------------------------- sweep */
typedef struct { i64 n, d; } Rat;

static int build_rats(int H, Rat *out, int cap)
{
    int n = 0;
    for (i64 den = 1; den <= H; den++)
        for (i64 num = -H; num <= H; num++) {
            if (num == 0) continue;                     /* r = 0 is I10 = 0 */
            i64 a = num < 0 ? -num : num, b = den;
            while (b) { i64 t = a % b; a = b; b = t; }
            if (a != 1) continue;
            if (n < cap) { out[n].n = num; out[n].d = den; }
            n++;
        }
    return n;
}

#define MAXPR 16

static void do_sweep(int H, int npr, u32 *prs, int nver, u32 *vprs,
                     int lo, int hi, i64 recbound,
                     Rat *xrats, int nxr)
{
    int ntot = npr + nver;
    if (npr < 2 || npr > 6 || ntot > MAXPR) {
        fprintf(stderr, "need 2 <= npr <= 6 (CRT overflow bound) and npr+nver <= %d\n", MAXPR);
        exit(1);
    }
    u32 allp[MAXPR];
    for (int t = 0; t < npr; t++) allp[t] = prs[t];
    for (int t = 0; t < nver; t++) allp[npr + t] = vprs[t];

    Rat *rats; int nr;
    if (xrats) { rats = xrats; nr = nxr; }
    else {
        rats = (Rat *)malloc(sizeof(Rat) * 4 * (H + 1) * (H + 1));
        nr = build_rats(H, rats, 4 * (H + 1) * (H + 1));
    }
    if (hi > nr || hi <= 0) hi = nr;
    if (lo < 0) lo = 0;
    printf("SWEEP H=%d #slices=%d npr=%d nver=%d shard [%d,%d) recbound=%lld\n",
           H, nr, npr, nver, lo, hi, (long long)recbound);
    fflush(stdout);

    EK M[MAXPR];
    for (int t = 0; t < ntot; t++) ek_init(&M[t], allp[t]);

    /* residues needed per prime */
    unsigned char **needed = (unsigned char **)malloc(sizeof(void *) * ntot);
    unsigned char **tab = (unsigned char **)malloc(sizeof(void *) * ntot);
    for (int t = 0; t < ntot; t++) {
        needed[t] = (unsigned char *)calloc(allp[t], 1);
        tab[t] = (unsigned char *)malloc((size_t)allp[t] * allp[t]);
        memset(tab[t], 0xFF, (size_t)allp[t] * allp[t]);   /* 0xFF = not computed */
    }
    for (int i = lo; i < hi; i++) {
        for (int t = 0; t < ntot; t++) {
            u32 p = allp[t];
            i64 dd = ((rats[i].d % p) + p) % p, nn = ((rats[i].n % p) + p) % p;
            if (dd == 0 || nn == 0) continue;
            needed[t][mmul(&M[t], (u32)nn, M[t].invtab[dd])] = 1;
        }
    }
    /* row work items */
    typedef struct { int t; u32 r0; } RowJob;
    RowJob *jobs = (RowJob *)malloc(sizeof(RowJob) * (size_t)(ntot) * 2200);
    int njobs = 0;
    for (int t = 0; t < ntot; t++)
        for (u32 r0 = 0; r0 < allp[t]; r0++)
            if (needed[t][r0]) { jobs[njobs].t = t; jobs[njobs].r0 = r0; njobs++; }
    printf("SWEEP fiber rows to compute: %d\n", njobs);
    fflush(stdout);
    double t0 = omp_get_wtime();
    u64 fibsum = 0, degsum = 0;
    int rowsdone = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:fibsum,degsum)
    for (int j = 0; j < njobs; j++) {
        int t = jobs[j].t; u32 r0 = jobs[j].r0, p = allp[t];
        unsigned char *row = tab[t] + (size_t)r0 * p;
        for (u32 s0 = 0; s0 < p; s0++) {
            u32 mask = ek8_mask(&M[t], r0, s0, NULL, NULL);
            row[s0] = (unsigned char)mask;
            if (mask == MASK_DEGEN) degsum++;
            else if (mask) fibsum++;
        }
#pragma omp critical
        {
            rowsdone++;
            if (rowsdone % 200 == 0) {
                printf("PROGRESS rows %d/%d  %.1f s\n", rowsdone, njobs,
                       omp_get_wtime() - t0);
                fflush(stdout);
            }
        }
    }
    printf("SWEEP rows done: %d  mean_fiber=%.1f mean_degen=%.1f  %.1f s\n",
           njobs, njobs ? (double)fibsum / njobs : 0.0,
           njobs ? (double)degsum / njobs : 0.0, omp_get_wtime() - t0);
    fflush(stdout);

    /* CRT matching stage.
     * Phase A: per usable slice gather the npr CRT fibers.
     * Phase B: flatten work items (slice, index-in-fiber-0) so that a single
     *          big slice (validation mode) still uses every core.
     * Phase C: enumerate levels 1..npr-1 with hoisted partial CRT. */
    double t1 = omp_get_wtime();
    u64 done = 0, skipped = 0, cand = 0, combos = 0;
    int nsl = hi - lo;
    typedef struct {
        int usable;
        unsigned char *rowv[MAXPR];       /* all ntot rows (verify incl.) */
        u32 *fib[MAXPR];                  /* npr fiber lists */
        int nf[MAXPR];
    } Slice;
    Slice *sl = (Slice *)calloc(nsl, sizeof(Slice));
    for (int i = lo; i < hi; i++) {
        Slice *S = &sl[i - lo];
        S->usable = 1;
        for (int t = 0; t < ntot; t++) {
            u32 p = allp[t];
            i64 dd = ((rats[i].d % p) + p) % p, nn = ((rats[i].n % p) + p) % p;
            if (dd == 0 || nn == 0) { if (t < npr) S->usable = 0; S->rowv[t] = NULL; continue; }
            S->rowv[t] = tab[t] + (size_t)mmul(&M[t], (u32)nn, M[t].invtab[dd]) * p;
        }
        if (!S->usable) { skipped++; continue; }
        for (int t = 0; t < npr; t++) {
            u32 p = allp[t];
            S->fib[t] = (u32 *)malloc(sizeof(u32) * p);
            S->nf[t] = 0;
            for (u32 s0 = 0; s0 < p; s0++)
                if (S->rowv[t][s0]) S->fib[t][S->nf[t]++] = s0;
            if (S->nf[t] == 0) S->usable = 0;
        }
        if (S->usable) done++;
    }
    /* work items: flatten one fiber level, or two when there are few slices
     * (validation mode) so all cores stay busy */
    int depth2 = (nsl <= 64 && npr >= 3);
    i64 nwork = 0;
    for (int i = 0; i < nsl; i++)
        if (sl[i].usable)
            nwork += depth2 ? (i64)sl[i].nf[0] * sl[i].nf[1] : sl[i].nf[0];
    int *wsl = (int *)malloc(sizeof(int) * (nwork ? nwork : 1));
    u32 *wj0 = (u32 *)malloc(sizeof(u32) * (nwork ? nwork : 1));
    i32 *wj1 = (i32 *)malloc(sizeof(i32) * (nwork ? nwork : 1));
    {
        i64 w = 0;
        for (int i = 0; i < nsl; i++) {
            if (!sl[i].usable) continue;
            for (int j = 0; j < sl[i].nf[0]; j++) {
                if (depth2)
                    for (int k = 0; k < sl[i].nf[1]; k++)
                        { wsl[w] = i; wj0[w] = (u32)j; wj1[w] = k; w++; }
                else { wsl[w] = i; wj0[w] = (u32)j; wj1[w] = -1; w++; }
            }
        }
    }
    printf("SWEEP match stage: %llu usable slices, %lld work items (depth2=%d)\n",
           (unsigned long long)done, (long long)nwork, depth2);
    fflush(stdout);
    /* per-level constant CRT data: modulus product and inverse of the
     * previous product mod the level prime */
    i64 lmm[MAXPR]; u32 lminv[MAXPR];
    lmm[0] = allp[0];
    for (int t = 1; t < npr; t++) {
        u32 p = allp[t];
        u32 mred = (u32)(lmm[t - 1] % p);
        lminv[t] = M[t].invtab[mred];         /* mred != 0: distinct primes */
        lmm[t] = lmm[t - 1] * allp[t];
    }
    i64 reported = 0;
#pragma omp parallel for schedule(dynamic, 8) reduction(+:cand,combos)
    for (i64 w = 0; w < nwork; w++) {
        Slice *S = &sl[wsl[w]];
        int i = wsl[w] + lo;
        i64 x[MAXPR];
        int idx[MAXPR];
        x[0] = S->fib[0][wj0[w]];
        idx[0] = (int)wj0[w];
        int lev = 1;
        if (wj1[w] >= 0) {
            idx[1] = wj1[w];
            u32 p1 = allp[1];
            i64 d1 = ((i64)S->fib[1][idx[1]] - x[0]) % p1; if (d1 < 0) d1 += p1;
            x[1] = x[0] + lmm[0] * (i64)mmul(&M[1], (u32)d1, lminv[1]);
            lev = 2;
        }
        idx[lev] = 0;
        int minlev = (wj1[w] >= 0) ? 2 : 1;
        while (lev >= minlev) {
            if (idx[lev] >= S->nf[lev]) { lev--; if (lev >= minlev) idx[lev]++; continue; }
            {
                u32 pl = allp[lev];
                i64 dl = ((i64)S->fib[lev][idx[lev]] - x[lev - 1]) % pl; if (dl < 0) dl += pl;
                x[lev] = x[lev - 1] + lmm[lev - 1] * (i64)mmul(&M[lev], (u32)dl, lminv[lev]);
            }
            if (lev == npr - 1) {
                combos++;
                i64 sn, sd;
                if (ratrecon(x[lev], lmm[lev], recbound, &sn, &sd) && sn != 0) {
                    int vok = 1;
                    for (int t = npr; t < ntot && vok; t++) {
                        if (!S->rowv[t]) continue;
                        u32 p = allp[t];
                        i64 dd = ((sd % p) + p) % p, nn2 = ((sn % p) + p) % p;
                        if (dd == 0) continue;          /* cannot check */
                        u32 s0 = mmul(&M[t], (u32)nn2, M[t].invtab[dd]);
                        if (!S->rowv[t][s0]) vok = 0;
                    }
                    if (vok) {
#pragma omp critical
                        {
                            printf("CAND r=%lld/%lld s=%lld/%lld masks=",
                                   (long long)rats[i].n, (long long)rats[i].d,
                                   (long long)sn, (long long)sd);
                            for (int t = 0; t < npr; t++)
                                printf("%u:%u ", allp[t],
                                       S->rowv[t][S->fib[t][idx[t]]]);
                            printf("\n");
                            fflush(stdout);
                        }
                        cand++;
                    }
                }
                idx[lev]++;
            } else {
                lev++;
                idx[lev] = 0;
            }
        }
#pragma omp critical
        {
            reported++;
            if (reported % 2000 == 0) {
                printf("PROGRESS work %lld/%lld cand=%llu  %.1f s\n",
                       (long long)reported, (long long)nwork, (unsigned long long)cand,
                       omp_get_wtime() - t1);
                fflush(stdout);
            }
        }
    }
    printf("SWEEP_DONE H=%d shard [%d,%d) slices=%llu skipped=%llu combos=%llu "
           "candidates=%llu  fiber_stage=%.1f s match_stage=%.1f s\n",
           H, lo, hi, (unsigned long long)done, (unsigned long long)skipped,
           (unsigned long long)combos, (unsigned long long)cand,
           t1 - t0, omp_get_wtime() - t1);
    for (int i = 0; i < nsl; i++)
        if (sl[i].fib[0])
            for (int t = 0; t < npr; t++) free(sl[i].fib[t]);
    free(sl); free(wsl); free(wj0); free(wj1);
    for (int t = 0; t < ntot; t++) { free(needed[t]); free(tab[t]); ek_free(&M[t]); }
    free(needed); free(tab); free(jobs);
    if (!xrats) free(rats);
}

/* --------------------------------------------------------------------- main */
int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s selftest | ref FILE | fiber P rnum rden | "
                "bench P [nrows] | sweep H npr P.. nver V.. lo hi recbound | "
                "sweepr rnum rden npr P.. nver V.. recbound\n", argv[0]);
        return 1;
    }
    if (!strcmp(argv[1], "selftest")) do_selftest();
    else if (!strcmp(argv[1], "ref")) do_ref(argv[2]);
    else if (!strcmp(argv[1], "fiber"))
        do_fiber((u32)atoi(argv[2]), atoll(argv[3]), atoll(argv[4]));
    else if (!strcmp(argv[1], "bench"))
        do_bench((u32)atoi(argv[2]), argc > 3 ? atoi(argv[3]) : 1);
    else if (!strcmp(argv[1], "sweep")) {
        int H = atoi(argv[2]);
        int npr = atoi(argv[3]);
        u32 prs[MAXPR], vprs[MAXPR];
        int a = 4;
        for (int t = 0; t < npr; t++) prs[t] = (u32)atoi(argv[a++]);
        int nver = atoi(argv[a++]);
        for (int t = 0; t < nver; t++) vprs[t] = (u32)atoi(argv[a++]);
        int lo = atoi(argv[a++]), hi = atoi(argv[a++]);
        i64 rb = atoll(argv[a++]);
        do_sweep(H, npr, prs, nver, vprs, lo, hi, rb, NULL, 0);
    }
    else if (!strcmp(argv[1], "sweepr")) {
        /* single-slice sweep at an explicit rational r (validation mode) */
        Rat rr; rr.n = atoll(argv[2]); rr.d = atoll(argv[3]);
        int npr = atoi(argv[4]);
        u32 prs[MAXPR], vprs[MAXPR];
        int a = 5;
        for (int t = 0; t < npr; t++) prs[t] = (u32)atoi(argv[a++]);
        int nver = atoi(argv[a++]);
        for (int t = 0; t < nver; t++) vprs[t] = (u32)atoi(argv[a++]);
        i64 rb = atoll(argv[a++]);
        do_sweep(0, npr, prs, nver, vprs, 0, 1, rb, &rr, 1);
    }
    else { fprintf(stderr, "unknown mode\n"); return 1; }
    return 0;
}
