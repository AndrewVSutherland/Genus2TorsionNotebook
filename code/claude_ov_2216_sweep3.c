/*
 *  claude_ov_2216_sweep3.c   -- lane OV/2216   (fast rewrite of claude_ov_2216_sweep.c)
 *
 *  SAME MATHEMATICS as claude_ov_2216_sweep.c (read that file's header for the
 *  derivation; the underlying x-T module code/claude_ov_2216_delta.m is
 *  validated 310/310 against exact torsion in
 *  results/claude_ov_2216_delta_validate.log).  Only the sieve engine changed.
 *
 *  ------------------------------------------------------------------
 *  RECAP OF THE CONDITIONS.  Odd model f_{u,v} on M_1(8,2,2), marked class Q
 *  of order 8, L = Q^3 x K.  Q+T is 2-divisible iff  d1*d2, d1*d3 are rational
 *  squares and  d1*dK  is a square in K; the purely rational necessary part of
 *  the last one is  N(dK) = u*v*(u+v+1) = square, the SAME for all eight T.
 *  Reduction of the 8 twists (code/claude_ov_2216_symbolic.m):
 *
 *    Y0  (T = 0 and T = B_quad = 4Q):   u(u-1)v(v-1) = sq, u(u-1)(u+v+1)(u+v+2) = sq
 *    Y1  (T = A_r2r3 = C_r1):          -u(u-1)v(2u+v+1) = sq, u(u-1)(u-v)(u+v+1) = sq
 *    (the two remaining condition sets  A_r1r2 = C_r3  and  A_r1r3 = C_r2  are
 *     the images of Y1 under the S_3 that permutes {u, v, w = -(u+v+1)};
 *     verified in results/claude_ov_2216_s3check.log)
 *
 *  For FIXED u = p/q the first condition is a conic through (v,y) = (0,0);
 *  the pencil y = k*v, k = a/b, parametrizes ALL of its rational points:
 *      Y0:  vn =  An b^2,            vd = An b^2 - a^2 Ad
 *      Y1:  vn = -An (2p+q) b^2,     vd = q (a^2 Ad + An b^2)
 *  with An = p(p-q), Ad = q^2, and then, square factors dropping out,
 *      S1 = (p+q) vd + q vn,   S2 = S1 + q vd,   Dn = p vd - q vn
 *      (2)  Y0: An*S1*S2 square      Y1: An*Dn*S1 square
 *      (3)      p*vn*S1 square
 *
 *  Closed forms used by the fast tier (exact integer identities):
 *      Y0:  S1 = An(p+2q) b^2 - Ad(p+q) a^2      S2 = An(p+3q) b^2 - Ad(p+2q) a^2
 *           c3 = (p*An) * S1 * b^2               c2 = An * S1 * S2
 *      Y1:  S1 = q[(p+q) Ad a^2 - p An b^2]      Dn = q[p Ad a^2 + An(3p+q) b^2]
 *           c3 = (-p*An*(2p+q)) * S1 * b^2       c2 = An * Dn * S1
 *  i.e. in both cases
 *      S1 = alpha1 a^2 + beta1 b^2,   S2|Dn = alpha2 a^2 + beta2 b^2,
 *      c3 = G3 * S1 * b^2,            c2 = An * S1 * (S2|Dn).
 *
 *  ------------------------------------------------------------------
 *  ENGINE.  For fixed (p,q) and fixed b, the pass/fail of a prime P is a
 *  function of a mod P alone (it depends on a^2 mod P), and it depends on b
 *  only through b^2 mod P -- so the periodic bit patterns are CACHED per
 *  (p,q) and indexed by b^2 mod P, which makes their construction free.  So:
 *
 *    tier A  for each of NSMALL small primes build the length-P accept table
 *            from Legendre tables (no multiplications: chi(c3) = chi(G3)
 *            chi(S1) chi(b^2), chi(c2) = chi(An) chi(S1) chi(S2)), lay it out
 *            periodically as a bit array over a in [0,Hk], AND the NSMALL of
 *            them into a master bitmap, and scan the (very few) set bits.
 *            Zero residues always ACCEPT, so no true solution is discarded.
 *    tier B  the survivors get the NBIG large primes by the same exact modular
 *            arithmetic as v1 (full vn, vd, S1, S2, Dn -- no factor dropped).
 *    tier C  exact 128-bit degeneracy screen and gcd(a,b) = 1, as in v1.
 *
 *  Effective filter density (1/4)^(NSMALL+NBIG) ~ 2^-64, so a survivor is
 *  essentially always a genuine rational point of the norm surface; every
 *  survivor is re-verified in exact rational arithmetic by
 *  code/claude_ov_2216_hits2uv.py and then exact-tested in Magma.
 *
 *  a runs over [0,Hk] and b over [1,Hk] WITHOUT a coprimality prefilter
 *  (gcd is checked on survivors only); k = a/b in lowest terms is therefore
 *  covered whenever height(k) <= Hk, exactly as in v1.
 *
 *  Usage:  ./claude_ov_2216_sweep3 <surface 0|1> <Hu> <Hk> <nthreads> [qlo] [qhi]
 *  Output: "HIT <surface> <p> <q> <a> <b>"   (u = p/q, k = a/b)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

/* ---- small primes: the bitmap tier ---- */
#define NSMALL 8
static int32_t SP[NSMALL] = { 31, 37, 41, 43, 47, 53, 59, 61 };

/* ---- large primes: the exact tier (identical list to v1) ---- */
#define NBIG 24
static int32_t BP[NBIG] = {
    211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,
    277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349
};

static int8_t  *LEG[NSMALL];    /* Legendre symbol, LEG[i][x] in {-1,0,1}   */
static int32_t *SQT[NSMALL];    /* SQT[i][x] = x^2 mod SP[i]                */
static uint8_t *QRB[NBIG];      /* QRB[i][x] = 1 iff x is 0 or a QR mod BP[i] */

static void build_tables(void)
{
    for (int i = 0; i < NSMALL; i++) {
        int32_t P = SP[i];
        LEG[i] = (int8_t *)malloc((size_t)P);
        SQT[i] = (int32_t *)malloc(sizeof(int32_t) * (size_t)P);
        for (int32_t x = 0; x < P; x++) { LEG[i][x] = -1; SQT[i][x] = (int32_t)(((int64_t)x*x) % P); }
        LEG[i][0] = 0;
        for (int32_t x = 1; x < P; x++) LEG[i][SQT[i][x]] = 1;
    }
    for (int i = 0; i < NBIG; i++) {
        int32_t P = BP[i];
        QRB[i] = (uint8_t *)calloc((size_t)P, 1);
        QRB[i][0] = 1;
        for (int32_t x = 1; x < P; x++) QRB[i][(int32_t)(((int64_t)x*x) % P)] = 1;
    }
}

static inline int64_t igcd(int64_t a, int64_t b)
{
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

static inline int32_t mmod(int64_t x, int32_t P) { int32_t r = (int32_t)(x % P); return r < 0 ? r + P : r; }

int main(int argc, char **argv)
{
    if (argc < 5) {
        fprintf(stderr, "usage: %s <surface 0|1> <Hu> <Hk> <nthreads> [qlo] [qhi]\n", argv[0]);
        return 1;
    }
    int surface  = atoi(argv[1]);
    int64_t Hu   = atoll(argv[2]);
    int64_t Hk   = atoll(argv[3]);
    int nthreads = atoi(argv[4]);
    int64_t qlo  = (argc > 5) ? atoll(argv[5]) : 1;
    int64_t qhi  = (argc > 6) ? atoll(argv[6]) : Hu;
    omp_set_num_threads(nthreads);
    build_tables();

    printf("# claude_ov_2216_sweep3 surface=%d Hu=%lld Hk=%lld qlo=%lld qhi=%lld threads=%d nsmall=%d nbig=%d\n",
           surface, (long long)Hu, (long long)Hk, (long long)qlo, (long long)qhi,
           nthreads, NSMALL, NBIG);
    fflush(stdout);

    /* flat list of (p,q), coprime, p != 0, p != q  --  for load balance */
    int64_t npq = 0;
    for (int64_t q = qlo; q <= qhi; q++)
        for (int64_t p = -Hu; p <= Hu; p++) {
            if (p == 0 || p == q) continue;
            if (igcd(p, q) != 1) continue;
            if (surface == 1 && 2*p + q == 0) continue;
            npq++;
        }
    int32_t *PL = (int32_t *)malloc(sizeof(int32_t) * (size_t)npq);
    int32_t *QL = (int32_t *)malloc(sizeof(int32_t) * (size_t)npq);
    {
        int64_t i = 0;
        for (int64_t q = qlo; q <= qhi; q++)
            for (int64_t p = -Hu; p <= Hu; p++) {
                if (p == 0 || p == q) continue;
                if (igcd(p, q) != 1) continue;
                if (surface == 1 && 2*p + q == 0) continue;
                PL[i] = (int32_t)p; QL[i] = (int32_t)q; i++;
            }
    }
    fprintf(stderr, "# (p,q) list size %lld,  a in [0,%lld], b in [1,%lld]\n",
            (long long)npq, (long long)Hk, (long long)Hk);

    int64_t WMASK = (Hk + 1 + 63) / 64;          /* words of the a-bitmap  */
    int64_t total_hits = 0, total_degen = 0, total_tierB = 0;
    int64_t done_pq = 0;

#pragma omp parallel reduction(+:total_hits,total_degen,total_tierB)
    {
        uint64_t *master = (uint64_t *)malloc(sizeof(uint64_t) * (size_t)(WMASK + 2));
        /* periodic patterns, CACHED per (p,q) and indexed by b^2 mod P */
        uint64_t *patc[NSMALL];
        uint8_t  *patv[NSMALL];
        int8_t   *acc[NSMALL];
        int64_t   PW[NSMALL];
        int32_t   STEP[NSMALL];
        for (int i = 0; i < NSMALL; i++) {
            PW[i]   = (SP[i] + 128) / 64 + 3;
            STEP[i] = 64 % SP[i];
            patc[i] = (uint64_t *)malloc(sizeof(uint64_t) * (size_t)(PW[i] * SP[i]));
            patv[i] = (uint8_t *)malloc((size_t)SP[i]);
            acc[i]  = (int8_t *)malloc((size_t)SP[i]);
        }

#pragma omp for schedule(dynamic,256)
        for (int64_t idx = 0; idx < npq; idx++) {
            int64_t p = PL[idx], q = QL[idx];
            int64_t An = p * (p - q);
            int64_t Ad = q * q;
            int64_t twopq = 2*p + q;
            for (int i = 0; i < NSMALL; i++) memset(patv[i], 0, (size_t)SP[i]);

            /* per-prime reduced constants */
            int32_t sAn[NSMALL], sA1[NSMALL], sB1[NSMALL], sA2[NSMALL], sB2[NSMALL], sG3[NSMALL];
            for (int i = 0; i < NSMALL; i++) {
                int32_t P = SP[i];
                int32_t rp = mmod(p, P), rq = mmod(q, P);
                int32_t rAn = mmod(An, P), rAd = mmod(Ad, P);
                sAn[i] = rAn;
                if (surface == 0) {
                    sA1[i] = mmod(-(int64_t)rAd * mmod(p + q,   P), P);
                    sB1[i] = mmod( (int64_t)rAn * mmod(p + 2*q, P), P);
                    sA2[i] = mmod(-(int64_t)rAd * mmod(p + 2*q, P), P);
                    sB2[i] = mmod( (int64_t)rAn * mmod(p + 3*q, P), P);
                    sG3[i] = mmod( (int64_t)rp  * rAn, P);
                } else {
                    sA1[i] = mmod( (int64_t)rq * mmod((int64_t)mmod(p + q, P) * rAd, P), P);
                    sB1[i] = mmod(-(int64_t)rq * mmod((int64_t)rp * rAn, P), P);
                    sA2[i] = mmod( (int64_t)rq * mmod((int64_t)rp * rAd, P), P);
                    sB2[i] = mmod( (int64_t)rq * mmod((int64_t)rAn * mmod(3*p + q, P), P), P);
                    sG3[i] = mmod(-(int64_t)rp * mmod((int64_t)rAn * mmod(twopq, P), P), P);
                }
            }

            for (int64_t b = 1; b <= Hk; b++) {
                /* ---- tier A: build the master bitmap over a in [0,Hk] ---- */
                int empty = 0;
                for (int64_t w = 0; w < WMASK; w++) master[w] = ~0ULL;
                /* clear the bits past Hk */
                {
                    int64_t nb = Hk + 1;
                    int64_t lastw = (nb - 1) >> 6; int32_t lastbit = (nb - 1) & 63;
                    if (lastbit != 63) master[lastw] &= (1ULL << (lastbit + 1)) - 1;
                    for (int64_t w = lastw + 1; w < WMASK; w++) master[w] = 0;
                }
                for (int i = 0; i < NSMALL && !empty; i++) {
                    int32_t P = SP[i];
                    int32_t b2 = SQT[i][(int32_t)(b % P)];
                    uint64_t *pt = patc[i] + (int64_t)b2 * PW[i];
                    if (!patv[i][b2]) {
                        int8_t lb2 = LEG[i][b2];
                        int8_t lAn = LEG[i][sAn[i]];
                        int8_t lG3 = LEG[i][sG3[i]];
                        /* accept table over t = a^2 mod P, then over x = a mod P */
                        int32_t S1 = mmod((int64_t)sB1[i] * b2, P);
                        int32_t S2 = mmod((int64_t)sB2[i] * b2, P);
                        int8_t *ac = acc[i];
                        for (int32_t t = 0; t < P; t++) {
                            int ok = (lG3 * LEG[i][S1] * lb2) >= 0 && (lAn * LEG[i][S1] * LEG[i][S2]) >= 0;
                            ac[t] = (int8_t)ok;
                            S1 += sA1[i]; if (S1 >= P) S1 -= P;
                            S2 += sA2[i]; if (S2 >= P) S2 -= P;
                        }
                        for (int64_t w = 0; w < PW[i]; w++) pt[w] = 0;
                        int32_t xm = 0;
                        for (int32_t x = 0; x < P + 64; x++) {
                            if (ac[SQT[i][xm]]) pt[x >> 6] |= 1ULL << (x & 63);
                            if (++xm == P) xm = 0;
                        }
                        patv[i][b2] = 1;
                    }
                    /* AND the periodic pattern into master */
                    int32_t s = 0; uint64_t any = 0; int32_t st = STEP[i];
                    for (int64_t w = 0; w < WMASK; w++) {
                        int64_t wi = s >> 6; int32_t off = s & 63;
                        uint64_t val = pt[wi] >> off;
                        if (off) val |= pt[wi + 1] << (64 - off);
                        uint64_t mw = master[w] & val;
                        master[w] = mw; any |= mw;
                        s += st; if (s >= P) s -= P;
                    }
                    if (!any) empty = 1;
                }
                if (empty) continue;

                /* ---- scan survivors ---- */
                for (int64_t w = 0; w < WMASK; w++) {
                    uint64_t m = master[w];
                    while (m) {
                        int64_t a = (w << 6) + __builtin_ctzll(m);
                        m &= m - 1;
                        total_tierB++;
                        /* ---- tier B: exact modular arithmetic, v1 formulas ---- */
                        int ok = 1;
                        for (int i = 0; i < NBIG && ok; i++) {
                            int32_t P = BP[i];
                            int64_t ar = a % P, br = b % P;
                            int64_t a2 = (ar * ar) % P, b2 = (br * br) % P;
                            int64_t rAn = mmod(An, P), rAd = mmod(Ad, P);
                            int64_t rp = mmod(p, P), rq = mmod(q, P);
                            int64_t vn, vd;
                            if (surface == 0) {
                                vn = (rAn * b2) % P;
                                vd = (vn - (a2 * rAd) % P + P) % P;
                            } else {
                                vn = (P - ((rAn * mmod(twopq, P)) % P) * b2 % P) % P;
                                int64_t t = ((a2 * rAd) % P + (rAn * b2) % P) % P;
                                vd = (rq * t) % P;
                            }
                            int64_t qvd = (rq * vd) % P, pvd = (rp * vd) % P, qvn = (rq * vn) % P;
                            int64_t S1 = (pvd + qvn + qvd) % P;
                            int64_t c3 = ((rp * vn) % P) * S1 % P;
                            if (!QRB[i][c3]) { ok = 0; break; }
                            int64_t c2;
                            if (surface == 0) {
                                int64_t S2 = (S1 + qvd) % P;
                                c2 = ((rAn * S1) % P) * S2 % P;
                            } else {
                                int64_t Dn = (pvd - qvn + P) % P;
                                c2 = ((rAn * Dn) % P) * S1 % P;
                            }
                            if (!QRB[i][c2]) { ok = 0; break; }
                        }
                        if (!ok) continue;
                        if (igcd(a, b) != 1) continue;      /* k not in lowest terms */

                        /* ---- tier C: exact 128-bit degeneracy screen ---- */
                        __int128 B2 = (__int128)b * b, A2 = (__int128)a * a;
                        __int128 vn, vd;
                        if (surface == 0) {
                            vn = (__int128)An * B2;
                            vd = vn - A2 * (__int128)Ad;
                        } else {
                            vn = -(__int128)An * twopq * B2;
                            vd = (__int128)q * (A2 * (__int128)Ad + (__int128)An * B2);
                        }
                        if (vd == 0 || vn == 0 || vn == vd) { total_degen++; continue; }
                        __int128 qvd = (__int128)q * vd, pvd = (__int128)p * vd, qvn = (__int128)q * vn;
                        __int128 S1 = pvd + qvn + qvd, S2 = S1 + qvd, Dn = pvd - qvn;
                        if (S1 == 0 || S2 == 0 || Dn == 0) { total_degen++; continue; }
                        total_hits++;
#pragma omp critical
                        {
                            printf("HIT %d %lld %lld %lld %lld\n", surface,
                                   (long long)p, (long long)q, (long long)a, (long long)b);
                            fflush(stdout);
                        }
                    }
                }
            }
#pragma omp atomic
            done_pq++;
            if ((idx & 8191) == 0) {
                int64_t d;
#pragma omp atomic read
                d = done_pq;
                fprintf(stderr, "PROGRESS pq_done=%lld / %lld\n", (long long)d, (long long)npq);
                fflush(stderr);
            }
        }
        free(master);
        for (int i = 0; i < NSMALL; i++) { free(patc[i]); free(patv[i]); free(acc[i]); }
    }

    /* pairs enumerated = npq * Hk * (Hk+1)   (a in [0,Hk], b in [1,Hk]) */
    printf("# TOTAL pq=%lld pairs=%lld tierB=%lld survivors=%lld degenerate_dropped=%lld\n",
           (long long)npq, (long long)(npq * Hk * (Hk + 1)),
           (long long)total_tierB, (long long)total_hits, (long long)total_degen);
    printf("SEARCH_DONE surface=%d Hu=%lld Hk=%lld qlo=%lld qhi=%lld nsmall=%d nbig=%d\n",
           surface, (long long)Hu, (long long)Hk, (long long)qlo, (long long)qhi, NSMALL, NBIG);
    return 0;
}
