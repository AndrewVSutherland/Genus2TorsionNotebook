/*
 *  claude_ov_2216_sweep.c   -- lane OV/2216
 *
 *  Conic-parametrized sweep of the norm surfaces attached to halving an
 *  order-8 class on the odd M_1(8,2,2) chart (target [2,2,16], order 64).
 *
 *  ------------------------------------------------------------------
 *  DERIVATION (code/claude_ov_2216_symbolic.m; the underlying x-T module
 *  code/claude_ov_2216_delta.m is validated 372/372 against Magma's
 *  IsDivisibleBy in results/claude_ov_2216_delta_validate.log)
 *
 *  Odd model f_{u,v}, marked class Q of order 8, L = Q^3 x K:
 *      delta'(Q) = ( u(u-1), v(v-1), (u+v+1)(u+v+2), theta+1 ).
 *  Q+T is 2-divisible  <=>  exists d in Q* with delta'_1,2,3 in d(Q*)^2
 *  and delta'_K in d(K*)^2  <=>  delta'_1 delta'_2, delta'_1 delta'_3
 *  squares in Q and delta'_1 delta'_K a square in K.
 *
 *  The EIGHT classes Q+T (T in J[2](Q)) give only FOUR distinct condition
 *  sets: 4Q is the quadratic-pair 2-torsion class and is ALWAYS
 *  2-divisible, so T and T+4Q give the same condition.  The three
 *  nontrivial ones form a single S_3 orbit (S_3 permutes the triple
 *  {u, v, w}, w = -(u+v+1), of non-normalized sextic roots).  So there
 *  are exactly TWO genuinely distinct norm surfaces:
 *
 *   Y0  (T = 0 and T = 4Q):
 *        (1)  u(u-1) * v(v-1)            = square
 *        (2)  u(u-1) * (u+v+1)(u+v+2)    = square
 *   Y1  (T = the class of the Weierstrass point u; v- and w- twists are
 *        the S_3 images):
 *        (1) -u(u-1) * v * (2u+v+1)      = square
 *        (2)  u(u-1) * (u-v) * (u+v+1)   = square
 *
 *  and in BOTH cases the field condition delta'_1 delta'_K in (K*)^2
 *  forces the purely rational NORM condition, which is the SAME for all
 *  eight classes:
 *        (3)  u * v * (u+v+1)            = square
 *  (this is exactly "the y-coordinate of the marked point is a square").
 *
 *  ------------------------------------------------------------------
 *  METHOD.  For FIXED u, put A = u(u-1).  Condition (1) is a conic in
 *  (v,y) through the rational point (0,0), so the pencil y = k v
 *  parametrizes ALL of its rational points:
 *        Y0:  v =  A / (A - k^2)
 *        Y1:  v = -A (2u+1) / (k^2 + A)
 *  Conditions (2) and (3) then remain.  With u = p/q, k = a/b in lowest
 *  terms, as integers:
 *        An = p(p-q),  Ad = q^2
 *   Y0:  vn =  An b^2,            vd = An b^2 - a^2 Ad
 *   Y1:  vn = -An (2p+q) b^2,     vd = q (a^2 Ad + An b^2)
 *        S1 = p*vd + q*vn + q*vd        (numerator of u+v+1 over q*vd)
 *        S2 = S1 + q*vd                 (numerator of u+v+2 over q*vd)
 *        Dn = p*vd - q*vn               (numerator of u-v   over q*vd)
 *  and, square factors q^2 and (q*vd)^2 dropping out,
 *        (2)  Y0: An*S1*S2 square      Y1: An*Dn*S1 square
 *        (3)      p*vn*S1 square
 *
 *  Stage 1 (this program) applies these as EXACT necessary conditions
 *  modulo NPRIME small primes -- all residues are computed from p,q,a,b
 *  by integer arithmetic mod the prime, so no overflow is possible and
 *  no true solution can ever be discarded.  Degenerate parameters are
 *  removed exactly in 128-bit arithmetic.  Survivors are verified in
 *  exact rational arithmetic downstream (claude_ov_2216_exact.py).
 *
 *  Since k ranges over ALL rationals of height <= Hk and covers every
 *  rational point of each conic fiber, the sweep is exhaustive for
 *  height(u) <= Hu and height(k) <= Hk.
 *
 *  Usage:  ./claude_ov_2216_sweep <surface 0|1> <Hu> <Hk> <nthreads>
 *  Output: "HIT <surface> <p> <q> <a> <b>"   (u = p/q, k = a/b)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>

#define NPRIME 24
static int32_t PR[NPRIME] = {
    211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,
    277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349
};
static uint8_t *QRTAB[NPRIME];

static void build_qr(void)
{
    for (int i = 0; i < NPRIME; i++) {
        int32_t p = PR[i];
        QRTAB[i] = (uint8_t *)calloc((size_t)p, 1);
        QRTAB[i][0] = 1;                       /* 0 is allowed */
        for (int32_t x = 1; x < p; x++)
            QRTAB[i][(int32_t)(((int64_t)x * x) % p)] = 1;
    }
}

static inline int64_t igcd(int64_t a, int64_t b)
{
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

int main(int argc, char **argv)
{
    if (argc < 5) {
        fprintf(stderr, "usage: %s <surface 0|1> <Hu> <Hk> <nthreads>\n", argv[0]);
        return 1;
    }
    int surface = atoi(argv[1]);
    int64_t Hu = atoll(argv[2]);
    int64_t Hk = atoll(argv[3]);
    int nthreads = atoi(argv[4]);
    omp_set_num_threads(nthreads);
    build_qr();

    printf("# claude_ov_2216_sweep surface=%d Hu=%lld Hk=%lld threads=%d nprime=%d\n",
           surface, (long long)Hu, (long long)Hk, nthreads, NPRIME);
    fflush(stdout);

    /* coprime (a,b) list, a>=0 (k and -k give the same v) */
    int64_t nk = 0;
    for (int64_t b = 1; b <= Hk; b++)
        for (int64_t a = 0; a <= Hk; a++)
            if (igcd(a, b) == 1) nk++;
    int32_t *KA = (int32_t *)malloc(sizeof(int32_t) * (size_t)nk);
    int32_t *KB = (int32_t *)malloc(sizeof(int32_t) * (size_t)nk);
    {
        int64_t i = 0;
        for (int64_t b = 1; b <= Hk; b++)
            for (int64_t a = 0; a <= Hk; a++)
                if (igcd(a, b) == 1) { KA[i] = (int32_t)a; KB[i] = (int32_t)b; i++; }
    }
    fprintf(stderr, "# k-list size %lld\n", (long long)nk);

    int64_t total_pairs = 0, total_hits = 0, total_degen = 0;

#pragma omp parallel for schedule(dynamic,1) reduction(+:total_pairs,total_hits,total_degen)
    for (int64_t q = 1; q <= Hu; q++) {
        int64_t loc_pairs = 0, loc_hits = 0, loc_degen = 0;
        for (int64_t p = -Hu; p <= Hu; p++) {
            if (p == 0 || p == q) continue;              /* u = 0 or 1 */
            if (igcd(p, q) != 1) continue;
            int64_t An = p * (p - q);
            int64_t Ad = q * q;
            int64_t twopq = 2 * p + q;
            if (surface == 1 && twopq == 0) continue;    /* v == 0 identically */

            int32_t Anr[NPRIME], Adr[NPRIME], pr_[NPRIME], qr_[NPRIME], tpq[NPRIME];
            for (int i = 0; i < NPRIME; i++) {
                int32_t P = PR[i];
                Anr[i] = (int32_t)(((An    % P) + P) % P);
                Adr[i] = (int32_t)(((Ad    % P) + P) % P);
                pr_[i] = (int32_t)(((p     % P) + P) % P);
                qr_[i] = (int32_t)(((q     % P) + P) % P);
                tpq[i] = (int32_t)(((twopq % P) + P) % P);
            }

            for (int64_t j = 0; j < nk; j++) {
                int64_t a = KA[j], b = KB[j];
                loc_pairs++;
                int ok = 1;
                for (int i = 0; i < NPRIME; i++) {
                    int32_t P = PR[i];
                    int64_t ar = a % P, br = b % P;
                    int64_t a2 = (ar * ar) % P;
                    int64_t b2 = (br * br) % P;
                    int64_t vn, vd;
                    if (surface == 0) {
                        vn = ((int64_t)Anr[i] * b2) % P;
                        vd = (vn - (a2 * Adr[i]) % P + P) % P;
                    } else {
                        vn = (P - (((int64_t)Anr[i] * tpq[i]) % P) * b2 % P) % P;
                        int64_t t = ((a2 * Adr[i]) % P + ((int64_t)Anr[i] * b2) % P) % P;
                        vd = ((int64_t)qr_[i] * t) % P;
                    }
                    int64_t qvd = ((int64_t)qr_[i] * vd) % P;
                    int64_t pvd = ((int64_t)pr_[i] * vd) % P;
                    int64_t qvn = ((int64_t)qr_[i] * vn) % P;
                    int64_t S1  = (pvd + qvn + qvd) % P;
                    /* condition (3): p * vn * S1 */
                    int64_t c3 = (((int64_t)pr_[i] * vn) % P) * S1 % P;
                    if (!QRTAB[i][c3]) { ok = 0; break; }
                    /* condition (2) */
                    int64_t c2;
                    if (surface == 0) {
                        int64_t S2 = (S1 + qvd) % P;
                        c2 = (((int64_t)Anr[i] * S1) % P) * S2 % P;
                    } else {
                        int64_t Dn = (pvd - qvn + P) % P;
                        c2 = (((int64_t)Anr[i] * Dn) % P) * S1 % P;
                    }
                    if (!QRTAB[i][c2]) { ok = 0; break; }
                }
                if (!ok) continue;

                /* exact degeneracy screen in 128-bit */
                __int128 b2 = (__int128)b * b, a2 = (__int128)a * a;
                __int128 vn, vd;
                if (surface == 0) {
                    vn = (__int128)An * b2;
                    vd = vn - a2 * (__int128)Ad;
                } else {
                    vn = -(__int128)An * twopq * b2;
                    vd = (__int128)q * (a2 * (__int128)Ad + (__int128)An * b2);
                }
                if (vd == 0 || vn == 0 || vn == vd) { loc_degen++; continue; }
                __int128 qvd = (__int128)q * vd;
                __int128 pvd = (__int128)p * vd;
                __int128 qvn = (__int128)q * vn;
                __int128 S1 = pvd + qvn + qvd;
                __int128 S2 = S1 + qvd;
                __int128 Dn = pvd - qvn;
                if (S1 == 0 || S2 == 0 || Dn == 0) { loc_degen++; continue; }

                loc_hits++;
#pragma omp critical
                {
                    printf("HIT %d %lld %lld %lld %lld\n", surface,
                           (long long)p, (long long)q, (long long)a, (long long)b);
                }
            }
        }
        total_pairs += loc_pairs;
        total_hits  += loc_hits;
        total_degen += loc_degen;
        if (q % 25 == 0)
            fprintf(stderr, "PROGRESS q=%lld pairs=%lld hits=%lld degen=%lld\n",
                    (long long)q, (long long)loc_pairs, (long long)loc_hits,
                    (long long)loc_degen);
    }

    printf("# TOTAL pairs=%lld survivors=%lld degenerate_dropped=%lld\n",
           (long long)total_pairs, (long long)total_hits, (long long)total_degen);
    printf("SEARCH_DONE surface=%d Hu=%lld Hk=%lld nprime=%d\n",
           surface, (long long)Hu, (long long)Hk, NPRIME);
    return 0;
}
