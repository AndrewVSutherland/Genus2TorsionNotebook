/* Lane 8 (overnight 2026-07-25): exhaustive genus-2 "local landscape" over F_p.
 *
 * Enumerates EVERY genus-2 curve y^2 = f(x) over F_p with deg f in {5,6} and
 * leading coefficient in {1, nu} (nu a fixed nonsquare) -- i.e. every curve up to
 * the quadratic twist/scaling y -> c y -- and histograms the Frobenius data
 * (s1, s2) where chi(T) = T^4 - s1 T^3 + s2 T^2 - p s1 T + p^2.
 *
 *   #J(F_p) = p^2 + 1 - (p+1) s1 + s2 .
 *
 * The histogram is postprocessed (claude_ov_lane8_landscape.gp) to answer, for a
 * target torsion order N: which isogeny classes over F_p have N | #J, and are any
 * of them (geometrically) SIMPLE.  This is the local step-0 for the HLP split
 * targets Z/63, Z/45, [7,7], [5,10], Z/35 -- a curve over Q with a rational
 * N-torsion point and good reduction at p must reduce into one of these classes.
 *
 * build: gcc -O3 -march=native -o landscape claude_ov_lane8_landscape.c
 * usage: ./landscape <p>          (writes "s1 s2 count" lines to stdout)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int32_t p;
static int32_t nu;              /* fixed nonsquare mod p */
static int8_t  chi[64];         /* quadratic character on F_p, chi[0]=0        */
static int8_t *Q2;              /* Q2[m*(2p-1) + A] = 1 + chi[(A*A - m) mod p] */
static int8_t *Q1;              /* Q1[A] = 1 + chi[A]  (A in 0..2p-2, folded)  */

static inline int32_t md(int32_t a) { a %= p; if (a < 0) a += p; return a; }

/* polynomial gcd degree over F_p; returns 1 iff f is squarefree */
static int squarefree(int32_t *c, int32_t deg)
{
    int32_t a[8], b[8], da, db, i, j;
    for (i = 0; i <= deg; i++) a[i] = c[i];
    da = deg;
    for (i = 0; i < deg; i++) b[i] = md((int64_t)(i + 1) * c[i + 1] % p);
    db = deg - 1;
    while (db >= 0 && b[db] == 0) db--;
    if (db < 0) return 0;                       /* f' == 0  =>  not squarefree */
    while (db >= 0) {
        /* a mod b */
        int32_t inv = 1, e = p - 2, base = b[db];
        while (e) { if (e & 1) inv = (int32_t)((int64_t)inv * base % p); base = (int32_t)((int64_t)base * base % p); e >>= 1; }
        while (da >= db) {
            int32_t co = (int32_t)((int64_t)a[da] * inv % p);
            if (co) for (j = 0; j <= db; j++)
                a[da - db + j] = md(a[da - db + j] - (int32_t)((int64_t)co * b[j] % p));
            a[da] = 0;
            da--;
            while (da >= 0 && a[da] == 0) da--;
            if (da < db) break;
        }
        /* swap */
        for (j = 0; j < 8; j++) { int32_t t = a[j]; a[j] = b[j]; b[j] = t; }
        { int32_t t = da; da = db; db = t; }
        while (db >= 0 && b[db] == 0) db--;
    }
    return (da <= 0);                            /* gcd constant => squarefree */
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s p\n", argv[0]); return 1; }
    p = atoi(argv[1]);
    int32_t p2 = p * p, i, j, k;

    for (i = 0; i < p; i++) chi[i] = 0;
    for (i = 1; i < p; i++) chi[(int32_t)((int64_t)i * i % p)] = 1;
    for (i = 1; i < p; i++) if (!chi[i]) chi[i] = -1;
    chi[0] = 0;
    for (nu = 2; nu < p; nu++) if (chi[nu] == -1) break;

    Q2 = malloc((size_t)p * (2 * p - 1));
    for (int32_t m = 0; m < p; m++)
        for (int32_t A = 0; A < 2 * p - 1; A++) {
            int32_t Am = A % p;
            Q2[m * (2 * p - 1) + A] = 1 + chi[md((int32_t)((int64_t)Am * Am % p) - m)];
        }
    Q1 = malloc(2 * p - 1);
    for (int32_t A = 0; A < 2 * p - 1; A++) Q1[A] = 1 + chi[A % p];

    /* enumerate F_{p^2} = F_p[i]/(i^2 - nu) as z = a + b i, index j = a*p + b */
    int32_t *za = malloc(p2 * sizeof(int32_t)), *zb = malloc(p2 * sizeof(int32_t));
    for (i = 0; i < p; i++) for (j = 0; j < p; j++) { za[i * p + j] = i; zb[i * p + j] = j; }

    /* histogram over (s1, s2): s1 in [-4sqrt p, 4sqrt p], s2 in [-2p*..,..] */
    int32_t S1 = 40, S2 = 200;   /* generous half-ranges, p <= 23 */
    int64_t *hist = calloc((size_t)(2 * S1 + 1) * (2 * S2 + 1), sizeof(int64_t));
    int64_t ncurves = 0, nsing = 0;

    int32_t gA[1024], gB[1024], hA[1024], nuB2[1024];
    int32_t c[8];

    for (int32_t topdeg = 6; topdeg >= 5; topdeg--) {
        int32_t lcs[2] = {1, nu};
        for (int32_t li = 0; li < 2; li++) {
            memset(c, 0, sizeof(c));
            c[topdeg] = lcs[li];
            int32_t nfree = topdeg - 1;       /* c[topdeg-1] .. c[1] outer, c[0] inner */
            int64_t ntup = 1; for (i = 0; i < nfree; i++) ntup *= p;
            for (int64_t t = 0; t < ntup; t++) {
                int64_t tt = t;
                for (i = 1; i <= nfree; i++) { c[i] = (int32_t)(tt % p); tt /= p; }
                /* g(z) = sum_{k>=1} c_k z^k  (c0 added in the inner loop) */
                for (j = 0; j < p2; j++) {
                    int32_t A = 0, B = 0;
                    for (k = topdeg; k >= 1; k--) {
                        int32_t nA = md((int32_t)(((int64_t)A * za[j] + (int64_t)nu * B % p * zb[j]) % p));
                        int32_t nB = md((int32_t)(((int64_t)A * zb[j] + (int64_t)B * za[j]) % p));
                        A = md(nA + c[k]); B = nB;
                    }
                    {   /* one more multiply by z: the Horner loop above stops at c_1 */
                        int32_t nA = md((int32_t)(((int64_t)A * za[j] + (int64_t)nu * B % p * zb[j]) % p));
                        int32_t nB = md((int32_t)(((int64_t)A * zb[j] + (int64_t)B * za[j]) % p));
                        A = nA; B = nB;
                    }
                    hA[j] = A;
                    nuB2[j] = (int32_t)((int64_t)nu * B % p * B % p);
                }
                for (int32_t c0 = 0; c0 < p; c0++) {
                    c[0] = c0;
                    if (!squarefree(c, topdeg)) { nsing++; continue; }
                    int32_t N1 = (topdeg == 6) ? (1 + chi[c[6]]) : 1;
                    int32_t N2 = (topdeg == 6) ? 2 : 1;
                    for (j = 0; j < p2; j++)
                        N2 += Q2[nuB2[j] * (2 * p - 1) + hA[j] + c0];
                    for (int32_t xx = 0; xx < p; xx++) {  /* b = 0 slice: j = xx*p */
                        N1 += Q1[hA[xx * p] + c0];
                    }
                    int32_t s1 = p + 1 - N1;
                    int32_t s2n = N2 - p2 - 1 + s1 * s1;
                    if (s2n & 1) { fprintf(stderr, "parity error\n"); return 2; }
                    int32_t s2 = s2n / 2;
                    if (s1 < -S1 || s1 > S1 || s2 < -S2 || s2 > S2) {
                        fprintf(stderr, "range overflow s1=%d s2=%d\n", s1, s2); return 3;
                    }
                    hist[(size_t)(s1 + S1) * (2 * S2 + 1) + (s2 + S2)]++;
                    ncurves++;
                }
            }
        }
    }
    fprintf(stderr, "p=%d nu=%d curves=%lld singular_skipped=%lld\n",
            p, nu, (long long)ncurves, (long long)nsing);
    printf("# p %d curves %lld\n", p, (long long)ncurves);
    for (int32_t s1 = -S1; s1 <= S1; s1++)
        for (int32_t s2 = -S2; s2 <= S2; s2++) {
            int64_t v = hist[(size_t)(s1 + S1) * (2 * S2 + 1) + (s2 + S2)];
            if (v) printf("%d %d %lld\n", s1, s2, (long long)v);
        }
    return 0;
}
