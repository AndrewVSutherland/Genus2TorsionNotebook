/* opus_abc_search.c  (2026-07-19)
 *
 * Fast search for rational points on Jen's (2,2,2,12) moduli surface
 *
 *   S:  y^2 = F1(A,B,C),  z^2 = F2(A,B,C)   in P^2_{A,B,C}
 *
 * using the factorizations (verified symbolically in opus_abc_verify.m)
 *
 *   D  = A^2 - B^2
 *   W  = C^2 + D                      (shared conic factor)
 *   P1 = (B^2-A^2)(C^2-B^2) = -D(C^2-B^2)
 *   Q  = C^4 + D*C^2 + A^2*D
 *   F1 = P1*W,   F2 = Q*W
 *
 * Both sextics are even in each of A,B,C, so WLOG A,B,C >= 0.
 *
 * Funnel per (A,B): build small-prime admissibility tables indexed by
 * C mod p (cost ~sum(p), negligible), then the inner C-loop is a chain
 * of table lookups with early exit, then medium primes, then exact
 * 128-bit square tests.
 *
 * Degenerate loci (excluded by Jen's discriminant condition) are
 * skipped or flagged: A=0, B=0, C=0, |A|=|B|, |C|=|B|, W=0.
 *
 * Build: cc -O3 -march=native -o opus_abc_search opus_abc_search.c -lm
 * Usage: opus_abc_search Alo Ahi H
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef __int128 i128;
typedef int64_t i64;

#define NS 8
static const int SP[NS] = {5, 7, 11, 13, 17, 19, 23, 29};
static signed char LEG[NS][32];
static unsigned char TAB[NS][32];

#define NM 3
static const int MP[NM] = {1009, 2003, 4003};
static signed char *LEGM[NM];

static void init_leg(void)
{
    for (int i = 0; i < NS; i++) {
        int p = SP[i];
        for (int v = 0; v < p; v++) LEG[i][v] = -1;
        LEG[i][0] = 0;
        for (int r = 1; r < p; r++) LEG[i][(r * r) % p] = 1;
    }
    for (int i = 0; i < NM; i++) {
        int p = MP[i];
        LEGM[i] = (signed char *)malloc(p);
        for (int v = 0; v < p; v++) LEGM[i][v] = -1;
        LEGM[i][0] = 0;
        for (int r = 1; r < p; r++) LEGM[i][(int)(((i64)r * r) % p)] = 1;
    }
}

static inline int issquare128(i128 x, i128 *rt)
{
    if (x < 0) return 0;
    if (x == 0) { *rt = 0; return 1; }
    i128 r = (i128)sqrtl((long double)x);
    while (r > 0 && r * r > x) r--;
    while ((r + 1) * (r + 1) <= x) r++;
    if (r * r == x) { *rt = r; return 1; }
    return 0;
}

/* Jen's full discriminant condition: the genus-2 curve degenerates iff any
 * of these 16 factors vanishes.  Returns 1 if the point is degenerate. */
static int degenerate(i64 A, i64 B, i64 C)
{
    i128 a = A, b = B, c = C;
    i128 a2 = a * a, b2 = b * b, c2 = c * c;
    if (!C || !B || !A) return 1;
    if (b == c || b == -c) return 1;
    if (a == b || a == -b) return 1;
    if (a == c || a == -c) return 1;
    if (-a2 + b2 - 2 * c2 == 0) return 1;
    if (-a2 + b2 - c2 == 0) return 1;
    if (-a2 + b2 - a * c - c2 == 0) return 1;
    if (-a2 + b2 + a * c - c2 == 0) return 1;
    if (-a2 * b2 + b2 * b2 + 2 * a2 * c2 - 2 * b2 * c2 + c2 * c2 == 0) return 1;
    if (-a2 * a2 + a2 * b2 - a2 * c2 + b2 * c2 - c2 * c2 == 0) return 1;
    if (a2 * a2 - 2 * a2 * b2 + b2 * b2 + 2 * a2 * c2 - b2 * c2 == 0) return 1;
    return 0;
}

static void print_i128(char *buf, i128 v)
{
    if (v == 0) { buf[0] = '0'; buf[1] = 0; return; }
    int neg = v < 0;
    if (neg) v = -v;
    char tmp[64]; int n = 0;
    while (v > 0) { tmp[n++] = '0' + (int)(v % 10); v /= 10; }
    int k = 0;
    if (neg) buf[k++] = '-';
    while (n > 0) buf[k++] = tmp[--n];
    buf[k] = 0;
}

int main(int argc, char **argv)
{
    if (argc < 4) { fprintf(stderr, "usage: %s Alo Ahi H\n", argv[0]); return 1; }
    int Alo = atoi(argv[1]), Ahi = atoi(argv[2]), H = atoi(argv[3]);
    init_leg();

    long long nexact = 0, nhit = 0, ndeg = 0;
    for (int A = Alo; A <= Ahi; A++) {
        i64 A2 = (i64)A * A;
        for (int B = 1; B <= H; B++) {
            i64 B2 = (i64)B * B;
            i64 D = A2 - B2;
            if (D == 0) continue;                 /* |A|=|B|: degenerate line */
            i64 A2D = A2 * D;

            for (int i = 0; i < NS; i++) {
                int p = SP[i];
                int dp  = (int)(((D % p) + p) % p);
                int b2p = (int)(B2 % p);
                int adp = (int)(((A2D % p) + p) % p);
                int mdp = (p - dp) % p;           /* -D mod p */
                for (int r = 0; r < p; r++) {
                    int c  = (r * r) % p;
                    int w  = (c + dp) % p;
                    int f1 = (int)(((i64)mdp * ((c - b2p + p) % p)) % p);
                    f1 = (int)(((i64)f1 * w) % p);
                    int q  = (int)((((i64)c * c) % p + (i64)dp * c % p + adp) % p);
                    int f2 = (int)(((i64)w * q) % p);
                    TAB[i][r] = (LEG[i][f1] >= 0 && LEG[i][f2] >= 0) ? 1 : 0;
                }
            }

            /* F1 and F2 are symmetric under (A,B,C) -> (C,B,A), so we may
             * restrict to A <= C and halve the work.  Counters are seeded
             * so that r_i = C mod SP[i] after the first increment. */
            int r0 = (A - 1) % SP[0], r1 = (A - 1) % SP[1];
            int r2 = (A - 1) % SP[2], r3 = (A - 1) % SP[3];
            int r4 = (A - 1) % SP[4], r5 = (A - 1) % SP[5];
            int r6 = (A - 1) % SP[6], r7 = (A - 1) % SP[7];
            for (int C = A; C <= H; C++) {
                /* advance residue counters first */
                if (++r0 == SP[0]) r0 = 0;
                if (++r1 == SP[1]) r1 = 0;
                if (++r2 == SP[2]) r2 = 0;
                if (++r3 == SP[3]) r3 = 0;
                if (++r4 == SP[4]) r4 = 0;
                if (++r5 == SP[5]) r5 = 0;
                if (++r6 == SP[6]) r6 = 0;
                if (++r7 == SP[7]) r7 = 0;

                if (!TAB[0][r0]) continue;
                if (!TAB[1][r1]) continue;
                if (!TAB[2][r2]) continue;
                if (!TAB[3][r3]) continue;
                if (!TAB[4][r4]) continue;
                if (!TAB[5][r5]) continue;
                if (!TAB[6][r6]) continue;
                if (!TAB[7][r7]) continue;

                i64 c = (i64)C * C;
                i64 W = c + D;
                if (W == 0) continue;             /* degenerate conic */
                if (c == B2) continue;            /* |C|=|B|: degenerate */

                /* medium-prime sieve */
                int dead = 0;
                for (int i = 0; i < NM && !dead; i++) {
                    int p = MP[i];
                    int cp  = (int)(c % p);
                    int dp  = (int)(((D % p) + p) % p);
                    int b2p = (int)(B2 % p);
                    int adp = (int)(((A2D % p) + p) % p);
                    int w   = (cp + dp) % p;
                    int f1  = (int)(((i64)((p - dp) % p) * ((cp - b2p + p) % p)) % p);
                    f1 = (int)(((i64)f1 * w) % p);
                    int q   = (int)((((i64)cp * cp) % p + (i64)dp * cp % p + adp) % p);
                    int f2  = (int)(((i64)w * q) % p);
                    if (LEGM[i][f1] < 0 || LEGM[i][f2] < 0) dead = 1;
                }
                if (dead) continue;

                nexact++;
                i128 F1 = (i128)(-D) * (i128)(c - B2) * (i128)W;
                i128 rt1;
                if (!issquare128(F1, &rt1)) continue;
                i128 Qv = (i128)c * c + (i128)D * c + (i128)A2D;
                i128 F2 = Qv * (i128)W;
                i128 rt2;
                if (!issquare128(F2, &rt2)) continue;

                if (degenerate(A, B, C)) { ndeg++; continue; }

                nhit++;
                char b1[64], b2[64];
                print_i128(b1, rt1);
                print_i128(b2, rt2);
                printf("HIT A=%d B=%d C=%d y=%s z=%s\n", A, B, C, b1, b2);
                fflush(stdout);
            }
        }
        if ((A % 50) == 0) {
            fprintf(stderr, "PROGRESS A=%d exact=%lld deg=%lld hits=%lld\n", A, nexact, ndeg, nhit);
            fflush(stderr);
        }
    }
    printf("DONE Alo=%d Ahi=%d H=%d exact_tests=%lld degenerate=%lld hits=%lld\n",
           Alo, Ahi, H, nexact, ndeg, nhit);
    return 0;
}
