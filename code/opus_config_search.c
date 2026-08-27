/* opus_config_search.c  (2026-07-19)
 *
 * CONFIGURATION search for (2,2,2,12) curves in Jen's family.
 *
 * Structure (discovered from the orbits of the two known curves, and
 * verified exactly for both): a curve corresponds to a configuration
 *
 *     (x1,x2,x3 ; b1,b2)   with   x1^2+x2^2+x3^2 = b1^2+b2^2 = N
 *
 * where each x_l is a rational point on the SINGLE elliptic curve
 *
 *     E_{b1,b2}:  w^2 = (x^2-b1^2)(x^2-b2^2)
 *                 <=>  Y^2 = X(X-b1^2)(X-b2^2),  X = x^2, Y = x w
 *
 * The 12 surface points (A,B,C) = (x_i, b_k, x_j), i!=j, k=1,2, are the
 * 12 known representations of the curve (S_3 x Z/2).
 *
 * KEY CONSEQUENCE: since X_l = x_l^2 > 0 and sum X_l = N, every x_l is
 * bounded by sqrt(N).  So for each (b1,b2) the x-search is FINITE, and
 * scanning b1<b2<=NB is complete for every configuration with
 * b1,b2 <= NB -- no unbounded Mordell-Weil enumeration is needed.
 *
 * Also reports the distribution of #{x : x on E, x < sqrt(N)}, which
 * drives the heuristic for how many such curves exist.
 *
 * Build: cc -O3 -o opus_config_search opus_config_search.c -lm
 * Usage: opus_config_search b2lo b2hi NB
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

typedef int64_t i64;

static unsigned char sq63[63], sq65[65], sq11[11];

static void init_filters(void)
{
    memset(sq63, 0, sizeof sq63);
    memset(sq65, 0, sizeof sq65);
    memset(sq11, 0, sizeof sq11);
    for (int r = 0; r < 63; r++) sq63[(r * r) % 63] = 1;
    for (int r = 0; r < 65; r++) sq65[(r * r) % 65] = 1;
    for (int r = 0; r < 11; r++) sq11[(r * r) % 11] = 1;
}

static inline int is_square(i64 x, i64 *root)
{
    if (x < 0) return 0;
    if (!((0x0202021202030213ULL >> (x & 63)) & 1)) return 0;
    if (!sq63[x % 63]) return 0;
    if (!sq65[x % 65]) return 0;
    if (!sq11[x % 11]) return 0;
    i64 r = (i64)sqrtl((long double)x);
    while (r > 0 && r * r > x) r--;
    while ((r + 1) * (r + 1) <= x) r++;
    if (r * r == x) { if (root) *root = r; return 1; }
    return 0;
}

#define MAXPTS 4096

int main(int argc, char **argv)
{
    if (argc < 4) { fprintf(stderr, "usage: %s b2lo b2hi NB\n", argv[0]); return 1; }
    int b2lo = atoi(argv[1]), b2hi = atoi(argv[2]), NB = atoi(argv[3]);
    (void)NB;
    init_filters();

    long long histo[16];
    memset(histo, 0, sizeof histo);
    long long npairs = 0, ntriples = 0;
    i64 xs[MAXPTS];

    for (int b2 = b2lo; b2 <= b2hi; b2++) {
        i64 B2 = (i64)b2 * b2;
        for (int b1 = 1; b1 < b2; b1++) {
            i64 B1 = (i64)b1 * b1;
            i64 N = B1 + B2;
            npairs++;
            int n = 0;
            /* x must satisfy x^2 < N since the three X_l are positive */
            int xmax = (int)sqrtl((long double)(N - 2));
            for (int x = 1; x <= xmax; x++) {
                i64 X = (i64)x * x;
                if (X == B1 || X == B2) continue;      /* degenerate */
                i64 g = (X - B1) * (X - B2);
                if (is_square(g, NULL)) {
                    if (n < MAXPTS) xs[n++] = X;
                }
            }
            histo[n < 15 ? n : 15]++;
            if (n < 3) continue;
            /* look for triples of X-values summing to N */
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++) {
                    i64 r = N - xs[i] - xs[j];
                    if (r <= 0) break;
                    for (int k = j + 1; k < n; k++) {
                        if (xs[k] == r) {
                            i64 r1, r2, r3;
                            is_square(xs[i], &r1);
                            is_square(xs[j], &r2);
                            is_square(xs[k], &r3);
                            ntriples++;
                            printf("CONFIG b1=%d b2=%d  x=(%lld,%lld,%lld)  N=%lld\n",
                                   b1, b2, (long long)r1, (long long)r2, (long long)r3,
                                   (long long)N);
                            fflush(stdout);
                        }
                    }
                }
        }
        if ((b2 % 200) == 0) {
            fprintf(stderr, "PROGRESS b2=%d pairs=%lld triples=%lld\n", b2, npairs, ntriples);
            fflush(stderr);
        }
    }
    printf("DONE b2 in [%d,%d] pairs=%lld triples=%lld\n", b2lo, b2hi, npairs, ntriples);
    printf("HISTO #points-on-E : count\n");
    for (int i = 0; i < 16; i++)
        if (histo[i]) printf("   %2d : %lld\n", i, histo[i]);
    return 0;
}
