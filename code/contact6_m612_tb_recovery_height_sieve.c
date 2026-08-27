/*
 * Conservative projective-height sieve on the rational recovery family.
 *
 * The family depends only on u^2, so it is enough to scan primitive
 * 0 <= numerator <= H, 1 <= denominator <= H.  Denominator residues zero
 * modulo p pass conservatively (u=infinity).  For finite u, the machine mask
 * contains the union of exact extra-3-allowed and bad-reduction residues.
 *
 * CRT primes 17,79,83 give 5*19*15 finite numerator residues for a generic
 * denominator and modulus 111469 > H when H=100000.
 */

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXP 256
#define MAXMASKS 64

typedef struct {
    int p;
    unsigned char pass[2][MAXP];
    int count[2];
} Mask;

static int64_t gcd64(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) {
        int64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static int invmod(int a, int p) {
    int t = 0, nt = 1, r = p, nr = a % p;
    while (nr) {
        int q = r / nr;
        int z = t - q * nt; t = nt; nt = z;
        z = r - q * nr; r = nr; nr = z;
    }
    if (r != 1) return -1;
    t %= p;
    if (t < 0) t += p;
    return t;
}

static int find_mask(Mask *masks, int nmasks, int p) {
    for (int i = 0; i < nmasks; ++i) if (masks[i].p == p) return i;
    return -1;
}

static int load_masks(const char *path, Mask *masks) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
        exit(2);
    }
    int nmasks = 0;
    for (;;) {
        int p, branch, n;
        int got = fscanf(fp, "%d %d %d", &p, &branch, &n);
        if (got == EOF) break;
        if (got != 3 || p >= MAXP || n < 0 || n > p) {
            fprintf(stderr, "bad mask header\n");
            exit(2);
        }
        int idx = find_mask(masks, nmasks, p);
        if (idx < 0) {
            if (nmasks >= MAXMASKS) exit(2);
            idx = nmasks++;
            memset(&masks[idx], 0, sizeof(Mask));
            masks[idx].p = p;
        }
        int bi = branch == -1 ? 0 : 1;
        masks[idx].count[bi] = n;
        for (int j = 0; j < n; ++j) {
            int z;
            if (fscanf(fp, "%d", &z) != 1 || z < 0 || z >= p) exit(2);
            masks[idx].pass[bi][z] = 1;
        }
    }
    fclose(fp);
    return nmasks;
}

int main(int argc, char **argv) {
    int H = argc > 1 ? atoi(argv[1]) : 100000;
    const char *mask_path = argc > 2
        ? argv[2] : "data/contact6_m612_tb_recovery_masks_p97.machine.txt";
    if (H <= 0) return 2;

    Mask masks[MAXMASKS];
    memset(masks, 0, sizeof(masks));
    int nmasks = load_masks(mask_path, masks);

    const int cp[3] = {17, 79, 83};
    int ci[3];
    int64_t mod = 1;
    for (int j = 0; j < 3; ++j) {
        ci[j] = find_mask(masks, nmasks, cp[j]);
        if (ci[j] < 0) return 2;
        mod *= cp[j];
    }
    int64_t crtcoef[3];
    for (int j = 0; j < 3; ++j) {
        int64_t mj = mod / cp[j];
        crtcoef[j] = mj * invmod((int)(mj % cp[j]), cp[j]);
    }

    printf("CONTACT6_M612_TB_RECOVERY_HEIGHT_SIEVE H=%d modulus=%" PRId64
           " masks=%d\n", H, mod, nmasks);

    for (int bi = 0; bi < 2; ++bi) {
        int branch = bi == 0 ? -1 : 1;
        uint64_t crt_tuples = 0, interval_candidates = 0;
        uint64_t primitive = 0, survivors = 0;

        for (int d = 1; d <= H; ++d) {
            int residues[3][MAXP];
            int nr[3] = {0, 0, 0};
            for (int j = 0; j < 3; ++j) {
                int p = cp[j];
                int dm = d % p;
                if (dm == 0) {
                    for (int z = 0; z < p; ++z) residues[j][nr[j]++] = z;
                } else {
                    for (int u = 0; u < p; ++u) {
                        if (masks[ci[j]].pass[bi][u])
                            residues[j][nr[j]++] = (dm * u) % p;
                    }
                }
            }

            for (int i0 = 0; i0 < nr[0]; ++i0)
            for (int i1 = 0; i1 < nr[1]; ++i1)
            for (int i2 = 0; i2 < nr[2]; ++i2) {
                ++crt_tuples;
                int64_t n0 = (residues[0][i0] * crtcoef[0]
                            + residues[1][i1] * crtcoef[1]
                            + residues[2][i2] * crtcoef[2]) % mod;
                /* The family is even in u, so scan only 0 <= n <= H. */
                for (int64_t n = n0; n <= H; n += mod) {
                    ++interval_candidates;
                    if (gcd64(n, d) != 1) continue;
                    ++primitive;
                    int ok = 1;
                    for (int mi = 0; mi < nmasks; ++mi) {
                        int p = masks[mi].p;
                        int dm = d % p;
                        if (dm == 0) continue; /* projective infinity */
                        int im = invmod(dm, p);
                        int nm = (int)(n % p);
                        int u = (nm * im) % p;
                        if (!masks[mi].pass[bi][u]) {
                            ok = 0;
                            break;
                        }
                    }
                    if (!ok) continue;
                    ++survivors;
                    printf("SURVIVOR branch=%d n=%" PRId64 " d=%d\n",
                           branch, n, d);
                }
            }
            if (d % 10000 == 0)
                fprintf(stderr, "branch %d d=%d survivors=%" PRIu64 "\n",
                        branch, d, survivors);
        }
        printf("BRANCH_DONE branch=%d crt=%" PRIu64
               " interval=%" PRIu64 " primitive=%" PRIu64
               " survivors=%" PRIu64 "\n",
               branch, crt_tuples, interval_candidates, primitive, survivors);
    }
    return 0;
}
