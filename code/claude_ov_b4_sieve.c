/* claude_ov_b4_sieve.c -- Lane 4 (route B4): exhaustive range sieve of the
 * three 2-rank-raise loci over the Flynn (6,2) stream.
 *
 * The stream's members are in bijection with E(Q)\{O} = { nG : n in Z\{0} },
 * E = 92.a1 : Y^2 = w^3+3w^2+2w+1, G = (0,1) (rank 1, trivial torsion).
 * A 2-rank raise of the Richelot codomain requires at least one of
 *     A1(nG), D1(nG), W(nG)   to be a rational square
 * (genus 3, 3 and 6 double covers of E respectively).  If x in Q* is a square
 * then chi_p(x) != -1 for every prime p where x is a p-adic unit, so a single
 * prime whose Legendre table says '0' at n mod N_p kills n outright.
 *
 * Input: the table file produced by code/claude_ov_b4_mktables.gp
 *        block = "p N_p" then 2 (or 3 with -w) lines of N_p characters.
 * Usage: ./b4sieve <tablefile> <NMAX> [nloci]
 *
 * cc -O3 -march=native -o b4sieve claude_ov_b4_sieve.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXP 100000
#define MAXN 8192

typedef struct { int32_t p; int32_t np; char *mask[3]; } Tab;

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s tablefile NMAX [nloci]\n", argv[0]); return 1; }
    char *fn = argv[1];
    int64_t NMAX = strtoll(argv[2], NULL, 10);
    int32_t nloci = (argc > 3) ? atoi(argv[3]) : 2;

    FILE *f = fopen(fn, "r");
    if (!f) { perror("open"); return 1; }
    static Tab tab[MAXP];
    int32_t ntab = 0;
    char line[1 << 22];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#') continue;
        int32_t p, np;
        if (sscanf(line, "%d %d", &p, &np) != 2) continue;
        tab[ntab].p = p; tab[ntab].np = np;
        for (int32_t k = 0; k < nloci; k++) {
            if (!fgets(line, sizeof line, f)) { fprintf(stderr, "truncated table\n"); return 1; }
            int32_t L = (int32_t)strlen(line);
            while (L > 0 && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
            if (L != np) { fprintf(stderr, "bad mask length %d != %d at p=%d\n", L, np, p); return 1; }
            tab[ntab].mask[k] = strdup(line);
        }
        ntab++;
        if (ntab >= MAXP) break;
    }
    fclose(f);
    printf("TABLES %d primes, nloci=%d, NMAX=%lld\n", ntab, nloci, (long long)NMAX);

    int64_t span = 2*NMAX + 1;                 /* n = -NMAX .. NMAX, index n+NMAX */
    char *alive = malloc(span);
    char *lname[3] = { "A1", "D1", "W" };

    for (int32_t k = 0; k < nloci; k++) {
        memset(alive, 1, span);
        alive[NMAX] = 0;                        /* n = 0 is not a member */
        int64_t remaining = span - 1;
        for (int32_t i = 0; i < ntab && remaining > 0; i++) {
            int32_t np = tab[i].np;
            char *m = tab[i].mask[k];
            /* n = -NMAX corresponds to residue r0 = ((-NMAX) mod np) */
            int64_t r0 = ((-NMAX) % np + np) % np;
            int64_t idx = 0;
            int64_t r = r0;
            while (idx < span) {
                if (m[r] == '0') alive[idx] = 0;
                idx++;
                if (++r == np) r = 0;
            }
            /* recount every few primes */
            if ((i & 7) == 7 || i == ntab - 1) {
                remaining = 0;
                for (int64_t j = 0; j < span; j++) remaining += alive[j];
                printf("  %s after p<=%d : survivors=%lld\n", lname[k], tab[i].p, (long long)remaining);
                fflush(stdout);
            }
        }
        int64_t nsurv = 0;
        for (int64_t j = 0; j < span; j++) if (alive[j]) {
            nsurv++;
            if (nsurv <= 200) printf("  SURVIVOR %s n=%lld\n", lname[k], (long long)(j - NMAX));
        }
        printf("LOCUS %s NMAX=%lld primes=%d SURVIVORS=%lld\n",
               lname[k], (long long)NMAX, ntab, (long long)nsurv);
        fflush(stdout);
    }
    printf("SIEVE_DONE\n");
    return 0;
}
