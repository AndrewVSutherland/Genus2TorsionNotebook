// t5surf.c — production sweep of the T5 hit-surface S of M(2,2,2,6) in the
// Pythagorean chart (u,g):  rho' = (2gu + (u-1)(g^2+1)) / (2g*q(u)),
// q(u) = 4u^2-6u+3.  ALL (2,2,2,12) T5-hits lie on S (product identity), and
// on S only THREE conditions are independent (derived+validated in surf_id.gp):
//   W1 = (u-1)*q*[(u-1)g^2+2ug+(u-1)]      square
//   W2 = q*[g^2+(6-8u)g+1]                 square
//   W4 = (2u-1)*q*[(2u-1)g^2+(6-4u)g+(2u-1)] square
// (then W3 == W1*W2*W4 mod squares is automatic, positivity included).
// u = p/q0 (gcd 1, q0>0, u != 0,1,1/2), g = a/b (gcd 1, b>0, |a|>=b: g~1/g).
// Integer-cleared (denominators are perfect squares), QQ = 4p^2-6pq0+3q0^2 > 0:
//   Y1 = (p-q0)*QQ*((p-q0)a^2 + 2p*ab + (p-q0)b^2)
//   Y2 = QQ*q0*(q0*(a^2+b^2) + (6q0-8p)ab)     [bracket denom q0*b^2: odd q0]
//   Y4 = (2p-q0)*QQ*((2p-q0)a^2 + (6q0-4p)ab + (2p-q0)b^2)
// HIT = all three positive squares.  NEAR2 = exactly two.  rho' = rn/rd with
//   rn0 = q0*(2abp + (p-q0)(a^2+b^2)), rd0 = 2ab*QQ, reduced, rd>0.
// Known hit must appear at (p,q0,a,b) = (-97,48,725,288) and (133,145,-725,288).
// Usage: t5surf N G [nthreads] [pqlo pqhi]   (optional q0 chunk [pqlo,pqhi])
//        t5surf --test p q0 a b
//        t5surf --selftest
// Compile: gcc -O3 -march=native -fopenmp -o t5surf t5surf.c -lm
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>

typedef __int128 i128;
static uint8_t sq256[256];
static uint8_t sq45045[45045];
static uint8_t sq97[97], sq101[101];
static void init_tabs(void) {
    for (int i = 0; i < 256; i++) sq256[(i * i) & 255] = 1;
    for (int64_t i = 0; i < 45045; i++) sq45045[(i * i) % 45045] = 1;
    for (int i = 0; i < 97; i++) sq97[(i * i) % 97] = 1;
    for (int i = 0; i < 101; i++) sq101[(i * i) % 101] = 1;
}
static int is_sq128(i128 x) {
    if (x <= 0) return 0;
    if (!sq256[(uint64_t)x & 255]) return 0;
    if (!sq45045[(uint64_t)(x % 45045)]) return 0;
    if (!sq97[(uint64_t)(x % 97)]) return 0;
    if (!sq101[(uint64_t)(x % 101)]) return 0;
    long double xd = (long double)x;
    int64_t r = (int64_t)sqrtl(xd);
    for (int64_t t = r - 2; t <= r + 2; t++)
        if (t >= 0 && (i128)t * t == x) return 1;
    return 0;
}
static int64_t gcd64(int64_t x, int64_t y) {
    if (x < 0) x = -x; if (y < 0) y = -y;
    while (y) { int64_t t = x % y; x = y; y = t; }
    return x;
}
// coprime bitmask table: row b (1..G), bit a (1..G) set iff gcd(a,b)==1
static uint64_t *coptab; static int64_t Gmax; static size_t words_per_row;
static void build_coptab(int64_t G) {
    Gmax = G; words_per_row = (size_t)(G + 64) / 64 + 1;
    coptab = calloc((size_t)(G + 1) * words_per_row, 8);
    if (!coptab) { fprintf(stderr, "coptab alloc failed\n"); exit(1); }
    for (int64_t b = 1; b <= G; b++) {
        uint64_t *row = coptab + (size_t)b * words_per_row;
        for (size_t w = 0; w < words_per_row; w++) row[w] = ~0ULL;
        int64_t m = b;
        for (int64_t pr = 2; pr * pr <= m; pr++) if (m % pr == 0) {
            while (m % pr == 0) m /= pr;
            for (int64_t x = pr; x <= G; x += pr) row[x >> 6] &= ~(1ULL << (x & 63));
        }
        if (m > 1) for (int64_t x = m; x <= G; x += m) row[x >> 6] &= ~(1ULL << (x & 63));
    }
}
static inline int cop(int64_t a, int64_t b) {
    return (coptab[(size_t)b * words_per_row + (a >> 6)] >> (a & 63)) & 1;
}
// evaluate the three conditions; returns bitmask Y1|Y2<<1|Y4<<2 of squareness
static inline int ymask(int64_t p, int64_t q0, int64_t a, int64_t b, int nmin) {
    int64_t QQ = 4*p*p - 6*p*q0 + 3*q0*q0;
    int64_t e1 = p - q0, e4 = 2*p - q0;
    int64_t ab = a * b, a2 = a * a, b2 = b * b;
    i128 t2 = (i128)q0 * (a2 + b2) + (i128)(6*q0 - 8*p) * ab;
    i128 Y2 = (i128)QQ * q0 * t2;
    int m = 0, cnt = 0;
    if (is_sq128(Y2)) { m |= 2; cnt++; }
    i128 t1 = (i128)e1 * a2 + (i128)(2*p) * ab + (i128)e1 * b2;
    i128 Y1 = (i128)e1 * QQ * t1;
    if (is_sq128(Y1)) { m |= 1; cnt++; }
    if (cnt >= nmin - 1) {   // Y4 can still matter (nmin=2: log NEAR2)
        i128 t4 = (i128)e4 * a2 + (i128)(6*q0 - 4*p) * ab + (i128)e4 * b2;
        i128 Y4 = (i128)e4 * QQ * t4;
        if (is_sq128(Y4)) m |= 4;
    }
    return m;
}
static void report(char *tag, int m, int64_t p, int64_t q0, int64_t a, int64_t b) {
    // rho' = rn0/rd0 reduced
    i128 rn0 = (i128)q0 * ((i128)2*a*b*p + (i128)(p - q0) * ((i128)a*a + (i128)b*b));
    i128 rd0 = (i128)2*a*b * (4*p*p - 6*p*q0 + 3*q0*q0);
    if (rd0 < 0) { rd0 = -rd0; rn0 = -rn0; }
    long long rn = (long long)rn0, rd = (long long)rd0;   // fits: ~6 N^2 G^2
    long long gg = gcd64(rn, rd); if (gg) { rn /= gg; rd /= gg; }
    printf("%s mask=%d p=%lld q=%lld a=%lld b=%lld rn=%lld rd=%lld\n",
           tag, m, (long long)p, (long long)q0, (long long)a, (long long)b, rn, rd);
    fflush(stdout);
}
// ---- reference X-forms from t5sweep.c (for --selftest cross-validation) ----
static int test_one_X(int64_t p, int64_t q, int64_t rn, int64_t rd) {
    int64_t QQ = 4*p*p - 6*p*q + 3*q*q;
    int64_t C3 = 16*p*p*p*p - 40*p*p*p*q + 40*p*p*q*q - 18*p*q*q*q + 3*q*q*q*q;
    int64_t C2 = -16*p*p*p*p + 32*p*p*p*q - 28*p*p*q*q + 10*p*q*q*q - q*q*q*q;
    int64_t C1 = (8*p*p*p - 12*p*p*q + 10*p*q*q - 3*q*q*q)*q;
    int64_t C0 = (-2*p + q)*q*q*q;
    i128 X1 = (i128)rn * ((i128)QQ*rn - (i128)(2*p-q)*q*rd);
    i128 X2 = (i128)QQ*rn*rn - (i128)(4*p*p-4*p*q+2*q*q)*rn*rd + (i128)(2*p-q)*q*rd*rd;
    i128 inner = ((i128)C3*rn*rn*rn) + ((i128)C2*rn*rn)*rd + ((i128)C1*rn)*rd*rd
                 + (i128)C0*rd*rd*rd;
    i128 X3 = (i128)rn * inner;
    i128 X4 = (i128)C3*rn*rn + (i128)(-16*p*p*p+28*p*p*q-18*p*q*q+4*q*q*q)*q*rn*rd
              + (i128)(2*p-q)*(2*p-q)*q*q*rd*rd;
    return (is_sq128(X1)?1:0) | (is_sq128(X2)?2:0) | (is_sq128(X3)?4:0) | (is_sq128(X4)?8:0);
}
int main(int argc, char *argv[]) {
    init_tabs();
    if (argc >= 6 && !strcmp(argv[1], "--test")) {
        int64_t p = atoll(argv[2]), q0 = atoll(argv[3]), a = atoll(argv[4]), b = atoll(argv[5]);
        int m = ymask(p, q0, a, b, 2);
        report(m == 7 ? "HIT" : "PT", m, p, q0, a, b);
        return 0;
    }
    if (argc >= 2 && !strcmp(argv[1], "--selftest")) {
        // cross-validate Y-masks vs X-masks on a grid of small points
        int64_t bad = 0, tot = 0;
        for (int64_t p = -12; p <= 12; p++) for (int64_t q0 = 1; q0 <= 9; q0++) {
            if (p == 0 || p == q0 || 2*p == q0 || gcd64(p, q0) != 1) continue;
            for (int64_t b = 1; b <= 7; b++) for (int64_t a = -13; a <= 13; a++) {
                if (a == 0 || labs(a) < b || gcd64(a, b) != 1) continue;
                int64_t QQ = 4*p*p - 6*p*q0 + 3*q0*q0;
                i128 rn0 = (i128)q0 * ((i128)2*a*b*p + (i128)(p - q0)*(a*a + b*b));
                i128 rd0 = (i128)2*a*b*QQ;
                if (rd0 < 0) { rd0 = -rd0; rn0 = -rn0; }
                int64_t rn = (int64_t)rn0, rd = (int64_t)rd0;
                int64_t gg = gcd64(rn, rd); if (gg) { rn /= gg; rd /= gg; }
                if (rn == 0 || rn == rd) continue;       // rho' = 0,1 degenerate
                int mx = test_one_X(p, q0, rn, rd);
                int my = ymask(p, q0, a, b, 2);
                // mapping: Y1~X1, Y2~X2, Y4~X4; Y-forms valid where nondegenerate
                int mxr = (mx & 1) | (mx & 2) | ((mx >> 1) & 4);
                tot++;
                if (mxr != (my & 7)) {
                    // Y4 may be untested (my bit 2 unset due to early exit): recheck
                    int fully = ymask(p, q0, a, b, 0);
                    if (mxr != (fully & 7)) {
                        bad++;
                        printf("MISMATCH p=%lld q=%lld a=%lld b=%lld X=%d Yfull=%d\n",
                               (long long)p, (long long)q0, (long long)a, (long long)b, mx, fully);
                    }
                }
            }
        }
        printf("selftest: %lld points, %lld mismatches\n", (long long)tot, (long long)bad);
        return 0;
    }
    if (argc >= 5 && !strcmp(argv[1], "--gfiber")) {
        // fixed g = a/b, sweep u = p/q0 with q0 <= N: 1-dim deep scan
        int64_t a = atoll(argv[2]), b = atoll(argv[3]), N = atoll(argv[4]);
        int nth = argc > 5 ? atoi(argv[5]) : 1;
        omp_set_num_threads(nth);
        int64_t near2 = 0, hits = 0;
        #pragma omp parallel for schedule(dynamic, 16) reduction(+:near2,hits)
        for (int64_t q0 = 1; q0 <= N; q0++)
            for (int64_t p = -N; p <= N; p++) {
                if (p == 0 || p == q0 || 2*p == q0 || gcd64(p, q0) != 1) continue;
                int m = ymask(p, q0, a, b, 2);
                int cnt = (m & 1) + ((m >> 1) & 1) + ((m >> 2) & 1);
                if (cnt >= 2) {
                    #pragma omp critical(out)
                    report(cnt == 3 ? "HIT" : "NEAR2", m, p, q0, a, b);
                    if (cnt == 3) hits++; else near2++;
                }
            }
        fprintf(stderr, "GFIBER DONE g=%lld/%lld N=%lld near2=%lld hits=%lld\n",
                (long long)a, (long long)b, (long long)N, (long long)near2, (long long)hits);
        return 0;
    }
    if (argc >= 5 && !strcmp(argv[1], "--ufiber")) {
        // fixed u = p/q0, sweep g = a/b, |a| >= b, b <= G: 1-dim deep scan
        int64_t p = atoll(argv[2]), q0 = atoll(argv[3]), G = atoll(argv[4]);
        int nth = argc > 5 ? atoi(argv[5]) : 1;
        omp_set_num_threads(nth);
        int64_t near2 = 0, hits = 0;
        #pragma omp parallel for schedule(dynamic, 16) reduction(+:near2,hits)
        for (int64_t b = 1; b <= G; b++)
            for (int64_t aa = b; aa <= G; aa++) {
                if (gcd64(aa, b) != 1) continue;
                for (int64_t as = -1; as <= 1; as += 2) {
                    int64_t a = as * aa;
                    int m = ymask(p, q0, a, b, 2);
                    int cnt = (m & 1) + ((m >> 1) & 1) + ((m >> 2) & 1);
                    if (cnt >= 2) {
                        #pragma omp critical(out)
                        report(cnt == 3 ? "HIT" : "NEAR2", m, p, q0, a, b);
                        if (cnt == 3) hits++; else near2++;
                    }
                }
            }
        fprintf(stderr, "UFIBER DONE u=%lld/%lld G=%lld near2=%lld hits=%lld\n",
                (long long)p, (long long)q0, (long long)G, (long long)near2, (long long)hits);
        return 0;
    }
    if (argc < 3) { puts("usage: t5surf N G [nthreads] [q0lo q0hi]"); return 1; }
    int64_t N = atoll(argv[1]), G = atoll(argv[2]);
    int nth = argc > 3 ? atoi(argv[3]) : 2;
    int64_t qlo = argc > 5 ? atoll(argv[4]) : 1, qhi = argc > 5 ? atoll(argv[5]) : N;
    omp_set_num_threads(nth);
    build_coptab(G);
    int64_t hits = 0, near2 = 0, pass1 = 0;
    double t0 = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:hits,near2,pass1)
    for (int64_t q0 = qlo; q0 <= qhi; q0++) {
        for (int64_t p = -N; p <= N; p++) {
            if (p == 0 || p == q0 || 2*p == q0 || gcd64(p, q0) != 1) continue;
            int64_t QQ = 4*p*p - 6*p*q0 + 3*q0*q0;
            int64_t e1 = p - q0, e4 = 2*p - q0;
            int64_t K2 = 6*q0 - 8*p, K1 = 2*p, K4 = 6*q0 - 4*p;
            for (int64_t b = 1; b <= G; b++) {
                uint64_t *row = coptab + (size_t)b * words_per_row;
                int64_t b2 = b * b;
                for (int64_t as = -1; as <= 1; as += 2) {
                    for (int64_t aa = b; aa <= G; aa++) {
                        if (!((row[aa >> 6] >> (aa & 63)) & 1)) continue;
                        int64_t a = as * aa;
                        int64_t a2 = aa * aa, ab = a * b;
                        // Y2 quick reject (cheapest)
                        i128 t2 = (i128)q0 * (a2 + b2) + (i128)K2 * ab;
                        i128 Y2 = (i128)QQ * q0 * t2;
                        int cnt = 0, m = 0;
                        if (Y2 > 0 && sq256[(uint64_t)Y2 & 255] && is_sq128(Y2)) { m |= 2; cnt++; }
                        i128 t1 = (i128)e1 * a2 + (i128)K1 * ab + (i128)e1 * b2;
                        i128 Y1 = (i128)e1 * QQ * t1;
                        if (Y1 > 0 && sq256[(uint64_t)Y1 & 255] && is_sq128(Y1)) { m |= 1; cnt++; }
                        if (!cnt) continue;
                        pass1++;
                        i128 t4 = (i128)e4 * a2 + (i128)K4 * ab + (i128)e4 * b2;
                        i128 Y4 = (i128)e4 * QQ * t4;
                        if (is_sq128(Y4)) { m |= 4; cnt++; }
                        if (cnt >= 2) {
                            #pragma omp critical(out)
                            report(cnt == 3 ? "HIT" : "NEAR2", m, p, q0, a, b);
                            if (cnt == 3) hits++; else near2++;
                        }
                    }
                }
            }
        }
        if (q0 % 25 == 0)
            fprintf(stderr, "  q0=%lld done (%.1fs)\n", (long long)q0, omp_get_wtime() - t0);
    }
    fprintf(stderr, "DONE N=%lld G=%lld q0=[%lld,%lld]: 1cond-passes=%lld near2=%lld hits=%lld time=%.1fs\n",
            (long long)N, (long long)G, (long long)qlo, (long long)qhi,
            (long long)pass1, (long long)near2, (long long)hits, omp_get_wtime() - t0);
    return 0;
}
