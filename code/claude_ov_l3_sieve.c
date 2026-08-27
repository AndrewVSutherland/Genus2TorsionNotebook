/* ===========================================================================
 * Lane 3 / route B2' : search the [1,1,2,2] + ord(D_inf)=11 locus.
 *
 * Chart (the factored realization of the brief's "the 2-structure is FREE"
 * observation -- each cubic acquires a rational root <=> f factors as):
 *
 *      y^2 = f(x) = x (x - w) (x^2 + a x + b) (x^2 + c x + d)
 *
 * monic sextic, factor type [1,1,2,2] => 2-rank 2 built in by construction.
 * The ONLY remaining condition is ord(D_inf) = 11, and then
 *      J(Q)_tors  >=  Z/2 x Z/2 x Z/11  =  [2,22].
 *
 * GATE.  If ord_Q(D_inf) = 11 then for EVERY good prime p the reduction has
 * ord_{F_p}(D_inf) = 11 EXACTLY (torsion injects).  That is a ~1.5% gate per
 * prime -- vastly sharper than "11 | #J(F_p)" (~18%).  We precompute, for
 * each small prime p, a p^5-bit table of the residue tuples (w,a,b,c,d) mod p
 * that pass; degenerate (non-squarefree mod p) residues are marked PASS so no
 * true global hit is ever lost.
 *
 * The CF order over F_p is the polynomial-continued-fraction quasi-period sum
 * (Platonov-Petrunin), skill `pell-cf-order`, validated f14->14, f18->18.
 *
 * POSITIVE CONTROLS (must be reported by --control):
 *   19044.h.2   (w,a,b,c,d) = (1, -3,   8, 4,  27)
 *   BLP C4corr  (w,a,b,c,d) = (15,-34, 160, 31,-50)
 *
 * build:  gcc -O3 -march=native -fopenmp -o claude_ov_l3_sieve claude_ov_l3_sieve.c
 * =========================================================================== */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef int32_t i32;
typedef int64_t i64;
typedef uint64_t u64;
typedef uint32_t u32;

#define PMAX 256
#define NPOLY 16

/* ---------------- polynomial arithmetic over F_p (tiny degrees) ---------- */

static inline int pdeg(i32 *A, int n) { int i; for (i = n; i >= 0; i--) if (A[i]) return i; return -1; }

/* Q = A div B, R = A mod B ; degA = da, degB = db >= 0, B[db]!=0 */
static void pdivmod(i32 *A, int da, i32 *B, int db, i32 *Q, int *dq, i32 *R, int *dr,
                    i32 p, const i32 *inv)
{
    i32 tmp[NPOLY];
    int i, k;
    for (i = 0; i < NPOLY; i++) { tmp[i] = 0; Q[i] = 0; }
    for (i = 0; i <= da; i++) tmp[i] = A[i];
    i32 ilc = inv[B[db]];
    for (k = da - db; k >= 0; k--) {
        i32 co = (i32)(((i64)tmp[k + db] * ilc) % p);
        Q[k] = co;
        if (co) {
            for (i = 0; i <= db; i++) {
                i32 v = tmp[k + i] - (i32)(((i64)co * B[i]) % p);
                v %= p; if (v < 0) v += p;
                tmp[k + i] = v;
            }
        }
    }
    for (i = 0; i < NPOLY; i++) R[i] = 0;
    for (i = 0; i < db; i++) R[i] = tmp[i];
    *dq = pdeg(Q, da - db < 0 ? 0 : da - db);
    *dr = pdeg(R, db - 1 < 0 ? 0 : db - 1);
}

/* C = A*B */
static void pmul(i32 *A, int da, i32 *B, int db, i32 *C, int *dc, i32 p)
{
    int i, j;
    for (i = 0; i < NPOLY; i++) C[i] = 0;
    if (da < 0 || db < 0) { *dc = -1; return; }
    for (i = 0; i <= da; i++) {
        if (!A[i]) continue;
        for (j = 0; j <= db; j++) {
            if (!B[j]) continue;
            C[i + j] = (i32)(((i64)C[i + j] + (i64)A[i] * B[j]) % p);
        }
    }
    *dc = pdeg(C, da + db);
}

/* squarefree test: deg gcd(f, f') == 0 */
static int is_squarefree(i32 *f, int df, i32 p, const i32 *inv)
{
    i32 A[NPOLY], B[NPOLY], Q[NPOLY], R[NPOLY];
    int da, db, dq, dr, i;
    for (i = 0; i < NPOLY; i++) { A[i] = 0; B[i] = 0; }
    for (i = 0; i <= df; i++) A[i] = f[i];
    da = pdeg(A, df);
    for (i = 1; i <= df; i++) B[i - 1] = (i32)(((i64)f[i] * (i % p)) % p);
    db = pdeg(B, df - 1);
    while (db >= 0) {
        pdivmod(A, da, B, db, Q, &dq, R, &dr, p, inv);
        for (i = 0; i < NPOLY; i++) A[i] = B[i];
        da = db;
        for (i = 0; i < NPOLY; i++) B[i] = R[i];
        db = dr;
    }
    return da <= 0;
}

/* exact order of D_inf over F_p via polynomial CF; 0 if none within budget.
   Early-abort as soon as the running total exceeds `cap`. */
static int cf_order_modp(i32 *f, i32 p, const i32 *inv, int cap, int maxsteps)
{
    i32 s[NPOLY], P[NPOLY], Q[NPOLY], ai[NPOLY], Pn[NPOLY], t[NPOLY], num[NPOLY], Qn[NPOLY], R[NPOLY];
    int i, dq, dr, dai, dPn, dt, dnum, dQn, dQ, dP, total = 0, step;
    i32 i2 = inv[2 % p];

    for (i = 0; i < NPOLY; i++) { s[i] = 0; P[i] = 0; Q[i] = 0; }
    /* s = x^3 + s2 x^2 + s1 x + s0 with f - s^2 of degree <= 2 (f monic sextic) */
    s[3] = 1;
    s[2] = (i32)(((i64)f[5] * i2) % p);
    { i64 v = (i64)f[4] - (i64)s[2] * s[2] % p; v %= p; if (v < 0) v += p;
      s[1] = (i32)((v * i2) % p); }
    { i64 v = (i64)f[3] - 2 * ((i64)s[1] * s[2] % p); v %= p; if (v < 0) v += p;
      s[0] = (i32)((v * i2) % p); }

    Q[0] = 1; dQ = 0; dP = -1;

    for (step = 0; step <= maxsteps; step++) {
        if (dQ < 0) return 0;
        /* ai = (P + s) div Q */
        for (i = 0; i < NPOLY; i++) t[i] = 0;
        for (i = 0; i <= 3; i++) t[i] = s[i];
        for (i = 0; i <= (dP < 0 ? -1 : dP); i++) { i32 v = t[i] + P[i]; if (v >= p) v -= p; t[i] = v; }
        dt = pdeg(t, 3);
        if (dt < dQ) { dai = -1; for (i = 0; i < NPOLY; i++) ai[i] = 0; }
        else pdivmod(t, dt, Q, dQ, ai, &dai, R, &dr, p, inv);
        if (dai < 0) return 0;
        total += dai;
        if (total > cap) return 0;
        /* Pn = ai*Q - P */
        pmul(ai, dai, Q, dQ, Pn, &dPn, p);
        for (i = 0; i <= (dP < 0 ? -1 : dP); i++) { i32 v = Pn[i] - P[i]; if (v < 0) v += p; Pn[i] = v; }
        dPn = pdeg(Pn, dPn > dP ? dPn : dP);
        /* num = f - Pn^2 */
        pmul(Pn, dPn, Pn, dPn, num, &dnum, p);
        for (i = 0; i < NPOLY; i++) num[i] = (i32)((p - num[i]) % p);
        for (i = 0; i <= 6; i++) { i32 v = num[i] + f[i]; if (v >= p) v -= p; num[i] = v; }
        dnum = pdeg(num, 6 > dnum ? 6 : dnum);
        if (dnum < 0) return 0;
        pdivmod(num, dnum, Q, dQ, Qn, &dQn, R, &dr, p, inv);
        if (dr >= 0) return 0;                 /* division must be exact */
        for (i = 0; i < NPOLY; i++) { P[i] = Pn[i]; Q[i] = Qn[i]; }
        dP = dPn; dQ = dQn;
        if (step >= 1 && dQ <= 0 && dQ >= 0) return total;
    }
    return 0;
}

static void build_f(i64 w, i64 a, i64 b, i64 c, i64 d, i32 p, i32 *f)
{
    i64 t4 = 1, t3 = a + c, t2 = b + d + a * c, t1 = a * d + b * c, t0 = b * d;
    i64 F[7];
    F[6] = t4;
    F[5] = t3 - w * t4;
    F[4] = t2 - w * t3;
    F[3] = t1 - w * t2;
    F[2] = t0 - w * t1;
    F[1] = -w * t0;
    F[0] = 0;
    for (int i = 0; i <= 6; i++) { i64 v = F[i] % p; if (v < 0) v += p; f[i] = (i32)v; }
}

/* 3-valued gate:  2 = DEGENERATE (f mod p not squarefree),
                   1 = OK   (ord_{F_p}(D_inf) == 11),
                   0 = FAIL (definitely not a global order-11 curve).          */
static int gate3_modp(i64 w, i64 a, i64 b, i64 c, i64 d, i32 p, const i32 *inv)
{
    i32 f[NPOLY];
    for (int i = 0; i < NPOLY; i++) f[i] = 0;
    build_f(w, a, b, c, d, p, f);
    if (pdeg(f, 6) != 6) return 2;
    if (!is_squarefree(f, 6, p, inv)) return 2;
    return cf_order_modp(f, p, inv, 11, 14) == 11 ? 1 : 0;
}

/* ------------------------------ tables ---------------------------------- */
static int NP;
static i32 PR[16];
static uint64_t *TOK[16], *TDG[16];      /* OK bitmap and DEGENERATE bitmap */
static i32 INV[16][PMAX];
static double RATE_OK[16], RATE_DG[16];

static void mkinv(i32 p, i32 *inv)
{
    inv[0] = 0; if (p > 1) inv[1] = 1;
    for (i32 i = 2; i < p; i++) inv[i] = (i32)((i64)(p - (p / i)) * inv[p % i] % p);
}

static void build_tables(void)
{
    for (int k = 0; k < NP; k++) {
        i32 p = PR[k];
        mkinv(p, INV[k]);
        u64 n = 1; for (int i = 0; i < 5; i++) n *= (u64)p;
        TOK[k] = (uint64_t *)calloc((n + 63) / 64, 8);
        TDG[k] = (uint64_t *)calloc((n + 63) / 64, 8);
        u64 cok = 0, cdg = 0;
#pragma omp parallel for schedule(static) reduction(+:cok) reduction(+:cdg)
        for (i64 w = 0; w < p; w++) {
            for (i64 a = 0; a < p; a++)
             for (i64 b = 0; b < p; b++)
              for (i64 c = 0; c < p; c++)
               for (i64 d = 0; d < p; d++) {
                   int g = gate3_modp(w, a, b, c, d, p, INV[k]);
                   if (!g) continue;
                   u64 idx = ((((u64)w * p + a) * p + b) * p + c) * p + d;
                   if (g == 1) { TOK[k][idx >> 6] |= (1ULL << (idx & 63)); cok++; }
                   else        { TDG[k][idx >> 6] |= (1ULL << (idx & 63)); cdg++; }
               }
        }
        RATE_OK[k] = (double)cok / (double)n;
        RATE_DG[k] = (double)cdg / (double)n;
        fprintf(stderr, "TABLE p=%d n=%llu  ok=%llu (%.6f)  deg=%llu (%.6f)\n",
                p, (unsigned long long)n, (unsigned long long)cok, RATE_OK[k],
                (unsigned long long)cdg, RATE_DG[k]);
        fflush(stderr);
    }
}

static inline int bitof(uint64_t *T, u64 idx) { return (T[idx >> 6] >> (idx & 63)) & 1ULL; }

/* extra on-the-fly primes for survivors */
static i32 XP[24]; static int NXP; static i32 XINV[24][PMAX];

int main(int argc, char **argv)
{
    i64 Ba = 60, Bb = 200;
    int control = 0, MAXDEG = 4, MAXDEG1 = 2;
    i64 wlo = 1, whi = 12;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--control")) control = 1;
        else if (!strncmp(argv[i], "Ba=", 3)) Ba = atoll(argv[i] + 3);
        else if (!strncmp(argv[i], "Bb=", 3)) Bb = atoll(argv[i] + 3);
        else if (!strncmp(argv[i], "wlo=", 4)) wlo = atoll(argv[i] + 4);
        else if (!strncmp(argv[i], "whi=", 4)) whi = atoll(argv[i] + 4);
        else if (!strncmp(argv[i], "maxdeg=", 7)) MAXDEG = atoi(argv[i] + 7);
        else if (!strncmp(argv[i], "maxdeg1=", 8)) MAXDEG1 = atoi(argv[i] + 8);
    }

    /* table primes, DESCENDING (best selectivity first for early exit) */
    NP = 5; PR[0] = 23; PR[1] = 19; PR[2] = 17; PR[3] = 13; PR[4] = 11;
    NXP = 0;
    { int ps[] = {29,31,37,41,43,47,53,59,61,67,71,73};
      for (unsigned i = 0; i < sizeof(ps)/sizeof(ps[0]); i++) { XP[NXP] = ps[i]; mkinv(XP[NXP], XINV[NXP]); NXP++; } }

    build_tables();

    if (control) {
        i64 ctl[2][5] = { {1, -3, 8, 4, 27}, {15, -34, 160, 31, -50} };
        const char *nm[2] = { "19044.h.2", "BLP-C4corr" };
        for (int t = 0; t < 2; t++) {
            printf("CONTROL %-11s (w,a,b,c,d)=(%lld,%lld,%lld,%lld,%lld):",
                   nm[t], (long long)ctl[t][0], (long long)ctl[t][1], (long long)ctl[t][2],
                   (long long)ctl[t][3], (long long)ctl[t][4]);
            int fail = 0, nd = 0;
            for (int k = 0; k < NP; k++) {
                i32 p = PR[k];
                i64 r[5]; for (int j = 0; j < 5; j++) { r[j] = ctl[t][j] % p; if (r[j] < 0) r[j] += p; }
                u64 idx = ((((u64)r[0] * p + r[1]) * p + r[2]) * p + r[3]) * p + r[4];
                int ok = bitof(TOK[k], idx), dg = bitof(TDG[k], idx);
                printf(" %d:%s", p, ok ? "OK" : (dg ? "DEG" : "FAIL"));
                if (!ok && !dg) fail = 1; if (dg) nd++;
            }
            for (int k = 0; k < NXP; k++) {
                int g = gate3_modp(ctl[t][0], ctl[t][1], ctl[t][2], ctl[t][3], ctl[t][4], XP[k], XINV[k]);
                printf(" %d:%s", XP[k], g == 1 ? "OK" : (g == 2 ? "DEG" : "FAIL"));
                if (!g) fail = 1; if (g == 2) nd++;
            }
            printf("  => %s (ndeg=%d)\n", (!fail && nd <= MAXDEG) ? "PASS" : "REJECT", nd);
        }
        fflush(stdout);
    }

    fprintf(stderr, "SWEEP w=[%lld..%lld] |a|,|c|<=%lld |b|,|d|<=%lld maxdeg1=%d maxdeg=%d\n",
            (long long)wlo, (long long)whi, (long long)Ba, (long long)Bb, MAXDEG1, MAXDEG);
    fflush(stderr);

    u64 total_tested = 0, hits = 0, stage1 = 0;
    double t0 = 0, t1 = 0;
#ifdef _OPENMP
    t0 = omp_get_wtime();
#endif

#pragma omp parallel for schedule(dynamic,1) collapse(2) reduction(+:total_tested) reduction(+:hits) reduction(+:stage1)
    for (i64 w = wlo; w <= whi; w++) {
      for (i64 a = -Ba; a <= Ba; a++) {
        i32 wm[16], am[16];
        for (int k = 0; k < NP; k++) { wm[k] = (i32)(w % PR[k]); if (wm[k] < 0) wm[k] += PR[k];
                                       am[k] = (i32)(a % PR[k]); if (am[k] < 0) am[k] += PR[k]; }
        for (i64 b = -Bb; b <= Bb; b++) {
            i32 bm[16];
            for (int k = 0; k < NP; k++) { bm[k] = (i32)(b % PR[k]); if (bm[k] < 0) bm[k] += PR[k]; }
            for (i64 c = a; c <= Ba; c++) {          /* dedupe quadratic swap: a <= c */
                u32 mok[16], mdg[16];
                for (int k = 0; k < NP; k++) {
                    i32 p = PR[k];
                    i32 cm = (i32)(c % p); if (cm < 0) cm += p;
                    u64 base = ((((u64)wm[k] * p + am[k]) * p + bm[k]) * p + cm) * p;
                    u32 mo = 0, md = 0;
                    for (i32 dd = 0; dd < p; dd++) {
                        if (bitof(TOK[k], base + dd)) mo |= (1u << dd);
                        else if (bitof(TDG[k], base + dd)) md |= (1u << dd);
                    }
                    mok[k] = mo; mdg[k] = md;
                }
                i64 dstart = (a == c) ? b : -Bb;     /* a==c: also need b <= d */
                i32 dm[16];
                for (int k = 0; k < NP; k++) { dm[k] = (i32)(dstart % PR[k]); if (dm[k] < 0) dm[k] += PR[k]; }
                for (i64 d = dstart; d <= Bb; d++) {
                    total_tested++;
                    int nd = 0, ok = 1;
                    for (int k = 0; k < NP; k++) {
                        u32 sh = (u32)dm[k];
                        if ((mok[k] >> sh) & 1u) continue;
                        if ((mdg[k] >> sh) & 1u) { if (++nd > MAXDEG1) { ok = 0; break; } continue; }
                        ok = 0; break;
                    }
                    if (ok) {
                        stage1++;
                        for (int k = 0; k < NXP && ok; k++) {
                            int g = gate3_modp(w, a, b, c, d, XP[k], XINV[k]);
                            if (g == 0) ok = 0;
                            else if (g == 2 && ++nd > MAXDEG) ok = 0;
                        }
                        if (ok) {
                            hits++;
#pragma omp critical
                            { printf("HIT %lld %lld %lld %lld %lld\n",
                                     (long long)w, (long long)a, (long long)b, (long long)c, (long long)d);
                              fflush(stdout); }
                        }
                    }
                    for (int k = 0; k < NP; k++) { dm[k]++; if (dm[k] == PR[k]) dm[k] = 0; }
                }
            }
        }
      }
    }
#ifdef _OPENMP
    t1 = omp_get_wtime();
#endif
    printf("SEARCH_DONE tested=%llu stage1=%llu hits=%llu wall=%.1fs w=[%lld..%lld] Ba=%lld Bb=%lld maxdeg1=%d maxdeg=%d\n",
           (unsigned long long)total_tested, (unsigned long long)stage1, (unsigned long long)hits, t1 - t0,
           (long long)wlo, (long long)whi, (long long)Ba, (long long)Bb, MAXDEG1, MAXDEG);
    return 0;
}
