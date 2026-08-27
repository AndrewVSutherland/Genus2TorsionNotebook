# Cyclic gaps {31,35,37,38} via CF/Pell exact-order scan: CALIBRATED COLD (2026-07-31)

Lane 2B of the 2026-07-31 campaign (primary lane: [8,8]).  Goal: find a monic
squarefree sextic whose class at infinity `D_inf` has exact order in
{31,35,37,38} — such a hit IS a rational point of that order, closing an open
cyclic gap directly.  Tool: `CFOrder` (polynomial Pell / Platonov–Petrunin),
exact and Jacobian-free; see `.claude/skills/pell-cf-order`.

## What was run

`code/cf_cyclic_gaps_scan.m` — new scan with a 3-stage funnel that makes wide
boxes affordable:

1. CF order mod p1 (first good prime of 101,103,107,109) — exact in J(F_p);
   a global order-N class reduces to exact order N at EVERY good prime, so no
   true hit can be lost;
2. CF order mod p2 must agree (kills mod-p coincidences);
3. exact CF order over Q for survivors.

Self-test (mandatory) passed: f14->14, f18->18, f28-trap->7, over Q and mod
101/103.  Boxes: monic sextics, c5 in [0..3] (translation/negation orbit
reps), |ci| <= H otherwise.

| Run | curves | 2-prime survivors | exact hits order 25..40 | log |
|---|---|---|---|---|
| H=3 | 67,228 | 1 | **0** (1 near-miss) | results/cf_cyclic_gaps_H3.log |
| H=5 | 644,204 | 3 | **0** (3 near-misses) | results/cf_cyclic_gaps_H5.log |

(near-miss = mod-p orders agree but no rational quasi-period within budget —
exactly the expected false-positive channel of the funnel.)

## Verdict: box scanning cannot see these orders — STOP

The calibration rule (plan §Lane 2B): if a box yields zero hits at the
REALIZED orders 25–30, it cannot be trusted to see 31+.  Both boxes yield
zero.  This independently reconfirms the recorded 2026-07-03 findings
(`notes/agent_a12_224_descent.md`): grid H<=8 on the structured 2-rank-2
family gave NO D_inf order >= 10, generic members have infinite-order D_inf,
and the PP order-14/18 examples were CONSTRUCTED, not scanned.  The order-N
CF locus has enormous height for N in the 20s already; for 31–38 the scan
route is dead at any accessible height.  The backward construction
(prescribed CF period) is a separately recorded dead end — do not restart
without a new structural idea.

## Z/31 twist probe: closed by literature, no computation needed

Elkies' N=31 example (g2_tors page): y^2 = 5x^6-4x^5+20x^4-2x^3+24x^2+20x+5
has a rational 31-SUBGROUP whose Galois action multiplies by {1,5,25} mod 31
— a character of ORDER 3 (cubic).  A quadratic twist can only fix a
quadratic character, so **no twist of this curve rationalizes a generator**.
(J appears RM-by-sqrt(2), isogenous to factor 245H of J0(245).)  The planned
`z31_elkies_subgroup_probe.m` is therefore unnecessary; its question is
answered negatively.  Alessandrì–Coppola's conjecture still PREDICTS 31 and
37 possible for GL_2-type surfaces — but realizing them needs a modular /
Eisenstein-congruence construction (find weight-2 newforms f with quadratic
Hecke field and an Eisenstein congruence mod 31/37, à la 245H mod 31), not a
height scan.  That is a genuine research lane, not attempted here.

## Where this leaves the cyclic gaps

- **Z/35**: only viable route remains A_1(7) (see
  notes/claude_prod_04_35_a17_chart.md, same date) — blocked on literature
  access, demoted per the session decision rule.
- **Z/31, Z/37**: modular/Eisenstein route only (above).  Z/38 = 2x19: no
  construction idea on record for a 19-part.
- The CF machinery itself remains validated and in place
  (code/cf_cyclic_gaps_scan.m self-tests in 2 s) for any future structured
  family that carries a plausible large-order D_inf.
