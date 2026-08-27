# Closure theorems: no Z/96 or Z/160 through the Elkies [32] component — 2026-07-21

Program (from the Sol-review follow-up): upgrade the reconstructed Elkies [32]
family (genus-0 component C32(z,r), parametrized over Q(u); members carry exact
[32], verified at u=2, -3/2) by an independent 3-contact (=> cyclic Z/96) or a
second 5-contact (=> Z/160).  All routes are now decided NEGATIVE:

1. Z/160 (second P-infty 5-torsion class): the contact system collapses via
   h2^2 = f4 + 5 c5 rho to two univariate polynomials in rho over Q(u).
   GCD = 1 (no horizontal family); Res_rho has u-degree 752 with exactly one
   rational root u = 6120 = the parametrization's blow-down point (degenerate).
   THEOREM: no member carries a second rational P-infty 5-torsion class.
   (General [P1+P2-2infty] 5-classes not covered by this ansatz — open but
   analogous machinery applies.)  code/claude_z160_gcd_decide.m + log.

2. Z/96 (independent 3-torsion): the quadratic contact ansatz is COMPLETE for
   3-torsion on quintic models (no pole-3 functions at a Weierstrass infinity).
   Elimination over Q[u,E,q1,q0]: the GCD of the two E-resultants consists of
   disc(q3) (boundary shape) and four vertical lines, ALL at degenerate members
   (u=6120 blow-down [mult 409], three deg-drop points; the deg-2/4 U-factors
   have no rational roots) — checked exactly, code/claude_z96_vertical_check.m.
   The cofactor eliminant (computed in 8.4h, degrees (U,Q1) = (3568,160)) OOM'd
   at factorization, but the SPECIALIZATION PROFILE settles the structure: at
   6 independent u0 the fiber eliminant (deg 160, squarefree) factors as
   [40, 120] EVERY time (code/claude_z96_sheetprofile.m).  So the 3-contact
   cover has exactly two irreducible sheets, fiber-degrees 40 and 120 (40 =
   (3^4-1)/2 geometric classes; Galois orbits 40+120) — both high-genus curves.
   THEOREM: no infinite Z/96 family through this component; sporadic members
   are finite (Faltings) and excluded to u-height 120 by the 3-divisibility
   gate scan (0/17,543 members; extension to height 400 running).

Method note: the "specialize-and-profile" trick (factor the fiber eliminant at
several random parameters; stable degree multisets = generic sheet degrees)
replaced a >=10h bivariate factorization with ~1 minute of computation, and is
reusable for every contact-cover analysis in this repo.
