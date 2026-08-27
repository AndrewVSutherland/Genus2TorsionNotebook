# Ari's (A,B,C) surface for (2,2,2,12): verification, THIRD CURVE, and the multigrade theorem

Session 2026-07-20. Input: Ari's email giving the moduli surface of the two known
(2,2,2,12) curves as a double-cover of P^2_{A,B,C} (two simultaneous square conditions),
with (s,m,n)-formulas back to the M(2,2,2,6) chart, and the "(2,2,2,12) scratch" section
of Rational_torsion_points_on_genus_2_Jacobians.tex (halving P_b = (-b,0)-infty via the
four square-conditions S1..S4; S3- and S2-conditions parametrized away; the two leftover
conditions = Ari's surface).

## 1. Verification (code/claude_ari_surface_verify.m, results/claude_ari_surface_verify.log)

- Email formulas s,m,n(A,B,C) agree with the document chain for n and m
  (the document's s-line has a typo or differs by normalization; the email versions are the
  correct ones, verified end-to-end below).
- S1(s,m,n(A,B,C)) / F = square in Q(A,B,C); S2, S3 become squares identically;
  S4/G = square TIMES F (so S4-square <=> G-square on the surface). All checked in Magma.
- Both sextics FACTOR (huge simplification vs "generic K3 double covers"):
    F = (B^2-A^2)(B^2-C^2)(B^2-A^2-C^2)                       [4 lines + Pythagorean conic]
    G = (B^2-A^2-C^2)(B^2(A^2+C^2) - (A^4+A^2C^2+C^4))        [same conic x quartic]
  The two branch sextics SHARE the conic B^2 = A^2+C^2; the 4 lines are tangent to it at
  (1:+-1:0),(0:1:+-1), and the quartic has double contact with it at the SAME 4 points.
  The quartic's only singular point is (0:1:0).
- Known curve #2 = (A:B:C) = (120:218:143) (t = 65/88 on the conic pencil below);
  known curve #1 = (408:437:143) (t = 13/187). Both: (s:m:n) match recorded values
  projectively, quintic model G2Invariants match the recorded minimal models exactly.

## 2. Geometry: pencil structure and the special fibers (both rank 0)

Pencil of conics through the four tangency points: M_t: B^2 = A^2 + tAC + C^2
(t = pencil parameter; every surface point lies on the member with t = (B^2-A^2-C^2)/(AC)).
Pulled back to X, the generic member carries a genus-5 curve (so each fiber has finitely
many rational points, Faltings). Special members:

- t = 0: the shared branch conic (fully degenerate).
- t = +-1: members B^2 = A^2 +- AC + C^2 lie ENTIRELY on X (y = AC(A+C), z = ACB) but are
  inside the discriminant (factor B^2-A^2-+AC-C^2) -- the "DEGEN flood" in the sieve.
- t = +-2: member degenerates to a line pair B = +-(A+-C); on the line B = A+C the fiber is
  the genus-1 curve E2: y'^2 = 2(2A+C)(A+2C), z'^2 = 2(2A^2+AC+2C^2).
- t = +-1/2: the z-cover splits identically (tA^2+AC+tC^2 = ((A+C)/2)^2 at t=1/2); fiber is
  the genus-1 curve Eh: quartic y''^2 = 2(A(u)+2C(u))(2A(u)+C(u)) over the member conic.

RESULT (Magma, unconditional RankBounds): E2 and Eh are BOTH the elliptic curve
[1,0,1,-19,26] of conductor 30, rank 0, torsion Z/2 x Z/6; all their rational points are
degenerate surface points. So the two natural genus-1 loci give NO infinite family.
Lines on X: the T_iT_j lines either degenerate or give E2 again; all pairs of line-arrangement
nodes lie on the arrangement lines themselves (degenerate). No other line-type curves.

## 3. Sieve over (A,B,C) (code/claude_abc_sieve.c, data/claude_abc_sieve_h20000.txt)

Positive coprime triples, A > C (A<->C is curve-preserving), B in the two sign-allowed
ranges; layered QR filters + exact 128-bit square tests. H = 20000 in 760 s on 30 threads.
Result: 650 GOOD hits = 13 primitive points = EXACTLY 3 curves:

- curve #2 orbit (6 reps): (133,109,60), (143,218,120), (143,241,120), (266,218,143),
  (266,241,120), (266,241,143)
- curve #1 orbit (6 reps): (408,437,143), (408,1013,143), (1015,437,143), (1015,437,408),
  (1015,1013,143), (1015,1013,408)
- NEW third curve (1 rep <= 20000): (5364, 19661, 4165), t = 320/21

This explains Ari's "up to automorphisms only the two": at height 500 there are 7 points
but they all belong to the two known orbits.

## 4. THIRD CURVE (code/claude_ari_new3_validate.m, results/claude_ari_new3_validate.log)

(A:B:C) = (5364:19661:4165), (s:m:n) = (436901823 : -11212618846 : 920719872).
Reduced minimal model:
  y^2 + (x^2+x)y = 3703062294195264x^6 - 360079374491052216x^5 + 8901721379573296848x^4
    - 5397945250386334945x^3 - 86737535708373850908x^2 + 36346694984390901540x
    + 43035470132681030400
- TorsionSubgroup EXACTLY [2,2,2,12] (order 96), 0.4 s on the minimized model.
- Simplicity certificates at p = 37, 127, 131, 179 (Frobenius charpoly irreducible AND
  12th-power transform irreducible of degree 4).
- G2-invariants distinct from curves #1,#2 (so not even a twist).
- Cross-check: integral orbit-mate (21456, 27593, 16660) gives the same curve.

## 5. THE MULTIGRADE THEOREM (the big structural discovery)

Empirically each curve's orbit is {three X-values} x {two B-values}; the (A,C) run over
2-subsets of {X1,X2,X3} and B over {B1,B2}. The algebraic characterization, PROVEN
symbolically (scratchpad idcheck; identities below):

  THEOREM. Let x_i = X_i^2, u = B1^2, v = B2^2. On the variety V defined by
      X1^2+X2^2+X3^2 = B1^2+B2^2   and   X1^4+X2^4+X3^4 = B1^4+B2^4
  (equivalently {u,v} are the roots of T^2 - e1(x)T + e2(x)), the identities
      F(X1,B1,X2)*(u-x3)^2 = x1^2*x2^2*x3      and      G(X1,B1,X2) = (x3-v)^2*v
  hold IDENTICALLY. Hence both square conditions are automatic and every point of V
  (off the degeneracy divisors) yields a (2,2,2,12) curve. Conversely a surface point
  determines (x3, v) rationally (linear solve), and X3 = y(u-x3)/(X1X2), B2 = z/(x3-v)
  are forced rational; so V is birational to Ari's surface X.

  In words: THIS COMPONENT OF M(2,2,2,12) IS the classical multigrade variety of
  "equal sums of 2nd and 4th powers, three terms vs two terms" -- a (2,4) complete
  intersection surface in P^4 with S3 x S2 x signs symmetry. The 6 representations
  = choices of 2-subset of X's and one B. Solutions always have X3 > X1+X2
  (disc = prod(X1 +- X2 +- X3) > 0), and the degenerate B1=B2 locus is exactly
  X3 = X1+X2 with B^2 = X1^2+X1X2+X2^2 (the t=-+1 conic family).

Known orbits:
  curve #2: X = {120, 143, 266},  B = {218, 241}     (both sums: 105605; 4th: 5631933137)
  curve #1: X = {143, 408, 1015}, B = {437, 1013}
  curve #3: X = {16660, 21456, 78793}, B = {78644, 27593}  (= {4165,5364,78793/4;19661,27593/4})
Note curve #3's primitive X-triple has a denominator 4 in the (A,B,C)-chart -- why the
(A,B,C)-sieve saw only one representative.

## 6. Orbit-level sieve (code/claude_sym_sieve.c, data/claude_sym_sieve_h50000.txt)

Since X3 > X1+X2, enumerate the two SMALL values X1 < X2 <= H2 and recover X3 from
divisor pairs of 4*x1*x2 (dd = (x3-x1-x2)^2 - 4x1x2 = w^2 gives (D1-w)(D1+w) = 4x1x2);
X3 is then essentially unbounded. H2 = 50000, X3 <= 2*10^6, ~4 min:
ONLY the three known orbits (plus scalings). The orbit census in this huge region is
complete: exactly 3 curves.

POST-REVIEW FIX (Codex P1 on PR #5, 2026-07-20): the original binary computed
N = 4*x1*x2 in uint64, which wraps once X1*X2 > 2^31 — about 1.04% of the H2=50000
pair region (13.1M pairs). Fixed by keeping N in unsigned __int128, capping divisor
generation at sqrt(N) = 2*X1*X2 exactly (lossless: a product > sqrt(N) only has larger
multiples), and taking the cofactor N/r in 128 bits with an early skip when it already
exceeds 2*X3MAX^2. The full H2=50000, X3MAX=2e6 rerun with the fixed binary produced
BYTE-IDENTICAL output (473 hits, same 3 primitive orbits; the old run also had zero
hits in the affected corner, so nothing spurious had been reported). The census claim
stands, now on sound arithmetic. Points on this moduli surface are extremely sparse --
consistent with a (2,4) CI in P^4 being of general type (K = O(1) on the smooth model
before singular corrections) and Bombieri-Lang behavior: rational points confined to
(rank-0 / degenerate) special curves plus finitely many sporadics.

## 7. Other components of M(2,2,2,12)

The document's component halves P_b. The divisibility criterion (calibrated against
curve #2; x-T convention with x0 = -e literal: {-e} cup {r-e : r other roots} all squares)
gives sibling surfaces for halving P_0, P_a, P_c, P_d. A scan of all five point-type
criteria over (s:m:n) height <= 60 finds ZERO points (code/claude_ari_other_components.m).
The 10 pair-type components ([(-e,0)+(-e',0)-2infty] divisible) are UNEXPLORED -- deriving
their multigrade-style models is the natural next step; they may be entirely different
surfaces with their own (possibly positive-rank) special curves.

## 8. Status of the infinitude question

- The two obvious genus-1 special fibers: rank 0 (unconditional). Dead.
- All fiber-type curves of the conic pencil otherwise have genus 5. Lines fully classified.
- Empirics (3 orbits in a range ~100x beyond the knowns) support general-type sparseness.
- Remaining routes: (i) known parametric families of the multigrade system in the classical
  literature (Gloden, Choudhry "Symmetric Diophantine systems"; the mixed-size (3,2) case
  of degrees (2,4) -- a literature dig for the team); (ii) rational/elliptic curves on V
  found through its rich singular locus / trivial-solution lines (X1=0, {X2,X3}={B1,B2});
  (iii) the 10 pair-type sibling components; (iv) elliptic fibrations with section on the
  two branch K3s individually, matched over P^2.
- A single nondegenerate parametric multigrade family = infinitude of (2,2,2,12) curves
  (up to the standard generic-simplicity argument; the family's moduli map is nonconstant
  because n = A/B varies).

## Files
- code/claude_ari_surface_verify.m, code/claude_ari_newpoints_check.m,
  code/claude_ari_identify.m, code/claude_ari_new3_validate.m,
  code/claude_ari_other_components.m, code/claude_ari_orbit_structure.m
- code/claude_abc_sieve.c (H=20000 run), code/claude_sym_sieve.c (orbit-level)
- results/claude_ari_*.log, data/claude_abc_sieve_h20000.txt, data/claude_sym_sieve_h50000.txt
- data/claude_ari_curve3.txt (third-curve record)
