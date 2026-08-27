// (3,12) production 2026-07-18: new points found on the S12 z=-5/4 fiber via
// its genus-2 quotient (see notes/claude_prod_01_312.md, scratchpad t312prod/t54).
//
// ============ A. FOURTH EXACT [3,12] REALIZATION (SPLIT) ============
// (z, r) = (-5/4, -32/65) on the M(2,12) chart — NOT in the previously known
// hit list {(-5/3,-3/5), (-5/3,-9/35), (-5/4,-32/35)}; found as the second
// r-value over x=3 on the quotient (the involution partner of the known hit
// r=-32/35, exactly parallel to the (5,+-144) pair on the z=-5/3 quotient).
// Extra 3-torsion contact triples (U, V, M), ALL with M a rational square:
//   (-91/132, 845/7744, 352836/28561),  352836/28561 = (594/169)^2
//   (-211/396, 4849/69696, 2916),       2916 = 54^2
//   (-113/396, 455/23232, 88209),       88209 = 297^2
// (contact identity f5 + M q^3 = perfect square verified for all three).
P<x> := PolynomialRing(Rationals());
C4th := HyperellipticCurve(52*x^6 + 156*x^5 - 1043*x^4 - 2346*x^3 - 629*x^2 + 570*x + 225);
// TorsionSubgroup(Jacobian(C4th)) = [3, 12]  (exact, Magma).
// G2-invariants differ from all three known split hits (new curve).
// SPLIT: automorphism group of order 4; elliptic quotients (minimal models):
//   E1: [1, 0, 0, -471900, 124722000],  conductor 4290, torsion Z/12
//   E2: [1, 0, 1, -28098, -1802744],    conductor 4290, torsion Z/6
// disc primes of C4th: [2, 3, 5, 11, 13]; conductor divisible by 10 as the
// universal (3,12) bad-reduction lemma requires.
//
// ============ B. SECOND SIMPLE NEAR-MISS (Q(zeta_3) 3-class) ============
// (z, r) = (-5/4, 811171/648100), from the quotient point pair with
// x = -79/4 (analogue of (-16,+-3615) at z=-5/3).
// W-root -1018841/4666320, U = -5094205/52533756,
// M = -6388415554323/2500 = -3 * (1459271/50)^2  => extra 3-class rational
// exactly over Q(sqrt(-3)) = Q(zeta_3);  V = 6482743551025/2759795519467536.
// TorsionSubgroup = Z/12 over Q; GEOMETRICALLY SIMPLE: L_p irreducible AND
// 12th-power transform irreducible at p = 29, 67, 71, 83, 97, 101, 103
// (7 primes).  disc primes: [2,3,5,11,1493,3793,6481,31081,132661,811171].
// (Integral model = IntegralModel of y^2 = f5(z0,r0); reproduce with
//  t54/split54.m.)
//
// ============ C. PATTERN / CONJECTURE ============
// Both non-split rational points of S12 found this session (z=-5/3 and
// z=-5/4 fibers) have M in -3 * (Q^*)^2: the extra 3-torsion class is
// rational exactly over Q(zeta_3).  CONJECTURE: for rational points of S12,
// the square class of M is always 1 or -3 (mu_3/cyclotomic constraint via
// the Weil-pairing isotropy defining S12).  If true, the (3,12) carrier
// double cover m^2 = M of S12 is obstructed exactly by the -3 twist class.
