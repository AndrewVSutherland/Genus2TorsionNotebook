// (3,12) production session 2026-07-18: z = -5/3 fiber of the S12 carrier
// surface (M(2,12) chart) — genus-4 curve over Q, its genus-2 quotient, and
// the new non-split rational point found via the quotient.
// Companion file: claude_prod_01_312_fiber_z53.m (the degree-12 model G(W;r),
// W = (r+1) U, verified mod 10007 (full match) and mod 31013 (6 points)).
//
// [1] Genus of the fiber over Q: 4 (Magma FunctionField/Genus, 376 s).
// [2] Unique nontrivial subfield of its degree-12 function field over Q(r):
//     degree 4 over Q(r), genus 2 => the genus-4 curve has a degree-3
//     (etale, by Riemann-Hurwitz 6 = 3*2 + 0) map to a genus-2 curve E.
Rq<r> := PolynomialRing(Rationals()); Pq<W> := PolynomialRing(Rq);
g4 := W^4 + (250000*r^2 - 534400*r + 141760)*W^3 + (-47616000000*r^5 - 23844960000*r^4 - 75102720000*r^3 + 152998656000*r^2 - 63305472000*r + 7247731200)*W^2 + (-762649600000000*r^8 - 12175372800000000*r^7 + 173627328000000*r^6 + 9137798896000000*r^5 + 13729258416000000*r^4 - 21413537235200000*r^3 + 9880603293440000*r^2 - 2263263772800000*r + 273913207168000)*W - 236421376000000000000*r^10 - 855324224000000000000*r^9 + 948588229299200000000*r^8 + 1026687349888000000000*r^7 - 702517992555520000000*r^6 - 1116044748712960000000*r^5 + 1245182645699584000000*r^4 - 520019902944256000000*r^3 + 151147240667750400000*r^2 - 39846584978022400000*r + 5396543556812800000;
// [3] Hyperelliptic models of E (Magma IsHyperelliptic on the plane model
//     g4(W;r) = 0, then minimization):
//     minimal model  Hm : y^2 + (x^2+1) y = x^5 + 61x^4 - 86x^3 - 624x^2 - 5657x + 18400
//     simplified     Hs : y^2 = 4x^5 + 245x^4 - 344x^3 - 2494x^2 - 22628x + 73601
//     bad primes of Hs-Jacobian: [2, 3, 5, 53629]
// [4] Known rational points of Hs up to search bound 1e5 (6 points):
//     infinity, (5, +-144), (11/4, 0)  [Weierstrass], (-16, +-3615)
//     images of the two split-hit triples: r=-3/5 -> (5,-144); r=-9/35 -> (5,144)
//     (each known split triple is exactly ONE fiber of the degree-3 map).
// [5] Mordell-Weil group of Jac(Hs) PROVED by Magma MordellWeilGroupGenus2:
//     MW = Z/2 x Z x Z  (rank exactly 2, generators proved; 143 s).
//     Height pairing of D1=[(5,144)-inf], D2=[(-16,3615)-inf]:
//     det = 9.456... (independent).  Rank = genus = 2: classical Chabauty
//     is NOT applicable; a Mordell-Weil sieve is the remaining route.
// [6] Correspondence between r (on the base P^1 of the fiber) and x (on Hs):
//     Phi(r, x) of bidegree (2, 4); over x=-16 the two places have
//     r = 10161/6025 and r = infinity (chart boundary); over x=11/4:
//     r = -81/800 (double); over x=5: r = -3/5, -9/35; over x=inf: r = inf.
// [7] NEW RATIONAL POINT of the genus-4 fiber (first non-split point of S12):
//     r0 = 10161/6025, W-root -1549/3615, i.e. U0 = -7745/48558;
//     contact data M = -16767142144/16875 (NOT a square => the extra
//     3-torsion class is only Q(sqrt(M))-rational), V = 227028025/37726069824.
// [8] The corresponding curve at (z, r) = (-5/3, 10161/6025), integral model:
P5<x5> := PolynomialRing(Rationals());
Cnew := HyperellipticCurve(325305196451086346597305440*x5^5 + 542872556586838090927651041*x5^4 - 228505537090735857283193730*x5^3 + 32330997896447299987534725*x5^2 - 1943818647740666708782500*x5 + 42038339338757762562500);
//     TorsionSubgroup = Z/12 (NOT (3,12): M nonsquare).
//     NOTE: M = -3 * (129488/225)^2, so the extra 3-torsion class is rational
//     exactly over Q(sqrt(-3)) = Q(zeta_3): over the cyclotomic field the
//     (geometrically simple) Jacobian has torsion containing Z/3 x Z/12.
//     GEOMETRICALLY SIMPLE: L_p irreducible AND 12th-power transform
//     irreducible at p = 29 AND p = 53 (certificate passes twice).
//     disc primes: [2, 3, 5, 19, 241, 379, 739, 1129, 8093].
//     G2-invariants differ from both split-hit curves (genuinely new curve).
// [9] STRUCTURAL: rational points of S12 do NOT automatically give [3,12]:
//     the extra 3-class is rational iff additionally M is a square.  The true
//     carrier is the double cover m^2 = M of S12.  Conditional on
//     Hs(Q) = the 6 points in [4], the z=-5/3 fiber has EXACTLY 7 rational
//     points: the 6 split ones (M square) and the new one in [7] (M nonsquare)
//     — i.e. NO new [3,12] curve on this fiber.
