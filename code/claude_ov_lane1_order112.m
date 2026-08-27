// claude_ov_lane1_order112.m -- Lane 1, overnight 2026-07-25.
//
// THE ORDER-112 ([2,2,2,14]) LOCUS ON THE u=-1/2 [2,2,14] FAMILY, DECIDED.
//
// On the family (code/claude_ov_lane1_family.m) parametrised by
//   E1 : w^2 = g(m) := -m^4 + 6m^2 + 4m = -m(m+2)(m^2-2m-2)   (rank 1)
// the residual quadratic of the contact-7 quintic has
//   Delta = B^2-4C = 4 * g(-m) / ( (m-1)^4 (m+2)^2 ),
// so the quintic SPLITS COMPLETELY (2-rank 4, torsion [2,2,2,14], ORDER 112,
// a new record) exactly when g(-m) is also a rational square.
//
// So the order-112 locus on this family is the genus-5 curve
//   C : w^2 = g(m),  z^2 = g(-m)
// i.e. "m and -m are BOTH x-coordinates of rational points of E1".
// C carries a (Z/2)^2 action over P^1_m with three intermediate quotients:
//   E1 : w^2 = g(m)             (rank 1)
//   E1': z^2 = g(-m) ~ E1       (rank 1)
//   H  : (wz/m)^2 = g(m)g(-m)/m^2 = (m^2-4)(m^4-8m^2+4)   -- genus 2
// and H is even in m, so Jac(H) ~ Ea x Eb splits.  This script computes
// H, its Jacobian's rank via the two elliptic quotients, does the point
// search / Chabauty, and pulls back to C.
//
// Usage: magma -b code/claude_ov_lane1_order112.m

SetColumns(0);
if not assigned MemGB then MemGB := 16; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned MBOUND then MBOUND := 2000; elif Type(MBOUND) eq MonStgElt then MBOUND := StringToInteger(MBOUND); end if;
Q := Rationals(); Z := Integers();
Pm<m> := PolynomialRing(Q);

g := -m^4 + 6*m^2 + 4*m;
gm := Evaluate(g, -m);
printf "g(m)  = %o = %o\n", g, Factorization(g);
printf "g(-m) = %o = %o\n", gm, Factorization(gm);
prod := g*gm;
printf "g(m)g(-m) = %o\n", Factorization(prod);
sext := prod div m^2;
printf "sextic (g*g(-)/m^2) = %o = %o\n", sext, Factorization(sext);
assert sext eq (m^2-4)*(m^4-8*m^2+4);

/* ---- 1. the genus-2 quotient H and its split Jacobian ---- */
H := HyperellipticCurve(sext);
printf "\nH : y^2 = %o  genus %o\n", sext, Genus(H);
JH := Jacobian(H);
printf "H is even in m -> Jac(H) splits.\n";
Pn<n> := PolynomialRing(Q);
fn := (n-4)*(n^2-8*n+4);
Ea := EllipticCurve(fn);                      // y^2 = f(n),   n = m^2
Cb := HyperellipticCurve(n*fn);               // v^2 = n f(n), quartic
Eb := EllipticCurve(Cb, Cb![0,0]);
printf "Ea : y^2 = %o    cond %o\n", fn, Conductor(Ea);
printf "Eb : v^2 = %o  ->  %o   cond %o\n", n*fn, Eb, Conductor(Eb);
la,ha := RankBounds(Ea); printf "rank(Ea) in [%o,%o]  torsion %o\n", la, ha, Invariants(TorsionSubgroup(Ea));
lb,hb := RankBounds(Eb); printf "rank(Eb) in [%o,%o]  torsion %o\n", lb, hb, Invariants(TorsionSubgroup(Eb));
printf "=> rank Jac(H) in [%o,%o]\n", la+lb, ha+hb;
try
  lj, hj := RankBounds(JH);
  printf "direct RankBounds(Jac H) = [%o,%o]\n", lj, hj;
catch e
  printf "direct RankBounds(Jac H) failed: %o\n", e`Object;
end try;

/* ---- 2. rational points on H ---- */
ptsH := Points(H : Bound := 10000);
printf "\nPoints(H : Bound:=10000) = %o\n", ptsH;

/* ---- 3. Chabauty on H if the rank permits ---- */
if ha + hb le 1 then
  printf "rank <= 1: attempting Chabauty on H\n";
  try
    JHm, mp := MordellWeilGroup(JH);
    printf "MW(Jac H) = %o\n", JHm;
    ptsC := Chabauty(mp(JHm.Ngens(JHm)));
    printf "Chabauty gives H(Q) = %o\n", ptsC;
  catch e
    printf "Chabauty attempt failed: %o\n", e`Object;
  end try;
end if;

/* ---- 4. direct search on C : both g(m) and g(-m) squares ---- */
printf "\n=== direct search: m = a/b, |a|,|b| <= %o, g(m) and g(-m) both squares ===\n", MBOUND;
nfound := 0; nhalf := 0;
for b in [1..MBOUND] do
  for a in [-MBOUND..MBOUND] do
    if GCD(Abs(a),b) ne 1 then continue; end if;
    if a eq 0 then continue; end if;
    mm := a/b;
    // g(m) = -m^4+6m^2+4m ; clear denominators: b^4 g = -a^4 + 6a^2b^2 + 4ab^3
    N1 := -a^4 + 6*a^2*b^2 + 4*a*b^3;
    if N1 lt 0 then continue; end if;
    if not IsSquare(N1) then continue; end if;
    nhalf +:= 1;
    N2 := -a^4 + 6*a^2*b^2 - 4*a*b^3;
    if N2 ge 0 and IsSquare(N2) then
      printf "ORDER112 CANDIDATE m = %o  (g=%o^2, g(-)=%o^2)\n", mm, Isqrt(N1), Isqrt(N2);
      nfound +:= 1;
    end if;
  end for;
  if b mod 200 eq 0 then printf "PROGRESS b=%o halfhits=%o fullhits=%o\n", b, nhalf, nfound; end if;
end for;
printf "SEARCH m-height <= %o : %o m with g(m) square, %o with BOTH square\n", MBOUND, nhalf, nfound;
printf "LANE1_ORDER112_DONE\n";
quit;
