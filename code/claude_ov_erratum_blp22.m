// claude_ov_erratum_blp22.m -- genus2red ERRATUM propagation (2026-07-25).
//
// The second genus2red conductor recorded in this repository (the first being
// the six/eleven contact-7 [2,2,14] curves): the BLP C4-CORRECTED [2,22]
// curve of notes/claude_generic_222_2214_plan_2026_07_23.md, recorded there as
// "conductor odd part 645^2".  gp confirms genus2red returns 416025 = 645^2
// WITH the [2,-1] sentinel (results/claude_ov_erratum_sentinel.log), so that
// is indeed only the odd part.  Recompute the full conductor with Magma.
//
// The curve is RM(sqrt 5) by the note's real-subfield-disc census, hence (if
// the RM is rational) of GL2-type, hence of SQUARE conductor -- so v_2(N)
// should come out even.  That is a free consistency check on the 2-part.
//
// Usage:   magma -b MemGB:=16 code/claude_ov_erratum_blp22.m
// Markers: BLP22 / ERRATUM_BLP22_DONE
SetColumns(0);
if not assigned MemGB then MemGB := 16; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();

F := x^6-18*x^5-4001*x^4-22524*x^3+859039*x^2-1926258*x-9043839;
assert F eq (x-9)*(x+21)*(x^2-80*x+439)*(x^2+50*x+109);
C := HyperellipticCurve(F);
Cm := ReducedMinimalWeierstrassModel(C);
Dm := Z!Discriminant(Cm);
v2 := Valuation(AbsoluteValue(Dm), 2);
printf "BLP22 minimal model = %o\n", Cm;
printf "BLP22 v_2(disc_min) = %o  ogg_guaranteed = %o\n", v2, v2 lt 12;
N := Conductor(Cm);
printf "BLP22 N = %o = %o\n", N, Factorization(N);
printf "BLP22 oddpart(N) = %o (should be 416025 = 645^2: %o)   2part(N) = 2^%o\n",
    N div 2^Valuation(N,2), (N div 2^Valuation(N,2)) eq 416025, Valuation(N,2);
issq, rt := IsSquare(N);
printf "BLP22 N is a perfect square: %o", issq;
if issq then printf "  (N = %o^2)", rt; end if;
printf "\n";
// exact torsion, as an independent re-check.  TorsionSubgroup REQUIRES a model
// y^2 = f with f integral, so pass through SimplifiedModel (y^2 = 4f + h^2).
Cs := SimplifiedModel(Cm);
printf "BLP22 simplified model = %o\n", Cs;
T := TorsionSubgroup(Jacobian(Cs));
printf "BLP22 torsion = %o (order %o)\n", Invariants(T), #T;
print "ERRATUM_BLP22_DONE";
quit;
