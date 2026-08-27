// claude_ov_lane7_conductor2.m -- Lane 7 (overnight 2026-07-25).
//
// CORRECTION RUN.  PARI's genus2red does NOT compute the conductor exponent at
// p = 2: it returns the sentinel exponent -1 in its factorisation matrix and
// then omits the 2-part from N entirely (checked directly --
// genus2red(x^5-x) returns N = 1, and for the RM curve it returns the matrix
// [2,-1; 3,2; 5,2; 13,2] with N = 38025).  Every conductor recorded in
// notes/claude_generic_2214_found_2026_07_23.md is therefore the ODD PART only,
// and the note's remark "all conductors are odd non-squares" is an artefact.
//
// This script recomputes the conductor of ONE contact-7 three-root curve with
// Magma (Ogg's formula; Magma itself warns when v_2(disc) >= 12, in which case
// the 2-part is not guaranteed -- we record the warning status via v_2(disc)).
// One curve per invocation, because Conductor() on the large-discriminant
// members exceeds the 3 GB lab cap when several are done in one session.
//
// Usage: magma -b IDX:=k code/claude_ov_lane7_conductor2.m
// Markers: COND2 / LANE7_COND2_DONE
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();
if not assigned IDX then IDX := 1; elif Type(IDX) eq MonStgElt then IDX := StringToInteger(IDX); end if;

Gfun := func< v | -(v^5 - v^3 - v^2/2)/(v+1)^2 >;
trips := [
  [-10, -10/7, -1/2], [-5, -15/8, -15/22], [-3, -3/4, -3/5],
  [-15/8, -15/19, -1/2], [-5/18, -10/49, 4/17], [-4/9, -4/25, 4/17],
  [-511/61, -511/625, -1/2], [-165/41, -33/16, -165/289],
  [-164/297, -1/2, 164/361], [-17/50, -34/189, 34/121], [-1/2, -13/49, 13/50]
];
T := trips[IDX]; s := T[1]; t := T[2]; u := T[3];
c4 := (Gfun(s)-Gfun(t))/(s^2-t^2); c0 := Gfun(s) - c4*s^2;
b := c4 - 2; a := 9/2 - c0 - c4;
h := 1 - (7/2)*x + a*x^2 + b*x^3;
f := (h^2 + (x-1)^7) div x^2;
for w in [s,t,u] do assert Evaluate(f, 1-w^2) eq 0; end for;
den := LCM([Denominator(co) : co in Coefficients(f)]);
F := P![ co*den^2 : co in Coefficients(f) ];
C := HyperellipticCurve(F);
Cm := ReducedMinimalWeierstrassModel(C);
Dm := Z!Discriminant(Cm);
v2 := Valuation(AbsoluteValue(Dm), 2);
N := Conductor(Cm);
printf "COND2 idx=%o (s,t,u)=(%o, %o, %o) v_2(disc_min)=%o oggreliable=%o\n", IDX, s, t, u, v2, v2 lt 12;
printf "COND2 idx=%o N = %o = %o\n", IDX, N, Factorization(N);
printf "COND2 idx=%o oddpart(N) = %o   2part(N) = 2^%o\n", IDX, N div 2^Valuation(N,2), Valuation(N,2);
print "LANE7_COND2_DONE";
quit;
