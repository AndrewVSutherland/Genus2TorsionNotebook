// claude_ov_erratum_cond5.m -- genus2red ERRATUM propagation (2026-07-25).
//
// The one gap in Lane 7's conductor table: contact-7 three-root curve #5,
// (s,t,u) = (-5/18, -10/49, 4/17).  Magma's Conductor OOM'd there at ~40 GB
// (results/claude_ov_lane7_condaudit_5.log) and again at the 3 GB lab cap
// (results/claude_ov_lane7_conductor2.log).  v_2(disc_min) = 30 for this
// curve, so even a successful Ogg run carries no correctness guarantee -- but
// the number is worth having, flagged.
//
// Usage:   magma -b MemGB:=48 code/claude_ov_erratum_cond5.m
// Markers: COND5 / ERRATUM_COND5_DONE
SetColumns(0);
if not assigned MemGB then MemGB := 48; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();

Gfun := func< v | -(v^5 - v^3 - v^2/2)/(v+1)^2 >;
s := -5/18; t := -10/49; u := 4/17;
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
printf "COND5 (s,t,u)=(%o, %o, %o)  v_2(disc_min)=%o  ogg_guaranteed=%o\n", s, t, u, v2, v2 lt 12;
printf "COND5 minimal model = %o\n", Cm;
N := Conductor(Cm);
printf "COND5 N = %o = %o\n", N, Factorization(N);
printf "COND5 oddpart(N) = %o   2part(N) = 2^%o\n", N div 2^Valuation(N,2), Valuation(N,2);
print "ERRATUM_COND5_DONE";
quit;
