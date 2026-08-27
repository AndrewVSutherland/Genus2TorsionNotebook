// verify_split11.m — independent verification of the first split-Jacobian
// realization of torsion [11] (session 2026-08-26), from integer
// coefficients alone.  Asserts:
//   (1) exact torsion [11] on the minimal model;
//   (2) point-count identities with the Weil restriction of the X1(11)
//       quadratic-point curve E/Q(sqrt 11) (s0 = 4/5): at split p,
//       #J(F_p) = #E(k_p1) * #E(k_p2); at inert p, #J(F_p) = #E(F_p^2);
//   (3) J is Q-simple: E is not isogenous to E^sigma (trace mismatch), and
//       E is not a Q-curve (j not rational), so no Q-isogeny factor exists;
//   (4) E is non-CM, E(K)_tors = [11].
// Run from product/code/:  magma -b verify_split11.m
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();
Px<x> := PolynomialRing(Q);
t0 := Cputime();

// the curve (minimal model of the descended twist-2 glue)
fmin := 1204142*x^6 - 5109634*x^5 + 31412066*x^4 - 65405928*x^3 + 99564424*x^2 + 94188680*x + 9471400;
C := HyperellipticCurve(fmin);
J := Jacobian(C);
I := Invariants(TorsionSubgroup(J));
assert I eq [11];
printf "(1) exact torsion [11] on y^2 = %o\n", fmin;

// the elliptic factor over K = Q(sqrt 11): X1(11) raw-model point s0 = 4/5
K<w> := QuadraticField(11);
OK := Integers(K);
RK<xk> := PolynomialRing(K);
s0 := K!(4/5);
r0 := Roots(xk^2 - (s0^3-3*s0^2+4*s0)*xk + s0)[1][1];
c0 := s0*(r0-1); b0 := r0*c0;
E := EllipticCurve([1-c0, -b0, -b0, 0, 0]);
sig := hom< K -> K | -w >;
Es := EllipticCurve([sig(a) : a in aInvariants(E)]);
assert Order(E![0,0]) eq 11;
assert Invariants(TorsionSubgroup(E)) eq [11];
assert not jInvariant(E) in Q;
assert not HasComplexMultiplication(E);
printf "(4) E/K: X1(11) point, E(K)_tors = [11], j irrational, non-CM\n";

// (2) count identities at all good p < 200
dsc := Integers()!Discriminant(fmin);
nsplit := 0; ninert := 0;
for p in PrimesInInterval(13, 200) do
    if p eq 11 or dsc mod p eq 0 then continue; end if;
    fp := PolynomialRing(GF(p))!fmin;
    if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
    nj := #Jacobian(HyperellipticCurve(fp));
    dec := Decomposition(OK, p);
    if #dec eq 2 then
        pred := #Reduction(E, dec[1][1]) * #Reduction(E, dec[2][1]);
        assert nj eq pred;
        nsplit +:= 1;
    else
        pred := #Reduction(E, dec[1][1]);
        assert nj eq pred;
        ninert +:= 1;
    end if;
    assert nj mod 11 eq 0;
end for;
assert nsplit ge 10 and ninert ge 10;
printf "(2) #J(F_p) = Weil-restriction counts at %o split + %o inert primes; 11 | #J everywhere\n", nsplit, ninert;

// (3) Q-simplicity: E not isogenous to E^sigma
mm := 0; nt := 0;
for p in PrimesInInterval(5, 300) do
    dec := Decomposition(OK, p);
    if #dec ne 2 then continue; end if;
    okr := true; t1 := 0; t2 := 0;
    try
        t1 := TraceOfFrobenius(Reduction(E, dec[1][1]));
        t2 := TraceOfFrobenius(Reduction(Es, dec[1][1]));
    catch e; okr := false; end try;
    if not okr then continue; end if;
    nt +:= 1;
    if t1 ne t2 then mm +:= 1; end if;
end for;
assert mm gt 0;
printf "(3) E not isogenous to E^sigma (%o/%o trace mismatches) and j not in Q: J is Q-simple\n", mm, nt;

printf "SPLIT [11] VERIFIED (%.1o s)\n", Cputime() - t0;
quit;
