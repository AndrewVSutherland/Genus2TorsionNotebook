// [2,22] hunt: fresh-session verifier for a confirmed sieve hit on
//   f = (x^2-1)(x^2+A x+B)(x^2+C x+D).
// Asserts: exact TorsionSubgroup contains (Z/2)^2 x Z/11; prints invariants,
// factor type, minimal model; strict Frobenius root-power simplicity
// certificates at >= 2 primes (NB an RM surface also passes these -- the
// non-RM claim is certified separately by code/verify_222_nonrm_hit.sage).
//
// Run: magma -b Av:=.. Bv:=.. Cv:=.. Dv:=.. code/verify_222_nonrm_hit.m

SetColumns(0);
SetSeed(1);
if not assigned Av then error "pass Av,Bv,Cv,Dv"; elif Type(Av) eq MonStgElt then Av := StringToInteger(Av); end if;
if assigned Bv and Type(Bv) eq MonStgElt then Bv := StringToInteger(Bv); end if;
if assigned Cv and Type(Cv) eq MonStgElt then Cv := StringToInteger(Cv); end if;
if assigned Dv and Type(Dv) eq MonStgElt then Dv := StringToInteger(Dv); end if;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals(); Z := Integers();
P<x> := PolynomialRing(Q);

f := (x^2-1)*(x^2+Av*x+Bv)*(x^2+Cv*x+Dv);
assert Discriminant(f) ne 0;
printf "TUPLE (A,B,C,D) = (%o,%o,%o,%o)\n", Av, Bv, Cv, Dv;
printf "F = %o\n", f;
printf "FACTOR_DEGREES = %o\n", {* Degree(fp[1])^^fp[2] : fp in Factorization(f) *};

C0 := HyperellipticCurve(f);
Cmin := ReducedMinimalWeierstrassModel(C0);
fmin, hmin := HyperellipticPolynomials(Cmin);
printf "REDUCED_MINIMAL = %o, %o\n", fmin, hmin;
disc := Z!Discriminant(Cmin);
printf "MINIMAL_DISC_FACTORIZATION = %o\n", Factorization(disc);

Cs := SimplifiedModel(Cmin);
J := Jacobian(Cs);
T := TorsionSubgroup(J);
invs := Invariants(T);
printf "TORSION_INVARIANTS = %o\n", invs;
printf "TORSION_ORDER = %o\n", #T;
assert #T mod 11 eq 0;
assert #invs ge 2 and invs[#invs-1] mod 2 eq 0;   // 2-rank >= 2
printf "CONTAINS_2_2_11 = true\n";

// strict geometric-simplicity witnesses
nc := 0;
for pp in PrimesInInterval(17, 500) do
    if disc mod pp eq 0 then continue; end if;
    Cp := ChangeRing(Cmin, GF(pp));
    chi := P!Reverse(Coefficients(LPolynomial(Cp)));
    if not IsIrreducible(chi) or Degree(chi) ne 4 then continue; end if;
    K<pi> := NumberField(chi);
    degs := [Degree(MinimalPolynomial(pi^e)) : e in [2..12]];
    if &and[dd eq 4 : dd in degs] then
        nc +:= 1;
        printf "SIMPLICITY_WITNESS p=%o chi=%o\n", pp, chi;
        if nc ge 3 then break; end if;
    end if;
end for;
printf "WITNESS_COUNT = %o\n", nc;
assert nc ge 2;
printf "VERIFIED_TORSION_AND_SIMPLICITY (End=Z certification: run the .sage companion)\n";
quit;
