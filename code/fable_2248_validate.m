//////////////////////////////////////////////////////////////////////
// fable_2248_validate.m  (2026-07-18, Fable session)
//
// Fresh-session validation of a candidate curve from the bank walk:
// exact TorsionSubgroup + two-prime Frobenius root-power geometric-
// simplicity certificate (chi_p irreducible over Q AND
// Degree(MinimalPolynomial(pi^n)) = 4 for n = 2..12).
//
// Run:  magma -b fstr:="<integer coefficient list c0,...,c6>" \
//              code/fable_2248_validate.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned fstr then
    error "pass fstr:=\"c0,c1,...\" (ascending integer coefficients)";
end if;
coeffs := [StringToInteger(s) : s in Split(fstr, ",")];
f := P!coeffs;
printf "CURVE f=%o\n", f;

C := HyperellipticCurve(f);
J := Jacobian(C);
time T, mp := TorsionSubgroup(J);
inv := Invariants(T);
printf "TORSION inv=%o order=%o\n", inv, #T;

function RootPowerCert(f, p)
    Zf := P![Z!c : c in Coefficients(f)];
    Fp := GF(p); Pp := PolynomialRing(Fp);
    fp := Pp![Fp!c : c in Coefficients(Zf)];
    if Degree(fp) ne Degree(Zf) or not IsSeparable(fp) then
        return false, "bad reduction";
    end if;
    Cp := HyperellipticCurve(fp);
    Lp := LPolynomial(Cp);
    // Frobenius characteristic polynomial: chi(T) = T^4 * Lp(1/T)
    cs := Coefficients(Lp);            // degree 4, cs[1] = 1
    chi := P!Reverse([Q!c : c in cs]);
    if not IsIrreducible(chi) then return false, Sprintf("chi_%o reducible", p); end if;
    K<pi> := NumberField(chi);
    for n in [2..12] do
        if Degree(MinimalPolynomial(pi^n)) ne 4 then
            return false, Sprintf("power %o drops degree at p=%o", n, p);
        end if;
    end for;
    return true, Sprintf("chi_%o = %o", p, chi);
end function;

ncert := 0;
for p in [13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73] do
    if ncert ge 2 then break; end if;
    ok, msg := RootPowerCert(f, p);
    if ok then
        ncert +:= 1;
        printf "CERT p=%o OK: %o\n", p, msg;
    end if;
end for;
if ncert ge 2 then
    printf "SIMPLICITY_CERTIFIED (two independent good primes)\n";
else
    printf "SIMPLICITY_NOT_CERTIFIED — run the Sage/Lombardo endomorphism test\n";
end if;
quit;
