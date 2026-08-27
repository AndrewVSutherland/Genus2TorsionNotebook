// hitverify.m — JACKPOT protocol for a T3/T5 pencil hit on M(2,2,2,6).
// Args (magma -b): cls:=3 or 5; pp, qq (u = pp/qq); rn, rd (rho = rn/rd).
// Builds the curve y^2 = prod(A_i + B_i x), computes TorsionSubgroup, then a
// multi-prime geometric-simplicity certificate (irreducible chi AND chi^12).
P<x> := PolynomialRing(Rationals());
cl := StringToInteger(cls);
u := StringToInteger(pp)/StringToInteger(qq);
r := StringToInteger(rn)/StringToInteger(rd);
if cl eq 3 then
    D := u - 1 + r; N := (2-4*r)*u - (1-4*r);
    s := u*D; m := -u*N; n := D;
else
    s := u*(r-1); m := 4*r*u*(u-1); n := r-1;
end if;
den := LCM([Denominator(z) : z in [s,m,n]]);
si := Integers()!(s*den); mi := Integers()!(m*den); ni := Integers()!(n*den);
g := GCD([si,mi,ni]); si div:= g; mi div:= g; ni div:= g;
printf "(s,m,n) = (%o,%o,%o)\n", si, mi, ni;
A := [1,1,1,2,2];
B := [2*si^2-si*ni, 2*si^2+si*mi-2*si*ni-mi*ni, 2*si^2+si*mi-si*ni-mi*ni,
      -mi*ni, 4*si^2-4*si*ni-mi*ni];
printf "B = %o\n", B;
f := &*[A[i] + B[i]*x : i in [1..5]];
if Degree(f) ne 5 or Discriminant(f) eq 0 then
    printf "DEGENERATE (deg %o) — not a valid curve\n", Degree(f); quit;
end if;
C := HyperellipticCurve(f);
J := Jacobian(C);
T := TorsionSubgroup(J);
inv := Invariants(T);
printf "TORSION INVARIANTS = %o\n", inv;
tgt := [2,2,2,12]; k := #inv;
cont := k ge 4 and &and[IsDivisibleBy(inv[k-i], tgt[4-i]) : i in [0..3]];
printf "contains (2,2,2,12): %o\n", cont;
// multi-prime simplicity scan
DD := Integers()!Discriminant(C);
ncert := 0;
for p in PrimesInInterval(17, 500) do
    if DD mod p ne 0 then
        Cp := HyperellipticCurve(PolynomialRing(GF(p))!f);
        Lp := LPolynomial(Cp);
        chi := Parent(Lp)!Reverse(Coefficients(Lp));
        if IsIrreducible(chi) and Degree(chi) eq 4 then
            M := CompanionMatrix(chi);
            chi12 := CharacteristicPolynomial(M^12);
            if IsIrreducible(chi12) and Degree(chi12) eq 4 then
                printf "SIMPLICITY CERT p=%o: chi=%o irred, chi12 irred deg 4\n", p, chi;
                ncert +:= 1;
                if ncert ge 4 then break; end if;
            end if;
        end if;
    end if;
end for;
printf "certificates found: %o (need >=1, prefer several)\n", ncert;
printf "G2INV = %o\n", G2Invariants(C);
quit;
