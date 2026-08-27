// batch5.m — JACKPOT protocol for the 5 new law-sweep hit members + known hit.
// For each: build curve, G2Invariants (dedupe), exact TorsionSubgroup,
// multi-prime simplicity certificates (chi irred deg 4 AND chi^12 irred deg 4).
P<x> := PolynomialRing(Rationals());
A := [1,1,1,2,2];
pts := [
  <-97, 48, -49, 240, "KNOWN">,
  <169, 48, -169, 34800, "A">,
  <-23, 75, -9025, 3519, "B">,
  <265, 169, -4225, 3519, "C">,
  <169, 290, -169, 5760, "D">,
  <-23, 338, -33124, 10557, "E">
];
G2s := [];
for pt in pts do
    u := pt[1]/pt[2]; r := pt[3]/pt[4];
    s := u*(r-1); m := 4*r*u*(u-1); n := r-1;
    den := LCM([Denominator(z) : z in [s,m,n]]);
    si := Integers()!(s*den); mi := Integers()!(m*den); ni := Integers()!(n*den);
    g := GCD([si,mi,ni]); si div:= g; mi div:= g; ni div:= g;
    B := [2*si^2-si*ni, 2*si^2+si*mi-2*si*ni-mi*ni, 2*si^2+si*mi-si*ni-mi*ni,
          -mi*ni, 4*si^2-4*si*ni-mi*ni];
    f := &*[A[i] + B[i]*x : i in [1..5]];
    printf "=== %o: u=%o/%o r=%o/%o ===\n", pt[5], pt[1], pt[2], pt[3], pt[4];
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        printf "DEGENERATE\n"; continue;
    end if;
    C := HyperellipticCurve(f);
    g2 := G2Invariants(C);
    Append(~G2s, <pt[5], g2>);
    printf "(s,m,n) = (%o,%o,%o)\n", si, mi, ni;
    J := Jacobian(C);
    T := TorsionSubgroup(J);
    printf "TORSION INVARIANTS = %o\n", Invariants(T);
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
                    printf "SIMPLICITY CERT p=%o: chi=%o\n", p, chi;
                    ncert +:= 1;
                    if ncert ge 3 then break; end if;
                end if;
            end if;
        end if;
    end for;
    printf "certificates: %o\n", ncert;
end for;
printf "=== G2 DEDUPE ===\n";
for i in [1..#G2s] do
    for j in [i+1..#G2s] do
        if G2s[i][2] eq G2s[j][2] then
            printf "SAME CURVE: %o == %o\n", G2s[i][1], G2s[j][1];
        end if;
    end for;
end for;
for t in G2s do printf "%o: %o\n", t[1], t[2]; end for;
quit;
