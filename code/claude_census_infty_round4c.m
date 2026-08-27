// Certify the Richelot-partner families: partner #5 of [2,2,4,4] (-> [4,4])
// and partner #4 of [2,2,2,8] (-> [2,8]): strict primes + invariants.
SetColumns(0);
SetSeed(1);
SetMemoryLimit(8*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);
function StrictPrime(C)
    fC, hC := HyperellipticPolynomials(C);
    g := 4*fC + hC^2;
    dsc := Integers()!Numerator(Discriminant(g)/16);
    for p in PrimesInInterval(11, 600) do
        if dsc mod p eq 0 then continue; end if;
        Pp := P!Reverse(Coefficients(LPolynomial(ChangeRing(SimplifiedModel(C), GF(p)))));
        if not IsIrreducible(Pp) then continue; end if;
        K<pi> := NumberField(Pp);
        ok := true;
        for n in [1..12] do
            if Degree(MinimalPolynomial(pi^n)) ne 4 then ok := false; break; end if;
        end for;
        if ok then return p, true; end if;
    end for;
    return 0, false;
end function;
procedure PartnerCert(f, idx, tag, target, ~invs)
    den := LCM([Denominator(c) : c in Coefficients(f)]);
    J := Jacobian(SimplifiedModel(HyperellipticCurve(den^2*f)));
    Rs := RichelotIsogenousSurfaces(J);
    s := Rs[idx];
    Cs := Curve(s);
    fs, hs := HyperellipticPolynomials(SimplifiedModel(Cs));
    dn := LCM([Denominator(co) : co in Coefficients(fs)]);
    Ci := HyperellipticCurve(dn^2*fs);
    Include(~invs, G2Invariants(Ci));
    T := Invariants(TorsionSubgroup(Jacobian(Ci)));
    printf "  %o partner%o torsion %o", tag, idx, T;
    if Sprint(T) eq target then
        sp, oks := StrictPrime(Ci);
        printf "  STRICT p=%o(%o)", sp, oks;
    end if;
    printf "\n";
end procedure;
printf "A) [4,4] via partner 5 of the [2,2,4,4]-family\n";
invA := {};
for s in [Q!3, Q!5, Q!7, Q!9, Q!11] do
    a := (s-2)*(s^2-s+1/2)*(s^2-3/2)*(s^4+s^2+9/4);
    b := (s^2-2*s+1/2)*(s^2-2*s+9/2)*(s^2-3/2)*(s^2+1);
    c := -(s-2)*(2*s^2-2*s+1)*(s^4+s^2+9/4);
    d := (s-2)*(2*s^2-2*s+1)*(2*s^2-3)*(s^2+1);
    f := x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
    PartnerCert(f, 5, Sprintf("s=%o",s), "[ 4, 4 ]", ~invA);
end for;
printf "A) distinct invariants: %o\n", #invA;
printf "B) [2,8] via partner 4 of the [2,2,2,8]-family\n";
invB := {};
for t in [Q!2, Q!3, Q!5, Q!7, Q!-2] do
    a := -4*t^2*(t+1)/(t^2+t+1)^2;
    b := -t/(t+1); c := Q!1; d := t;
    f := x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
    PartnerCert(f, 4, Sprintf("t=%o",t), "[ 2, 8 ]", ~invB);
end for;
printf "B) distinct invariants: %o\n", #invB;
printf "ROUND4C_DONE\n";
quit;
