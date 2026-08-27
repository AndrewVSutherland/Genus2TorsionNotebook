// Round 5: partner tabulation for [2,10]/[2,2,8] hunts.
// A) partners of the [2,20]-family (z-parametrized contact-5 locus)
// B) second-generation: partners of the [2,8]-partner family
SetColumns(0);
SetSeed(1);
SetMemoryLimit(8*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);
procedure Partners(f, tag)
    den := LCM([Denominator(c) : c in Coefficients(f)]);
    C := HyperellipticCurve(den^2*f);
    J := Jacobian(SimplifiedModel(C));
    T0 := Invariants(TorsionSubgroup(J));
    printf "FIBER %o source torsion %o\n", tag, T0;
    try
        Rs := RichelotIsogenousSurfaces(J);
        for i in [1..#Rs] do
            s := Rs[i];
            if Type(s) eq SetCart then continue; end if;
            Cs := 0;
            try Cs := Curve(s); catch e ; end try;
            if Type(Cs) eq RngIntElt then continue; end if;
            try
                fs, hs := HyperellipticPolynomials(SimplifiedModel(Cs));
                dn := LCM([Denominator(co) : co in Coefficients(fs)]);
                Ti := Invariants(TorsionSubgroup(Jacobian(HyperellipticCurve(dn^2*fs))));
                printf "  partner %o: %o\n", i, Ti;
            catch e2 ; end try;
        end for;
    catch e printf "  RERR\n"; end try;
end procedure;
// A) [2,20]: contact-5 with t(z) from the paper
printf "A) [2,20]-family partners\n";
for z in [Q!2, Q!3, Q!4] do
    t := -(z^4+4*z+4)/(z^4+4*z^3+8*z^2+8*z+4);
    b := (t^2-1)/2;
    f := (1+t*x+b*x^2)^2 - ((t+1)^4/4)*x^5;
    Partners(f, Sprintf("z=%o", z));
end for;
// B) second-generation from [2,8]-partners
printf "B) partners of the [2,8]-partner family\n";
for t in [Q!2, Q!3] do
    a := -4*t^2*(t+1)/(t^2+t+1)^2;
    b := -t/(t+1); c := Q!1; d := t;
    f := x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
    den := LCM([Denominator(co) : co in Coefficients(f)]);
    J := Jacobian(SimplifiedModel(HyperellipticCurve(den^2*f)));
    Rs := RichelotIsogenousSurfaces(J);
    Cs := Curve(Rs[4]);
    fs, hs := HyperellipticPolynomials(SimplifiedModel(Cs));
    dn := LCM([Denominator(co) : co in Coefficients(fs)]);
    Partners(dn^2*fs, Sprintf("gen2 t=%o", t));
end for;
printf "ROUND5_DONE\n";
quit;
