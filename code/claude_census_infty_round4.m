// Round 4 (user's isogeny idea): Richelot partners of the certified
// [2,2,4,4]- and [2,2,2,8]-families; tabulate partner torsions per fiber.
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
            if Type(s) eq SetCart then printf "  partner %o: E x E product\n", i; continue; end if;
            Cs := 0;
            try Cs := Curve(s); catch e ; end try;
            if Type(Cs) eq RngIntElt then printf "  partner %o: (no curve)\n", i; continue; end if;
            try
                fs, hs := HyperellipticPolynomials(SimplifiedModel(Cs));
                dn := LCM([Denominator(co) : co in Coefficients(fs)]);
                Ci := HyperellipticCurve(dn^2*fs);
                Ti := Invariants(TorsionSubgroup(Jacobian(Ci)));
                printf "  partner %o: torsion %o\n", i, Ti;
            catch e2
                printf "  partner %o: torsfail %o\n", i, Substring(Sprint(e2`Object),1,60);
            end try;
        end for;
    catch e
        printf "  RICHELOT_ERR %o\n", Substring(Sprint(e`Object),1,80);
    end try;
end procedure;

// [2,2,4,4]-family fibers (paper's s-family)
for s in [Q!3, Q!5, Q!7] do
    a := (s-2)*(s^2-s+1/2)*(s^2-3/2)*(s^4+s^2+9/4);
    b := (s^2-2*s+1/2)*(s^2-2*s+9/2)*(s^2-3/2)*(s^2+1);
    c := -(s-2)*(2*s^2-2*s+1)*(s^4+s^2+9/4);
    d := (s-2)*(2*s^2-2*s+1)*(2*s^2-3)*(s^2+1);
    f := x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
    Partners(f, Sprintf("2244 s=%o", s));
end for;

// [2,2,2,8]-family fibers (K3 rational curve)
for t in [Q!2, Q!3, Q!5] do
    a := -4*t^2*(t+1)/(t^2+t+1)^2;
    b := -t/(t+1); c := Q!1; d := t;
    f := x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
    Partners(f, Sprintf("2228 t=%o", t));
end for;
printf "ROUND4_DONE\n";
quit;
