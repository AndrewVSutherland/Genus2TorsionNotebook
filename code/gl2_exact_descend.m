// 37-hunt endgame: exact descent test of the two reconstructed curves.
// Curve 1 / K3 = Q[t]/(t^3+3t-54)   (star-stable type-(1,2) lattice)
// Curve 2 / K4 = Q[t]/(t^4-36t^2-282t-108)  (star-stable type-(1,3) lattice)
// For each: exact G2 invariants -> minimal polynomials over Q.  If rational,
// Mestre (HyperellipticCurveFromG2Invariants) over Q, twist-match to the
// LMFDB traces of 2190.2.a.v, TorsionSubgroup, 37-check.
//
// Run: magma -b code/gl2_exact_descend.m > results/gl2_exact_descend.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(16*10^9);

Q := Rationals();
trtargets := [<7,2>, <11,-6>, <13,-2>, <17,6>, <19,10>];

procedure HandleCurve(C, tag)
    printf "== %o over %o\n", tag, DefiningPolynomial(BaseRing(C));
    g2 := G2Invariants(C);
    degs := [];
    for k in [1..3] do
        mp := MinimalPolynomial(g2[k], Q);
        printf "%o G2[%o] mindeg=%o\n", tag, k, Degree(mp);
        Append(~degs, Degree(mp));
    end for;
    if Max(degs) ne 1 then
        printf "%o NOT_Q_MODULI degs=%o\n", tag, degs;
        // print the minimal polynomials for the record
        for k in [1..3] do
            printf "%o G2[%o] minpoly %o\n", tag, k, MinimalPolynomial(g2[k], Q);
        end for;
        return;
    end if;
    g2Q := [ Q!g2[k] : k in [1..3] ];
    printf "%o G2_RATIONAL %o\n", tag, g2Q;
    CQ := HyperellipticCurveFromG2Invariants(g2Q);
    if Type(BaseRing(CQ)) ne FldRat then
        printf "%o MESTRE_FIELD %o (obstruction?)\n", tag, BaseRing(CQ);
        return;
    end if;
    CQ := ReducedMinimalWeierstrassModel(CQ);
    printf "%o QMODEL %o\n", tag, CQ;
    fC, hC := HyperellipticPolynomials(CQ);
    gC := 4*fC + hC^2;
    for dt in [1,-1,2,-2,3,-3,5,-5,6,-6,10,-10,15,-15,30,-30,73,-73,146,-146,219,-219,365,-365,438,-438,730,-730,1095,-1095,2190,-2190] do
        Cd := HyperellipticCurve(dt*gC);
        match := true;
        for tt in trtargets do
            p := tt[1];
            if Integers()!Discriminant(Cd) mod p eq 0 then continue; end if;
            Cp := ChangeRing(Cd, GF(p));
            chi := Reverse(Coefficients(LPolynomial(Cp)));
            if -Integers()!chi[2] ne tt[2] then match := false; break; end if;
        end for;
        if match then
            printf "%o TWIST_MATCH d=%o\n", tag, dt;
            Cmin := Cd;
            try Cmin := ReducedMinimalWeierstrassModel(Cd); catch e ; end try;
            Cs := SimplifiedModel(Cmin);
            printf "%o MATCHED_CURVE %o\n", tag, Cs;
            Tt := TorsionSubgroup(Jacobian(Cs));
            printf "%o TORSION %o (order %o)\n", tag, Invariants(Tt), #Tt;
            if #Tt mod 37 eq 0 then
                printf "*** THEOREM: genus-2 Jacobian over Q with a rational point of order 37 ***\n";
                printf "CURVE %o\n", Cmin;
            end if;
            break;
        end if;
    end for;
end procedure;

// Curve 1
K3<t> := NumberField(Polynomial(Q, [-54, 3, 0, 1]));
P3<x> := PolynomialRing(K3);
f1 := 1/70761*(304360*t^2 + 194256*t - 5052516)*x^6
    + 1/23587*(-413320*t^2 - 1162368*t + 4982718)*x^5
    + 1/94348*(315736*t^2 - 397440*t - 17223893)*x^4
    + 1/141522*(-1137241*t^2 + 801168*t - 5469876)*x^3
    + 1/94348*(786865*t^2 + 2043498*t + 17488258)*x^2
    + 1/23587*(-50136*t^2 - 165408*t - 1491258)*x
    + 1/94348*(-12331*t^2 + 184258*t + 411795);
C1 := HyperellipticCurve(f1);
HandleCurve(C1, "C1");

// Curve 2
K4<u> := NumberField(Polynomial(Q, [-108, -282, -36, 0, 1]));
P4<xx> := PolynomialRing(K4);
f2 := 1/1044523*(-38249867*u^3 + 16092346*u^2 + 1354458784*u + 10415745162)*xx^6
    + 1/1044523*(-46415533*u^3 + 14269734*u^2 + 1635962256*u + 13064166648)*xx^5
    + 1/12534276*(-690913172*u^3 + 220668552*u^2 + 24448195152*u + 193092189417)*xx^4
    + 1/10445230*(-354580313*u^3 + 165432492*u^2 + 12524675850*u + 96565503951)*xx^3
    + 1/20890460*(-408985824*u^3 + 287648556*u^2 + 14744286420*u + 109620654143)*xx^2
    + 1/15667845*(-84430879*u^3 + 55331541*u^2 + 3129252165*u + 22701390888)*xx
    + 1/47003535*(-69415217*u^3 + 32969403*u^2 + 2542072275*u + 18447136674);
C2 := HyperellipticCurve(f2);
HandleCurve(C2, "C2");

printf "GL2_DESCEND_DONE\n";
quit;
