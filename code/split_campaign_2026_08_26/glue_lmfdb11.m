// test the five LMFDB quadratic 11-torsion pairs for sigma-2-congruence and
// glue any match (same validated pipeline as lane_qglue.m).
SetColumns(0);
SetMemoryLimit(12*10^9);
Q := Rationals();
Px<x> := PolynomialRing(Q);
QT<T> := PolynomialRing(Q);

// label, field polynomial, ainvs (pairs c0,c1 w.r.t. [1, K.1])
DATA := [*
  <"2.2.8.1-46.2-b2",   T^2-2,   [[1,1],[-1,-1],[1,0],[-3,1],[3,-2]]>,
  <"2.2.17.1-172.2-f1", T^2-T-4, [[1,0],[-1,-1],[1,1],[-27,11],[455,-179]]>,
  <"2.0.7.1-268.3-b1",  T^2-T+2, [[1,0],[-1,1],[1,1],[-2,-2],[7,-5]]>,
  <"2.2.13.1-828.2-d1", T^2-T-3, [[1,0],[0,0],[0,0],[-104,-297],[-29541,16497]]>,
  <"2.0.8.1-4338.4-b2", T^2+2,   [[1,1],[-1,1],[1,0],[122,62],[-483,286]]>
*];

for rec in DATA do
    lbl := rec[1];
    K<wg> := NumberField(rec[2]);
    // nontrivial automorphism
    orts := [ rr[1] : rr in Roots(PolynomialRing(K)!rec[2]) ];
    conj := [ rr : rr in orts | rr ne wg ][1];
    sig := hom< K -> K | conj >;
    ai := [ K!(c[1] + c[2]*wg) : c in rec[3] ];
    E := EllipticCurve(ai);
    Es := EllipticCurve([ sig(a) : a in ai ]);
    tors := Invariants(TorsionSubgroup(E));
    printf "== %o: tors %o\n", lbl, tors;
    if tors ne [11] then printf "   (skip: torsion)\n"; continue; end if;
    // 2-division cubic match over K
    f1 := HyperellipticPolynomials(WeierstrassModel(E));
    g1 := HyperellipticPolynomials(WeierstrassModel(Es));
    RK := PolynomialRing(K);
    if not IsIrreducible(RK!f1) then printf "   reducible 2-cubic?!\n"; continue; end if;
    LC := ext< K | RK!f1 >;
    match := #Roots(PolynomialRing(LC)!(RK!g1)) gt 0;
    printf "   sigma-2-congruent: %o\n", match;
    if not match then continue; end if;
    // ---- glue + descent (validated pipeline) ----
    L := SplittingField(f1);
    RL := PolynomialRing(L);
    al := [ rr[1] : rr in Roots(RL!f1) ];
    gm := [ rr[1] : rr in Roots(RL!g1) ];
    if #al ne 3 or #gm ne 3 then printf "   ROOTFAIL\n"; continue; end if;
    DfK := Discriminant(f1); DgK := Discriminant(g1);
    GA := Automorphisms(L);
    GK := [ gA : gA in GA | gA(L!wg) eq L!wg ];
    Gaction := func< rs0 | [[ Index(rs0, gA(rr)) : rr in rs0 ] : gA in GK ] >;
    aK := [[ al[i]-al[j] : j in [1..3]] : i in [1..3]];
    RtL := PolynomialRing(L);
    for sgm in SymmetricGroup(3) do
        beta := [ gm[i^sgm] : i in [1..3] ];
        bb := [[ beta[i]-beta[j] : j in [1..3]] : i in [1..3]];
        if al[1]*bb[3][2] + al[2]*bb[1][3] + al[3]*bb[2][1] eq 0 then continue; end if;
        if Gaction(al) ne Gaction(beta) then continue; end if;
        a1g := aK[3][2]^2/bb[3][2]+aK[2][1]^2/bb[2][1]+aK[1][3]^2/bb[1][3];
        a2g := al[1]*bb[3][2]+al[2]*bb[1][3]+al[3]*bb[2][1];
        b1g := bb[3][2]^2/aK[3][2]+bb[2][1]^2/aK[2][1]+bb[1][3]^2/aK[1][3];
        b2g := beta[1]*aK[3][2]+beta[2]*aK[1][3]+beta[3]*aK[2][1];
        A := DgK*a1g/a2g;  B := DfK*b1g/b2g;
        tL := RtL.1;
        sext := -(A*aK[2][1]*aK[1][3]*tL^2+B*bb[2][1]*bb[1][3])
                *(A*aK[3][2]*aK[2][1]*tL^2+B*bb[3][2]*bb[2][1])
                *(A*aK[1][3]*aK[3][2]*tL^2+B*bb[1][3]*bb[3][2]);
        okS := true; SX := RK!0;
        try SX := RK![ K!c : c in Coefficients(sext) ]; catch e; okS := false; end try;
        if not okS then continue; end if;
        printf "   glue sextic over K (matching %o)\n", sgm;
        // descent
        a6 := Coefficient(SX,6); a4 := Coefficient(SX,4);
        a2e := Coefficient(SX,2); a0 := Coefficient(SX,0);
        if a6 eq 0 or a0 eq 0 or a2e eq 0 or a4 eq 0 then printf "   DEGEN\n"; continue; end if;
        MKK := 0; haveM := false;
        SXs := RK![ sig(cc) : cc in Coefficients(SX) ];
        lamP := sig(a6)/a6;
        if SXs eq lamP*SX then
            MKK := Matrix(K,2,2,[1,0,0,1]); haveM := true;
        else
            mu1 := sig(a6)/a0;
            k2 := sig(a4)/(mu1*a2e);
            if sig(a2e) eq mu1*a4*k2^2 and sig(a0) eq mu1*a6*k2^3 then
                okk, kap := IsSquare(k2);
                MKK := Matrix(K,2,2,[0, okk select kap else k2, 1,0]); haveM := true;
            else
                mu2 := sig(a6)/a6;
                l2 := sig(a4)/(mu2*a4);
                if sig(a2e) eq mu2*a2e*l2^2 and sig(a0) eq mu2*a0*l2^3 then
                    okk, lm := IsSquare(l2);
                    if okk then MKK := Matrix(K,2,2,[lm,0,0,1]); haveM := true; end if;
                end if;
            end if;
        end if;
        if not haveM then printf "   NODESCENTSHAPE\n"; continue; end if;
        Mks := Matrix(K,2,2,[ sig(e) : e in Eltseq(MKK) ]);
        CCm := Mks*MKK;
        if not (CCm[1,2] eq 0 and CCm[2,1] eq 0 and CCm[1,1] eq CCm[2,2]) then
            printf "   COCYCLEFAIL\n"; continue;
        end if;
        FQ := Px!0;
        for Aseq in [[K|1,0,0,1],[K|1,1,0,1],[K|1,wg,0,1],[K|1,0,wg,1],[K|2,wg,1,1],[K|1,2,wg,1]] do
            AM := Matrix(K,2,2,Aseq);
            AMs := Matrix(K,2,2,[sig(e) : e in Eltseq(AM)]);
            NM := AM + MKK*AMs;
            if Determinant(NM) eq 0 then continue; end if;
            nu := NM[1,1]*RK.1 + NM[1,2]; de := NM[2,1]*RK.1 + NM[2,2];
            T6 := &+[ Coefficient(SX,i)*nu^i*de^(6-i) : i in [0..6] ];
            degs := [ i : i in [0..6] | Coefficient(T6,i) ne 0 ];
            if #degs eq 0 then continue; end if;
            dg := Max(degs);
            T6s := RK![ sig(co) : co in Coefficients(T6) ];
            lam2 := Coefficient(T6s,dg)/Coefficient(T6,dg);
            if T6s ne lam2*T6 then continue; end if;
            for ccx in [K| 1, wg, 1+wg, 2+wg, 1-2*wg, 3+wg, 1+2*wg ] do
                mu := ccx + lam2*sig(ccx);
                if mu eq 0 then continue; end if;
                TT := mu*T6;
                if RK![ sig(co) : co in Coefficients(TT) ] eq TT then
                    FQ := Px![ Q!co : co in Coefficients(TT) ];
                    break;
                end if;
            end for;
            if FQ ne 0 then break; end if;
        end for;
        if FQ eq 0 then printf "   DESCENTFAIL\n"; continue; end if;
        // integralize + twist sieve + exact torsion
        denl := LCM([ Denominator(co) : co in Coefficients(FQ) ]);
        fZ := Px!(FQ*denl^2);
        cint := GCD([ Integers()!co : co in Coefficients(fZ) ]);
        sq := 1;
        for pr in Factorization(cint) do sq *:= pr[1]^(2*(pr[2] div 2)); end for;
        fZ := Px![ co div sq : co in [Integers()!c : c in Coefficients(fZ)] ];
        if Discriminant(fZ) eq 0 then printf "   SINGULAR\n"; continue; end if;
        dsc := Integers()!Discriminant(fZ);
        goodp := [ p : p in PrimesInInterval(13, 300) | p ne 11 and dsc mod p ne 0 ];
        dK := Integers()!Discriminant(MaximalOrder(K));
        tsupport := { pf[1] : pf in TrialDivision(AbsoluteValue(dK), 100) } join {2, 3, 5, 7, 11, 23};
        for dsupp in Subsets(tsupport) do
            d0 := &*[ Integers() | p : p in dsupp ];
            for d in [d0, -d0] do
                ok := true; ntest := 0;
                for p in goodp do
                    if ntest ge 12 then break; end if;
                    if d mod p eq 0 then continue; end if;
                    fp := PolynomialRing(GF(p))!(d*fZ);
                    if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
                    ntest +:= 1;
                    if #Jacobian(HyperellipticCurve(fp)) mod 11 ne 0 then ok := false; break; end if;
                end for;
                if ok and ntest ge 9 then
                    fT := Px!(d*fZ);
                    JT := Jacobian(HyperellipticCurve(fT));
                    I := Invariants(TorsionSubgroup(JT));
                    printf "   TORSION d=%o invs=%o\n", d, I;
                    if 11 in I or (#I gt 0 and I[#I] mod 11 eq 0) then
                        Cmin := ReducedMinimalWeierstrassModel(HyperellipticCurve(fT));
                        printf "HIT %o: invs=%o MINMODEL %o\n", lbl, I, HyperellipticPolynomials(Cmin);
                    end if;
                end if;
            end for;
        end for;
    end for;
end for;
printf "GLUE_LMFDB11_DONE\n";
quit;
