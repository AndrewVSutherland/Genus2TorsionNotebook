// lane_lmfdbglue.m — systematic sigma-glue sweep over ALL LMFDB elliptic
// curves E/K, K quadratic, with torsion order >= 5 (AVS directive
// 2026-08-26): for each curve, test whether a Jacobian can be cooked up in
// the isogeny class of Res_{K/Q}(E):
//   quadratic-j rows: sigma-2-congruence (generalized etale-algebra match of
//     the 2-division cubics) -> on match, the full validated glue + descent
//     + twist + exact-torsion pipeline (HIT when invariants not in KNOWN);
//     sigma-3-congruence (trace mod 3 then quartic module match) -> SURV3
//     (constructor to be built on demand);
//   rational-j rows (twist-type pairs; 2-glue is Kani-degenerate): the
//     chi_d-twisted 3-congruence necessary condition -> SURVT3.
// Data: ../data/lmfdb_quad_tors5.txt (label|D|a-coeffs|torsion|jrat).
// Usage: magma -b Part:=1 NParts:=1 lane_lmfdbglue.m > ../logs/lmfdbglue_p1.log
SetColumns(0);
if not assigned MemGB then MemGB := 10; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned Part then Part := 1; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;

load "split_lab.m";  // KNOWN (incl. [11] 2026-08-26) + Drew's spec + genus2.m

Q := Rationals();
Px<x> := PolynomialRing(Q);
QT<T> := PolynomialRing(Q);

lines := Split(Read("../data/lmfdb_quad_tors5.txt"), "\n");
printf "rows: %o\n", #lines;

// field cache by signed discriminant
FC := AssociativeArray();
function GetField(D)
    // returns K, sigma for the quadratic field of discriminant D (LMFDB
    // canonical: D = 1 mod 4 -> x^2-x-(D-1)/4; D = 0 mod 4 -> x^2 - D/4)
    if D mod 4 eq 1 then
        pol := T^2 - T - (D-1) div 4;
    else
        pol := T^2 - D div 4;
    end if;
    K := NumberField(pol);
    conj := [ rr[1] : rr in Roots(PolynomialRing(K)!pol) | rr[1] ne K.1 ][1];
    sig := hom< K -> K | conj >;
    return K, sig;
end function;

nrow := 0; nsurv2 := 0; nsurv3 := 0; nsurvt3 := 0; nhit := 0; nwarn := 0;
t0c := Cputime();
for li in [1..#lines] do
    if li mod NParts ne (Part - 1) then continue; end if;
    ln := lines[li];
    if #ln eq 0 then continue; end if;
    fld := Split(ln, "|");
    lbl := fld[1];
    Dfld := StringToInteger(fld[2]);
    cs := [ StringToInteger(c) : c in Split(fld[3], ",") ];
    torsdb := [ StringToInteger(c) : c in Split(fld[4], ",") ];
    jrat := StringToInteger(fld[5]) eq 1;
    if not IsDefined(FC, Dfld) then
        K0, s0m := GetField(Dfld);
        FC[Dfld] := <K0, s0m>;
    end if;
    K := FC[Dfld][1]; sig := FC[Dfld][2];
    wg := K.1;
    RK := PolynomialRing(K);
    ai := [ K | cs[2*i-1] + cs[2*i]*wg : i in [1..5] ];
    E := 0; okE := true;
    try E := EllipticCurve(ai); catch e; okE := false; end try;
    if not okE then printf "BADROW %o\n", lbl; continue; end if;
    nrow +:= 1;
    if nrow le 20 or nrow mod 500 eq 0 then
        // periodic data-integrity check
        IT := Invariants(TorsionSubgroup(E));
        if IT ne torsdb then printf "DATAWARN %o db=%o computed=%o\n", lbl, torsdb, IT; nwarn +:= 1; end if;
    end if;
    Es := EllipticCurve([ sig(a) : a in ai ]);
    OK := Integers(K);

    if jrat then
        // twist-type: chi_d-twisted 3-congruence necessary condition
        istw, dtw := IsQuadraticTwist(E, Es);
        if not istw then continue; end if;
        oksq, _ := IsSquare(dtw);
        if oksq then continue; end if;   // base change: E^sigma = E
        alive := true; ninert := 0;
        for p in PrimesInInterval(5, 500) do
            if ninert ge 10 then break; end if;
            for dd in Decomposition(OK, p) do
                pr := dd[1];
                okr := true; ap := 0; chid := 0;
                try
                    kp, red := ResidueClassField(pr);
                    dr := red(dtw);
                    if dr eq 0 then okr := false; end if;
                    if okr then
                        chid := IsSquare(dr) select 1 else -1;
                        ap := TraceOfFrobenius(Reduction(E, pr));
                    end if;
                catch e; okr := false; end try;
                if not okr then continue; end if;
                if chid eq -1 then
                    ninert +:= 1;
                    if ap mod 3 ne 0 then alive := false; break; end if;
                end if;
            end for;
            if not alive then break; end if;
        end for;
        if alive and ninert ge 8 then
            nsurvt3 +:= 1;
            printf "SURVT3 %o tors=%o D=%o\n", lbl, torsdb, Dfld;
        end if;
        continue;
    end if;

    // quadratic-j: sigma-2-congruence (generalized etale match)
    f1 := HyperellipticPolynomials(WeierstrassModel(E));
    g1 := HyperellipticPolynomials(WeierstrassModel(Es));
    fac1 := Factorization(RK!f1);
    fac2 := Factorization(RK!g1);
    sh1 := Sort([ Degree(ff[1]) : ff in fac1 ]);
    sh2 := Sort([ Degree(ff[1]) : ff in fac2 ]);
    m2 := false;
    if sh1 eq sh2 then
        if sh1 eq [1,1,1] then
            m2 := true;
        elif sh1 eq [1,2] then
            q1 := [ ff[1] : ff in fac1 | Degree(ff[1]) eq 2 ][1];
            q2 := [ ff[1] : ff in fac2 | Degree(ff[1]) eq 2 ][1];
            m2 := IsSquare(Discriminant(q1)*Discriminant(q2));
        else
            LC := ext< K | fac1[1][1] >;
            m2 := #Roots(PolynomialRing(LC)!(RK!g1)) gt 0;
        end if;
    end if;
    // sigma-3-congruence: quick trace test, then quartic module match
    m3 := false;
    okt3 := true; ntr := 0;
    for p in PrimesInInterval(5, 200) do
        if ntr ge 20 or not okt3 then break; end if;
        for dd in Decomposition(OK, p) do
            pr := dd[1];
            okr := true; t1 := 0; t2 := 0;
            try
                t1 := TraceOfFrobenius(Reduction(E, pr));
                t2 := TraceOfFrobenius(Reduction(Es, pr));
            catch e; okr := false; end try;
            if not okr then continue; end if;
            ntr +:= 1;
            if (t1 - t2) mod 3 ne 0 then okt3 := false; break; end if;
        end for;
    end for;
    if okt3 and ntr ge 15 then
        p3 := RK!DivisionPolynomial(E, 3);
        p3s := RK!DivisionPolynomial(Es, 3);
        ff3 := Factorization(p3); ff3s := Factorization(p3s);
        d3 := Sort([Degree(ff[1]) : ff in ff3]); d3s := Sort([Degree(ff[1]) : ff in ff3s]);
        if d3 eq d3s then
            if Degree(ff3[#ff3][1]) eq 1 then
                m3 := true;
            else
                L3 := ext< K | ff3[#ff3][1] >;
                m3 := exists{ ffx : ffx in ff3s | Degree(ffx[1]) eq Degree(ff3[#ff3][1])
                              and #Roots(PolynomialRing(L3)!ffx[1]) gt 0 };
            end if;
        end if;
    end if;
    if m3 then
        nsurv3 +:= 1;
        printf "SURV3 %o tors=%o D=%o\n", lbl, torsdb, Dfld;
    end if;
    if not m2 then continue; end if;
    nsurv2 +:= 1;
    printf "SURV2 %o tors=%o D=%o\n", lbl, torsdb, Dfld;
    // glue only rows whose odd part could yield a NEW group (>= 11): smaller
    // odd parts give [3]/[5]/[7]/[9]-type J-torsion, all realized.
    TORD := &*torsdb div 2^Valuation(&*torsdb, 2);   // odd part injects
    if TORD lt 11 then continue; end if;
    // ---- glue over K ----
    Es := EllipticCurve([sig(a) : a in aInvariants(E)]);
    f1 := HyperellipticPolynomials(WeierstrassModel(E));
    g1 := HyperellipticPolynomials(WeierstrassModel(Es));
    okL := true; L := 0;
    try L := SplittingField(f1); catch e; okL := false; end try;
    if not okL then printf "SPLITFAIL %o\n", lbl; continue; end if;
    RL := PolynomialRing(L);
    al := [ rr[1] : rr in Roots(RL!f1) ];
    gm := [ rr[1] : rr in Roots(RL!g1) ];
    if #al ne 3 or #gm ne 3 then printf "ROOTFAIL %o\n", lbl; continue; end if;
    DfK := Discriminant(f1); DgK := Discriminant(g1);
    // twist-pair shortcut (X1(13)-type: hyperelliptic involution = diamond, so
    // E^sigma is a quadratic twist of E): the equivariant matching is beta_i =
    // dtw * alpha_i directly -- no Galois-group computation needed.
    twistpair := jInvariant(E) eq jInvariant(Es);
    sgms := [];
    if twistpair then
        // root scaling beta = sc * alpha, sc in K (twist factor mod squares)
        for j0 in [1..3] do
            if al[1] eq 0 then break; end if;
            sc := gm[j0]/al[1];
            if not sc in K then continue; end if;
            prm := [ Index(gm, sc*al[i]) : i in [1..3] ];
            if 0 notin prm then
                Append(~sgms, SymmetricGroup(3)!prm);
                break;
            end if;
        end for;
    end if;
    if #sgms eq 0 then
        GA := Automorphisms(L);
        GK := [ gA : gA in GA | gA(L!wg) eq L!wg ];
        Gaction := func< rs0 | [[ Index(rs0, gA(rr)) : rr in rs0 ] : gA in GK ] >;
        for sgm0 in SymmetricGroup(3) do
            beta0 := [ gm[i^sgm0] : i in [1..3] ];
            if Gaction(al) eq Gaction(beta0) then Append(~sgms, sgm0); end if;
        end for;
    end if;
    aK := [[ al[i]-al[j] : j in [1..3]] : i in [1..3]];
    RtL := PolynomialRing(L);
    for sgm in sgms do
        beta := [ gm[i^sgm] : i in [1..3] ];
        bb := [[ beta[i]-beta[j] : j in [1..3]] : i in [1..3]];
        if al[1]*bb[3][2] + al[2]*bb[1][3] + al[3]*bb[2][1] eq 0 then continue; end if;
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
        // ---- descent (even sextic) ----
        a6 := Coefficient(SX,6); a4 := Coefficient(SX,4);
        a2e := Coefficient(SX,2); a0 := Coefficient(SX,0);
        if a6 eq 0 or a0 eq 0 or a2e eq 0 or a4 eq 0 then printf "DEGEN %o\n", lbl; continue; end if;
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
        if not haveM then printf "NODESCENTSHAPE %o\n", lbl; continue; end if;
        Mks := Matrix(K,2,2,[ sig(e) : e in Eltseq(MKK) ]);
        CCm := Mks*MKK;
        if not (CCm[1,2] eq 0 and CCm[2,1] eq 0 and CCm[1,1] eq CCm[2,2]) then
            printf "COCYCLEFAIL %o\n", lbl; continue;
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
        if FQ eq 0 then printf "DESCENTFAIL %o\n", lbl; continue; end if;
        // ---- integralize + twist sieve + exact torsion ----
        denl := LCM([ Denominator(co) : co in Coefficients(FQ) ]);
        fZ := Px!(FQ*denl^2);
        cint := GCD([ Integers()!co : co in Coefficients(fZ) ]);
        sq := 1;
        for pr in Factorization(cint) do sq *:= pr[1]^(2*(pr[2] div 2)); end for;
        fZ := Px![ co div sq : co in [Integers()!c : c in Coefficients(fZ)] ];
        if Discriminant(fZ) eq 0 then printf "SINGULAR %o\n", lbl; continue; end if;
        dsc := Integers()!Discriminant(fZ);
        goodp := [ p : p in PrimesInInterval(13, 300) | dsc mod p ne 0 ];
        tsupport := { pf[1] : pf in TrialDivision(AbsoluteValue(Dfld), 100) }
                    join {2, 3, 5, 7, 11, 13};
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
                    if #Jacobian(HyperellipticCurve(fp)) mod TORD ne 0 then ok := false; break; end if;
                end for;
                if ok and ntest ge 9 then
                    fT := Px!(d*fZ);
                    JT := Jacobian(HyperellipticCurve(fT));
                    I := Invariants(TorsionSubgroup(JT));
                    printf "TORSION %o d=%o invs=%o g=%o\n", lbl, d, I, fT;
                    if not I in KNOWN and #I gt 0 then
                        nhit +:= 1;
                        printf "HIT lmfdbglue %o invs=%o (NOT IN KNOWN)\n", lbl, I;
                    end if;
                end if;
            end for;
        end for;
    end for;
    if nrow mod 500 eq 0 then
        printf "PROGRESS rows=%o surv2=%o surv3=%o survt3=%o hit=%o %os\n",
            nrow, nsurv2, nsurv3, nsurvt3, nhit, Cputime()-t0c;
    end if;
end for;
printf "SEARCH_DONE lmfdbglue part %o/%o rows=%o surv2=%o surv3=%o survt3=%o hits=%o warns=%o %.1o s\n",
    Part, NParts, nrow, nsurv2, nsurv3, nsurvt3, nhit, nwarn, Cputime()-t0c;
quit;
