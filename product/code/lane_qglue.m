// lane_qglue.m — Q-simple geometrically-split realizations from QUADRATIC
// points of X1(N), N in {11, 13} (AVS relaxation 2026-08-26: factors over
// quadratic fields allowed).  Pipeline (validated on the [11] hit, s0=4/5):
//
//   1. derive the raw Kubert model F_N(r,s) = 0 at runtime as the relevant
//      factor of psi_N(0) for E(b,c), (b,c) = (r*s*(r-1), s*(r-1)), P=(0,0);
//      keep factors of r-degree 2 (rational s-fibers = quadratic points);
//   2. for rational s0 with nonsquare r-discriminant D(s0): E over
//      K = Q(sqrt D) with E(K) >= Z/N;
//   3. sigma-congruence test: E[2] ~ E^sigma[2] over K iff the conjugate
//      2-division cubic has a root in K[x]/(cubic) (both irreducible: a
//      K-rational 2-torsion point would give Z/2N, dead over quadratic);
//   4. on match: HLP/BHLS (2,2)-glue over the splitting field, sextic lands
//      in K[x] for the Galois-equivariant matching(s);
//   5. Weil descent: the glue output is an EVEN sextic; the conjugation acts
//      as x -> kap/x or x -> lam*x on coefficients; matrix Hilbert-90
//      N = A + M*A^sigma, substitute, scalar-H90 rescale -> model over Q;
//   6. twist sieve (N | #J_d(F_p) at >= 10 good primes) then exact
//      TorsionSubgroup.  HIT lines carry the full model.
//
// Usage (from product/code/):
//   magma -b N:=11 SH:=30 Part:=1 NParts:=1 lane_qglue.m > ../logs/qglue11_p1.log
SetColumns(0);
if not assigned MemGB then MemGB := 12; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned N then N := 11; elif Type(N) eq MonStgElt then N := StringToInteger(N); end if;
if not assigned SH then SH := 30; elif Type(SH) eq MonStgElt then SH := StringToInteger(SH); end if;
if not assigned Part then Part := 1; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;

Q := Rationals();
Px<x> := PolynomialRing(Q);

// ---- 1. raw model ----
Qrs<r,s> := FunctionField(Q, 2);
cS := s*(r-1); bS := r*cS;
Egen := EllipticCurve([1-cS, -bS, -bS, 0, 0]);
psiN := DivisionPolynomial(Egen, N);
F := Numerator(Evaluate(psiN, 0));
P2<R,S> := PolynomialRing(Q, 2);
hh := hom< Parent(F) -> P2 | [R, S] >;
FP := hh(F);
cands := [ fc[1] : fc in Factorization(FP) | Degree(fc[1], R) eq 2 and Degree(fc[1], S) ge 1 ];
printf "raw-model candidate factors (deg_R = 2): %o\n", #cands;
// branch A (genus-0-in-r): factor with deg_R = 2 whose quadratic points give
// order-N curves.  branch B (genus 2, e.g. N = 13): hyperelliptic fibers.
FN := P2!0;
for fc in cands do
    okf := 0;
    for sv in [Q| 4/5, 7/3, -2/7 ] do
        fr := UnivariatePolynomial(Evaluate(fc, [P2.1, P2!sv]));
        if Degree(fr) ne 2 then continue; end if;
        dd := Discriminant(fr);
        if dd eq 0 or IsSquare(dd) then continue; end if;
        dn := Numerator(dd)*Denominator(dd);
        K0 := QuadraticField(Squarefree(dn));
        rts := Roots(PolynomialRing(K0)!fr);
        if #rts eq 0 then continue; end if;
        rv := rts[1][1];
        cv := sv*(rv-1); bv := rv*cv;
        okE := true; E0 := 0;
        try E0 := EllipticCurve([1-cv, -bv, -bv, 0, 0]); catch e; okE := false; end try;
        if not okE then continue; end if;
        okO := false;
        try okO := Order(E0![0,0]) eq N; catch e; end try;
        if okO then okf +:= 1; end if;
    end for;
    if okf ge 2 then FN := fc; break; end if;
end for;
hyper := FN eq 0;
XH := 0; mph := 0; hpol := Px!0;
F13a := P2!0; XPa := P2!0; ZPa := P2!0;
if hyper then
    printf "no deg_R=2 factor; trying hyperelliptic branch\n";
    bigf := [ fc[1] : fc in Factorization(FP) | Degree(fc[1],R) ge 2 and Degree(fc[1],S) ge 2 ];
    error if #bigf eq 0, "no candidate raw-model factor for N =", N;
    A2q<uu,vv> := AffineSpace(Q,2);
    hcr := hom< P2 -> CoordinateRing(A2q) | [uu,vv] >;
    CX := ProjectiveClosure(Curve(A2q, hcr(bigf[1])));
    error if Genus(CX) ne 2, "raw curve genus", Genus(CX), "unhandled";
    okh, XH, mph := IsHyperelliptic(CX);
    error if not okh, "genus-2 raw curve not seen as hyperelliptic";
    hpol := Px!HyperellipticPolynomials(XH);
    printf "hyperelliptic branch: y^2 = %o\n", hpol;
    // the x-coordinate of the hyperelliptic map as (r,s)-polynomials (z=1)
    DPm := DefiningPolynomials(mph);
    RA := CoordinateRing(Ambient(CX));
    hz1 := hom< RA -> P2 | [P2.1, P2.2, 1] >;
    XPa := hz1(DPm[1]); ZPa := hz1(DPm[3]);
    F13a := 0;
    for fc in Factorization(FP) do
        if Degree(fc[1],R) ge 2 and Degree(fc[1],S) ge 2 then F13a := fc[1]; break; end if;
    end for;
else
    printf "raw model verified: %o\n", FN;
end if;

// produce a quadratic point (K, b0, c0) from the scan parameter s0
function GetPoint(s0)
    if not hyper then
        fr := UnivariatePolynomial(Evaluate(FN, [P2.1, P2!s0]));
        if Degree(fr) ne 2 then return false, 0, 0, 0, 0; end if;
        D := Discriminant(fr);
        if D eq 0 or IsSquare(D) then return false, 0, 0, 0, 0; end if;
        dn := Numerator(D)*Denominator(D);
        K := QuadraticField(Squarefree(dn));
        rts := Roots(PolynomialRing(K)!fr);
        if #rts eq 0 then return false, 0, 0, 0, 0; end if;
        r0 := rts[1][1];
        c0 := s0*(r0-1); b0 := r0*c0;
        return true, K, b0, c0, Squarefree(dn);
    else
        hv := Evaluate(hpol, s0);
        if hv eq 0 or IsSquare(hv) then return false, 0, 0, 0, 0; end if;
        dn := Numerator(hv)*Denominator(hv);
        K := QuadraticField(Squarefree(dn));
        oky, yv := IsSquare(K!hv);
        if not oky then return false, 0, 0, 0, 0; end if;
        // fiber system {F13 = 0, X(r,s) = s0 * Z(r,s)} solved over K
        R2K := PolynomialRing(K, 2);
        h2 := hom< P2 -> R2K | [R2K.1, R2K.2] >;
        V := [];
        try
            IK := ideal< R2K | h2(F13a), h2(XPa) - s0*h2(ZPa) >;
            V := Variety(IK);
        catch e;
            return false, 0, 0, 0, 0;
        end try;
        for pv in V do
            r0 := pv[1]; sr := pv[2];
            if r0 in Q and sr in Q then continue; end if;   // cusps
            c0 := sr*(r0-1); b0 := r0*c0;
            if b0 eq 0 then continue; end if;
            okE := true; E0 := 0;
            try E0 := EllipticCurve([1-c0, -b0, -b0, 0, 0]); catch e; okE := false; end try;
            if not okE then continue; end if;
            okO := false;
            try okO := Order(E0![0,0]) eq N; catch e; end try;
            if okO then return true, K, b0, c0, Squarefree(dn); end if;
        end for;
        return false, 0, 0, 0, 0;
    end if;
end function;

function HeightRats(H)
    S0 := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S0, Q!a/b); Append(~S0, Q!-a/b); end if;
    end for;
    return S0;
end function;

svals := HeightRats(SH);
npts := 0; nmatch := 0; nhit := 0;
t0c := Cputime();

for si in [1..#svals] do
    if si mod NParts ne (Part - 1) then continue; end if;
    s0 := svals[si];
    okpt, K, b0, c0, dK := GetPoint(s0);
    if not okpt or b0 eq 0 then continue; end if;
    w := K.1;
    RK := PolynomialRing(K);
    E := 0; okE := true;
    try E := EllipticCurve([1-c0, -b0, -b0, 0, 0]); catch e; okE := false; end try;
    if not okE then continue; end if;
    okO := false;
    try okO := Order(E![0,0]) eq N; catch e; end try;
    if not okO then continue; end if;
    npts +:= 1;
    sig := hom< K -> K | -w >;
    // 2-division cubic match
    a1 := 1-c0;
    b2c := a1^2 - 4*b0; b4c := -a1*b0; b6c := b0^2;
    cub := RK.1^3 + b2c*RK.1^2 + 8*b4c*RK.1 + 16*b6c;
    cubs := RK![ sig(Coefficient(cub,i)) : i in [0..3] ];
    if not IsIrreducible(cub) then printf "REDUCIBLE2 s0=%o (unexpected)\n", s0; continue; end if;
    LC := ext< K | cub >;
    if #Roots(PolynomialRing(LC)!cubs) eq 0 then continue; end if;
    nmatch +:= 1;
    printf "MATCH s0=%o D=%o\n", s0, dK;
    // ---- glue over K ----
    Es := EllipticCurve([sig(a) : a in aInvariants(E)]);
    f1 := HyperellipticPolynomials(WeierstrassModel(E));
    g1 := HyperellipticPolynomials(WeierstrassModel(Es));
    okL := true; L := 0;
    try L := SplittingField(f1); catch e; okL := false; end try;
    if not okL then printf "SPLITFAIL s0=%o\n", s0; continue; end if;
    RL := PolynomialRing(L);
    al := [ rr[1] : rr in Roots(RL!f1) ];
    gm := [ rr[1] : rr in Roots(RL!g1) ];
    if #al ne 3 or #gm ne 3 then printf "ROOTFAIL s0=%o\n", s0; continue; end if;
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
        GK := [ gA : gA in GA | gA(L!w) eq L!w ];
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
        if a6 eq 0 or a0 eq 0 or a2e eq 0 or a4 eq 0 then printf "DEGEN s0=%o\n", s0; continue; end if;
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
        if not haveM then printf "NODESCENTSHAPE s0=%o\n", s0; continue; end if;
        Mks := Matrix(K,2,2,[ sig(e) : e in Eltseq(MKK) ]);
        CCm := Mks*MKK;
        if not (CCm[1,2] eq 0 and CCm[2,1] eq 0 and CCm[1,1] eq CCm[2,2]) then
            printf "COCYCLEFAIL s0=%o\n", s0; continue;
        end if;
        FQ := Px!0;
        for Aseq in [[K|1,0,0,1],[K|1,1,0,1],[K|1,w,0,1],[K|1,0,w,1],[K|2,w,1,1],[K|1,2,w,1]] do
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
            for ccx in [K| 1, w, 1+w, 2+w, 1-2*w, 3+w, 1+2*w ] do
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
        if FQ eq 0 then printf "DESCENTFAIL s0=%o\n", s0; continue; end if;
        // ---- integralize + twist sieve + exact torsion ----
        denl := LCM([ Denominator(co) : co in Coefficients(FQ) ]);
        fZ := Px!(FQ*denl^2);
        cint := GCD([ Integers()!co : co in Coefficients(fZ) ]);
        sq := 1;
        for pr in Factorization(cint) do sq *:= pr[1]^(2*(pr[2] div 2)); end for;
        fZ := Px![ co div sq : co in [Integers()!c : c in Coefficients(fZ)] ];
        if Discriminant(fZ) eq 0 then printf "SINGULAR s0=%o\n", s0; continue; end if;
        dsc := Integers()!Discriminant(fZ);
        goodp := [ p : p in PrimesInInterval(13, 300) | p ne N and dsc mod p ne 0 ];
        tsupport := { pf[1] : pf in TrialDivision(AbsoluteValue(dK), 100) }
                    join {2, 5, N};
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
                    if #Jacobian(HyperellipticCurve(fp)) mod N ne 0 then ok := false; break; end if;
                end for;
                if ok and ntest ge 9 then
                    fT := Px!(d*fZ);
                    JT := Jacobian(HyperellipticCurve(fT));
                    I := Invariants(TorsionSubgroup(JT));
                    printf "TORSION s0=%o d=%o invs=%o g=%o\n", s0, d, I, fT;
                    if N in I or (#I gt 0 and I[#I] mod N eq 0) then
                        nhit +:= 1;
                        printf "HIT qglue N=%o s0=%o K=Q(sqrt %o) invs=%o\n", N, s0, dK, I;
                    end if;
                end if;
            end for;
        end for;
    end for;
    if npts mod 100 eq 0 then
        printf "PROGRESS pts=%o match=%o hit=%o %os\n", npts, nmatch, nhit, Cputime()-t0c;
    end if;
end for;
printf "SEARCH_DONE qglue N=%o SH=%o part %o/%o pts=%o matches=%o hits=%o %.1o s\n",
    N, SH, Part, NParts, npts, nmatch, nhit, Cputime()-t0c;
quit;
