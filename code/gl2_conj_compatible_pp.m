// 37-hunt v14: conjugation-COMPATIBLE principal polarizations.
// A PP member (T, E) is defined over R only if star preserves the lattice
// AND acts antisymplectically on E:  S E S^t = -E.  Compute the exact star
// action on the NS span <E1,E2>, restrict to its (-1)-eigenray, and only
// principalize those forms -- on L itself (Bcontrol) and on the three
// Eisenstein 37-quotient lattices.  Numeric invariants must then be REAL;
// recognize with PowerRelation; Mestre + twist + torsion on rational hits.
//
// Run: magma -b Prec:=100 code/gl2_conj_compatible_pp.m > results/gl2_conj_compatible_pp.log

SetColumns(0);
SetSeed(1);
if not assigned Prec then Prec := 100; elif Type(Prec) eq MonStgElt then Prec := StringToInteger(Prec); end if;
SetMemoryLimit(24*10^9);

Attach("~/.claude/jobs/a1db5dd4/tmp/polredabs_shim.m");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/endomorphisms/endomorphisms/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/quartic/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/curve_reconstruction/magma/spec");

SetDefaultRealField(RealField(Prec + 10));
FX := RationalsExtra(Prec);
CC := FX`CC;
Q := Rationals();

Lv := 2190;
PP := 37;
trtargets := [<7,2>, <11,-6>, <13,-2>, <17,6>, <19,10>];
badtargets := [<2,1>, <3,1>, <5,1>, <73,73>];

M := ModularSymbols(Lv, 2, 0);
S := CuspidalSubspace(M);
NS_ := NewSubspace(S);
D := NewformDecomposition(NS_);
target := 0;
for i in [1..#D] do
    if Dimension(D[i]) ne 4 then continue; end if;
    ok := true;
    for tt in trtargets do
        if Trace(HeckeOperator(D[i], tt[1])) ne 2*tt[2] then ok := false; break; end if;
    end for;
    if ok then target := i; break; end if;
end for;
error if target eq 0, "piece not found";
Af := D[target];
printf "PIECE %o\n", target;

BAf := ChangeRing(BasisMatrix(VectorSpace(Af)), Q);
Lb := ChangeRing(BasisMatrix(Lattice(Af)), Q);
Cm := Solution(BAf, Lb);
SIm := ChangeRing(StarInvolution(Af), Q);
SIL := Cm * SIm * Cm^-1;
error if not &and[Denominator(v) eq 1 : v in Eltseq(SIL)], "SIL not integral";

// Eisenstein kernel
F := GF(PP);
V4 := VectorSpace(F, 4);
W := V4;
ells := [ l : l in [7,11,13,17,19,23,29,31,41,43,47,53] | Lv mod l ne 0 ];
for l in ells cat [bt[1] : bt in badtargets] do
    if Lv mod l eq 0 then
        av := [ bt[2] : bt in badtargets | bt[1] eq l ][1];
    else
        av := 1 + l;
    end if;
    Tl := HeckeOperator(Af, l);
    TlL := Cm * ChangeRing(Tl, Q) * Cm^-1;
    TlF := ChangeRing(Matrix(Integers(),4,4,[Integers()!x : x in Eltseq(TlL)]), F);
    W := W meet Kernel(TlF - ScalarMatrix(F, 4, F!av));
end for;
printf "EIS_KERNEL_DIM %o\n", Dimension(W);
error if Dimension(W) ne 2, "unexpected kernel dim";
SI37 := ChangeRing(Matrix(Integers(),4,4,[Integers()!x : x in Eltseq(SIL)]), F);
wplus := 0; wminus := 0;
Bw := Basis(W);
for co in CartesianPower([0..PP-1], 2) do
    wv := co[1]*Bw[1] + co[2]*Bw[2];
    if wv eq V4!0 then continue; end if;
    im := wv*SI37;
    if im eq wv and Type(wplus) eq RngIntElt then wplus := wv; end if;
    if im eq -wv and Type(wminus) eq RngIntElt then wminus := wv; end if;
end for;
printf "STAR_EIGVECS plus=%o minus=%o\n", wplus, wminus;

prs := [<1,2>,<1,3>,<1,4>,<2,3>,<2,4>,<3,4>];
f1c := [ 0, 1, -1, -1, 2, 0 ];
f2c := [ -3, 4, -2, 2, 0, 3 ];
function FMat(f)
    FM := ZeroMatrix(Q, 4, 4);
    for k in [1..6] do
        FM[prs[k][1]][prs[k][2]] := f[k];
        FM[prs[k][2]][prs[k][1]] := -f[k];
    end for;
    return FM;
end function;
F1 := FMat(f1c); F2 := FMat(f2c);
function Pf(FM) return FM[1][2]*FM[3][4] - FM[1][3]*FM[2][4] + FM[1][4]*FM[2][3]; end function;
E1 := Pf(F1)*F1^-1; E2 := Pf(F2)*F2^-1;

// exact star action on the E-span: S Ei S^t = A2[i][1] E1 + A2[i][2] E2
function SpanCoords(X, E1, E2)
    // solve X = a E1 + b E2 (16 eqns, 2 unknowns), exact
    v1 := Vector(Q, Eltseq(E1)); v2 := Vector(Q, Eltseq(E2)); vx := Vector(Q, Eltseq(X));
    Mv := Matrix(Q, 2, 16, Eltseq(v1) cat Eltseq(v2));
    sol, ker := Solution(Mv, vx);
    error if Dimension(ker) ne 0, "span solve not unique";
    // verify
    error if sol[1]*E1 + sol[2]*E2 ne X, "not in span";
    return sol;
end function;
SE1 := SIL*ChangeRing(E1,Q)*Transpose(SIL);
SE2 := SIL*ChangeRing(E2,Q)*Transpose(SIL);
r1 := SpanCoords(SE1, ChangeRing(E1,Q), ChangeRing(E2,Q));
r2 := SpanCoords(SE2, ChangeRing(E1,Q), ChangeRing(E2,Q));
A2 := Matrix(Q, 2, 2, [r1[1], r1[2], r2[1], r2[2]]);
printf "STAR_NS_ACTION %o\n", Eltseq(A2);
// (-1)-eigenvectors: rows v with v*A2 = -v  (v = (a,b) coefficient rows)
Km := KernelMatrix(A2 + IdentityMatrix(Q,2));
printf "MINUS_EIGENSPACE dim %o : %o\n", Nrows(Km), Nrows(Km) gt 0 select Eltseq(Km[1]) else [];
Kp := KernelMatrix(A2 - IdentityMatrix(Q,2));
printf "PLUS_EIGENSPACE dim %o : %o\n", Nrows(Kp), Nrows(Kp) gt 0 select Eltseq(Kp[1]) else [];
error if Nrows(Km) eq 0, "no conjugation-compatible polarization ray";
// primitive integer generator of the compatible ray
den := LCM([Denominator(x) : x in Eltseq(Km[1])]);
ab0 := [ Integers()!(den*Km[1][1]), Integers()!(den*Km[1][2]) ];
g0 := GCD(ab0);
ab0 := [ ab0[1] div g0, ab0[2] div g0 ];
printf "COMPATIBLE_RAY (a,b)=%o\n", ab0;

nterms := Ceiling(Lv * Prec * Log(10)/(2*Pi(RealField(20)))) + 200;
printf "PERIOD_TERMS %o\n", nterms;
Pv := Periods(Af, nterms);
PMb := Matrix(CC, 4, 2, &cat[ [CC!Pv[i][1], CC!Pv[i][2]] : i in [1..4] ]);
printf "PERIODS ok\n";
RCC<xc> := PolynomialRing(CC);

procedure TryCompatible(Mb, SILB, tagbase, PMb, F1, F2, E1, E2, ab0)
    CC := BaseRing(PMb);
    Q := Rationals();
    PB := ChangeRing(Mb, CC) * PMb;
    P0old := Transpose(PMb);
    P0oldc := Matrix(CC,2,4,[Conjugate(x) : x in Eltseq(P0old)]);
    // compatible E-forms integral on this lattice: multiples k*(a0 E1 + b0 E2)
    E0 := ab0[1]*ChangeRing(E1,Q) + ab0[2]*ChangeRing(E2,Q);
    EB0 := Mb * E0 * Transpose(Mb);
    dmin := LCM([Denominator(x) : x in Eltseq(EB0)]);
    EBZ := Matrix(Integers(),4,4,[Integers()!(dmin*x) : x in Eltseq(EB0)]);
    gE := GCD([x : x in Eltseq(EBZ) | x ne 0]);
    EBZ := EBZ div gE;
    // positivity via F-side on old frame
    FE := ChangeRing(ab0[1]*F1 + ab0[2]*F2, Q);
    Hm := CC.1 * ChangeRing(P0old,CC) * ChangeRing(FE,CC) * Transpose(P0oldc);
    sgnE := 1;
    if not (Re(Hm[1][1]) gt 0 and Re(Hm[1][1]*Hm[2][2]-Hm[1][2]*Hm[2][1]) gt 0) then
        sgnE := -1;
        Hm := -Hm;
        error if not (Re(Hm[1][1]) gt 0 and Re(Hm[1][1]*Hm[2][2]-Hm[1][2]*Hm[2][1]) gt 0),
            "compatible form not definite either sign";
    end if;
    EBu := sgnE*EBZ;
    // antisymplectic sanity on this lattice
    anti := SILB*ChangeRing(EBu,Q)*Transpose(SILB) eq -ChangeRing(EBu,Q);
    edv := ElementaryDivisors(EBu);
    printf "%o COMPAT_FORM edivs=%o antisympl=%o\n", tagbase, edv, anti;
    d := edv[3] div edv[1];
    if edv[1] ne edv[2] or edv[3] ne edv[4] then
        printf "%o NONSTANDARD_TYPE %o (skip)\n", tagbase, edv;
        return;
    end if;
    latlist := [];
    if d eq 1 then
        latlist := [ IdentityMatrix(Q,4) ];
    else
        Zd := Integers(d);
        Einvp := ChangeRing(EBu, Q)^-1;
        Gg := [ Vector(Zd, [ Zd!(Integers()!(d*x) mod d) : x in Eltseq(Einvp[i]) ]) : i in [1..4] ];
        Gm := sub< RSpace(Zd,4) | Gg >;
        EZd := ChangeRing(EBu, Zd);
        SId := ChangeRing(Matrix(Integers(),4,4,[Integers()!x : x in Eltseq(SILB)]), Zd);
        seensub := {};
        for x in Gm do
            if x eq Gm!0 then continue; end if;
            xord := d div GCD([d] cat [Integers()!cc : cc in Eltseq(x)]);
            if xord ne d then continue; end if;
            Hx := sub< Gm | x >;
            key := { Eltseq(t) : t in Hx };
            if key in seensub then continue; end if;
            Include(~seensub, key);
            if not &and[ IsZero(InnerProduct(u*EZd, v)) : u, v in Generators(Hx) ] then continue; end if;
            if not &and[ u*SId in Hx : u in Generators(Hx) ] then continue; end if;
            Hgens := [ xx : xx in Generators(Hx) ];
            lifts := [ Vector(Q, [ (Integers()!cc)/d : cc in Eltseq(u) ]) : u in Hgens ];
            B := VerticalJoin(ChangeRing(IdentityMatrix(Integers(),4), Q),
                              Matrix(Q, #lifts, 4, &cat[Eltseq(v) : v in lifts]));
            dn := LCM([Denominator(x) : x in Eltseq(B)]);
            BZ := Matrix(Integers(), Nrows(B), 4, [Integers()!(dn*x) : x in Eltseq(B)]);
            Hn := HermiteForm(BZ);
            Bp := Matrix(Q, 4, 4, [ Hn[i][j]/dn : j in [1..4], i in [1..4] ]);
            Ep := Bp * ChangeRing(EBu, Q) * Transpose(Bp);
            if LCM([Denominator(x) : x in Eltseq(Ep)]) ne 1 then continue; end if;
            EpZ := Matrix(Integers(),4,4,[Integers()!x : x in Eltseq(Ep)]);
            edpp := ElementaryDivisors(EpZ);
            if edpp[1] ne edpp[4] then continue; end if;
            // FINAL antisymplectic check in the new frame
            SILnew := Bp * SILB * Bp^-1;
            if not &and[Denominator(v) eq 1 : v in Eltseq(SILnew)] then continue; end if;
            if SILnew*ChangeRing(EpZ,Q)*Transpose(SILnew) ne -ChangeRing(EpZ,Q) then
                printf "%o d%o: isotropic ok but NOT antisymplectic (skip)\n", tagbase, d;
                continue;
            end if;
            Append(~latlist, Bp);
        end for;
    end if;
    printf "%o N_COMPAT_PRINCIPAL %o (type d=%o)\n", tagbase, #latlist, d;
    nn := 0;
    for Bp in latlist do
        nn +:= 1;
        Ep := Bp * ChangeRing(EBu, Q) * Transpose(Bp);
        EpZ := Matrix(Integers(),4,4,[Integers()!x : x in Eltseq(Ep)]);
        edpp := ElementaryDivisors(EpZ);
        EpN := EpZ div edpp[1];
        tag := Sprintf("%o_d%o_n%o", tagbase, d, nn);
        donebp := false;
        for sgn in [1,-1] do
            if donebp then break; end if;
            _, T := FrobeniusFormAlternating(sgn*EpN);
            Om := ChangeRing(ChangeRing(T, Q) * Bp, CC) * PB;
            for blk in [1,2] do
                if blk eq 1 then
                    bigP := Matrix(CC, 2, 4, [ Om[i][j] : i in [1..4], j in [1..2] ]);
                else
                    ordr := [3,4,1,2];
                    bigP := Matrix(CC, 2, 4, [ Om[ordr[i]][j] : i in [1..4], j in [1..2] ]);
                end if;
                isb := false;
                try isb := IsBigPeriodMatrix(bigP); catch e ; end try;
                if not isb then continue; end if;
                printf "%o BIGP_OK sgn=%o blk=%o\n", tag, sgn, blk;
                try
                    tau := SmallPeriodMatrix(bigP);
                    taunew, gamma := ReduceSmallPeriodMatrix(tau);
                    Ag := Submatrix(gamma, 1,1, 2,2);
                    Bg := Submatrix(gamma, 1,3, 2,2);
                    Cg := Submatrix(gamma, 3,1, 2,2);
                    Dg := Submatrix(gamma, 3,3, 2,2);
                    Pnew := bigP * Transpose(BlockMatrix([[Ag, Bg], [Cg, Dg]]));
                    P2new := Submatrix(Pnew, 1,3, 2,2);
                    P2inew := P2new^(-1);
                    wsv := [ [[0,1],[0,1]], [[0,1],[1,1]], [[1,0],[1,0]],
                             [[1,0],[1,1]], [[1,1],[0,1]], [[1,1],[1,0]] ];
                    rats := [ ];
                    for wv in wsv do
                        w := (1/2)*taunew*Transpose(Matrix(CC, [wv[1]]))
                             + (1/2)*Transpose(Matrix(CC, [wv[2]]));
                        tds := ThetaDerivatives(taunew, w);
                        Hh := Matrix(CC, [ tds ]) * P2inew;
                        seq := Eltseq(Hh);
                        add := true;
                        if Abs(seq[2]) lt Abs(seq[1]) then
                            if Abs(seq[2]/seq[1])^2 lt CC`epscomp then add := false; end if;
                        end if;
                        if add then Append(~rats, -seq[1]/seq[2]); end if;
                    end for;
                    fCC := &*[ xc - rat : rat in rats ];
                    JJ := IgusaInvariants(fCC);
                    J2 := JJ[1]; J4 := JJ[2]; J6 := JJ[3]; J10 := JJ[5];
                    g2s := [ J2^5/J10, J2^3*J4/J10, J2^2*J6/J10 ];
                    maxim := Max([ Abs(Im(z)) / Max(Abs(z), CC!1) : z in g2s ]);
                    printf "%o INV_REALITY relIm=%o\n", tag, RealField(6)!maxim;
                    allQ := true;
                    g2Q := [ Q | ];
                    for k in [1..3] do
                        printf "%o G2n_%o %o\n", tag, k, ComplexField(15)!g2s[k];
                        reck := false;
                        for dg in [1,2,3,4] do
                            okp := false;
                            try
                                pol := PowerRelation(g2s[k], dg : Precision := Prec - 25);
                                okp := true;
                            catch e ; end try;
                            if not okp then continue; end if;
                            ht := Max([Abs(cc2) : cc2 in Coefficients(pol)]);
                            ev := Abs(Evaluate(ChangeRing(pol, CC), g2s[k]));
                            rel := ev / Max(Abs(g2s[k]), CC!1)^dg / ht;
                            if rel lt 10^(-12) and ht lt 10^((Prec-25) div 2) then
                                printf "%o INV_%o deg=%o pol=%o\n", tag, k, dg, pol;
                                if dg eq 1 then
                                    Append(~g2Q, -Coefficient(pol,0)/Coefficient(pol,1));
                                else
                                    allQ := false;
                                end if;
                                reck := true;
                                break;
                            end if;
                        end for;
                        if not reck then
                            printf "%o INV_%o UNREC\n", tag, k;
                            allQ := false;
                        end if;
                    end for;
                    if allQ and #g2Q eq 3 then
                        printf "%o *** RATIONAL_INVARIANTS %o ***\n", tag, g2Q;
                        CQ := HyperellipticCurveFromG2Invariants(g2Q);
                        if Type(BaseRing(CQ)) eq FldRat then
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
                        else
                            printf "%o MESTRE_OBSTRUCTION %o\n", tag, BaseRing(CQ);
                        end if;
                    end if;
                catch e
                    ms := Sprint(e`Object);
                    printf "%o NUMINV_ERR %o\n", tag, #ms gt 120 select Substring(ms,1,120) else ms;
                end try;
                donebp := true;
                break;
            end for;
        end for;
    end for;
end procedure;

choices := [* *];
Append(~choices, <"Bcontrol", []>);
if Type(wplus) ne RngIntElt then Append(~choices, <"Bplus", [wplus]>); end if;
if Type(wminus) ne RngIntElt then Append(~choices, <"Bminus", [wminus]>); end if;
Append(~choices, <"Bfull", [V4!b : b in Basis(W)]>);

for ch in choices do
    tagbase := ch[1];
    wws := ch[2];
    rows := [ ];
    for i in [1..4] do Append(~rows, [ Q | i eq j select 1 else 0 : j in [1..4] ]); end for;
    for ww in wws do
        Append(~rows, [ Q!(Integers()!ww[j])/PP : j in [1..4] ]);
    end for;
    B0 := Matrix(Q, #rows, 4, &cat rows);
    dn := LCM([Denominator(x) : x in Eltseq(B0)]);
    BZ := Matrix(Integers(), Nrows(B0), 4, [Integers()!(dn*x) : x in Eltseq(B0)]);
    Hn := HermiteForm(BZ);
    Mb := Matrix(Q, 4, 4, [ Hn[i][j]/dn : j in [1..4], i in [1..4] ]);
    printf "%o LATTICE_INDEX 1/%o\n", tagbase, Determinant(Mb)^-1;
    SILB := Mb * SIL * Mb^-1;
    if not &and[Denominator(v) eq 1 : v in Eltseq(SILB)] then
        printf "%o STAR_NOT_STABLE (skip)\n", tagbase;
        continue;
    end if;
    TryCompatible(Mb, SILB, tagbase, PMb, F1, F2, E1, E2, ab0);
end for;
printf "GL2_COMPAT_DONE\n";
quit;
