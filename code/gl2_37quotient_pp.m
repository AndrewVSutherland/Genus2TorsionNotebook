// 37-hunt v12: principal polarizations on the EISENSTEIN QUOTIENT side.
// The rational 37-point lives on B = A_f/C (C = rational Eisenstein
// 37-subgroup), so enumerate principal lattices of L_B = L + (1/37)lift(C)
// (and the star companions), not of L.  For each: derive the integral NS
// forms on L_B from F1,F2, find positive types, star-stable principal
// members, reconstruct in free-field mode, take EXACT invariants, test Q.
//
// Run: magma -b Prec:=100 code/gl2_37quotient_pp.m > results/gl2_37quotient_pp.log

SetColumns(0);
SetSeed(1);
if not assigned Prec then Prec := 100; elif Type(Prec) eq MonStgElt then Prec := StringToInteger(Prec); end if;
SetMemoryLimit(24*10^9);

Attach("~/.claude/jobs/a1db5dd4/tmp/polredabs_shim.m");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/endomorphisms/endomorphisms/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/quartic/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/curve_reconstruction/magma/spec");
SetVerbose("CurveRec", 1);

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

// ---- Eisenstein kernel W in L/37L ----
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
    error if not &and[Denominator(v) eq 1 : v in Eltseq(TlL)], "TlL not integral";
    TlF := ChangeRing(Matrix(Integers(),4,4,[Integers()!x : x in Eltseq(TlL)]), F);
    W := W meet Kernel(TlF - ScalarMatrix(F, 4, F!av));
end for;
printf "EIS_KERNEL_DIM %o\n", Dimension(W);
error if Dimension(W) ne 2, "unexpected kernel dim";

// star eigenvectors in W
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

nterms := Ceiling(Lv * Prec * Log(10)/(2*Pi(RealField(20)))) + 200;
printf "PERIOD_TERMS %o\n", nterms;
Pv := Periods(Af, nterms);
PMb := Matrix(CC, 4, 2, &cat[ [CC!Pv[i][1], CC!Pv[i][2]] : i in [1..4] ]);
printf "PERIODS ok\n";

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

// E-side (intersection) forms: Pfaffian adjugate E = Pf(F)*F^-1, LINEAR in F
function Pf(F)
    return F[1][2]*F[3][4] - F[1][3]*F[2][4] + F[1][4]*F[2][3];
end function;
function PfAdj(F) return Pf(F)*F^-1; end function;
E1 := PfAdj(F1); E2 := PfAdj(F2);
// linearity guard
error if PfAdj(F1+F2) ne E1+E2, "Pfaffian adjugate not linear -- convention bug";
printf "E1 %o\n", Eltseq(E1);
printf "E2 %o\n", Eltseq(E2);

RCC<xc> := PolynomialRing(CC);

procedure ExamineCurve(Yf, tag)
    printf "%o CURVE over %o\n%o\n", tag, DefiningPolynomial(BaseRing(Yf)), Yf;
    if Type(BaseRing(Yf)) eq FldRat then
        printf "%o OVER_Q_DIRECT\n", tag;
    end if;
    g2 := G2Invariants(Yf);
    degs := [ Degree(MinimalPolynomial(g2[k], Q)) : k in [1..3] ];
    printf "%o G2_DEGREES %o\n", tag, degs;
    if Max(degs) eq 1 then
        g2Q := [ Q!(-Coefficient(MinimalPolynomial(g2[k], Q),0)) : k in [1..3] ];
        printf "%o G2_RATIONAL %o\n", tag, g2Q;
        CQ := HyperellipticCurveFromG2Invariants(g2Q);
        if Type(BaseRing(CQ)) ne FldRat then
            printf "%o MESTRE_OBSTRUCTION field %o\n", tag, BaseRing(CQ);
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
    end if;
end procedure;

procedure TryLattice(Mb, SILB, tagbase, PMb, F1, F2, E1, E2)
    CC := BaseRing(PMb);
    Q := Rationals();
    PB := ChangeRing(Mb, CC) * PMb;   // 4x2 periods of the new lattice
    P0old := Transpose(PMb);          // 2x4, OLD frame (positivity is frame-independent)
    P0oldc := Matrix(CC,2,4,[Conjugate(x) : x in Eltseq(P0old)]);
    // E-side integrality lattice: {(a,b) : Mb (a E1 + b E2) Mb^t integral}
    G1 := Mb * ChangeRing(E1,Q) * Transpose(Mb);
    G2m := Mb * ChangeRing(E2,Q) * Transpose(Mb);
    dd2 := LCM([Denominator(x) : x in Eltseq(G1)] cat [Denominator(x) : x in Eltseq(G2m)]);
    H1 := Matrix(Integers(), 1, 16, [Integers()!(dd2*x) : x in Eltseq(G1)]);
    H2 := Matrix(Integers(), 1, 16, [Integers()!(dd2*x) : x in Eltseq(G2m)]);
    Mcong := VerticalJoin(H1, H2);
    Kbig := KernelMatrix(VerticalJoin(Mcong, dd2*IdentityMatrix(Integers(),16)));
    Kab := Submatrix(Kbig, 1, 1, Nrows(Kbig), 2);
    Lab := Lattice(Kab);
    Bab := Basis(LLL(Lab));
    printf "%o AB_BASIS %o %o (dd2=%o)\n", tagbase, Eltseq(Bab[1]), Eltseq(Bab[2]), dd2;
    cands := [];
    seenE := {};
    for c1 in [-8..8] do for c2 in [-8..8] do
        if c1 eq 0 and c2 eq 0 then continue; end if;
        vab := c1*Bab[1] + c2*Bab[2];
        a := vab[1]; b := vab[2];
        EB := a*G1 + b*G2m;
        if LCM([Denominator(x) : x in Eltseq(EB)]) ne 1 then continue; end if;
        if Determinant(EB) eq 0 then continue; end if;
        EBZ := Matrix(Integers(),4,4,[Integers()!x : x in Eltseq(EB)]);
        gE := GCD([x : x in Eltseq(EBZ) | x ne 0]);
        EBZ2 := EBZ div gE;
        // check the DIVIDED form is still integral-lattice-consistent: it is E-side on L_B by construction only for the full EB; keep primitive rep only if (a,b)/gE still in Lab
        if gE gt 1 then
            va := Vector(Rationals(), [a/gE, b/gE]);
            if Denominator(va[1]) eq 1 and Denominator(va[2]) eq 1 then
                if Vector(Integers(), [Integers()!va[1], Integers()!va[2]]) in Lab then
                    continue; // primitive version appears elsewhere in the scan
                end if;
            end if;
            EBZ := EBZ2; a := a/gE; b := b/gE; // primitive E integral on L_B only via gcd division: E entries still integers
        end if;
        if Eltseq(EBZ) in seenE then continue; end if;
        // positivity via F-side in the OLD frame: F(a,b) = a F1 + b F2 (Pf-adjugate linearity)
        FE := ChangeRing(a*F1 + b*F2, Q);
        Hm := CC.1 * ChangeRing(P0old,CC) * ChangeRing(FE,CC) * Transpose(P0oldc);
        h11 := Re(Hm[1][1]); det2 := Re(Hm[1][1]*Hm[2][2] - Hm[1][2]*Hm[2][1]);
        sgnE := 0;
        if h11 gt 0 and det2 gt 0 then sgnE := 1;
        else
            Hm := -Hm;
            if Re(Hm[1][1]) gt 0 and Re(Hm[1][1]*Hm[2][2]-Hm[1][2]*Hm[2][1]) gt 0 then sgnE := -1; end if;
        end if;
        if sgnE eq 0 then continue; end if;
        EBu := sgnE*EBZ;
        Include(~seenE, Eltseq(EBZ));
        edv := ElementaryDivisors(EBu);
        Append(~cands, <edv[3] div edv[1], a, b, EBu>);
    end for; end for;
    Sort(~cands, func<u,v | u[1] - v[1]>);
    printf "%o N_POSFORMS %o types %o\n", tagbase, #cands, Sort(Setseq({c[1] : c in cands}));
    ntried := 0;
    for c in cands do
        d := c[1]; EuZ := c[4];
        if d gt 500 then continue; end if;
        if ntried ge 8 then break; end if;
        latlist := [];
        if d eq 1 then
            latlist := [ IdentityMatrix(Q,4) ];
        else
            Zd := Integers(d);
            Einvp := ChangeRing(EuZ, Q)^-1;
            Gg := [ Vector(Zd, [ Zd!(Integers()!(d*x) mod d) : x in Eltseq(Einvp[i]) ]) : i in [1..4] ];
            Gm := sub< RSpace(Zd,4) | Gg >;
            EZd := ChangeRing(EuZ, Zd);
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
                Ep := Bp * ChangeRing(EuZ, Q) * Transpose(Bp);
                if LCM([Denominator(x) : x in Eltseq(Ep)]) ne 1 then continue; end if;
                EpZ := Matrix(Integers(),4,4,[Integers()!x : x in Eltseq(Ep)]);
                edpp := ElementaryDivisors(EpZ);
                if edpp[1] ne edpp[4] then continue; end if;
                Append(~latlist, Bp);
            end for;
        end if;
        for Bp in latlist do
            ntried +:= 1;
            Ep := Bp * ChangeRing(EuZ, Q) * Transpose(Bp);
            EpZ := Matrix(Integers(),4,4,[Integers()!x : x in Eltseq(Ep)]);
            edpp := ElementaryDivisors(EpZ);
            EpN := EpZ div edpp[1];
            tag := Sprintf("%o_d%o_a%o_b%o_n%o", tagbase, d, c[2], c[3], ntried);
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
                        Yf, hf, okf := ReconstructCurve(bigP, FX);
                        if Type(Yf) in {CrvHyp, Crv} then
                            ExamineCurve(Yf, tag);
                        else
                            printf "%o FREEFIELD_FAIL\n", tag;
                        end if;
                    catch e
                        ms := Sprint(e`Object);
                        printf "%o FREEFIELD_ERR %o\n", tag, #ms gt 120 select Substring(ms,1,120) else ms;
                    end try;
                    donebp := true;
                    break;
                end for;
            end for;
        end for;
    end for;
end procedure;

// build quotient lattices for w+, w-, and full W
choices := [* *];
// control validated (types [2,3,11]; d2 members = the cubic curve) -- removed
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
    TryLattice(Mb, SILB, tagbase, PMb, F1, F2, E1, E2);
end for;
printf "GL2_37Q_DONE\n";
quit;
