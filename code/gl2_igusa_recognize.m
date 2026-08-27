// 37-hunt v10: numeric G2/Igusa invariants of the star-stable principal
// lattices + LLL (PowerRelation) recognition of their minimal polynomials.
// This reads off the true field of moduli, bypassing ReconstructCurve's
// Base:=true coefficient algebraization (which fails whenever the natural
// sextic model needs a bigger field than K, even if the curve descends).
// Also tries ReconstructCurve WITHOUT Base (NumberFieldExtra mode).
//
// Run: magma -b Prec:=100 code/gl2_igusa_recognize.m > results/gl2_igusa_recognize.log

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

Lv := 2190;
trtargets := [<7,2>, <11,-6>, <13,-2>, <17,6>, <19,10>];

M := ModularSymbols(Lv, 2, 0);
S := CuspidalSubspace(M);
NS := NewSubspace(S);
D := NewformDecomposition(NS);
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

SIm := ChangeRing(StarInvolution(Af), Rationals());
BAf := ChangeRing(BasisMatrix(VectorSpace(Af)), Rationals());
Lb := ChangeRing(BasisMatrix(Lattice(Af)), Rationals());
Cm := Solution(BAf, Lb);
SIL := Cm * SIm * Cm^-1;
error if not &and[Denominator(v) eq 1 : v in Eltseq(SIL)], "SIL not integral";

nterms := Ceiling(Lv * Prec * Log(10)/(2*Pi(RealField(20)))) + 200;
printf "PERIOD_TERMS %o\n", nterms;
Pv := Periods(Af, nterms);
PMb := Matrix(CC, 4, 2, &cat[ [CC!Pv[i][1], CC!Pv[i][2]] : i in [1..4] ]);
printf "PERIODS ok\n";

prs := [<1,2>,<1,3>,<1,4>,<2,3>,<2,4>,<3,4>];
f1 := [ 0, 1, -1, -1, 2, 0 ];
f2 := [ -3, 4, -2, 2, 0, 3 ];
function FMat(f)
    FM := ZeroMatrix(Rationals(), 4, 4);
    for k in [1..6] do
        FM[prs[k][1]][prs[k][2]] := f[k];
        FM[prs[k][2]][prs[k][1]] := -f[k];
    end for;
    return FM;
end function;
F1 := FMat(f1); F2 := FMat(f2);
P0 := Transpose(PMb);
Ecands := [];
for a in [-12..12] do for b in [-12..12] do
    if a eq 0 and b eq 0 then continue; end if;
    FF := a*F1 + b*F2;
    if Determinant(FF) eq 0 then continue; end if;
    Hm := CC.1 * ChangeRing(P0,CC) * ChangeRing(FF,CC)
          * Transpose(Matrix(CC,2,4,[Conjugate(x) : x in Eltseq(P0)]));
    h11 := Re(Hm[1][1]); det2 := Re(Hm[1][1]*Hm[2][2] - Hm[1][2]*Hm[2][1]);
    if h11 le 0 or det2 le 0 then continue; end if;
    EE := FF^-1;
    dE := LCM([Denominator(x) : x in Eltseq(EE)]);
    EEZ := Matrix(Integers(),4,4,[Integers()!(dE*x) : x in Eltseq(EE)]);
    gE := GCD([x : x in Eltseq(EEZ) | x ne 0]);
    EEZ := EEZ div gE;
    Append(~Ecands, <a, b, EEZ>);
end for; end for;
printf "N_POSDEF %o\n", #Ecands;

RCC<xc> := PolynomialRing(CC);
procedure RecognizeInv(z, tag)
    for dg in [1,2,3,4,6] do
        okp := false;
        pol := RCC!0;
        try
            pol := PowerRelation(z, dg : Precision := Prec - 20);
            okp := true;
        catch e ; end try;
        if not okp then continue; end if;
        // plausibility: small height relative to precision, and actual vanishing
        ht := Max([Abs(c) : c in Coefficients(pol)]);
        ev := Abs(Evaluate(ChangeRing(pol, CC), z));
        rel := ev / Max(Abs(z), CC!1)^dg / ht;
        if rel lt 10^(-15) and ht lt 10^(Prec div 2) then
            printf "  %o deg=%o RECOGNIZED htdigits=%o pol=%o\n", tag, dg, Ceiling(Log(10, ht+1)), pol;
            return;
        end if;
    end for;
    printf "  %o UNRECOGNIZED (deg<=6)\n", tag;
end procedure;

seenlat := {};
for ec in Ecands do
    EZp := ec[3];
    edp := ElementaryDivisors(EZp);
    d := edp[3];
    if d gt 5 then continue; end if;
    Zd := Integers(d);
    Einvp := ChangeRing(EZp, Rationals())^-1;
    Gg := [ Vector(Zd, [ Zd!(Integers()!(d*x) mod d) : x in Eltseq(Einvp[i]) ]) : i in [1..4] ];
    Gm := sub< RSpace(Zd,4) | Gg >;
    EZd := ChangeRing(EZp, Zd);
    subs := [];
    seensub := {};
    for x in Gm do
        if x eq Gm!0 then continue; end if;
        xord := d div GCD([d] cat [Integers()!c : c in Eltseq(x)]);
        if xord ne d then continue; end if;
        Hx := sub< Gm | x >;
        key := { Eltseq(t) : t in Hx };
        if key in seensub then continue; end if;
        Include(~seensub, key);
        if &and[ IsZero(InnerProduct(u*EZd, v)) : u, v in Generators(Hx) ] then
            SId := ChangeRing(Matrix(Integers(),4,4,[Integers()!x : x in Eltseq(SIL)]), Zd);
            if &and[ u*SId in Hx : u in Generators(Hx) ] then
                Append(~subs, Hx);
            end if;
        end if;
    end for;
    for H in subs do
        Hgens := [ x : x in Generators(H) ];
        lifts := [ Vector(Rationals(), [ (Integers()!c)/d : c in Eltseq(u) ]) : u in Hgens ];
        B := VerticalJoin(ChangeRing(IdentityMatrix(Integers(),4), Rationals()),
                          Matrix(Rationals(), #lifts, 4, &cat[Eltseq(v) : v in lifts]));
        dn := LCM([Denominator(x) : x in Eltseq(B)]);
        BZ := Matrix(Integers(), Nrows(B), 4, [Integers()!(dn*x) : x in Eltseq(B)]);
        Hn := HermiteForm(BZ);
        Bp := Matrix(Rationals(), 4, 4, [ Hn[i][j]/dn : j in [1..4], i in [1..4] ]);
        Ep := Bp * ChangeRing(EZp, Rationals()) * Transpose(Bp);
        if LCM([Denominator(x) : x in Eltseq(Ep)]) ne 1 then continue; end if;
        EpZ := Matrix(Integers(),4,4,[Integers()!x : x in Eltseq(Ep)]);
        edpp := ElementaryDivisors(EpZ);
        if edpp[1] ne edpp[4] then continue; end if;
        latkey := <d, Eltseq(Bp)>;
        if latkey in seenlat then continue; end if;
        Include(~seenlat, latkey);
        EpN := EpZ div edpp[1];
        done := false;
        for sgn in [1,-1] do
            if done then break; end if;
            _, T := FrobeniusFormAlternating(sgn*EpN);
            Om := ChangeRing(ChangeRing(T, Rationals()) * Bp, CC) * PMb;
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
                printf "LATTICE a=%o b=%o type=(1,%o) sgn=%o blk=%o\n", ec[1], ec[2], d, sgn, blk;
                // numeric branch points via theta derivatives at the six
                // odd 2-torsion points (the ReconstructCurveG2 route), then
                // TWIST-INDEPENDENT absolute Igusa/G2 invariants.
                okros := false;
                fCC := RCC!0;
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
                    okros := true;
                catch e printf "  THETAFAIL %o\n", e`Object; end try;
                if okros then
                    JJ := IgusaInvariants(fCC);
                    J2 := JJ[1]; J4 := JJ[2]; J6 := JJ[3]; J10 := JJ[5];
                    g2s := [ J2^5/J10, J2^3*J4/J10, J2^2*J6/J10, J4^5/J10^2, J6^5/J10^3 ];
                    for k in [1..5] do
                        printf "  G2_%o abs=%o\n", k, ComplexField(15)!g2s[k];
                        RecognizeInv(g2s[k], Sprintf("G2_%o", k));
                    end for;
                end if;
                // free-field reconstruction (NumberFieldExtra mode)
                try
                    Yf, hf, okf := ReconstructCurve(bigP, FX);
                    if Type(Yf) in {CrvHyp, Crv} then
                        printf "  FREEFIELD_CURVE over %o : %o\n", BaseRing(Yf), Yf;
                    else
                        printf "  FREEFIELD_FAIL\n";
                    end if;
                catch e
                    ms := Sprint(e`Object);
                    printf "  FREEFIELD_ERR %o\n", #ms gt 120 select Substring(ms,1,120) else ms;
                end try;
                done := true;
                break;
            end for;
        end for;
    end for;
end for;
printf "GL2_IGUSA_DONE\n";
quit;
