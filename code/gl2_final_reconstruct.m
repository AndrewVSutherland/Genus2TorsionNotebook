// 37-hunt FINAL reconstruction: the raw Periods basis of the 2190.2.a.v
// piece carries an integral UNIMODULAR alternating Riemann form
//   E_true = [0,0,-2,-1 | 0,0,-1,-1 | 2,1,0,0 | 1,1,0,0]
// (derived from the periods themselves; results/gl2_reconstruct_37j.log),
// i.e. the period lattice is already principally polarized.  Reconstruct
// the corresponding genus-2 curve over Q, match the quadratic twist to the
// LMFDB traces of f, and compute the torsion.
//
// Run: magma -b Prec:=60 code/gl2_final_reconstruct.m > results/gl2_final_reconstruct.log

SetColumns(0);
SetSeed(1);
if not assigned Prec then Prec := 60; elif Type(Prec) eq MonStgElt then Prec := StringToInteger(Prec); end if;
SetMemoryLimit(24*10^9);

Attach("~/.claude/jobs/a1db5dd4/tmp/polredabs_shim.m");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/endomorphisms/endomorphisms/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/quartic/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/curve_reconstruction/magma/spec");
SetVerbose("CurveRec", 1);
SetVerbose("EndoFind", 1);

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

// star involution transported to the LATTICE (= Periods) frame:
// SIL = Cm * SIm * Cm^-1, Cm = Solution(BAf, Lb).  Verified integral and
// different from the Basis-frame matrix (results/gl2_star_frame_check.log).
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

// the two integral Riemann forms on the raw Periods basis (frame 1),
// from results/gl2_reconstruct_37j.log:
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
P0 := Transpose(PMb);   // 2x4
// find positive-definite principal combinations
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
    edE := ElementaryDivisors(EEZ);
    printf "POSDEF a=%o b=%o edivs=%o\n", a, b, edE;
    Append(~Ecands, <a, b, EEZ>);
end for; end for;
printf "N_POSDEF %o\n", #Ecands;
error if #Ecands eq 0, "no positive combination";

// take the smallest-type positive forms and principalize by isotropic
// modification: for E of type (1,d) (edivs [1,1,d,d]), the dual quotient
// L^#/L = (Z/d)^2; each maximal isotropic (cyclic order-d) subgroup H
// gives L' = L + lift(H)/d with E principal on L'.
RQ<xq> := PolynomialRing(Rationals());
FX3 := NumberFieldExtra(xq^2 - 3);   // Q(sqrt 3) with recognition attributes
function TryReconstruct(bigP)
    // attempt Q first, then Q(sqrt3) (natural field for RM-sqrt3 members)
    try
        Crec, hKL, okb := ReconstructCurve(bigP, FX : Base := true);
        if Type(Crec) in {CrvHyp, Crv} then return true, Crec; end if;
    catch e ; end try;
    try
        Crec3, hKL3, okb3 := ReconstructCurve(bigP, FX3 : Base := true);
        if Type(Crec3) in {CrvHyp, Crv} then
            printf "  SQRT3_CURVE over %o : %o\n", BaseRing(Crec3), Crec3;
        end if;
    catch e ; end try;
    return false, false;
end function;

seen := {};
for ec in Ecands do
    EZp := ec[3];
    edp := ElementaryDivisors(EZp);
    d := edp[3];
    if d gt 5 then continue; end if;   // smallest types
    printf "TRY_FORM a=%o b=%o type=(1,%o)\n", ec[1], ec[2], d;
    // dual quotient generators
    Zd := Integers(d);
    Einvp := ChangeRing(EZp, Rationals())^-1;
    Gg := [ Vector(Zd, [ Zd!(Integers()!(d*x) mod d) : x in Eltseq(Einvp[i]) ]) : i in [1..4] ];
    Gm := sub< RSpace(Zd,4) | Gg >;
    EZd := ChangeRing(EZp, Zd);
    // cyclic order-d isotropic subgroups
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
            // conjugation (star) stability in the periods frame
            SId := ChangeRing(Matrix(Integers(),4,4,[Integers()!x : x in Eltseq(SIL)]), Zd);
            if &and[ u*SId in Hx : u in Generators(Hx) ] then
                Append(~subs, Hx);
            end if;
        end if;
    end for;
    printf "  isotropics: %o\n", #subs;
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
        printf "  PRINCIPAL_LATTICE found (index %o)\n", d;
        EpN := EpZ div edpp[1];
        for sgn in [1,-1] do
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
                printf "  bigP valid (sgn=%o blk=%o)\n", sgn, blk;
                okR, Crec := TryReconstruct(bigP);
                if not okR then printf "  recfail\n"; continue; end if;
                printf "REC_TYPE %o\n", Type(Crec);
                if Type(Crec) notin {CrvHyp, Crv} then
                    printf "REC_OBJECT %o\n", Crec; continue;
                end if;
                printf "RECONSTRUCTED over %o:\n%o\n", BaseRing(Crec), Crec;
                if Type(BaseRing(Crec)) ne FldRat then continue; end if;
                CQ := ChangeRing(Crec, Rationals());
                ii := IgusaInvariants(CQ);
                if ii in seen then continue; end if;
                Include(~seen, ii);
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
                        printf "TWIST_MATCH d=%o\n", dt;
                        Cmin := Cd;
                        try Cmin := ReducedMinimalWeierstrassModel(Cd); catch e ; end try;
                        Cs := SimplifiedModel(Cmin);
                        Tt := TorsionSubgroup(Jacobian(Cs));
                        printf "TORSION %o (order %o)\n", Invariants(Tt), #Tt;
                        if #Tt mod 37 eq 0 then
                            printf "*** THEOREM: genus-2 Jacobian over Q with a rational point of order 37 ***\n";
                            printf "CURVE %o\n", Cmin;
                        end if;
                        break;
                    end if;
                end for;
            end for;
        end for;
    end for;
end for;
printf "GL2_FINAL_DONE\n";
quit;
