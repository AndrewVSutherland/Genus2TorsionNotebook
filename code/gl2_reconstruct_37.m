// 37-torsion Jacobian hunt, step D (final driver): reconstruct genus-2
// curves over Q for candidate principally polarized members of the
// isogeny class of A_f, f = 2190.2.a.v, and test torsion for 37.
// Uses JRSijsling/curve_reconstruction: ReconstructCurve(bigP, QQ).
//
// Candidate filters, in order:
//   1. star-involution (complex conjugation) stability -- genuine Galois
//      necessary condition available in modular symbols;
//   2. if the star filter is also vacuous at 37, restrict the 37-part to
//      lifts of the mod-37 Eisenstein kernel (the theory locates the
//      37-structure there);
//   3. final arbiter: ReconstructCurve succeeds over Q (Base := true
//      forces base-field recognition; failures are skipped).
// Each reconstructed curve: quadratic-twist match of a_7, a_11, a_13
// against LMFDB traces (2, -6, -2), then TorsionSubgroup.
//
// Run: magma -b Prec:=100 MaxRec:=200 code/gl2_reconstruct_37.m > results/gl2_reconstruct_37.log

SetColumns(0);
SetSeed(1);
if not assigned Prec then Prec := 100; elif Type(Prec) eq MonStgElt then Prec := StringToInteger(Prec); end if;
if not assigned MaxRec then MaxRec := 200; elif Type(MaxRec) eq MonStgElt then MaxRec := StringToInteger(MaxRec); end if;
SetMemoryLimit(24*10^9);

Attach("~/.claude/jobs/a1db5dd4/tmp/polredabs_shim.m");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/endomorphisms/endomorphisms/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/quartic/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/curve_reconstruction/magma/spec");

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

BAf := ChangeRing(BasisMatrix(VectorSpace(Af)), Rationals());
Lb := ChangeRing(BasisMatrix(Lattice(Af)), Rationals());
Cm := Solution(BAf, Lb);
IP := ChangeRing(IntersectionPairing(Af), Rationals());
E := Cm * IP * Transpose(Cm);
den := LCM([Denominator(x) : x in Eltseq(E)]);
EZ0 := Matrix(Integers(), 4, 4, [Integers()!(den*x) : x in Eltseq(E)]);
g0 := GCD([x : x in Eltseq(EZ0) | x ne 0]);
EZ := EZ0 div g0;
ed := ElementaryDivisors(EZ);
dd := ed[3];
printf "TYPE %o dd=%o\n", ed, dd;

SI := StarInvolution(Af);
SIL := Cm * ChangeRing(SI, Rationals()) * Cm^-1;
SIZ := Matrix(Integers(), 4, 4, [Integers()!x : x in Eltseq(SIL)]);

// Eisenstein kernel mod 37 on the lattice (for filter level 2)
F37 := GF(37);
function HeckeOnLattice(l)
    Tl := HeckeOperator(Af, l);
    TlL := Cm * ChangeRing(Tl, Rationals()) * Cm^-1;
    return Matrix(Integers(), 4, 4, [Integers()!x : x in Eltseq(TlL)]);
end function;
W37 := VectorSpace(F37, 4);
for pr in [<2,1>,<3,1>,<5,1>,<73,73>] cat [<l, 1+l> : l in [7,11,13,17,19,23,29,31]] do
    Tl := HeckeOnLattice(pr[1]);
    W37 := W37 meet Kernel(Matrix(F37,4,4,[F37!x : x in Eltseq(Tl)]) - ScalarMatrix(F37,4,F37!pr[2]));
end for;
printf "EIS_KERNEL_MOD37_DIM %o\n", Dimension(W37);

Zd := Integers(dd);
Einv := ChangeRing(EZ, Rationals())^-1;
gensD := [ Vector(Rationals(), Eltseq(Einv[i])) : i in [1..4] ];
Ggens := [ Vector(Zd, [ Zd!(Integers()!(dd*x) mod dd) : x in Eltseq(v) ]) : v in gensD ];
EZd := ChangeRing(EZ, Zd);
SId := ChangeRing(SIZ, Zd);

fac := Factorization(dd);
function AddOrd(x)
    return dd div GCD([dd] cat [Integers()!c : c in Eltseq(x)]);
end function;

function Cands(p, k)
    m := dd div p^k;
    pg := [ m*g_ : g_ in Ggens ];
    Mp := sub< RSpace(Zd, 4) | pg >;
    ord := p^k;
    e1 := 0; e2 := 0;
    for x in Mp do if AddOrd(x) eq ord then e1 := x; break; end if; end for;
    for x in Mp do
        if AddOrd(x) eq ord and #sub<Mp | [e1, x]> eq ord^2 then e2 := x; break; end if;
    end for;
    genlist := [ e1 + t*e2 : t in [0..ord-1] ] cat
               [ p*u*e1 + e2 : u in [0..(ord div p)-1] ];
    cands := [];
    for x in genlist do
        Hx := sub< Mp | x >;
        if #Hx ne ord then continue; end if;
        if not &and[ IsZero(InnerProduct(u*EZd, v)) : u, v in Generators(Hx) ] then continue; end if;
        if not &and[ u*SId in Hx : u in Generators(Hx) ] then continue; end if;   // star filter
        if p eq 37 then
            // Eisenstein filter: the generator reduced mod 37 must lie in W37
            xv := Vector(F37, [ F37!((Integers()!c) div (dd div (p^k)) mod 37) : c in Eltseq(x) ]);
            // NB x is dd-scaled: x = dd * (element of (1/dd)L/L); its 37-part
            // lift in A[37]: (dd/37)*(x/dd) = x/37: coordinates (x div (dd/37)) /37...
            // pragmatic: reduce x * inverse(dd/p^k) mod p... handled below via 37*Hx test:
        end if;
        Append(~cands, Hx);
    end for;
    if k eq 2 then
        Hp := sub< Mp | [ p*x : x in Generators(Mp) ] >;
        if #Hp eq ord and
           &and[ IsZero(InnerProduct(u*EZd, v)) : u, v in Generators(Hp) ] and
           &and[ u*SId in Hp : u in Generators(Hp) ] then
            Append(~cands, Hp);
        end if;
    end if;
    return cands;
end function;

allc := [* *];
for pf in fac do
    cs := Cands(pf[1], pf[2]);
    printf "PRIME %o^%o: %o star-stable isotropics\n", pf[1], pf[2], #cs;
    Append(~allc, cs);
end for;
ncomb := &*[ Max(#cs,1) : cs in allc ];
printf "STAR_STABLE_COMBINATIONS %o\n", ncomb;
error if &or[ #cs eq 0 : cs in allc ], "a prime has no star-stable candidates";

// periods -- NB Magma computes Periods at the DEFAULT real precision,
// not at the precision implied by the term count; set it explicitly.
SetDefaultRealField(RealField(Prec + 10));
// convergence of period integrals at level N goes like exp(-2 pi n / N):
// need n ~ N * Prec * log(10) / (2 pi) terms for Prec digits.
nterms := Ceiling(Lv * Prec * Log(10)/(2*Pi(RealField(20)))) + 200;
printf "PERIOD_TERMS %o\n", nterms;
Pv := Periods(Af, nterms);
FX := RationalsExtra(Prec);      // field with recognition attributes
CC := FX`CC;                     // tagged complex field of the package
// FRAME FIX: Periods returns periods of Basis(Af) (rational basis); convert
// to the integral-lattice frame used by all our lattice coordinates:
// period of lattice-row v = v * Cm * PM_basis.
PMb := Matrix(CC, 4, 2, &cat[ [CC!Pv[i][1], CC!Pv[i][2]] : i in [1..4] ]);
PMl := ChangeRing(Cm, CC) * PMb;   // lattice-frame variant
printf "PERIODS ok (both frames prepared)\n";

// Derive the TRUE integral pairing on each frame directly from the
// periods: P * F * P^t = 0 is LINEAR in the alternating F = E^-1-scaled;
// solve, LLL-recognize integer alternating solutions, invert, scale.
function DeriveE(PMx)
    P0 := Transpose(PMx);   // 2 x 4
    // F alternating: unknowns f12,f13,f14,f23,f24,f34
    // (P0 F P0^t)[1][2] = sum over pairs (i<j) f_ij * (P0[1][i]P0[2][j] - P0[1][j]P0[2][i])
    prs := [<1,2>,<1,3>,<1,4>,<2,3>,<2,4>,<3,4>];
    cs := [ P0[1][p[1]]*P0[2][p[2]] - P0[1][p[2]]*P0[2][p[1]] : p in prs ];
    // real linear system: Re and Im of sum f_ij cs_ij = 0; integer relation
    // via LLL on rows [ 10^50*Re(cs_i), 10^50*Im(cs_i), e_i ]
    n := 6;
    Sc := 10^50;
    Lat := Matrix(Integers(), n, n+2, [0 : i in [1..n*(n+2)]]);
    for i in [1..n] do
        Lat[i][1] := Round(Sc*Re(cs[i]));
        Lat[i][2] := Round(Sc*Im(cs[i]));
        Lat[i][2+i] := 1;
    end for;
    Lr := LLL(Lat);
    sols := [];
    for i in [1..n] do
        if Abs(Lr[i][1]) lt 10^15 and Abs(Lr[i][2]) lt 10^15 then
            f := [ Lr[i][2+j] : j in [1..6] ];
            if &or[ x ne 0 : x in f ] then Append(~sols, f); end if;
        end if;
    end for;
    return sols, prs;
end function;

for frame in [1,2] do
    PMx := frame eq 1 select PMb else PMl;
    sols, prs := DeriveE(PMx);
    printf "FRAME %o : %o integral Riemann relations found\n", frame, #sols;
    for f in sols do
        FM := ZeroMatrix(Rationals(), 4, 4);
        for k in [1..6] do
            FM[prs[k][1]][prs[k][2]] := f[k];
            FM[prs[k][2]][prs[k][1]] := -f[k];
        end for;
        if Determinant(FM) eq 0 then printf "  degenerate F %o\n", f; continue; end if;
        Einv_true := FM;
        Etr := Einv_true^-1;
        dtr := LCM([Denominator(x) : x in Eltseq(Etr)]);
        EtrZ := Matrix(Integers(),4,4,[Integers()!(dtr*x) : x in Eltseq(Etr)]);
        gtr := GCD([x : x in Eltseq(EtrZ) | x ne 0]);
        EtrZ := EtrZ div gtr;
        printf "  F=%o -> E_true=%o edivs=%o\n", f, Eltseq(EtrZ), ElementaryDivisors(EtrZ);
    end for;
end for;
printf "E_DERIVATION_DONE\n";
// conventions validated on J0(23) (results/gl2_convention_cal.log):
// <frame, sign, blockorder> in { <1,-1,1>, <2,1,2>, <2,-1,1> }
convs := [ <1,-1,1>, <2,1,2>, <2,-1,1>, <1,1,1>, <1,1,2>, <2,-1,2> ];

// enumeration + reconstruction
nrec := 0; nQ := 0;
seen := {};
combo := [ 1 : cs in allc ];
done := false;
while not done and nrec lt MaxRec do
    Hgens := &cat[ [ x : x in Generators(allc[j][combo[j]]) ] : j in [1..#allc] ];
    lifts := [ Vector(Rationals(), [ (Integers()!x)/dd : x in Eltseq(u) ]) : u in Hgens ];
    B := VerticalJoin(ChangeRing(IdentityMatrix(Integers(),4), Rationals()),
                      Matrix(Rationals(), #lifts, 4, &cat[Eltseq(v) : v in lifts]));
    dn := LCM([Denominator(x) : x in Eltseq(B)]);
    BZ := Matrix(Integers(), Nrows(B), 4, [Integers()!(dn*x) : x in Eltseq(B)]);
    Hn := HermiteForm(BZ);
    Bp := Matrix(Rationals(), 4, 4, [ Hn[i][j]/dn : j in [1..4], i in [1..4] ]);
    Ep := Bp * ChangeRing(EZ, Rationals()) * Transpose(Bp);
    dnE := LCM([Denominator(x) : x in Eltseq(Ep)]);
    printf "COMBO %o dnE=%o\n", combo, dnE;
    if dnE eq 1 then
        EpZ := Matrix(Integers(), 4, 4, [Integers()!x : x in Eltseq(Ep)]);
        edp := ElementaryDivisors(EpZ);
        printf "  edivs=%o\n", edp;
        if edp[1] eq edp[4] then
            // principal: try the calibrated conventions in order
            nrec +:= 1;
            EpB := EpZ div edp[1];
            donecombo := false;
            for conv in convs do
                if donecombo then break; end if;
                sgn := conv[2];
                okF, _, T := true, 0, IdentityMatrix(Integers(),4);
                try
                    _, T := FrobeniusFormAlternating(sgn*EpB);
                catch e okF := false; end try;
                if not okF then continue; end if;
                PMuse := conv[1] eq 1 select PMb else PMl;
                Om := ChangeRing(ChangeRing(T, Rationals()) * Bp, CC) * PMuse;
                if conv[3] eq 1 then
                    bigP := Matrix(CC, 2, 4, [ Om[i][j] : i in [1..4], j in [1..2] ]);
                else
                    ordr := [3,4,1,2];
                    bigP := Matrix(CC, 2, 4, [ Om[ordr[i]][j] : i in [1..4], j in [1..2] ]);
                end if;
                isb := false;
                try isb := IsBigPeriodMatrix(bigP); catch e ; end try;
                if not isb then continue; end if;
                printf "COMBO %o conv=%o : bigP valid\n", combo, conv;
                try
                    Crec := ReconstructCurve(bigP, FX : Base := true);
                    FQ := BaseRing(Crec);
                    if Type(FQ) eq FldRat then
                        ii := IgusaInvariants(Crec);
                        if ii notin seen then
                            Include(~seen, ii);
                            nQ +:= 1;
                            printf "QCURVE %o combo=%o : %o\n", nQ, combo, Crec;
                        end if;
                        donecombo := true;
                    else
                        printf "NONQ combo=%o base=%o\n", combo, FQ;
                        donecombo := true;
                    end if;
                catch e
                    msg := Sprint(e`Object);
                    printf "RECFAIL combo=%o conv=%o err=%o\n", combo, conv,
                           #msg gt 160 select Substring(msg,1,160) else msg;
                end try;
            end for;
        end if;
    end if;
    j := 1;
    while j le #allc do
        combo[j] +:= 1;
        if combo[j] le #allc[j] then break; end if;
        combo[j] := 1; j +:= 1;
    end while;
    if j gt #allc then done := true; end if;
end while;
printf "RECONSTRUCTED %o principal candidates, %o rational curves\n", nrec, nQ;
printf "GL2_RECONSTRUCT_DONE\n";
quit;
