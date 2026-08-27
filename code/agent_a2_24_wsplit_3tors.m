
//////////////////////////////////////////////////////////////////////
//  [2,24] DENSE search: W-split (2-rank 2 free) + 3-torsion on A(8).
//
//  Composite route 24 = 8 x 3: A(8) gives order-8 D8 free.  Instead of
//  a blind (r,p,t) scan (2-rank 2 is ~3e-5 rare), parametrize the
//  2-rank-2 sublocus directly: W = Q^2+q splits into two rational
//  quadratics.  The split condition
//    E(p) = beta^2*(2d+a+beta^2)^2 - 4*(d^2+c)*beta^2 - b^2 = 0
//  is QUADRATIC in p for fixed (r,t,beta) (from the [4,16] work), so
//  every (r,t,beta) with a rational root gives a 2-rank-2 A(8) curve.
//  Then only the (rare) 3-torsion remains: strong 3|#J(F_q) prefilter,
//  then exact torsion.  order-8 + 2-rank-2 + 3-torsion = [2,24].
//
//  Usage: magma -b RH:=20 TH:=20 BH:=20 NParts:=3 Part:=0 agent_a2_24_wsplit_3tors.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned RH then RH := 16; elif Type(RH) eq MonStgElt then RH := StringToInteger(RH); end if;
if not assigned TH then TH := 16; elif Type(TH) eq MonStgElt then TH := StringToInteger(TH); end if;
if not assigned BH then BH := 16; elif Type(BH) eq MonStgElt then BH := StringToInteger(BH); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned progress then progress := 100000;
elif Type(progress) eq MonStgElt then progress := StringToInteger(progress); end if;

PreP := [7,11,13,17,19,23,29,31,37,41,43,47,53,59];

function A8f(rv, pv, tv)
    e := tv^2 - 2*pv*tv/rv; d := e + 2*pv - rv^2; lambda := rv/tv;
    u := pv + rv*tv - 2*rv;
    v := e + rv^2 - rv*pv - rv^2*tv + 3*pv*tv - rv*tv^2;
    a := rv^2 - lambda;
    b := 2*rv*pv - 2*lambda*(pv + rv*tv) + 2*rv*lambda;
    c := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
         - lambda*(rv^2 + e) + 2*lambda*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
    Qpoly := x^2 + d; q := a*x^2 + b*x + c; g8 := x^2 + u*x + v;
    f := q*(Qpoly^2 + q);
    L := rv*x + (pv - rv^2);
    ellBase := -(q + Qpoly*L);
    return f, g8, ellBase;
end function;

// W-split condition as a quadratic in p; return rational p-roots
function SplitWPvals(rv, tv, bv)
    Pp<pp> := PolynomialRing(Q);
    e := tv^2 - 2*pp*tv/rv; d := e + 2*pp - rv^2; lambda := rv/tv;
    a := rv^2 - lambda;
    b := 2*rv*pp - 2*lambda*(pp + rv*tv) + 2*rv*lambda;
    c := pp^2 + 2*pp*rv^2 - rv^4 - rv^3*tv - rv*pp^2/tv
         - lambda*(rv^2 + e) + 2*lambda*(rv*pp + rv^2*tv - 3*pp*tv + rv*tv^2);
    E := bv^2*(2*d + a + bv^2)^2 - 4*(d^2 + c)*bv^2 - b^2;
    En := Numerator(E);
    if En eq 0 or Degree(En) lt 1 then return []; end if;
    return [rt[1] : rt in Roots(En)];
end function;

function IntModel(fv)
    L := 1;
    for i in [0..Degree(fv)] do L := LCM(L, Denominator(Coefficient(fv, i))); end for;
    return P!(L^2*fv), L;
end function;

function TwoRank(fp)
    degs := [Degree(g[1]) : g in Factorization(fp)];
    k := #degs; even := 0;
    for mask in [0..2^k-1] do
        ss := 0; for i in [1..k] do if (mask div 2^(i-1)) mod 2 eq 1 then ss +:= degs[i]; end if; end for;
        if ss mod 2 eq 0 then even +:= 1; end if;
    end for;
    return Ilog2(even) - 1;
end function;

function CountCurve(fp)
    Fq := BaseRing(Parent(fp)); cnt := 0;
    for xx in Fq do vv := Evaluate(fp, xx);
        if vv eq 0 then cnt +:= 1; elif IsSquare(vv) then cnt +:= 2; end if; end for;
    if IsSquare(LeadingCoefficient(fp)) then cnt +:= 2; end if;
    return cnt;
end function;

function SimplicityCertificate(fInt)
    RT := PolynomialRing(Q); T := RT.1;
    dsc := Discriminant(fInt);
    for pp in [3,5,7,11,13,17,19,23,29,31,37,41,43,47] do
        if (Z!LeadingCoefficient(fInt)) mod pp eq 0 then continue; end if;
        if (Z!Numerator(dsc)) mod pp eq 0 then continue; end if;
        PF := PolynomialRing(GF(pp));
        fp := PF![GF(pp)!co : co in Coefficients(fInt)];
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        PF2 := PolynomialRing(GF(pp^2));
        fp2 := PF2![GF(pp^2)!co : co in Coefficients(fInt)];
        a1 := pp + 1 - CountCurve(fp);
        a2 := (CountCurve(fp2) - pp^2 - 1 + a1^2) div 2;
        chi := T^4 - a1*T^3 + a2*T^2 - a1*pp*T + pp^2;
        if not IsIrreducible(chi) then continue; end if;
        K := NumberField(chi); pi := K.1; drop := false;
        for nn in [2..12] do
            if Degree(MinimalPolynomial(pi^nn)) lt 4 then drop := true; break; end if;
        end for;
        if not drop then return true, pp, chi; end if;
    end for;
    return false, 0, RT!0;
end function;

function ThreeTorsionPrefilter(fInt)
    for q in PreP do
        PF := PolynomialRing(GF(q));
        fp := PF![GF(q)!co : co in Coefficients(fInt)];
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        if #Jacobian(HyperellipticCurve(fp)) mod 3 ne 0 then return false; end if;
    end for;
    return true;
end function;

function HeightRationals(HH)
    vals := [];
    for den in [1..HH] do for num in [-HH..HH] do
        if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
    end for; end for;
    return Sort(Setseq(Seqset(vals)));
end function;

printf "WSPLIT+3TORS [2,24] search RH=%o TH=%o BH=%o Part=%o/%o\n", RH, TH, BH, Part, NParts;
rvals := HeightRationals(RH);
tvals := HeightRationals(TH);
bvals := HeightRationals(BH);
tested := 0; pcnt := 0; rank2 := 0; pre3 := 0; tors := 0; hits24 := 0; hits224 := 0;
torsHist := AssociativeArray();
ridx := 0;
for rv in rvals do
    ridx +:= 1;
    if (ridx mod NParts) ne Part then continue; end if;
    if rv eq 0 or rv eq 1 then continue; end if;
    for tv in tvals do
        if tv eq 0 or rv*tv eq 1 then continue; end if;
        for bv in bvals do
            if bv eq 0 then continue; end if;
            tested +:= 1;
            if tested mod progress eq 0 then
                printf "PROGRESS tested=%o p=%o rank2=%o pre3=%o tors=%o hits24=%o hits224=%o\n",
                    tested, pcnt, rank2, pre3, tors, hits24, hits224;
            end if;
            pvs := [];
            try pvs := SplitWPvals(rv, tv, bv); catch ee continue; end try;
            for pv in pvs do
                pcnt +:= 1;
                f, g8, ellBase := A8f(rv, pv, tv);
                if Degree(f) lt 5 then continue; end if;
                if Discriminant(f) eq 0 then continue; end if;
                // skip even (split)
                if &and[Coefficient(f,i) eq 0 : i in [1,3,5]] then continue; end if;
                fInt, L := IntModel(f);
                if TwoRank(fInt) lt 2 then continue; end if;   // confirm 2-rank 2
                rank2 +:= 1;
                if not ThreeTorsionPrefilter(fInt) then continue; end if;
                pre3 +:= 1;
                J := Jacobian(HyperellipticCurve(fInt));
                O := J!0;
                v8 := (-L*ellBase) mod g8;
                ok8 := true;
                try
                    D8 := J![g8, v8];
                    if 8*D8 ne O or 4*D8 eq O then ok8 := false; end if;
                catch ee ok8 := false; end try;
                if not ok8 then continue; end if;
                tors +:= 1;
                inv := Invariants(TorsionSubgroup(J));
                key := Sprint(inv);
                if IsDefined(torsHist,key) then torsHist[key] +:= 1; else torsHist[key] := 1; end if;
                has24 := #[n : n in inv | n mod 24 eq 0] ge 1;
                has2 := #[n : n in inv | n mod 2 eq 0] ge 2;
                if has24 then
                    hits24 +:= 1;
                    issimple, pp, chi := SimplicityCertificate(fInt);
                    printf "HIT24 r=%o t=%o b=%o p=%o torsion=%o simple=%o%o\n",
                        rv, tv, bv, pv, inv, issimple,
                        issimple select Sprintf(" (q=%o chi=%o)", pp, chi) else "";
                    printf "  f = %o\n", fInt;
                    if has2 then
                        hits224 +:= 1;
                        printf "TARGET_2_24 r=%o t=%o b=%o p=%o torsion=%o simple=%o\n",
                            rv, tv, bv, pv, inv, issimple;
                    end if;
                end if;
            end for;
        end for;
    end for;
end for;
printf "SEARCH_DONE tested=%o p=%o rank2=%o pre3=%o tors=%o hits24=%o hits224=%o\n",
    tested, pcnt, rank2, pre3, tors, hits24, hits224;
print "TORSION_HISTOGRAM";
for k in Sort([kk : kk in Keys(torsHist)]) do printf "  %o : %o\n", k, torsHist[k]; end for;
print "DONE";
quit;
