
//////////////////////////////////////////////////////////////////////
//  d=0 slice of the A(8) chart: a fast 2-parameter (r,t) probe.
//
//  d = e + 2p - r^2 = 0  <=>  p = r(r+t)/2   (clean!).
//  Then Q = x^2, f = q*(x^4 + q), and x=0 is a FREE rational point
//  (f(0) = q(0)^2).  The known simple Z/24 curve B (r=1/3,p=-1/9,t=-1)
//  lives here.  Because p is determined, this is a 2-dim scan (~1000x
//  faster than the full 3-param chart), so we can go to large height and
//  harvest Z/24 samples quickly, and directly catch any 2-rank-2 = [2,24].
//
//  Fast #J mod 3 prefilter via direct point counting (no Jacobian object).
//
//  Usage: magma -b H:=40 NParts:=3 Part:=0 agent_a2_24_d0slice.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned H then H := 40; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned progress then progress := 200000;
elif Type(progress) eq MonStgElt then progress := StringToInteger(progress); end if;

PreP := [7,13,19,31,37,43];   // primes = 1 mod 3 (richer 3-torsion signal)

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

// direct #C(F_q) for y^2=fp (deg 6, two points at infinity)
function CountC(fp)
    Fq := BaseRing(Parent(fp)); cnt := 0;
    for xx in Fq do vv := Evaluate(fp, xx);
        if vv eq 0 then cnt +:= 1; elif IsSquare(vv) then cnt +:= 2; end if; end for;
    if IsSquare(LeadingCoefficient(fp)) then cnt +:= 2; end if;
    return cnt;
end function;

// #J(F_q) mod 3 == 0 ?  via char poly from #C(F_q), #C(F_q^2)
function Div3JacFF(fInt, q)
    PF := PolynomialRing(GF(q));
    fp := PF![GF(q)!co : co in Coefficients(fInt)];
    if Degree(fp) lt 5 or not IsSquarefree(fp) then return true; end if;  // no info => pass
    PF2 := PolynomialRing(GF(q^2));
    fp2 := PF2![GF(q^2)!co : co in Coefficients(fInt)];
    n1 := CountC(fp); n2 := CountC(fp2);
    a1 := q + 1 - n1;
    a2 := (n2 - q^2 - 1 + a1^2) div 2;
    nJ := 1 - a1 + a2 - a1*q + q^2;    // = P(1)
    return nJ mod 3 eq 0;
end function;

function ThreeTorsionPrefilter(fInt)
    for q in PreP do
        if not Div3JacFF(fInt, q) then return false; end if;
    end for;
    return true;
end function;

function CountCurve(fp)
    return CountC(fp);
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

function HeightRationals(HH)
    vals := [];
    for den in [1..HH] do for num in [-HH..HH] do
        if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
    end for; end for;
    return Sort(Setseq(Seqset(vals)));
end function;

printf "d=0 SLICE (p=r(r+t)/2) H=%o Part=%o/%o\n", H, Part, NParts;
vals := HeightRationals(H);
tested := 0; pre3 := 0; tors := 0; hits24 := 0; hits224 := 0;
twoRankHist := AssociativeArray();
ridx := 0;
for rv in vals do
    ridx +:= 1;
    if (ridx mod NParts) ne Part then continue; end if;
    if rv eq 0 then continue; end if;
    for tv in vals do
        if tv eq 0 or rv*tv eq 1 then continue; end if;
        pv := rv*(rv+tv)/2;      // d=0
        tested +:= 1;
        if tested mod progress eq 0 then
            printf "PROGRESS tested=%o pre3=%o tors=%o hits24=%o hits224=%o\n",
                tested, pre3, tors, hits24, hits224;
        end if;
        f, g8, ellBase := A8f(rv, pv, tv);
        if Degree(f) lt 5 then continue; end if;
        if Discriminant(f) eq 0 then continue; end if;
        if Coefficient(f,1) eq 0 and Coefficient(f,3) eq 0 and Coefficient(f,5) eq 0 then continue; end if;
        fInt, Lden := IntModel(f);
        if not ThreeTorsionPrefilter(fInt) then continue; end if;
        pre3 +:= 1;
        J := Jacobian(HyperellipticCurve(fInt));
        O := J!0;
        v8 := (-Lden*ellBase) mod g8;
        ok8 := true;
        try
            D8 := J![g8, v8];
            if 8*D8 ne O or 4*D8 eq O then ok8 := false; end if;
        catch ee ok8 := false; end try;
        if not ok8 then continue; end if;
        tors +:= 1;
        inv := Invariants(TorsionSubgroup(J));
        has24 := #[n : n in inv | n mod 24 eq 0] ge 1;
        if not has24 then continue; end if;
        hits24 +:= 1;
        r2 := TwoRank(fInt);
        if IsDefined(twoRankHist, r2) then twoRankHist[r2] +:= 1; else twoRankHist[r2] := 1; end if;
        has2 := #[n : n in inv | n mod 2 eq 0] ge 2;
        printf "HIT24 r=%o t=%o p=%o torsion=%o 2rank=%o\n", rv, tv, pv, inv, r2;
        printf "  f=%o\n", fInt;
        if has2 then
            hits224 +:= 1;
            issimple, pp, chi := SimplicityCertificate(fInt);
            printf "TARGET_2_24 r=%o t=%o p=%o torsion=%o 2rank=%o simple=%o%o\n",
                rv, tv, pv, inv, r2, issimple,
                issimple select Sprintf(" (q=%o chi=%o)", pp, chi) else "";
        end if;
    end for;
end for;
printf "SEARCH_DONE tested=%o pre3=%o tors=%o hits24=%o hits224=%o\n",
    tested, pre3, tors, hits24, hits224;
print "TWO_RANK_HISTOGRAM (among Z/24 hits)";
for k in Sort([kk : kk in Keys(twoRankHist)]) do printf "  2rank %o : %o\n", k, twoRankHist[k]; end for;
print "DONE";
quit;
