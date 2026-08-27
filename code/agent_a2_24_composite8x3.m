
//////////////////////////////////////////////////////////////////////
//  [24] and [2,24] via coprime composition 24 = 8 x 3 on the A(8) chart.
//
//  KEY IDEA (avoids both the halving wall and the Pell wall): a curve
//  with a rational order-8 class (visible D8 on the A(8) chart) AND any
//  rational 3-torsion class has Z/8 + Z/3 = Z/24 in its torsion
//  automatically -- coprime orders always compose.
//
//  Funnel: sweep the A(8) chart (r,p,t); D8 gives 8 | #J(F_q) for free,
//  so prefilter 3 | #J(F_q) at several good primes (necessary for
//  rational 3-torsion); survivors get exact TorsionSubgroup; flag
//  [24]-containing torsion, and [2,24] if 2-rank >= 2.
//
//  Usage:
//    magma -b H:=6 NParts:=3 Part:=0 agent_a2_24_composite8x3.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned H then H := 6; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned progress then progress := 100000;
elif Type(progress) eq MonStgElt then progress := StringToInteger(progress); end if;

PreP := [7,11,13,17,19,23,29,31,37,41];   // strong 3-divisibility prefilter

// A(8) chart (m=1 gauge, from notes.tex / search_A8_* scripts)
function A8f(rv, pv, tv)
    e := tv^2 - 2*pv*tv/rv; d := e + 2*pv - rv^2; lambda := rv/tv;
    u := pv + rv*tv - 2*rv;
    v := e + rv^2 - rv*pv - rv^2*tv + 3*pv*tv - rv*tv^2;
    a := rv^2 - lambda;
    b := 2*rv*pv - 2*lambda*(pv + rv*tv) + 2*rv*lambda;
    c := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
         - lambda*(rv^2 + e) + 2*lambda*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
    Qpoly := x^2 + d; q := a*x^2 + b*x + c;
    g8 := x^2 + u*x + v;
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

// prefilter: 3 | #J(F_q) at all good primes in PreP (necessary for
// rational 3-torsion); skip bad-reduction primes (no info)
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

printf "COMPOSITE 8x3 SEARCH on A(8): H=%o Part=%o/%o prefilter 3|#J at %o\n",
    H, Part, NParts, PreP;
vals := HeightRationals(H);
tested := 0; smooth := 0; rank2 := 0; pre3 := 0; tors := 0; hits24 := 0; hits224 := 0;
torsHist := AssociativeArray();
ridx := 0;
for rv in vals do
    ridx +:= 1;
    if (ridx mod NParts) ne Part then continue; end if;
    if rv eq 0 then continue; end if;
    for tv in vals do
        if tv eq 0 or rv*tv eq 1 then continue; end if;   // lambda pole / a=0 handled below
        for pv in vals do
            tested +:= 1;
            if tested mod progress eq 0 then
                printf "PROGRESS tested=%o smooth=%o rank2=%o pre3=%o tors=%o hits224=%o\n",
                    tested, smooth, rank2, pre3, tors, hits224;
            end if;
            f, g8, ellBase := A8f(rv, pv, tv);
            if Degree(f) lt 5 then continue; end if;
            if Discriminant(f) eq 0 then continue; end if;
            // skip even sextics (bielliptic => split, never simple)
            if &and[Coefficient(f,i) eq 0 : i in [1,3,5]] then continue; end if;
            smooth +:= 1;
            fInt, L := IntModel(f);
            // need 2-rank >= 2 (the extra [2] factor) for [2,24]; cheap factor filter
            if TwoRank(fInt) lt 2 then continue; end if;
            rank2 +:= 1;
            if not ThreeTorsionPrefilter(fInt) then continue; end if;
            pre3 +:= 1;
            // verify D8 really has order 8 here, then exact torsion
            J := Jacobian(HyperellipticCurve(fInt));
            O := J!0;
            v8 := (-L*ellBase) mod g8;
            ok8 := true;
            try
                D8 := J![g8, v8];
                if 8*D8 ne O or 4*D8 eq O then ok8 := false; end if;
            catch e ok8 := false; end try;
            if not ok8 then continue; end if;
            tors +:= 1;
            inv := Invariants(TorsionSubgroup(J));
            key := Sprint(inv);
            if IsDefined(torsHist,key) then torsHist[key] +:= 1; else torsHist[key] := 1; end if;
            has24 := #[n : n in inv | n mod 24 eq 0] ge 1;
            if has24 then
                hits24 +:= 1;
                r2 := TwoRank(f);
                issimple, pp, chi := SimplicityCertificate(fInt);
                printf "HIT24 r=%o p=%o t=%o torsion=%o 2rank=%o simple=%o%o\n",
                    rv, pv, tv, inv, r2, issimple,
                    issimple select Sprintf(" (q=%o chi=%o)", pp, chi) else "";
                printf "  f = %o\n", fInt;
                has2 := #[n : n in inv | n mod 2 eq 0] ge 2;
                if has2 then
                    hits224 +:= 1;
                    printf "TARGET_2_24 r=%o p=%o t=%o torsion=%o simple=%o\n",
                        rv, pv, tv, inv, issimple;
                end if;
            end if;
        end for;
    end for;
end for;
printf "SEARCH_DONE tested=%o smooth=%o rank2=%o pre3=%o tors=%o hits24=%o hits224=%o\n",
    tested, smooth, rank2, pre3, tors, hits24, hits224;
print "TORSION_HISTOGRAM";
for k in Sort([kk : kk in Keys(torsHist)]) do printf "  %o : %o\n", k, torsHist[k]; end for;
print "DONE";
quit;
