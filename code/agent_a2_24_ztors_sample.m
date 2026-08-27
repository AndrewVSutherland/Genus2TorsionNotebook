
//////////////////////////////////////////////////////////////////////
//  Harvest the rational-3-torsion locus on the A(8) chart.
//
//  We have order-8 FREE on A(8).  Scan (r,p,t), prefilter 3|#J(F_q) at
//  many primes (necessary for rational 3-torsion), exact TorsionSubgroup,
//  and LOG EVERY curve with 24 | an invariant (i.e. a rational order-24
//  point => Z/24 or Z/2xZ/24).  NO 2-rank filter: we want the whole
//  Z/24 family to (a) catch any that happen to be 2-rank 2 = [2,24], and
//  (b) collect (r,p,t) samples so the 3-torsion locus can be identified
//  (fit a curve/surface, then intersect with W-split for [2,24]).
//
//  Outputs, for each hit: r,p,t, full torsion, 2-rank, simplicity.
//
//  Usage: magma -b H:=18 NParts:=3 Part:=0 agent_a2_24_ztors_sample.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned H then H := 18; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned progress then progress := 100000;
elif Type(progress) eq MonStgElt then progress := StringToInteger(progress); end if;

PreP := [7,11,13,17,19,23,29,31,37,41];

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
    return f, g8, ellBase, L;
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

printf "ZTORS SAMPLE (harvest 3-torsion locus) H=%o Part=%o/%o\n", H, Part, NParts;
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
        for pv in vals do
            tested +:= 1;
            if tested mod progress eq 0 then
                printf "PROGRESS tested=%o pre3=%o tors=%o hits24=%o hits224=%o\n",
                    tested, pre3, tors, hits24, hits224;
            end if;
            f, g8, ellBase, Lform := A8f(rv, pv, tv);
            if Degree(f) lt 5 then continue; end if;
            if Discriminant(f) eq 0 then continue; end if;
            if Coefficient(f,1) eq 0 and Coefficient(f,3) eq 0 and Coefficient(f,5) eq 0 then continue; end if;
            fInt, Lden := IntModel(f);
            if not ThreeTorsionPrefilter(fInt) then continue; end if;
            pre3 +:= 1;
            J := Jacobian(HyperellipticCurve(fInt));
            O := J!0;
            v8 := (-Lden*ellBase) mod g8;   // scaled model: v8 = Lden*(unscaled v8)
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
            issimple, pp, chi := SimplicityCertificate(fInt);
            has2 := #[n : n in inv | n mod 2 eq 0] ge 2;
            printf "HIT24 r=%o p=%o t=%o torsion=%o 2rank=%o simple=%o\n",
                rv, pv, tv, inv, r2, issimple;
            printf "  f=%o\n", fInt;
            if has2 then
                hits224 +:= 1;
                printf "TARGET_2_24 r=%o p=%o t=%o torsion=%o 2rank=%o simple=%o%o\n",
                    rv, pv, tv, inv, r2, issimple,
                    issimple select Sprintf(" (q=%o chi=%o)", pp, chi) else "";
            end if;
        end for;
    end for;
end for;
printf "SEARCH_DONE tested=%o pre3=%o tors=%o hits24=%o hits224=%o\n",
    tested, pre3, tors, hits24, hits224;
print "TWO_RANK_HISTOGRAM (among Z/24 hits)";
for k in Sort([kk : kk in Keys(twoRankHist)]) do printf "  2rank %o : %o\n", k, twoRankHist[k]; end for;
print "DONE";
quit;
