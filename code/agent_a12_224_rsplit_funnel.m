
//////////////////////////////////////////////////////////////////////
//  Z/2 x Z/24 focused rank-2 search on the A(12) chart.
//
//  Land directly on the 2-rank-2 locus (necessary for [2,24]):
//   * R splits for free:  disc(R) = 4p^2 + 4pr + 1 = w^2, so
//        r = (w^2 - 4p^2 - 1)/(4p),  giving f = R*F with R -> [1,1];
//   * require F (quartic) to have a rational root  =>  f -> [1,1,1,3]
//        (or richer), 2-rank >= 2.
//  On such curves P12 = P4 + P6 has order 12; a rational half of some
//  2-torsion translate P12 + T gives an order-24 point, and with the
//  extra rational 2-torsion the torsion contains Z/2 x Z/24.
//
//  Per point:  R-split param -> build f -> require Roots(F) nonempty
//  (cheap) -> mod-p 24-divisibility prefilter -> exact IsDivisibleBy
//  over the (small) rational 2-torsion translate group -> torsion +
//  simplicity, flag TARGET_2_24.
//
//  Usage:
//    magma -b Validate:=true agent_a12_224_rsplit_funnel.m
//    magma -b pH:=40 zH:=40 wH:=40 NParts:=3 Part:=0 agent_a12_224_rsplit_funnel.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned pH then pH := 20; elif Type(pH) eq MonStgElt then pH := StringToInteger(pH); end if;
if not assigned zH then zH := 20; elif Type(zH) eq MonStgElt then zH := StringToInteger(zH); end if;
if not assigned wH then wH := 20; elif Type(wH) eq MonStgElt then wH := StringToInteger(wH); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned Validate then Validate := false;
elif Type(Validate) eq MonStgElt then Validate := Validate in {"true","True","1","yes"}; end if;
if not assigned progress then progress := 50000;
elif Type(progress) eq MonStgElt then progress := StringToInteger(progress); end if;
PreP := 13;   // cheap 24-divisibility prefilter prime

function ChartData(pv, zv, rv)
    if pv eq 0 or zv eq 0 then return false, _, _, _, _, _; end if;
    sv := (zv^2 - 4*pv^2 + 1)/(2*zv);
    if sv^2 eq 1 then return false, _, _, _, _, _; end if;
    tv := (zv^2 + 4*pv^2 - 1)^2/(8*pv^2*zv);
    muv := ((sv^2 - 1)*(2*pv*rv + 1) - pv^2*(2*sv*tv - 4))/(4*pv^3);
    lav := (4 - muv^2)*pv^2/(sv^2 - 1);
    if lav eq 0 then return false, _, _, _, _, _; end if;
    T1v := pv*x + rv;
    Rv := (T1v^2 + x - 1)/lav;
    ellv := sv*x + tv;
    Qv := 2*T1v + muv*Rv;
    Fv := Rv*x^2 + 4*(Rv + x - 1)*(Rv - 1);
    fv := Rv*Fv;
    if Degree(fv) ne 6 or Discriminant(fv) eq 0 then return false, _, _, _, _, _; end if;
    return true, fv, Rv, Qv, ellv, Fv;
end function;

function IntModel(fv)
    L := 1;
    for i in [0..Degree(fv)] do L := LCM(L, Denominator(Coefficient(fv, i))); end for;
    return P!(L^2*fv), L;
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

function Has224(inv)
    return #[n : n in inv | n mod 24 eq 0] ge 1 and #[n : n in inv | n mod 2 eq 0] ge 2;
end function;

// cheap prefilter: 24 | #J(F_PreP)
function PrefilterOK(fInt)
    PF := PolynomialRing(GF(PreP));
    fp := PF![GF(PreP)!co : co in Coefficients(fInt)];
    if Degree(fp) lt 5 or not IsSquarefree(fp) then return true; end if;   // bad red: don't reject
    return #Jacobian(HyperellipticCurve(fp)) mod 24 eq 0;
end function;

// process a rank-2 (F-root) point: exact halving over 2-torsion translates
function ProcessRank2(fv, Rv, Qv, ellv)
    fInt, L := IntModel(fv);
    if not PrefilterOK(fInt) then return "prefilter", _; end if;
    J := Jacobian(HyperellipticCurve(fInt));
    O := J!0;
    u4 := Qv/LeadingCoefficient(Qv);
    u6 := (Rv + x - 1)/LeadingCoefficient(Rv + x - 1);
    try
        P4 := J![u4, (L*Rv*ellv) mod u4];
        P6 := J![u6, (L*x*Rv) mod u6];
    catch e return "badpt", _; end try;
    P12 := P4 + P6;
    if not (12*P12 eq O and &and[nn*P12 ne O : nn in [1..11]]) then return "not12", _; end if;
    // rational 2-torsion translate group
    T2, mp2 := TwoTorsionSubgroup(J);
    for tg in T2 do
        D := P12 + mp2(tg);
        if not (12*D eq O) then continue; end if;
        div2, half := IsDivisibleBy(D, 2);
        if div2 then
            inv := Invariants(TorsionSubgroup(J));
            return "HALVABLE", <inv, Order(half)>;
        end if;
    end for;
    return "rank2-no-half", _;
end function;

function HeightRationals(H)
    vals := [];
    for den in [1..H] do for num in [-H..H] do
        if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
    end for; end for;
    return Sort(Setseq(Seqset(vals)));
end function;

if Validate then
    print "== VALIDATION: R-split + F-root => 2-rank>=2, and disc(R)=w^2 ==";
    nchart := 0; nrsplit := 0; nfroot := 0; nrank2 := 0; nhalv := 0;
    for pv in HeightRationals(4) do
        if pv eq 0 then continue; end if;
        for zv in HeightRationals(3) do
            if zv eq 0 then continue; end if;
            for wv in HeightRationals(4) do
                rv := (wv^2 - 4*pv^2 - 1)/(4*pv);
                ok, fv, Rv, Qv, ellv, Fv := ChartData(pv, zv, rv);
                if not ok then continue; end if;
                nchart +:= 1;
                // disc(R) sanity: should be w^2 up to lambda scaling
                if #Roots(Rv) eq 2 then nrsplit +:= 1; else continue; end if;
                if #Roots(Fv) ge 1 then nfroot +:= 1; else continue; end if;
                // confirm 2-rank >= 2
                degs := [Degree(g[1]) : g in Factorization(fv)];
                k := #degs; even := 0;
                for mask in [0..2^k-1] do
                    ss := 0; for i in [1..k] do if (mask div 2^(i-1)) mod 2 eq 1 then ss +:= degs[i]; end if; end for;
                    if ss mod 2 eq 0 then even +:= 1; end if;
                end for;
                r2 := Ilog2(even) - 1;
                if r2 ge 2 then nrank2 +:= 1; end if;
                status, data := ProcessRank2(fv, Rv, Qv, ellv);
                if status eq "HALVABLE" then
                    nhalv +:= 1;
                    printf "  VALIDATION HALVABLE p=%o z=%o w=%o torsion=%o\n", pv, zv, wv, data[1];
                end if;
            end for;
        end for;
    end for;
    printf "validation: chart=%o R-split=%o F-root=%o rank>=2=%o halvable=%o\n",
        nchart, nrsplit, nfroot, nrank2, nhalv;
    printf "(R-split points that also have F-root should be 2-rank>=2: F-root=%o vs rank2=%o)\n",
        nfroot, nrank2;
    print "DONE";
    quit;
end if;

// ---- search ----
printf "RSPLIT SEARCH pH=%o zH=%o wH=%o Part=%o/%o (r=(w^2-4p^2-1)/(4p))\n",
    pH, zH, wH, Part, NParts;
pvals := HeightRationals(pH);
zvals := HeightRationals(zH);
wvals := HeightRationals(wH);
tested := 0; rsplit := 0; froot := 0; prefiltered := 0; halvable := 0;
pidx := 0;
for pv in pvals do
    pidx +:= 1;
    if (pidx mod NParts) ne Part then continue; end if;
    if pv eq 0 then continue; end if;
    for wv in wvals do
        rv := (wv^2 - 4*pv^2 - 1)/(4*pv);
        for zv in zvals do
            if zv eq 0 then continue; end if;
            tested +:= 1;
            if tested mod progress eq 0 then
                printf "PROGRESS tested=%o rsplit=%o froot=%o halvable=%o\n",
                    tested, rsplit, froot, halvable;
            end if;
            ok, fv, Rv, Qv, ellv, Fv := ChartData(pv, zv, rv);
            if not ok then continue; end if;
            rsplit +:= 1;
            if #Roots(Fv) eq 0 then continue; end if;   // need F-root for rank 2
            froot +:= 1;
            status, data := ProcessRank2(fv, Rv, Qv, ellv);
            if status eq "HALVABLE" then
                halvable +:= 1;
                inv := data[1];
                printf "HALVABLE p=%o z=%o w=%o r=%o torsion=%o half_order=%o\n",
                    pv, zv, wv, rv, inv, data[2];
                fInt := IntModel(fv);
                printf "  f = %o\n", fInt;
                if Has224(inv) then
                    issimple, pp, chi := SimplicityCertificate(fInt);
                    printf "TARGET_2_24 p=%o z=%o w=%o torsion=%o simple=%o%o\n",
                        pv, zv, wv, inv,
                        issimple, issimple select Sprintf(" (p=%o chi=%o)", pp, chi) else "";
                end if;
            end if;
        end for;
    end for;
end for;
printf "SEARCH_DONE tested=%o rsplit=%o froot=%o halvable=%o\n",
    tested, rsplit, froot, halvable;
print "DONE";
quit;
