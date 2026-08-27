
//////////////////////////////////////////////////////////////////////
//  Z/2 x Z/24 via exact CF order on a 2-rank-2 family.
//
//  Family:  f = (x^2 - 1)*(x^2 + a*x + b)*(x^2 + c*x + d)
//  monic deg 6, factor type [1,1,2,2] => 2-rank 2 (independent rational
//  2-torsion built in).  Compute the EXACT order of D_infty via the
//  polynomial continued fraction (agent_a2_24_cf.m); order 24 + 2-rank 2
//  gives torsion containing Z/2 x Z/24 (once D_infty's 2-part is
//  independent of the extra 2-torsion, checked by exact torsion).
//
//  The CF order is a cheap EXACT filter, so we hunt order 24 directly.
//  Survivors get simplicity + exact TorsionSubgroup.
//
//  Usage:
//    magma -b Validate:=true agent_a2_24_cf_search.m
//    magma -b H:=15 NParts:=3 Part:=0 agent_a2_24_cf_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned H then H := 12; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned Validate then Validate := false;
elif Type(Validate) eq MonStgElt then Validate := Validate in {"true","True","1","yes"}; end if;
if not assigned progress then progress := 200000;
elif Type(progress) eq MonStgElt then progress := StringToInteger(progress); end if;
if not assigned TrackFrom then TrackFrom := 16; elif Type(TrackFrom) eq MonStgElt then TrackFrom := StringToInteger(TrackFrom); end if;

function SqrtPolyPart(f)
    s := x^3;
    for k in [1..3] do
        d := f - s^2;
        if Degree(d) le 2 then break; end if;
        s := s + (Coefficient(d, 6-k)/(2*Coefficient(s,3)))*x^(3-k);
    end for;
    return s;
end function;

// exact D_infty order via CF; returns 0 if not periodic within maxsteps
function CFOrder(f, maxsteps)
    s := SqrtPolyPart(f);
    Pi := P!0; Qi := P!1; total := 0;
    for i in [0..maxsteps] do
        if Qi eq 0 then return 0; end if;
        ai := (Pi + s) div Qi;
        total +:= Degree(ai);
        Pn := ai*Qi - Pi;
        if (f - Pn^2) mod Qi ne 0 then return 0; end if;
        Qn := (f - Pn^2) div Qi;
        Pi := Pn; Qi := Qn;
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then return total; end if;
    end for;
    return 0;
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

if Validate then
    // f14 = (x^2+1)(x^4+5x^2+4x+4) is NOT in our family, but test CFOrder here:
    f14 := (x^2+1)*(x^4+5*x^2+4*x+4);
    printf "CFOrder(f14) = %o (expect 14)\n", CFOrder(f14, 60);
    // a family member: pick random a,b,c,d and print its order + torsion
    for tr in [[1,1,-1,2],[2,-3,1,1],[0,2,1,-1]] do
        f := (x^2-1)*(x^2+tr[1]*x+tr[2])*(x^2+tr[3]*x+tr[4]);
        if not IsSquarefree(f) then continue; end if;
        ord := CFOrder(f, 60);
        inv := Invariants(TorsionSubgroup(Jacobian(HyperellipticCurve(f))));
        printf "  a,b,c,d=%o: CFOrder=%o torsion=%o 2rank=%o\n", tr, ord, inv, TwoRank(f);
    end for;
    print "DONE"; quit;
end if;

printf "CF-ORDER RANK2 SEARCH H=%o Part=%o/%o family=(x^2-1)(x^2+ax+b)(x^2+cx+d)\n", H, Part, NParts;
tested := 0; sqfree := 0; ord24 := 0; hits := 0;
maxord := 0; maxinfo := "";
aidx := 0;
for a in [-H..H] do
  aidx +:= 1;
  if (aidx mod NParts) ne Part then continue; end if;
  for b in [-H..H] do for c in [-H..H] do for d in [-H..H] do
    tested +:= 1;
    if tested mod progress eq 0 then
      printf "PROGRESS tested=%o sqfree=%o ord24=%o hits=%o maxord=%o\n",
          tested, sqfree, ord24, hits, maxord;
    end if;
    f := (x^2-1)*(x^2+a*x+b)*(x^2+c*x+d);
    if Discriminant(f) eq 0 then continue; end if;
    sqfree +:= 1;
    ord := CFOrder(f, 40);
    if ord eq 0 then continue; end if;
    if ord gt maxord then maxord := ord; maxinfo := Sprintf("a,b,c,d=%o,%o,%o,%o", a,b,c,d); end if;
    if ord ge TrackFrom then
        printf "BIGORDER a,b,c,d=%o,%o,%o,%o CFOrder=%o\n", a, b, c, d, ord;
    end if;
    if ord mod 24 ne 0 then continue; end if;
    ord24 +:= 1;
    // confirm torsion + 2-rank + simplicity
    J := Jacobian(HyperellipticCurve(f));
    inv := Invariants(TorsionSubgroup(J));
    has24 := #[n : n in inv | n mod 24 eq 0] ge 1;
    has2 := #[n : n in inv | n mod 2 eq 0] ge 2;
    printf "ORDER24 a,b,c,d=%o,%o,%o,%o CFOrder=%o torsion=%o 2rank=%o\n",
        a, b, c, d, ord, inv, TwoRank(f);
    if has24 and has2 then
        hits +:= 1;
        issimple, pp, chi := SimplicityCertificate(f);
        printf "TARGET_2_24 a,b,c,d=%o,%o,%o,%o torsion=%o simple=%o%o\n  f=%o\n",
            a, b, c, d, inv, issimple, issimple select Sprintf(" (p=%o chi=%o)",pp,chi) else "", f;
    end if;
  end for; end for; end for;
end for;
printf "SEARCH_DONE tested=%o sqfree=%o ord24=%o hits=%o maxord=%o (%o)\n",
    tested, sqfree, ord24, hits, maxord, maxinfo;
print "DONE";
quit;
