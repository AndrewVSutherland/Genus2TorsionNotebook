
//////////////////////////////////////////////////////////////////////
//  [2,24] backward construction, step 1: mod-p seeds.
//
//  Family f = (x^2-1)(x^2+a x+b)(x^2+c x+d): 2-rank 2 built in; if
//  D_infty has order 24 then torsion >= Z/2 x Z/24 automatically
//  (12*D_infty is one 2-torsion class; the rank-2 group has another).
//
//  The order-24 locus is codim 2 in (a,b,c,d): a SURFACE.  Plan:
//  find F_p points (this script), Newton-lift p-adically along the CF
//  closure conditions, rationally reconstruct.
//
//  This script scans (a,b,c,d) in F_p^4, computes the exact order of
//  D_infty via the continued fraction over F_p, and records order-24
//  seeds together with their CF degree pattern (generic cell
//  [3,1,1,...,1] preferred for the 2-condition Newton lift).
//
//  Usage: magma -b p:=11 agent_a2_24_construct_seeds.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 2; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned p then p := 11; elif Type(p) eq MonStgElt then p := StringToInteger(p); end if;
if not assigned MaxSeeds then MaxSeeds := 40; elif Type(MaxSeeds) eq MonStgElt then MaxSeeds := StringToInteger(MaxSeeds); end if;

Fp := GF(p);
P<x> := PolynomialRing(Fp);

function SqrtPolyPart(f)
    s := x^3;
    for k in [1..3] do
        d := f - s^2;
        if Degree(d) le 2 then break; end if;
        s := s + (Coefficient(d,6-k)/(2*Coefficient(s,3)))*x^(3-k);
    end for;
    return s;
end function;

// CF order over F_p with degree pattern
function CFOrderPattern(f, ms)
    s := SqrtPolyPart(f);
    Pi := P!0; Qi := P!1; tot := 0; degs := [];
    for i in [0..ms] do
        if Qi eq 0 then return 0, degs; end if;
        if not IsUnit(LeadingCoefficient(Qi)) then return 0, degs; end if;
        ai := (Pi + s) div Qi;
        tot +:= Degree(ai);
        Append(~degs, Degree(ai));
        Pn := ai*Qi - Pi;
        if (f - Pn^2) mod Qi ne 0 then return 0, degs; end if;
        Qi := (f - Pn^2) div Qi;
        Pi := Pn;
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then return tot, degs; end if;
    end for;
    return 0, degs;
end function;

printf "SEED SCAN p=%o family=(x^2-1)(x^2+ax+b)(x^2+cx+d)\n", p;
nsf := 0; found := 0;
ordHist := AssociativeArray();
for a in Fp do for b in Fp do for c in Fp do for d in Fp do
    f := (x^2-1)*(x^2+a*x+b)*(x^2+c*x+d);
    if not IsSquarefree(f) then continue; end if;
    nsf +:= 1;
    ord, degs := CFOrderPattern(f, 40);
    if ord eq 0 then continue; end if;
    if IsDefined(ordHist, ord) then ordHist[ord] +:= 1; else ordHist[ord] := 1; end if;
    if ord eq 24 then
        found +:= 1;
        generic := (#degs eq 22) and degs[1] eq 3 and &and[degs[i] eq 1 : i in [2..22]];
        if found le MaxSeeds then
            printf "SEED24 a,b,c,d=%o,%o,%o,%o pattern=%o generic=%o\n",
                a, b, c, d, degs, generic;
        end if;
    end if;
end for; end for; end for; end for;
printf "squarefree=%o order-24 seeds=%o\n", nsf, found;
print "order histogram (finite orders over F_p):";
for k in Sort([kk : kk in Keys(ordHist)]) do
    printf "  %o: %o\n", k, ordHist[k];
end for;
print "DONE";
quit;
