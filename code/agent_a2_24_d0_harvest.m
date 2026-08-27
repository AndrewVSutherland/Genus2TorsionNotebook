
//////////////////////////////////////////////////////////////////////
//  Harvest simple Z/24 curves from the d=0 family E: y^2=x^3+x^2-4x
//  (rank 1).  Each rational point -> (r,t) on G -> genus-2 curve.
//  Verify torsion, simplicity, 2-rank.  Usage: magma -b agent_a2_24_d0_harvest.m
//////////////////////////////////////////////////////////////////////
SetColumns(0);
if not assigned MemGB then MemGB := 6; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals(); Z := Integers();
P<x> := PolynomialRing(Q);
R2<r,t> := PolynomialRing(Q, 2);

G := 3*r^3*t^3 + 2*r^3*t^2 - r^3*t - 2*r^2*t^4 + 8*r^2*t^3 - 22*r^2*t^2 + 4*r^2*t
     - r*t^5 + 6*r*t^4 - 13*r*t^3 + 12*r*t^2 + 8*r*t - 4;
PC := ProjectiveClosure(Curve(AffineSpace(R2), G));
ptB := PC ! [1/3, -1];
E, mp := EllipticCurve(PC, ptB);
gens := Generators(E);
printf "E = %o, rank gens = %o, torsion = %o\n", E, #gens, Invariants(TorsionSubgroup(E));

function A8f(rv, pv, tv)
    e := tv^2 - 2*pv*tv/rv; d := e + 2*pv - rv^2; lambda := rv/tv;
    a := rv^2 - lambda; b := 2*rv*pv - 2*lambda*(pv + rv*tv) + 2*rv*lambda;
    c := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
         - lambda*(rv^2 + e) + 2*lambda*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
    Qpoly := x^2 + d; q := a*x^2 + b*x + c; return q*(Qpoly^2 + q);
end function;
function IntModel(fv)
    L := 1; for i in [0..Degree(fv)] do L := LCM(L, Denominator(Coefficient(fv, i))); end for;
    return P!(L^2*fv); end function;
function TwoRank(fp)
    degs := [Degree(g[1]) : g in Factorization(fp)]; k := #degs; even := 0;
    for mask in [0..2^k-1] do ss := 0;
        for i in [1..k] do if (mask div 2^(i-1)) mod 2 eq 1 then ss +:= degs[i]; end if; end for;
        if ss mod 2 eq 0 then even +:= 1; end if; end for;
    return Ilog2(even) - 1; end function;
function CountCurve(fp)
    Fq := BaseRing(Parent(fp)); cnt := 0;
    for xx in Fq do vv := Evaluate(fp, xx);
        if vv eq 0 then cnt +:= 1; elif IsSquare(vv) then cnt +:= 2; end if; end for;
    if IsSquare(LeadingCoefficient(fp)) then cnt +:= 2; end if; return cnt; end function;
function IsSimple(fInt)
    RT := PolynomialRing(Q); T := RT.1; dsc := Discriminant(fInt);
    for pp in [3,5,7,11,13,17,19,23,29,31,37,41,43,47] do
        if (Z!LeadingCoefficient(fInt)) mod pp eq 0 then continue; end if;
        if (Z!Numerator(dsc)) mod pp eq 0 then continue; end if;
        PF := PolynomialRing(GF(pp)); fp := PF![GF(pp)!co : co in Coefficients(fInt)];
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        PF2 := PolynomialRing(GF(pp^2)); fp2 := PF2![GF(pp^2)!co : co in Coefficients(fInt)];
        a1 := pp + 1 - CountCurve(fp); a2 := (CountCurve(fp2) - pp^2 - 1 + a1^2) div 2;
        chi := T^4 - a1*T^3 + a2*T^2 - a1*pp*T + pp^2;
        if not IsIrreducible(chi) then continue; end if;
        K := NumberField(chi); pi := K.1; drop := false;
        for nn in [2..12] do if Degree(MinimalPolynomial(pi^nn)) lt 4 then drop := true; break; end if; end for;
        if not drop then return true, pp, chi; end if;
    end for; return false, 0, RT!0;
end function;

// build the list of E-points to pull back: multiples of the FREE
// generator, plus their translates by the 2-torsion point (all E-points)
E0 := E ! 0;
freegens := [g : g in gens | Order(g) eq 0];
torgens  := [g : g in gens | Order(g) ne 0];
P0 := freegens[1];
Tp := #torgens gt 0 select torgens[1] else E0;
pts := [];
for k in [-12..12] do
    Append(~pts, k*P0);
    Append(~pts, k*P0 + Tp);
end for;

printf "==== Z/24 curves from family points ====\n";
found := 0; simplecount := 0; r2hist := AssociativeArray(); seen := {};
for Ept in pts do
    if Ept eq E0 then continue; end if;
    clus := Ept @@ mp;                 // preimage on PC
    rp := RationalPoints(clus);
    for pcp in rp do
        cor := Coordinates(pcp);
        if cor[3] eq 0 then continue; end if;
        rv := cor[1]/cor[3]; tv := cor[2]/cor[3];
        if <rv,tv> in seen then continue; end if; Include(~seen, <rv,tv>);
        if rv eq 0 or tv eq 0 or rv*tv eq 1 then continue; end if;
        pv := rv*(rv+tv)/2;
        f := A8f(rv, pv, tv);
        if Degree(f) lt 5 or Discriminant(f) eq 0 then continue; end if;
        fInt := IntModel(f);
        J := Jacobian(HyperellipticCurve(fInt));
        inv := Invariants(TorsionSubgroup(J));
        has24 := #[nn : nn in inv | nn mod 24 eq 0] ge 1;
        if not has24 then
            printf "  (r,t)=(%o,%o) torsion=%o  [not 24]\n", rv, tv, inv; continue;
        end if;
        found +:= 1; r2 := TwoRank(fInt);
        if IsDefined(r2hist, r2) then r2hist[r2] +:= 1; else r2hist[r2] := 1; end if;
        iss, pp, chi := IsSimple(fInt);
        if iss then simplecount +:= 1; end if;
        printf "  (r,t)=(%o,%o) torsion=%o 2rank=%o simple=%o\n", rv, tv, inv, r2, iss;
    end for;
end for;
printf "found Z/24: %o (simple: %o)\n", found, simplecount;
print "2-rank histogram:";
for k in Sort([kk : kk in Keys(r2hist)]) do printf "  2rank %o : %o\n", k, r2hist[k]; end for;
print "DONE";
quit;
