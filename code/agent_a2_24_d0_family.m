
//////////////////////////////////////////////////////////////////////
//  The d=0 Z/24 family: rank of the genus-1 component + point harvest.
//
//  From agent_a2_24_d0_derive.m: on d=0 (p=r(r+t)/2), the "split-q3
//  through x=0" Z/24 locus factors; the relevant component G(r,t)=0 is a
//  degree-6 genus-1 plane curve through curve B=(1/3,-1).  If its
//  elliptic model has positive rank we get INFINITELY many simple Z/24
//  curves; each point -> (r,t) -> a genus-2 curve; we log torsion,
//  simplicity, 2-rank (hunting any 2-rank 2 = [2,24] on this slice).
//
//  Usage: magma -b agent_a2_24_d0_family.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 6; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);
R2<r,t> := PolynomialRing(Q, 2);

// degree-6 component (x3 to clear denominators)
G := 3*r^3*t^3 + 2*r^3*t^2 - r^3*t - 2*r^2*t^4 + 8*r^2*t^3 - 22*r^2*t^2 + 4*r^2*t
     - r*t^5 + 6*r*t^4 - 13*r*t^3 + 12*r*t^2 + 8*r*t - 4;
printf "G(1/3,-1) = %o (expect 0)\n", Evaluate(G, [1/3, -1]);

AS := AffineSpace(R2);
Cg := Curve(AS, G);
PC := ProjectiveClosure(Cg);
printf "genus = %o\n", Genus(PC);

// elliptic model using the rational point B=(1/3,-1)
ptB := PC ! [1/3, -1];
E, mp := EllipticCurve(PC, ptB);
E := MinimalModel(E);
printf "Elliptic model: %o\n", E;
printf "Conductor = %o\n", Conductor(E);
tt, T := TorsionSubgroup(E);
printf "Torsion(E) = %o\n", Invariants(tt);
rk, rkbound := Rank(E);
printf "Rank(E) = %o (proven up to bound %o)\n", rk, rkbound;
gens := Generators(E);
printf "#Generators = %o\n", #gens;

// A(8) f on d=0
function A8f(rv, pv, tv)
    e := tv^2 - 2*pv*tv/rv; d := e + 2*pv - rv^2; lambda := rv/tv;
    a := rv^2 - lambda;
    b := 2*rv*pv - 2*lambda*(pv + rv*tv) + 2*rv*lambda;
    c := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
         - lambda*(rv^2 + e) + 2*lambda*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
    Qpoly := x^2 + d; q := a*x^2 + b*x + c;
    return q*(Qpoly^2 + q);
end function;
function IntModel(fv)
    L := 1; for i in [0..Degree(fv)] do L := LCM(L, Denominator(Coefficient(fv, i))); end for;
    return P!(L^2*fv);
end function;
function TwoRank(fp)
    degs := [Degree(g[1]) : g in Factorization(fp)]; k := #degs; even := 0;
    for mask in [0..2^k-1] do ss := 0;
        for i in [1..k] do if (mask div 2^(i-1)) mod 2 eq 1 then ss +:= degs[i]; end if; end for;
        if ss mod 2 eq 0 then even +:= 1; end if; end for;
    return Ilog2(even) - 1;
end function;

// harvest points on E, pull back to (r,t), build curve, check torsion+2rank
Einv := Inverse(mp);
printf "==== Z/24 curves from points on the family ====\n";
found := 0; r2hist := AssociativeArray();
pts := [];
if #gens gt 0 then
    g1 := gens[1];
    for k in [-8..8] do Append(~pts, k*g1 + (#gens ge 2 select gens[2] else E!0)); end for;
end if;
for tp in tt do Append(~pts, T(tp)); end for;   // torsion points too
for Ept in pts do
    try
        Ppc := Einv(Ept);
        cor := Coordinates(Ppc);
        rv := cor[1]/cor[3]; tv := cor[2]/cor[3];
    catch ee continue; end try;
    if rv eq 0 or tv eq 0 or rv*tv eq 1 then continue; end if;
    pv := rv*(rv+tv)/2;
    f := A8f(rv, pv, tv);
    if Degree(f) lt 5 or Discriminant(f) eq 0 then continue; end if;
    fInt := IntModel(f);
    J := Jacobian(HyperellipticCurve(fInt));
    inv := Invariants(TorsionSubgroup(J));
    has24 := #[nn : nn in inv | nn mod 24 eq 0] ge 1;
    if not has24 then continue; end if;
    found +:= 1;
    r2 := TwoRank(fInt);
    if IsDefined(r2hist, r2) then r2hist[r2] +:= 1; else r2hist[r2] := 1; end if;
    printf "  (r,t)=(%o,%o)  torsion=%o  2rank=%o\n", rv, tv, inv, r2;
end for;
printf "Z/24 curves found from family points: %o\n", found;
print "2-rank histogram:";
for k in Sort([kk : kk in Keys(r2hist)]) do printf "  2rank %o : %o\n", k, r2hist[k]; end for;
print "DONE";
quit;
