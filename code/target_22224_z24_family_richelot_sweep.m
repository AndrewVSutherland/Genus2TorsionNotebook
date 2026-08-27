//////////////////////////////////////////////////////////////////////
// Low-height members of the positive-rank d=0 Z/24 family, followed by
// rational Richelot-neighbor tests.  This avoids beginning at huge negative
// multiples (which makes exact torsion unnecessarily expensive).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(4*10^9);
SetLogFile("results/target_22224_z24_family_richelot_sweep.log" : Overwrite := true);
Q := Rationals(); P<x> := PolynomialRing(Q); R2<r,t> := PolynomialRing(Q,2);
G := 3*r^3*t^3 + 2*r^3*t^2 - r^3*t - 2*r^2*t^4 + 8*r^2*t^3
     - 22*r^2*t^2 + 4*r^2*t - r*t^5 + 6*r*t^4 - 13*r*t^3
     + 12*r*t^2 + 8*r*t - 4;
PC := ProjectiveClosure(Curve(AffineSpace(R2),G));
ptB := PC![1/3,-1];
E,mp := EllipticCurve(PC,ptB); E := MinimalModel(E); Einv := Inverse(mp);
gens := Generators(E);
free := [T : T in gens | Order(T) eq 0];
assert #free ge 1;

function A8f(rv,pv,tv)
    e:=tv^2-2*pv*tv/rv; d:=e+2*pv-rv^2; lam:=rv/tv;
    aa:=rv^2-lam;
    bb:=2*rv*pv-2*lam*(pv+rv*tv)+2*rv*lam;
    cc:=pv^2+2*pv*rv^2-rv^4-rv^3*tv-rv*pv^2/tv
        -lam*(rv^2+e)+2*lam*(rv*pv+rv^2*tv-3*pv*tv+rv*tv^2);
    q:=aa*x^2+bb*x+cc; Qq:=x^2+d;
    return q*(Qq^2+q);
end function;

function Normalize(f)
    den:=LCM([Denominator(Coefficient(f,i)):i in [0..Degree(f)]]);
    return SimplifiedModel(HyperellipticCurve(P!(den^2*f)));
end function;

print "TARGET_22224_Z24_FAMILY_RICHELOT_START","E",E;
seen := {};
for n in [0,1,-1,2,-2,3,-3,4,-4,5,-5] do
    ep := n*free[1];
    try pp := Einv(ep); catch err continue; end try;
    co := Coordinates(pp); if co[3] eq 0 then continue; end if;
    rv:=co[1]/co[3]; tv:=co[2]/co[3];
    if rv eq 0 or tv eq 0 or rv*tv eq 1 then continue; end if;
    pv:=rv*(rv+tv)/2; f:=A8f(rv,pv,tv);
    if Degree(f) lt 5 or Discriminant(f) eq 0 then continue; end if;
    gi := <rv,tv>; if gi in seen then continue; end if; Include(~seen,gi);
    facdeg:=Sort([Degree(z[1]):z in Factorization(f)]);
    print "Z24_FAMILY_SOURCE","n",n,"r",rv,"t",tv,"factor_degrees",facdeg;
    C:=Normalize(f); Rs:=RichelotIsogenousSurfaces(Jacobian(C));
    print "Z24_FAMILY_RICHELOTS","n",n,"count",#Rs;
    // Full rational 2-torsion at the source would itself be decisive, so
    // exact-test only that rare case.
    if facdeg eq [1,1,1,1,1,1] then
        inv:=Invariants(TorsionSubgroup(Jacobian(C)));
        print "Z24_FAMILY_FULL2_SOURCE","n",n,"torsion",inv;
        if inv eq [2,2,2,24] then print "TARGET_22224_HIT_SOURCE",n,C; end if;
    end if;
    for j in [1..#Rs] do
        if Type(Rs[j]) ne JacHyp then continue; end if;
        ff,hh:=HyperellipticPolynomials(Curve(Rs[j]));
        D:=Normalize(hh eq 0 select P!ff else P!(hh^2+4*ff));
        inv:=Invariants(TorsionSubgroup(Jacobian(D)));
        print "Z24_FAMILY_CODOMAIN","n",n,"edge",j,"torsion",inv;
        if inv eq [2,2,2,24] then print "TARGET_22224_HIT_CODOMAIN",n,j,D; end if;
    end for;
end for;
print "TARGET_22224_Z24_FAMILY_RICHELOT_DONE","sources",#seen;
quit;
