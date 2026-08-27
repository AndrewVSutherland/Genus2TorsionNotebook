//////////////////////////////////////////////////////////////////////
//  Geometry of the exact point-contact-5 cover on the Q4-root
//  rational surface inside M(12), in the compact (k,v) chart.
//
//  Put D=1-k, q=2*k-1, d=q^2/D, and
//
//      a=k*q^2/(D*v^2),  r=v-q^2/D.
//
//  With X=v/(x-v), a square-equivalent odd quintic is F=P*Q3, where
//
//      P  = 4*k*D^2*X^2+(2*k*q^2-v*D)*X+k*q^2,
//      R  = -q^2*D*X^2+(2*k*q^2-v*D)*X+k*q^2,
//      M  = q^2*X+v*D,
//      Q3 = M^2*P+4*k*q^2*D*(X+1)^2*R.
//
//  The X^4 coefficient in Q3 cancels.  The open chart is
//
//      k*(1-k)*(2*k-1)*v != 0.
//
//  For a contact abscissa u, let Ai=F^(i)(u)/i!.  This script forms
//
//      E3=8*A0^2*A3-A1*(4*A0*A2-A1^2),
//      E4=64*A0^3*A4-(4*A0*A2-A1^2)^2,
//
//  removes global/boundary factors, eliminates u, factors the projected
//  contact curve in Q[k,v], and recovers rational u and the final square
//  condition F(u)=c^2 on low-degree components.
//
//  It also analyzes two distinguished slices.  At u=-1, the contact
//  point collides with the marked order-12 point and should contribute
//  only singular/boundary solutions.  At u=0,
//
//      F(0)=(k*q^2)^2*D*(v^2*D+4*k*q^2),
//
//  so the missing square condition is the elementary conic
//
//      z^2=D*(v^2*D+4*k*q^2).
//
//  Intended later geometry run (not part of the finite/exact search):
//      magma -b PointBound:=100 MaxGeometryDegree:=14 \
//          code/m12_q4root_contact5_kv_geometry.m
//
//  Audited search baseline (2026-07-10):
//
//    * direct k-height 20: 511 values, 261121 pairs, 47845 survivors
//      through p=7,11,13,17, 47753 smooth, no contact point;
//    * branch-free q-height 50, q=2*k-1: 3095 q-values and v-values,
//      9579025 pairs, 12672 raw mask survivors through
//      p=7,11,13,17,19,23,29,31,37,41,43,47.  Removing
//      q=-1,0,1 and v=0 leaves 295 exact states: 247 singular,
//      8 degree-drop, and 40 smooth with no common E3/E4 root.
//
//  Thus there is no contact through this meaningful low-height range, but
//  the finite masks have smooth contact residues at every tested prime.
//  After the separate a=-1/4 boundary elimination, the saturated affine
//  resultant/component analysis below is the next step; another blind
//  height increase should come only after that geometry is understood.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned MemGB then
    MemGB := 8;
elif Type(MemGB) eq MonStgElt then
    MemGB := StringToInteger(MemGB);
end if;
if not assigned PointBound then
    PointBound := 100;
elif Type(PointBound) eq MonStgElt then
    PointBound := StringToInteger(PointBound);
end if;
if not assigned MaxGeometryDegree then
    MaxGeometryDegree := 14;
elif Type(MaxGeometryDegree) eq MonStgElt then
    MaxGeometryDegree := StringToInteger(MaxGeometryDegree);
end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Z := Integers();
R3<k,v,u> := PolynomialRing(Q,3,"grevlex");
PX<X> := PolynomialRing(R3);

D := 1-k;
q := 2*k-1;
P2 := 4*k*D^2*X^2+(2*k*q^2-v*D)*X+k*q^2;
R2 := -q^2*D*X^2+(2*k*q^2-v*D)*X+k*q^2;
M := q^2*X+v*D;
Q3 := M^2*P2+4*k*q^2*D*(1+X)^2*R2;
F := P2*Q3;

assert Degree(Q3) eq 3;
assert Degree(F) eq 5;

A0 := R3!Evaluate(F,u);
A1 := R3!Evaluate(Derivative(F),u);
A2 := R3!Evaluate(Derivative(F,2),u)/2;
A3 := R3!Evaluate(Derivative(F,3),u)/6;
A4 := R3!Evaluate(Derivative(F,4),u)/24;
Q2 := 4*A0*A2-A1^2;
E3raw := 8*A0^2*A3-A1*Q2;
E4raw := 64*A0^3*A4-Q2^2;

function PrimitiveR3Polynomial(g)
    if g eq 0 then return g; end if;
    den := LCM([Denominator(c) : c in Coefficients(g)]);
    gz := [Z!(den*c) : c in Coefficients(g)];
    cont := GCD(gz);
    if cont eq 0 then cont := 1; end if;
    return R3!((Q!den/Q!Abs(cont))*g);
end function;

E3 := PrimitiveR3Polynomial(E3raw);
E4 := PrimitiveR3Polynomial(E4raw);
boundary := k*D*q*v;

print "M12 Q4-ROOT CONTACT5 COMPACT (k,v) GEOMETRY";
print "F degree",Degree(F),"terms",#Terms(F),
      "leading_coefficient",LeadingCoefficient(F),
      "constant_coefficient",Coefficient(F,0);
print "P2",P2;
print "Q3 degree",Degree(Q3),"terms",#Terms(Q3),Q3;
print "A0 degree",TotalDegree(A0),"terms",#Terms(A0);
print "E3 total_degree",TotalDegree(E3),"degree_u",Degree(E3,u),
      "terms",#Terms(E3);
print "E4 total_degree",TotalDegree(E4),"degree_u",Degree(E4,u),
      "terms",#Terms(E4);

print "FACTOR E3";
fac3 := Factorization(E3);
for fe in fac3 do
    print fe[2],"degree",TotalDegree(fe[1]),"degree_u",Degree(fe[1],u),
          "terms",#Terms(fe[1]),fe[1];
end for;
print "FACTOR E4";
fac4 := Factorization(E4);
for fe in fac4 do
    print fe[2],"degree",TotalDegree(fe[1]),"degree_u",Degree(fe[1],u),
          "terms",#Terms(fe[1]),fe[1];
end for;

common := GCD(E3,E4);
print "COMMON_GCD degree",TotalDegree(common),"degree_u",Degree(common,u),
      "terms",#Terms(common),common;
E3open := ExactQuotient(E3,common);
E4open := ExactQuotient(E4,common);

procedure StripKnownBoundary(~g,label)
    for h in [k,D,q,v] do
        multiplicity := 0;
        while g mod h eq 0 do
            g := ExactQuotient(g,h);
            multiplicity +:= 1;
        end while;
        if multiplicity gt 0 then
            print "STRIPPED",label,h,"multiplicity",multiplicity;
        end if;
    end for;
end procedure;

StripKnownBoundary(~E3open,"E3");
StripKnownBoundary(~E4open,"E4");
print "OPEN E3 degree",TotalDegree(E3open),"degree_u",Degree(E3open,u),
      "terms",#Terms(E3open);
print "OPEN E4 degree",TotalDegree(E4open),"degree_u",Degree(E4open,u),
      "terms",#Terms(E4open);

//////////////////////////////////////////////////////////////////////
// Distinguished contact-coordinate slices.
//////////////////////////////////////////////////////////////////////

print "SLICE u=0";
F0 := R3!Evaluate(A0,[k,v,R3!0]);
F0expected := (k*q^2)^2*D*(v^2*D+4*k*q^2);
print "F(0)",F0,"expected",F0expected,"identity",F0 eq F0expected;
E30 := R3!Evaluate(E3open,[k,v,R3!0]);
E40 := R3!Evaluate(E4open,[k,v,R3!0]);
print "E3(u=0) factors";
for fe in Factorization(E30) do print fe; end for;
print "E4(u=0) factors";
for fe in Factorization(E40) do print fe; end for;
g0 := GCD(E30,E40);
print "u=0 gcd degree",TotalDegree(g0),"terms",#Terms(g0),g0;
print "u=0 residual square cover z^2 =",D*(v^2*D+4*k*q^2);

print "SLICE u=-1 (marked D12 collision)";
Fm1 := R3!Evaluate(A0,[k,v,R3!-1]);
ym1 := (v*D-q^2)*(v*D-q^2+D);
print "F(-1)",Fm1,"expected",ym1^2,"identity",Fm1 eq ym1^2;
E3m1 := R3!Evaluate(E3open,[k,v,R3!-1]);
E4m1 := R3!Evaluate(E4open,[k,v,R3!-1]);
print "E3(u=-1) factors";
for fe in Factorization(E3m1) do print fe; end for;
print "E4(u=-1) factors";
for fe in Factorization(E4m1) do print fe; end for;
gm1 := GCD(E3m1,E4m1);
print "u=-1 gcd degree",TotalDegree(gm1),"terms",#Terms(gm1),gm1;

print "DISCRIMINANT F";
discF := R3!Discriminant(F);
for fe in Factorization(discF) do
    print "DISC_FACTOR",fe[2],"degree",TotalDegree(fe[1]),
          "terms",#Terms(fe[1]),fe[1];
end for;

print "RESULTANT eliminating u";
time ResRaw := Resultant(E3open,E4open,u);
Res := PrimitiveR3Polynomial(ResRaw);
print "RESULTANT total_degree",TotalDegree(Res),"terms",#Terms(Res);
resfac := Factorization(Res);
print "RESULTANT factors",#resfac;
for fe in resfac do
    isKnownBoundary := GCD(fe[1],boundary) ne 1;
    isDiscBoundary := GCD(fe[1],discF) ne 1;
    print "RES_FACTOR","multiplicity",fe[2],
          "degree",TotalDegree(fe[1]),"terms",#Terms(fe[1]),
          "known_boundary",isKnownBoundary,
          "discriminant_boundary",isDiscBoundary;
    if #Terms(fe[1]) le 120 then print fe[1]; end if;
end for;

//////////////////////////////////////////////////////////////////////
// Low-degree projected components and rational-point recovery.
//////////////////////////////////////////////////////////////////////

B<K,V> := PolynomialRing(Q,2);
toB := hom<R3 -> B | K,V,B!0>;

function PrimitiveBPolynomial(g)
    if g eq 0 then return g; end if;
    den := LCM([Denominator(c) : c in Coefficients(g)]);
    gz := [Z!(den*c) : c in Coefficients(g)];
    cont := GCD(gz);
    if cont eq 0 then cont := 1; end if;
    return B!((Q!den/Q!Abs(cont))*g);
end function;

function SpecializeToU(g,kv,vv)
    QU<U> := PolynomialRing(Q);
    mp := hom<R3 -> QU | QU!kv,QU!vv,U>;
    return QU!mp(g);
end function;

baseFactors := [];
for fe in resfac do
    g := PrimitiveBPolynomial(B!toB(fe[1]));
    if g ne 0 then Append(~baseFactors,<g,fe[2]>); end if;
end for;

print "PROJECTED COMPONENTS";
for i in [1..#baseFactors] do
    g := baseFactors[i][1];
    mult := baseFactors[i][2];
    knownBoundary := GCD(g,K*(1-K)*(2*K-1)*V) ne 1;
    print "BASE_COMPONENT",i,"multiplicity",mult,
          "degree",TotalDegree(g),"terms",#Terms(g),
          "known_boundary",knownBoundary;
    if knownBoundary or TotalDegree(g) gt MaxGeometryDegree then
        if TotalDegree(g) gt MaxGeometryDegree then
            print " geometry_skipped_degree_limit",MaxGeometryDegree;
        end if;
        continue;
    end if;
    try
        Caff := Curve(AffineSpace(B),g);
        Cp := ProjectiveClosure(Caff);
        print " projective_degree",Degree(Cp),
              "normalized_genus",Genus(Cp),
              "nonsingular",IsNonsingular(Cp);
        pts := Points(Cp : Bound := PointBound);
        print " points_bound",PointBound,"count",#pts;
        if #pts le 40 then print " points",pts; end if;
        for pt in pts do
            cc := Eltseq(pt);
            if cc[3] eq 0 then continue; end if;
            kv := Q!cc[1]/Q!cc[3];
            vv := Q!cc[2]/Q!cc[3];
            if kv*(1-kv)*(2*kv-1)*vv eq 0 then continue; end if;
            e3u := SpecializeToU(E3open,kv,vv);
            e4u := SpecializeToU(E4open,kv,vv);
            gu := GCD(e3u,e4u);
            for rt in Roots(gu) do
                uv := Q!rt[1];
                f0 := Q!Evaluate(A0,[kv,vv,uv]);
                square := f0 ne 0 and IsSquare(f0);
                print " RECOVERED","k",kv,"v",vv,"u",uv,
                      "f0",f0,"square",square;
            end for;
        end for;
    catch eg
        print " geometry_or_points_failed",eg`Object;
    end try;
end for;

quit;
