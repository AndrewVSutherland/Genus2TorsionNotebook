//////////////////////////////////////////////////////////////////////
//  Geometry of the exact point-contact-5 cover of the split M(12)
//  surface in the normalized rational-root chart.
//
//  Starting from T=a*x^2-x+r and a rational root w of T+1, put
//
//      b=a*w,  z=2*b-1,  r=w*(1-b)-1.
//
//  After xi=w*X and eta=w^2*Y, the odd quintic is
//
//      F = L * (L*A^2 + 4*b*(1+xi)^2*(w*L-xi^2)),
//      L = b+(2*b-1)*xi,
//      A = xi+w*(1+b*xi).
//
//  The marked order-12 point is xi=-1.  The open parameter chart has
//
//      b*w*(1-b)*(2*b-1) != 0.
//
//  For a contact abscissa u, let Ai=F^(i)(u)/i!.  The exact point-contact
//  equations are
//
//      E3 = 8*A0^2*A3-A1*(4*A0*A2-A1^2),
//      E4 = 64*A0^3*A4-(4*A0*A2-A1^2)^2.
//
//  This script factors E3,E4 and their common/boundary pieces, eliminates
//  u, analyzes low-degree projected components, tests the root-swap
//  involution b -> 1-b, w -> w*(1-b)/b, and recovers rational u and the
//  final square condition c^2=F(u) at small rational points.
//
//  Intended short run (from torsion_jac):
//      magma -b PointBound:=100 MaxGeometryDegree:=12 \
//          code/m12_contact5_bw_geometry.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned MemGB then
    MemGB := 6;
elif Type(MemGB) eq MonStgElt then
    MemGB := StringToInteger(MemGB);
end if;
if not assigned PointBound then
    PointBound := 100;
elif Type(PointBound) eq MonStgElt then
    PointBound := StringToInteger(PointBound);
end if;
if not assigned MaxGeometryDegree then
    MaxGeometryDegree := 12;
elif Type(MaxGeometryDegree) eq MonStgElt then
    MaxGeometryDegree := StringToInteger(MaxGeometryDegree);
end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Z := Integers();
R<b,w,u> := PolynomialRing(Q,3,"grevlex");
PX<xi> := PolynomialRing(R);

L := b+(2*b-1)*xi;
A := xi+w*(1+b*xi);
F := L*(L*A^2+4*b*(1+xi)^2*(w*L-xi^2));

A0 := R!Evaluate(F,u);
A1 := R!Evaluate(Derivative(F),u);
A2 := R!Evaluate(Derivative(F,2),u)/2;
A3 := R!Evaluate(Derivative(F,3),u)/6;
A4 := R!Evaluate(Derivative(F,4),u)/24;
Q2 := 4*A0*A2-A1^2;
E3raw := 8*A0^2*A3-A1*Q2;
E4raw := 64*A0^3*A4-Q2^2;

function PrimitiveQPolynomial(g)
    if g eq 0 then
        return g;
    end if;
    den := LCM([Denominator(c) : c in Coefficients(g)]);
    gz := [Z!(den*c) : c in Coefficients(g)];
    cont := GCD(gz);
    if cont eq 0 then
        cont := 1;
    end if;
    return R!((Q!den/Q!Abs(cont))*g);
end function;

E3 := PrimitiveQPolynomial(E3raw);
E4 := PrimitiveQPolynomial(E4raw);

boundary := b*w*(b-1)*(2*b-1);

print "M12 CONTACT5 NORMALIZED ROOT CHART";
print "F degree",Degree(F),"terms",#Terms(F),
      "leading_coefficient",LeadingCoefficient(F),
      "constant_coefficient",Coefficient(F,0);
print "F factorization";
for fe in Factorization(F) do
    print fe[2],"degree",Degree(fe[1]),"terms",#Terms(fe[1]),fe[1];
end for;
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
print "COMMON_GCD factors";
for fe in Factorization(common) do print fe; end for;

E3open := ExactQuotient(E3,common);
E4open := ExactQuotient(E4,common);

// Remove any remaining factors supported entirely on a known parameter
// boundary.  Irreducible factors meeting a boundary in codimension one are
// retained; only factors dividing the boundary polynomial are removed.
procedure StripKnownBoundary(~g, label)
    for h in [b,w,b-1,2*b-1] do
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
// Two distinguished contact-coordinate slices.
//////////////////////////////////////////////////////////////////////

print "SLICE u=0";
F0slice := R!Evaluate(A0,[b,w,R!0]);
E30 := R!Evaluate(E3open,[b,w,R!0]);
E40 := R!Evaluate(E4open,[b,w,R!0]);
print "F(0)",F0slice,"expected",b^2*w*(w+4*b),
      "identity",F0slice eq b^2*w*(w+4*b);
print "E3(u=0) factors";
for fe in Factorization(E30) do print fe; end for;
print "E4(u=0) factors";
for fe in Factorization(E40) do print fe; end for;
g0 := GCD(E30,E40);
print "u=0 gcd degree",TotalDegree(g0),"terms",#Terms(g0),g0;

// On the open b*w != 0 chart, c^2=F(0) is the conic
// v^2=w(w+4b).  Put v=t*w, hence b=w*(t^2-1)/4.
S<Wv,tv> := PolynomialRing(Q,2);
u0square := hom<R -> S | Wv*(tv^2-1)/4,Wv,S!0>;
E30sq := S!u0square(E30);
E40sq := S!u0square(E40);
print "u=0 square-param E3 factors";
for fe in Factorization(E30sq) do print fe; end for;
print "u=0 square-param E4 factors";
for fe in Factorization(E40sq) do print fe; end for;
g0sq := GCD(E30sq,E40sq);
print "u=0 square-param gcd degree",TotalDegree(g0sq),
      "terms",#Terms(g0sq),g0sq;

print "SLICE u=-1 (collision with marked D12 abscissa)";
Fm1slice := R!Evaluate(A0,[b,w,R!-1]);
print "F(-1)",Fm1slice,
      "expected",(1-b)^2*(w*(1-b)-1)^2,
      "identity",Fm1slice eq (1-b)^2*(w*(1-b)-1)^2;
E3m1 := R!Evaluate(E3open,[b,w,R!-1]);
E4m1 := R!Evaluate(E4open,[b,w,R!-1]);
print "E3(u=-1) factors";
for fe in Factorization(E3m1) do print fe; end for;
print "E4(u=-1) factors";
for fe in Factorization(E4m1) do print fe; end for;
gm1 := GCD(E3m1,E4m1);
print "u=-1 gcd degree",TotalDegree(gm1),"terms",#Terms(gm1),gm1;

// The discriminant gives the complete singular-fiber boundary in this
// denominator-free chart and helps classify resultant factors.
print "DISCRIMINANT F";
discF := R!Discriminant(F);
for fe in Factorization(discF) do
    print "DISC_FACTOR",fe[2],"degree",TotalDegree(fe[1]),
          "terms",#Terms(fe[1]),fe[1];
end for;

print "RESULTANT eliminating u";
time ResRaw := Resultant(E3open,E4open,u);
Res := PrimitiveQPolynomial(ResRaw);
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
    if #Terms(fe[1]) le 120 then
        print fe[1];
    end if;
end for;

//////////////////////////////////////////////////////////////////////
// Project to Q[b,w], check low-degree geometry and the root-swap map.
//////////////////////////////////////////////////////////////////////

B<br,wr> := PolynomialRing(Q,2);
toB := hom<R -> B | br,wr,B!0>;
BF := FieldOfFractions(B);

function PrimitiveBPolynomial(g)
    if g eq 0 then return g; end if;
    den := LCM([Denominator(c) : c in Coefficients(g)]);
    gz := [Z!(den*c) : c in Coefficients(g)];
    cont := GCD(gz);
    if cont eq 0 then cont := 1; end if;
    return B!((Q!den/Q!Abs(cont))*g);
end function;

function RootSwapPolynomial(g)
    val := Evaluate(BF!g,[BF!(1-br),BF!(wr*(1-br)/br)]);
    return PrimitiveBPolynomial(B!Numerator(val));
end function;

function SpecializeToU(g,bv,wv)
    QU<U> := PolynomialRing(Q);
    mp := hom<R -> QU | QU!bv,QU!wv,U>;
    return QU!mp(g);
end function;

print "PROJECTED COMPONENTS";
baseFactors := [];
for fe in resfac do
    g := PrimitiveBPolynomial(B!toB(fe[1]));
    if g eq 0 then continue; end if;
    Append(~baseFactors,<g,fe[2]>);
end for;

for i in [1..#baseFactors] do
    g := baseFactors[i][1];
    mult := baseFactors[i][2];
    knownBoundary := GCD(g,br*wr*(br-1)*(2*br-1)) ne 1;
    print "BASE_COMPONENT",i,"multiplicity",mult,
          "degree",TotalDegree(g),"terms",#Terms(g),
          "known_boundary",knownBoundary;

    swapped := RootSwapPolynomial(g);
    matches := [];
    for j in [1..#baseFactors] do
        h := baseFactors[j][1];
        gg := GCD(swapped,h);
        if TotalDegree(gg) gt 0 then
            Append(~matches,<j,TotalDegree(gg)>);
        end if;
    end for;
    print " root_swap_degree",TotalDegree(swapped),"matches",matches;

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

        // Recover rational contact abscissas from affine rational points and
        // impose the missing double-cover condition c^2=F(u).
        for pt in pts do
            cc := Eltseq(pt);
            if cc[3] eq 0 then continue; end if;
            bv := Q!cc[1]/Q!cc[3];
            wv := Q!cc[2]/Q!cc[3];
            if bv*wv*(bv-1)*(2*bv-1) eq 0 then continue; end if;
            e3u := SpecializeToU(E3open,bv,wv);
            e4u := SpecializeToU(E4open,bv,wv);
            gu := GCD(e3u,e4u);
            for rt in Roots(gu) do
                uv := Q!rt[1];
                f0 := Q!Evaluate(A0,[bv,wv,uv]);
                square := f0 ne 0 and IsSquare(f0);
                print " RECOVERED","b",bv,"w",wv,"u",uv,
                      "f0",f0,"square",square;
            end for;
        end for;
    catch eg
        print " geometry_or_points_failed",eg`Object;
    end try;
end for;

quit;
