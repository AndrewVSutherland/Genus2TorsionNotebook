//////////////////////////////////////////////////////////////////////
//  Exact projected P55 component of the M(12)+point-contact-5 cover.
//
//  Normalized root chart:
//    F=L*(L*A^2+4*b*(1+x)^2*(w*L-x^2)),
//    L=b+(2b-1)x, A=x+w*(1+b*x).
//
//  E3=E4=0 are the exact contact covariants.  After removing the common
//  b-boundary and eliminating w, the resultant has one irreducible factor
//  P55(b,u), bidegree (26,29), total degree 55.  This script reconstructs
//  it from first principles, records every excluded boundary, audits its
//  projective finite-field points, and optionally asks Magma for the
//  normalized geometric genus.
//
//  Modes: summary (default), local, genus.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
if mode notin {"summary","local","genus"} then
    error "mode must be summary, local, or genus";
end if;
if not assigned MemGB then MemGB := 8;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned DoF3Places then DoF3Places := false; end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals(); Z := Integers();
R<b,w,u> := PolynomialRing(Q,3,"grevlex");
PX<x> := PolynomialRing(R);
Ln := b+(2*b-1)*x;
An := x+w*(1+b*x);
Fn := Ln*(Ln*An^2+4*b*(1+x)^2*(w*Ln-x^2));

A0 := R!Evaluate(Fn,u);
A1 := R!Evaluate(Derivative(Fn),u);
A2 := R!Evaluate(Derivative(Fn,2),u)/2;
A3 := R!Evaluate(Derivative(Fn,3),u)/6;
A4 := R!Evaluate(Derivative(Fn,4),u)/24;
Q2 := 4*A0*A2-A1^2;
E3raw := 8*A0^2*A3-A1*Q2;
E4raw := 64*A0^3*A4-Q2^2;

function Primitive(poly)
    if poly eq 0 then return poly; end if;
    den := LCM([Denominator(c) : c in Coefficients(poly)]);
    vals := [Z!(den*c) : c in Coefficients(poly)];
    cont := GCD(vals); if cont eq 0 then cont := 1; end if;
    return Parent(poly)!((Q!den/Q!Abs(cont))*poly);
end function;

E3 := Primitive(E3raw);
E4 := Primitive(E4raw);
assert IsDivisibleBy(E3,b) and IsDivisibleBy(E4,b);
E3 := ExactQuotient(E3,b);
E4 := ExactQuotient(E4,b);
assert GCD(E3,E4) eq 1;
print "P55_SELF_TEST_PASS","E3",<TotalDegree(E3),Degree(E3,w),Degree(E3,u),#Terms(E3)>,
      "E4",<TotalDegree(E4),Degree(E4,w),Degree(E4,u),#Terms(E4)>;

print "P55_RESULTANT_W_BEGIN";
time Res3 := Resultant(E3,E4,w);
assert Degree(Res3,w) eq 0;
B2<B,U> := PolynomialRing(Q,2,"grevlex");
toB2 := hom<R -> B2 | B,B2!0,U>;
Res := Primitive(B2!toB2(Res3));
for h in [B,B-1,2*B-1] do
    while IsDivisibleBy(Res,h) do Res := ExactQuotient(Res,h); end while;
end for;
fac := Factorization(Res);
p55s := [fe[1] : fe in fac | TotalDegree(fe[1]) eq 55];
assert #p55s eq 1;
P55 := Primitive(p55s[1]);
assert IsIrreducible(P55);
print "P55_EXACT","degree",TotalDegree(P55),"degree_b",Degree(P55,B),
      "degree_u",Degree(P55,U),"terms",#Terms(P55),
      "multiplicity",[fe[2] : fe in fac | fe[1] eq p55s[1]][1];

// Full open boundary before the missing square cover c^2=F(u).
// Parameter boundary: b*w*(b-1)*(2b-1).
// Smoothness: discriminant factors (b*w-w+1) and D12(b,w).
// Contact boundary: F(u)=0, u=-1, and L(u)=0.
disc := Primitive(R!Discriminant(Fn));
for h in [b,w,b-1,2*b-1] do
    while IsDivisibleBy(disc,h) do disc := ExactQuotient(disc,h); end while;
end for;
discfac := Factorization(disc);
print "P55_OPEN_DISCRIMINANT_FACTORS";
for fe in discfac do
    print " DISC",fe[2],TotalDegree(fe[1]),#Terms(fe[1]),fe[1];
end for;
print "P55_CONTACT_BOUNDARY","u+1",u+1,"L(u)",b+(2*b-1)*u,
      "square_equation c^2=F(u)",A0;

// Exact specializations of the projected plane component.  These identify
// all visible affine/projective boundary intersections before normalization.
QB<T> := PolynomialRing(Q);
function SpecB(bv)
    mp := hom<B2 -> QB | QB!bv,T>;
    return QB!mp(P55);
end function;
function SpecU(uv)
    mp := hom<B2 -> QB | T,QB!uv>;
    return QB!mp(P55);
end function;
for bv in [Q!0,Q!1,Q!1/2] do
    g := SpecB(bv);
    print "P55_B_SPECIALIZATION",bv,"degree",Degree(g),
          "factorization",Factorization(g);
end for;
for uv in [Q!0,Q!-1] do
    g := SpecU(uv);
    print "P55_U_SPECIALIZATION",uv,"degree",Degree(g),
          "factorization",Factorization(g);
end for;

P2<Xb,Xu,Xz> := ProjectiveSpace(Q,2);
P55h := Homogenization(Evaluate(P55,[Xb,Xu]),Xz);
assert IsHomogeneous(P55h) and Degree(P55h) eq 55;
QInf<Ti> := PolynomialRing(Q);
inf1 := QInf!Evaluate(P55h,[QInf!1,Ti,QInf!0]);
infVertical := Q!Evaluate(P55h,[Q!0,Q!1,Q!0]);
print "P55_INFINITY_CHART degree",Degree(inf1),
      "factorization",Factorization(inf1),
      "vertical_point_value",infVertical;

if mode eq "summary" then
    print "P55_SUMMARY_DONE";
    quit;
end if;

if mode eq "local" then
    plist := [3,5,7,11,13,17,19,23,29,31,37,41,43];
    for p in plist do
        Fp := GF(p);
        Rp<bp,up,zp> := PolynomialRing(Fp,3,"grevlex");
        red := hom<Parent(P55h) -> Rp | bp,up,zp>;
        hp := red(P55h);
        db := Derivative(hp,bp); du := Derivative(hp,up); dz := Derivative(hp,zp);
        total := 0; smooth := 0; affine := 0; infinity := 0;
        procedure PrintLocalTangent(loc,label)
            nz := [m : m in Monomials(loc) | MonomialCoefficient(loc,m) ne 0];
            if #nz eq 0 then
                print "P55_LOCAL_ZERO_POLYNOMIAL",label;
                return;
            end if;
            md := Min([TotalDegree(m) : m in nz]);
            tc := &+[MonomialCoefficient(loc,m)*m : m in nz | TotalDegree(m) eq md];
            print "P55_TANGENT_CONE",label,"multiplicity",md,
                  "factorization",Factorization(tc);
        end procedure;
        for bv in Fp do for uv in Fp do
            pt := [bv,uv,Fp!1];
            if Evaluate(hp,pt) eq 0 then
                total +:= 1; affine +:= 1;
                if &or[Evaluate(d,pt) ne 0 : d in [db,du,dz]] then smooth +:= 1; end if;
                if p eq 3 then
                    print "P55_F3_PLANE_POINT",pt;
                    PrintLocalTangent(Evaluate(hp,[bp+bv,up+uv,Fp!1]),<"z=1",bv,uv>);
                end if;
            end if;
        end for; end for;
        for tv in Fp do
            pt := [Fp!1,tv,Fp!0];
            if Evaluate(hp,pt) eq 0 then
                total +:= 1; infinity +:= 1;
                if &or[Evaluate(d,pt) ne 0 : d in [db,du,dz]] then smooth +:= 1; end if;
                if p eq 3 then
                    print "P55_F3_PLANE_POINT",pt;
                    PrintLocalTangent(Evaluate(hp,[Fp!1,up+tv,zp]),<"b=1",tv>);
                end if;
            end if;
        end for;
        pt := [Fp!0,Fp!1,Fp!0];
        if Evaluate(hp,pt) eq 0 then
            total +:= 1; infinity +:= 1;
            if &or[Evaluate(d,pt) ne 0 : d in [db,du,dz]] then smooth +:= 1; end if;
            if p eq 3 then
                print "P55_F3_PLANE_POINT",pt;
                PrintLocalTangent(Evaluate(hp,[bp,Fp!1,zp]),<"u=1">);
            end if;
        end if;
        print "P55_LOCAL",p,"plane_points",total,"affine",affine,
              "infinity",infinity,"nonsingular_plane_points",smooth;
        if p in {3,5,7} then
            time ffp := Factorization(hp);
            print "P55_REDUCTION_FACTORS",p,
                  [<TotalDegree(fe[1]),fe[2],#Terms(fe[1])> : fe in ffp];
            if DoF3Places and p eq 3 and #ffp eq 1 and ffp[1][2] eq 1 then
                P2p<Yb,Yu,Yz> := ProjectiveSpace(Fp,2);
                hcurve := Evaluate(hp,[Yb,Yu,Yz]);
                C3 := Curve(P2p,hcurve);
                print "P55_F3_NORMALIZATION_BEGIN";
                time pl1 := Places(C3,1);
                print "P55_F3_DEGREE1_PLACES",#pl1;
                for pl in pl1 do print " F3_PLACE",pl; end for;
            end if;
        end if;
    end for;
    print "P55_LOCAL_DONE";
    quit;
end if;

print "P55_GENUS_BEGIN";
C55 := Curve(P2,P55h);
print "P55_PROJECTIVE irreducible",IsIrreducible(C55),
      "degree",Degree(C55),"arithmetic_genus",(55-1)*(55-2) div 2;
time sing := SingularPoints(C55);
print "P55_RATIONAL_SINGULAR_POINTS",#sing;
for pp in sing do print " P55_SING",pp,"mult",Multiplicity(C55,pp); end for;
time gg := Genus(C55);
print "P55_NORMALIZED_GENUS",gg;
print "P55_GENUS_DONE";
quit;
