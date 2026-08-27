//////////////////////////////////////////////////////////////////////
// Cubic-contact geometry on the removable s=1 branch of the contact-7
// two-root surface.
//
// Here r=1-s^2=0.  The limiting root condition is a=35/8, and a second
// root w=1-t^2 is imposed by
//
//   b=(t^7-1+(7/2)w-(35/8)w^2)/w^3.
//
// Thus f=(h^2+(x-1)^7)/x^2 has roots 0,w.  For cubic contact write
//
//   q=x^2+U*x+V, H=x^3+A*x^2+B*x+E,
//   H^2-q^3=M*f.
//
// High-to-low coefficient comparison reconstructs A,B,E.  Since f(0)=0,
// the remaining equations are G2=G1=0 and E^2=V^3.  Their open curve in
// (t,U,V,M), saturated by M, is the algebraic quotient of the exact
// rational-3 cover; the exact cover further adjoins L^2=M.
//
// Modes: summary, factor, eliminate.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
if mode notin {"summary","factor","eliminate"} then
    error "mode must be summary, factor, or eliminate";
end if;
if not assigned MemGB then MemGB := 8;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals(); Z := Integers();
R<t,U,V,M> := PolynomialRing(Q,4,"grevlex");
K := FieldOfFractions(R);
KX<x> := PolynomialRing(K);

w := 1-t^2;
a := Q!35/8;
b := K!((t^7-1+(Q!7/2)*w-a*w^2)/w^3);
h := 1-(Q!7/2)*x+a*x^2+b*x^3;
f := ExactQuotient(h^2+(x-1)^7,x^2);
c := [K!Coefficient(f,i) : i in [0..5]];
assert c[1] eq 0 and Evaluate(f,K!w) eq 0;

A := (M*c[6]+3*U)/2;
B := (M*c[5]+3*(U^2+V)-A^2)/2;
E := (M*c[4]+U^3+6*U*V-2*A*B)/2;
q := x^2+U*x+V;
H := x^3+A*x^2+B*x+E;
ident := H^2-q^3-M*f;
assert &and[Numerator(Coefficient(ident,i)) eq 0 : i in [3..6]];

function Primitive(f0)
    if f0 eq 0 then return R!0; end if;
    p := R!f0;
    den := LCM([Denominator(cc) : cc in Coefficients(p)]);
    ints := [Z!(den*cc) : cc in Coefficients(p)];
    cont := GCD(ints); if cont eq 0 then cont := 1; end if;
    return R!((Q!den/Q!Abs(cont))*p);
end function;

G2 := Primitive(Numerator(Coefficient(ident,2)));
G1 := Primitive(Numerator(Coefficient(ident,1)));
G0 := Primitive(Numerator(Coefficient(ident,0)));

print "CONTACT7_TWO_ROOT_PLUS3_S1",mode;
print "b numerator/denominator",<Degree(Numerator(b),1),Degree(Denominator(b),1),
      #Terms(Numerator(b)),#Terms(Denominator(b))>;
print "f coefficient num/den degrees",
      [<Degree(Numerator(ci),1),Degree(Denominator(ci),1)> : ci in c];
print "equations",[<TotalDegree(g),Degree(g,t),Degree(g,U),Degree(g,V),
      Degree(g,M),#Terms(g)> : g in [G2,G1,G0]];
disc := K!Discriminant(f);
print "discriminant num/den degrees",
      <Degree(Numerator(disc)),Degree(Denominator(disc))>;
print "discriminant numerator factors",
      [<Degree(fe[1]),fe[2],fe[1]> : fe in Factorization(Numerator(disc))];
print "discriminant denominator factors",
      [<Degree(fe[1]),fe[2],fe[1]> : fe in Factorization(Denominator(disc))];

if mode eq "summary" then
    print "b",b;
    print "f",f;
    quit;
end if;

if mode eq "factor" then
    for i in [1..3] do
        time fac := Factorization([G2,G1,G0][i]);
        print "FACTOR",i,[<fe[2],TotalDegree(fe[1]),Degree(fe[1],t),
              Degree(fe[1],U),Degree(fe[1],V),Degree(fe[1],M),#Terms(fe[1])>
              : fe in fac];
    end for;
    I := ideal<R|G2,G1,G0>;
    time dim := Dimension(I);
    print "RAW_DIMENSION",dim;
    quit;
end if;

I := ideal<R|G2,G1,G0>;
print "SATURATE_M_BEGIN";
time Isat := Saturation(I,ideal<R|M>);
print "SATURATE_M_DONE",Dimension(Isat),#Basis(Isat),
      [TotalDegree(g) : g in Basis(Isat)];

// Eliminate U,M to a plane model in (t,V).
RL<UL,ML,TL,VL> := PolynomialRing(Q,4,"lex");
mp := hom<R->RL|TL,UL,VL,ML>;
J := ideal<RL|[mp(g):g in Basis(Isat)]>;
print "LEX_BEGIN";
time gb := GroebnerBasis(J);
print "LEX_DONE",#gb,[TotalDegree(g):g in gb];
elim := [g:g in gb|Degree(g,UL) eq 0 and Degree(g,ML) eq 0];
print "TV_ELIM",#elim;
if #elim eq 0 then quit; end if;
P<TP,VP> := PolynomialRing(Q,2);
toP := hom<RL->P|P!0,P!0,TP,VP>;
core := P!toP(elim[1]);
for i in [2..#elim] do core := GCD(core,P!toP(elim[i])); end for;
for bd in [VP,TP-1,TP+1] do
    while core ne 0 and IsDivisibleBy(core,bd) do core:=ExactQuotient(core,bd); end while;
end for;
print "TV_CORE",<TotalDegree(core),Degree(core,TP),Degree(core,VP),#Terms(core)>;
fac := Factorization(core);
print "TV_FACTORS",#fac;
for i in [1..#fac] do
    pc:=fac[i][1];
    print " FACTOR",i,fac[i][2],<TotalDegree(pc),Degree(pc,TP),Degree(pc,VP),#Terms(pc)>;
    try
        Ca:=Curve(AffineSpace(P),pc); Cp:=ProjectiveClosure(Ca);
        print "  PROJECTIVE",Degree(Cp),Genus(Cp),IsNonsingular(Cp);
        Cn,nu:=Normalization(Cp);
        print "  NORMALIZATION",Genus(Cn),Degree(Cn),IsNonsingular(Cn);
        if Genus(Cn) ge 2 and Genus(Cn) le 30 then
            try
                aut,ma:=AutomorphismGroup(Cn);
                print "  AUTOMORPHISM_ORDER",#aut;
            catch err
                print "  AUTOMORPHISM_FAILED",err`Object;
            end try;
        end if;
    catch err
        print "  GEOMETRY_FAILED",err`Object;
    end try;
end for;
quit;
