//////////////////////////////////////////////////////////////////////
// Geometry of the contact-7 family with two rational Weierstrass
// points, pulled back to the cubic-contact cover for rational 3-torsion.
//
// Signs are absorbed in the signed parameters s,t:
//
//   r=1-s^2, w=1-t^2, h(r)=s^7, h(w)=t^7.
//
// These two linear conditions determine a,b in
//
//   h=1-(7/2)x+a*x^2+b*x^3,
//   f=(h^2+(x-1)^7)/x^2.
//
// We translate r to zero via x=r+z.  Thus g(z)=f(r+z) has roots
// 0 and w-r.  Keeping the second root unscaled makes the fixed-s
// elimination substantially smaller than normalizing both roots to 0,1.
//
// A cubic contact is
//
//   H^2-q^3=M*g,  q=z^2+U*z+V, H=z^3+A*z^2+B*z+E.
//
// Since g(0)=0, the constant equation is E^2=V^3.  On V != 0 its
// normalization is V=j^2,E=j^3 (the other sign is j -> -j).  The
// coefficients z^5 and z^1 reconstruct A and B directly.  The remaining
// z^4,z^3,z^2 equations are substantially smaller than the usual
// high-to-low recursive equations.  Rational 3-torsion additionally
// requires M=L^2.
//
// Modes:
//   summary  : symbolic construction/degrees and exact assertions
//   factor   : fix s=s0 and factor the three cleared contact equations
//   fiber    : fix s=s0, eliminate the triangular variables, and study
//              quotient and square-cover projections
//
// Examples:
//   magma -b mode:=summary code/contact7_two_root_plus3_geometry.m
//   magma -b mode:=fiber s0:=2 code/contact7_two_root_plus3_geometry.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
if mode notin {"summary","factor","fiber"} then
    error "mode must be summary, factor, or fiber";
end if;
if not assigned s0 then s0 := 2;
elif Type(s0) eq MonStgElt then s0 := StringToInteger(s0); end if;
if not assigned MemGB then MemGB := 8;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();

function Primitive(S,f)
    if f eq 0 then return S!0; end if;
    p := S!f;
    den := LCM([Denominator(c) : c in Coefficients(p)]);
    ints := [Integers()!(den*c) : c in Coefficients(p)];
    cont := GCD(ints);
    if cont eq 0 then cont := 1; end if;
    return S!((Q!den/Q!Abs(cont))*p);
end function;

procedure Strip(~f,h,label)
    n := 0;
    while f ne 0 and IsDivisibleBy(f,h) do
        f := ExactQuotient(f,h); n +:= 1;
    end while;
    if n gt 0 then print "STRIPPED",label,n; end if;
end procedure;

if mode eq "summary" then
    R<s,t,U,j,M> := PolynomialRing(Q,5,"grevlex");
    K := FieldOfFractions(R);
    KZ<z> := PolynomialRing(K);

    r := 1-s^2; w := 1-t^2; d := w-r;
    Ar := s^7-1+(Q!7/2)*r;
    Aw := t^7-1+(Q!7/2)*w;
    b := K!((Aw/w^2-Ar/r^2)/(w-r));
    a := K!(Ar/r^2-b*r);
    h := 1-(Q!7/2)*z+a*z^2+b*z^3;
    f := ExactQuotient(h^2+(z-1)^7,z^2);
    g := Evaluate(f,r+z);
    c := [K!Coefficient(g,i) : i in [0..5]];
    assert c[1] eq 0;
    assert Evaluate(g,K!0) eq 0 and Evaluate(g,d) eq 0;

    A := (M*c[6]+3*U)/2;
    E := j^3;
    B := (3*U*j^4+M*c[2])/(2*j^3);
    q := z^2+U*z+j^2;
    H := z^3+A*z^2+B*z+E;
    ident := H^2-q^3-M*g;
    assert &and[Numerator(Coefficient(ident,i)) eq 0 : i in [0,1,5,6]];
    F4 := Primitive(R,Numerator(Coefficient(ident,4)));
    F3 := Primitive(R,Numerator(Coefficient(ident,3)));
    F2 := Primitive(R,Numerator(Coefficient(ident,2)));

    print "CONTACT7_TWO_ROOT_PLUS3_SUMMARY";
    print "signed signs absorbed: (eps*s)^7 and (delta*t)^7";
    print "r",r,"w",w,"d",d;
    print "a numerator/denominator degrees",
          <TotalDegree(Numerator(a)),TotalDegree(Denominator(a))>,
          "terms",<#Terms(Numerator(a)),#Terms(Denominator(a))>;
    print "b numerator/denominator degrees",
          <TotalDegree(Numerator(b)),TotalDegree(Denominator(b))>,
          "terms",<#Terms(Numerator(b)),#Terms(Denominator(b))>;
    print "g coefficient num/den degrees",
          [<TotalDegree(Numerator(ci)),TotalDegree(Denominator(ci))> : ci in c];
    print "contact equations degrees",[TotalDegree(F) : F in [F4,F3,F2]],
          "terms",[#Terms(F) : F in [F4,F3,F2]];
    print "BOUNDARY r*w*d*j*M*disc(q)*disc(g) != 0";
    quit;
end if;

//////////////////////////////////////////////////////////////////////
// Fixed-s fiber.  Work directly over Q[t,U,j,M].
//////////////////////////////////////////////////////////////////////

R<t,U,j,M> := PolynomialRing(Q,4,"grevlex");
K := FieldOfFractions(R);
KZ<z> := PolynomialRing(K);
s := Q!s0;
r := 1-s^2; w := 1-t^2; d := w-r;
if r eq 0 then error "s0=+-1 is outside this chart"; end if;
Ar := s^7-1+(Q!7/2)*r;
Aw := t^7-1+(Q!7/2)*w;
b := K!((Aw/w^2-Ar/r^2)/(w-r));
a := K!(Ar/r^2-b*r);
h := 1-(Q!7/2)*z+a*z^2+b*z^3;
f := ExactQuotient(h^2+(z-1)^7,z^2);
g := Evaluate(f,r+z);
c := [K!Coefficient(g,i) : i in [0..5]];
assert c[1] eq 0;
assert Evaluate(g,K!0) eq 0 and Evaluate(g,d) eq 0;

A := (M*c[6]+3*U)/2;
E := j^3;
B := (3*U*j^4+M*c[2])/(2*j^3);
q := z^2+U*z+j^2;
H := z^3+A*z^2+B*z+E;
ident := H^2-q^3-M*g;
assert &and[Numerator(Coefficient(ident,i)) eq 0 : i in [0,1,5,6]];

F4 := Primitive(R,Numerator(Coefficient(ident,4)));
F3 := Primitive(R,Numerator(Coefficient(ident,3)));
F2 := Primitive(R,Numerator(Coefficient(ident,2)));
if mode eq "factor" then
    print "CONTACT7_TWO_ROOT_PLUS3_FACTOR s0",s;
    for k in [1..3] do
        FF := [F4,F3,F2][k];
        print "EQUATION",k,<TotalDegree(FF),#Terms(FF)>;
        time facFF := Factorization(FF);
        print " FACTORS",[
            <fe[2],TotalDegree(fe[1]),Degree(fe[1],t),Degree(fe[1],U),
             Degree(fe[1],j),Degree(fe[1],M),#Terms(fe[1])> : fe in facFF];
    end for;
    quit;
end if;
Iraw := ideal<R | F4,F3,F2>;

print "CONTACT7_TWO_ROOT_PLUS3_FIBER s0",s;
print "equation degrees",[TotalDegree(F) : F in [F4,F3,F2]],
      "terms",[#Terms(F) : F in [F4,F3,F2]];
print "raw dimension",Dimension(Iraw),"basis",#Basis(Iraw),
      "basis degrees",[TotalDegree(F) : F in Basis(Iraw)];

// Saturate only simple chart boundaries before elimination.  Family
// discriminant and resultant(q,g) are audited on projected components.
discq := U^2-4*j^2;
print "SATURATION_J_BEGIN";
time Isat := Saturation(Iraw,ideal<R | j>);
print "SATURATION_M_BEGIN";
time Isat := Saturation(Isat,ideal<R | M>);
print "SATURATION_DONE dimension",Dimension(Isat),"basis",#Basis(Isat),
      "degrees",[TotalDegree(F) : F in Basis(Isat)],
      "terms",[#Terms(F) : F in Basis(Isat)];

// Lexicographic elimination to a plane curve in (t,j).  U,M are first.
RL<UL,ML,TL,JL> := PolynomialRing(Q,4,"lex");
mp := hom<R -> RL | TL,UL,JL,ML>;
Jlex := ideal<RL | [mp(F) : F in Basis(Isat)]>;
print "LEX_GB_BEGIN";
time gb := GroebnerBasis(Jlex);
print "LEX_GB_DONE size",#gb,"degrees",[TotalDegree(F) : F in gb],
      "terms",[#Terms(F) : F in gb];
elim := [F : F in gb | Degree(F,1) eq 0 and Degree(F,2) eq 0];
print "TJ_ELIM_GENERATORS",#elim;
if #elim eq 0 then print "NO_TJ_ELIMINATION"; quit; end if;

P<TP,JP> := PolynomialRing(Q,2);
toP := hom<RL -> P | P!0,P!0,TP,JP>;
core := P!toP(elim[1]);
for i in [2..#elim] do core := GCD(core,P!toP(elim[i])); end for;
Strip(~core,JP,"j=0");
Strip(~core,TP-1,"t=1");
Strip(~core,TP+1,"t=-1");
Strip(~core,TP-s,"t=s");
Strip(~core,TP+s,"t=-s");
print "TJ_CORE",<TotalDegree(core),Degree(core,TP),Degree(core,JP),#Terms(core)>;
fac := Factorization(core);
print "TJ_FACTORS",#fac;
for ie in [1..#fac] do
    pc := fac[ie][1];
    print " TJ_FACTOR",ie,"mult",fac[ie][2],
          <TotalDegree(pc),Degree(pc,TP),Degree(pc,JP),#Terms(pc)>;
    if TotalDegree(pc) le 8 then print pc; end if;
    try
        Caff := Curve(AffineSpace(P),pc);
        Cp := ProjectiveClosure(Caff);
        print "  projective degree",Degree(Cp),"arith_genus",Genus(Cp),
              "nonsingular",IsNonsingular(Cp);
        Cn,nu := Normalization(Cp);
        print "  normalization genus",Genus(Cn),"degree",Degree(Cn),
              "nonsingular",IsNonsingular(Cn);
        if Genus(Cn) ge 2 and Genus(Cn) le 20 then
            try
                aut,mpaut := AutomorphismGroup(Cn);
                print "  automorphism_order",#aut;
            catch err
                print "  automorphism_failed",err`Object;
            end try;
        end if;
    catch err
        print "  geometry_failed",err`Object;
    end try;
end for;

print "CONTACT7_TWO_ROOT_PLUS3_FIBER_DONE";
quit;
