//////////////////////////////////////////////////////////////////////
// Modular generic component profile for the independent-[3,3] cover
// pulled back to the rational P8 R3-halving family.
//
// Work over F_p(u), put b=0 and a=1/e(u), eliminate N and R from the
// cubic-contact equations, and retain L explicitly through M=L^2.
// Saturation removes the contact-boundary factors
//
//   L*v*(U^2-4*v^2).
//
// A minimal polynomial of a generic linear form in the resulting finite
// affine algebra provides a bounded component-degree diagnostic.
//
// Usage:
//   magma -b p:=7 code/contact6_m612_p8_core_modular_profile_root.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned p then p:=7;
elif Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
if not assigned do_full_saturation then do_full_saturation:=false;
elif Type(do_full_saturation) eq MonStgElt then
    do_full_saturation:=do_full_saturation in {"true","True","1","yes"};
end if;
if not assigned do_decompose then do_decompose:=true;
elif Type(do_decompose) eq MonStgElt then
    do_decompose:=do_decompose in {"true","True","1","yes"};
end if;
if not assigned do_minpoly then do_minpoly:=false;
elif Type(do_minpoly) eq MonStgElt then
    do_minpoly:=do_minpoly in {"true","True","1","yes"};
end if;
assert IsPrime(p) and p notin {2,3,5};

k:=GF(p); K<u>:=FunctionField(k);
kk:=func<n|K!(k!n)>;
t:=kk(4)*(u^2+u-kk(6))/(u^2+kk(6));
ee:=-(kk(25)/kk(3))*t^2/
    (t^4-kk(25)*t^2+kk(1250)/kk(3));
assert ee ne 0;

R<L,U,v>:=PolynomialRing(K,3,"grevlex");
M:=L^2;
N:=(kk(3)*U+kk(6)*M)/kk(2);
Rc:=(kk(3)*U^2+kk(3)*v^2+(kk(2)/ee-kk(15))*M-N^2)/kk(2);
F3:=kk(2)*v^3+kk(2)*N*Rc-U^3-kk(6)*U*v^2-kk(22)*M;
F2:=Rc^2+kk(2)*N*v^3-kk(3)*U^2*v^2-kk(3)*v^4
    -(kk(1)/ee^2-kk(15))*M;
F1:=kk(2)*Rc*v^3-kk(3)*U*v^4-(kk(2)/ee+kk(6))*M;
eqs:=[R!Numerator(f):f in [F3,F2,F1]];

print "CONTACT6_M612_P8_CORE_MODULAR_PROFILE_ROOT";
print "PRIME",p;
print "E_NUM_DEN_DEGREES",Degree(Numerator(ee)),Degree(Denominator(ee));
print "EQUATION_SHAPES",[<TotalDegree(f),#Terms(f)>:f in eqs];

Iraw:=ideal<R|eqs>;
print "RAW_DIMENSION",Dimension(Iraw);
print "SATURATION_START";
t0:=Cputime();
I:=Saturation(Iraw,ideal<R|L>);
print " L_SATURATION_SECONDS",Cputime(t0),"DIMENSION",Dimension(I);
t1:=Cputime();
I:=Saturation(I,ideal<R|v>);
print " V_SATURATION_SECONDS",Cputime(t1),"DIMENSION",Dimension(I);
disc:=U^2-kk(4)*v^2;
Apre,qpre:=quo<R|I>;
disc_is_unit:=IsUnit(qpre(disc));
print " PRE_DISCRIMINANT_AFFINE_ALGEBRA_DIMENSION",Dimension(Apre);
print " DISCRIMINANT_IS_UNIT_AFTER_L_V_SATURATION",disc_is_unit;
Jdisc:=I+ideal<R|disc>;
dimdisc:=Dimension(Jdisc);
print " DISCRIMINANT_BOUNDARY_IDEAL_DIMENSION",dimdisc;
if dimdisc eq 0 then
    Adisc,qdisc:=quo<R|Jdisc>;
    print " DISCRIMINANT_BOUNDARY_LENGTH",Dimension(Adisc);
end if;
if not disc_is_unit and do_full_saturation then
    t2:=Cputime();
    I:=Saturation(I,ideal<R|disc>);
    print " DISCRIMINANT_SATURATION_SECONDS",Cputime(t2),
          "DIMENSION",Dimension(I);
end if;
print "SATURATION_SECONDS",Cputime(t0),"DIMENSION",Dimension(I);
assert Dimension(I) eq 0;

A,qmap:=quo<R|I>;
d:=Dimension(A);
print "AFFINE_ALGEBRA_DIMENSION",d;
assert d le 120;

if do_decompose then
    print "PRIMARY_DECOMPOSITION_START";
    td:=Cputime();
    comps:=PrimaryDecomposition(I);
    print "PRIMARY_DECOMPOSITION_SECONDS",Cputime(td),
          "COMPONENT_COUNT",#comps;
    for j in [1..#comps] do
        C:=comps[j]; AC,qC:=quo<R|C>;
        print " COMPONENT",j,"LENGTH",Dimension(AC),
              "PRIME",IsPrime(C),
              "DISCRIMINANT_UNIT",IsUnit(qC(disc)),
              "BASIS_SHAPES",
              [<TotalDegree(g),#Terms(g)>:g in Basis(C)];
    end for;
end if;

if do_minpoly then for c in [1..4] do
    zpoly:=L+kk(c)*U+kk(c*c+1)*v;
    z:=qmap(zpoly);
    mp:=MinimalPolynomial(z);
    fac:=Factorization(mp);
    print "LINEAR_FORM",c,"MINPOLY_DEGREE",Degree(mp),
          "SQUAREFREE",GCD(mp,Derivative(mp)) eq 1,
          "FACTOR_DEGREES_MULTIPLICITIES",
          [<Degree(fe[1]),fe[2]>:fe in fac];
    if dimdisc eq 0 then
        mpdisc:=MinimalPolynomial(qdisc(zpoly));
        print " BOUNDARY_MINPOLY_DEGREE",Degree(mpdisc),
              "BOUNDARY_FACTOR_DEGREES_MULTIPLICITIES",
              [<Degree(fe[1]),fe[2]>:fe in Factorization(mpdisc)];
    end if;
    if Degree(mp) eq d and GCD(mp,Derivative(mp)) eq 1 then
        print "PRIMITIVE_LINEAR_FORM",c;
        break;
    end if;
end for;
end if;

print "CONTACT6_M612_P8_CORE_MODULAR_PROFILE_ROOT_DONE";
quit;
