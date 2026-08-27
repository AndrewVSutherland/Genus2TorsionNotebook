//////////////////////////////////////////////////////////////////////
// Modular relative-3 resolvent after pullback to the rational P8
// parameter u.  This is the finite-characteristic counterpart of
// contact6_m612_p8_relative3_generic.m and is deliberately bounded.
//
// Usage: magma -b p:=7 code/contact6_m612_p8_relative3_modular.m
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetSeed(1);
if not assigned p then p:=7;
elif Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
if not assigned do_genus then do_genus:=false;
elif Type(do_genus) eq MonStgElt then
    do_genus:=do_genus in {"true","True","1","yes"};
end if;
require IsPrime(p) and p notin {2,3,5}: "use p away from 2,3,5";

k:=GF(p); K<u>:=FunctionField(k);
kk:=func<n|K!(k!n)>;
t:=kk(4)*(u^2+u-kk(6))/(u^2+kk(6));
e:=-(kk(25)/kk(3))*t^2/(t^4-kk(25)*t^2+kk(1250)/kk(3));

R<M,U,v>:=PolynomialRing(K,3,"grevlex");
N:=(kk(3)*U+kk(6)*M)/kk(2);
Rc:=(kk(3)*U^2+kk(3)*v^2+(kk(2)/e-kk(15))*M-N^2)/kk(2);
F3:=kk(2)*v^3+kk(2)*N*Rc-U^3-kk(6)*U*v^2-kk(22)*M;
F2:=Rc^2+kk(2)*N*v^3-kk(3)*U^2*v^2-kk(3)*v^4
    -(kk(1)/e^2-kk(15))*M;
F1:=kk(2)*Rc*v^3-kk(3)*U*v^4-(kk(2)/e+kk(6))*M;
eqs:=[R!Numerator(f):f in [F3,F2,F1]];

print "CONTACT6_M612_P8_RELATIVE3_MODULAR";
print "PRIME",p,"E_DEGREES",Degree(Numerator(e)),Degree(Denominator(e));
Iraw:=ideal<R|eqs>;
t0:=Cputime(); I:=Saturation(Iraw,ideal<R|M>);
print "M_SATURATION_SECONDS",Cputime(t0),"DIMENSION",Dimension(I);
assert Dimension(I) eq 0;
A,q:=quo<R|I>; d:=Dimension(A);
print "AFFINE_LENGTH",d;
assert d eq 40;

t1:=Cputime(); IL:=ChangeOrder(I,"lex"); GL:=GroebnerBasis(IL);
print "LEX_SECONDS",Cputime(t1),"BASIS_SIZE",#GL;
print "LEX_BASIS_SHAPES",
      [<Degree(g,1),Degree(g,2),Degree(g,3),#Terms(g)>:g in GL];
univ:=[g:g in GL|Degree(g,1) eq 0 and Degree(g,2) eq 0 and
                  Degree(g,3) gt 0];
print "V_UNIVARIATE_COUNT",#univ;
if #univ gt 0 then
    rv:=UnivariatePolynomial(univ[#univ]); rv/:=LeadingCoefficient(rv);
    fac:=Factorization(rv);
    print "V_RESOLVENT_DEGREE",Degree(rv),"SQUAREFREE",
          GCD(rv,Derivative(rv)) eq 1,
          "FACTOR_DEGREES",[<Degree(fe[1]),fe[2]>:fe in fac];
    f1:=[fe[1]:fe in fac|Degree(fe[1]) eq 1][1];
    f12:=[fe[1]:fe in fac|Degree(fe[1]) eq 12][1];
    print "BUILTIN_LINEAR_FACTOR",f1;
    print "RELATIVE12_COEFFICIENT_SHAPES",
          [<Degree(Numerator(Coefficient(f12,i))),
             Degree(Denominator(Coefficient(f12,i)))>:i in [0..12]];
    assert IsIrreducible(f12);

    // The first lex basis polynomial is linear in M and contains only
    // M and v.  Recover M in the degree-12 support field and decide
    // whether adjoining L with L^2=M splits or doubles this component.
    gm:=[g:g in GL|Degree(g,1) eq 1 and Degree(g,2) eq 0][1];
    lcM:=Coefficient(gm,1,1);
    assert Degree(lcM,1) eq 0 and Degree(lcM,2) eq 0 and
           Degree(lcM,3) eq 0;
    hv:=Evaluate(gm,[Parent(gm)!0,Parent(gm)!0,Parent(gm).3]);
    mpoly:=-UnivariatePolynomial(hv)/(K!lcM);
    mpoly mod:= f12;
    F12<w>:=ext<K|f12>;
    m12:=Evaluate(mpoly,w);
    msquare,mroot:=IsSquare(m12);
    print "RELATIVE12_M_IS_SQUARE",msquare;
    if msquare then
        print "RELATIVE12_M_SQRT_SHAPE",
              <Degree(Numerator(mroot)),Degree(Denominator(mroot))>;
    end if;
    KT<T>:=PolynomialRing(K); KTV<VV>:=PolynomialRing(KT);
    fV:=&+[(KT!Coefficient(f12,i))*VV^i:i in [0..12]];
    mV:=&+[(KT!Coefficient(mpoly,i))*VV^i:i in [0..Degree(mpoly)]];
    p24:=Resultant(fV,T^2-mV); p24/:=LeadingCoefficient(p24);
    fac24:=Factorization(p24);
    print "RELATIVE24_L_RESOLVENT_DEGREE",Degree(p24),
          "FACTOR_DEGREES",[<Degree(fe[1]),fe[2]>:fe in fac24];
    assert Degree(p24) eq 24 and #fac24 eq 1 and
           Degree(fac24[1][1]) eq 24 and fac24[1][2] eq 1;
    if do_genus then
        print "RELATIVE12_GENUS_START";
        tg12:=Cputime(); g12:=Genus(F12);
        print "RELATIVE12_GENUS",g12,"SECONDS",Cputime(tg12);
        if not msquare then
            PF12<z>:=PolynomialRing(F12);
            F24<ell>:=ext<F12|z^2-m12>;
            print "RELATIVE24_GENUS_START";
            tg24:=Cputime(); g24:=Genus(F24);
            print "RELATIVE24_GENUS",g24,"SECONDS",Cputime(tg24);
        end if;
    end if;
end if;
print "CONTACT6_M612_P8_RELATIVE3_MODULAR_DONE";
quit;
