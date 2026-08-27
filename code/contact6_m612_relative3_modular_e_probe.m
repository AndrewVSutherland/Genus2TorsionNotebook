//////////////////////////////////////////////////////////////////////
// Fast modular probe of the b=0 relative-3 support cover over F_p(e).
// This determines coefficient degree bounds for reconstructing the exact
// degree-12 orthogonal factor over Q(e) from fixed rational fibers.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned p then p:=7;
elif Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
assert IsPrime(p) and p notin {2,3,5};

k:=GF(p); K<e>:=FunctionField(k); kk:=func<n|K!(k!n)>;
R<M,U,v>:=PolynomialRing(K,3,"grevlex");
N:=(kk(3)*U+kk(6)*M)/kk(2);
Rc:=(kk(3)*U^2+kk(3)*v^2+(kk(2)/e-kk(15))*M-N^2)/kk(2);
F3:=kk(2)*v^3+kk(2)*N*Rc-U^3-kk(6)*U*v^2-kk(22)*M;
F2:=Rc^2+kk(2)*N*v^3-kk(3)*U^2*v^2-kk(3)*v^4
    -(kk(1)/e^2-kk(15))*M;
F1:=kk(2)*Rc*v^3-kk(3)*U*v^4-(kk(2)/e+kk(6))*M;
I:=ideal<R|[R!Numerator(f):f in [F3,F2,F1]]>;
I:=Saturation(I,ideal<R|M>);
assert Dimension(I) eq 0 and Dimension(quo<R|I>) eq 40;

print "CONTACT6_M612_RELATIVE3_MODULAR_E_PROBE","PRIME",p;
t0:=Cputime(); IL:=ChangeOrder(I,"lex"); GL:=GroebnerBasis(IL);
print "LEX_SECONDS",Cputime(t0),"BASIS_SIZE",#GL;
univ:=[g:g in GL|Degree(g,1) eq 0 and Degree(g,2) eq 0 and
                  Degree(g,3) gt 0];
assert #univ eq 1;
rv:=UnivariatePolynomial(univ[1]); rv/:=LeadingCoefficient(rv);
fac:=Factorization(rv);
print "FACTOR_DEGREES",[<Degree(fe[1]),fe[2]>:fe in fac];
print "COEFFICIENT_SHAPES_BY_FACTOR",
      [<Degree(fe[1]),
         [<Degree(Numerator(Coefficient(fe[1],i))),
            Degree(Denominator(Coefficient(fe[1],i)))>
          :i in [0..Degree(fe[1])]]>:fe in fac];
f12:=[fe[1]:fe in fac|Degree(fe[1]) eq 12][1];
gm:=[g:g in GL|Degree(g,1) eq 1 and Degree(g,2) eq 0][1];
gu:=[g:g in GL|Degree(g,1) eq 0 and Degree(g,2) eq 1][1];
cm:=K!Coefficient(gm,1,1); cu:=K!Coefficient(gu,2,1);
mrec:=(-UnivariatePolynomial(Evaluate(gm,[R!0,R!0,v]))/cm) mod f12;
urec:=(-UnivariatePolynomial(Evaluate(gu,[R!0,R!0,v]))/cu) mod f12;
print "M_RECOVERY_COEFFICIENT_SHAPES",
      [<Degree(Numerator(Coefficient(mrec,i))),
         Degree(Denominator(Coefficient(mrec,i)))>:i in [0..11]];
print "U_RECOVERY_COEFFICIENT_SHAPES",
      [<Degree(Numerator(Coefficient(urec,i))),
         Degree(Denominator(Coefficient(urec,i)))>:i in [0..11]];
KT<T>:=PolynomialRing(K); KTV<VV>:=PolynomialRing(KT);
fV:=&+[(KT!Coefficient(f12,i))*VV^i:i in [0..12]];
mV:=&+[(KT!Coefficient(mrec,i))*VV^i:i in [0..Degree(mrec)]];
pL:=Resultant(fV,T^2-mV); pL/:=LeadingCoefficient(pL);
print "L_RESOLVENT_DEGREE",Degree(pL),"L_RESOLVENT_FACTORS",
      [<Degree(fe[1]),fe[2]>:fe in Factorization(pL)];
print "CONTACT6_M612_RELATIVE3_MODULAR_E_PROBE_DONE";
quit;
