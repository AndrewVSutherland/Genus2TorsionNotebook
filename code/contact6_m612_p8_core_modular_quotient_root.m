//////////////////////////////////////////////////////////////////////
// Modular generic profile of the L -> -L quotient of the P8-pulled-back
// independent-3 cover.  Here M=L^2 is retained as a coordinate.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned p then p:=13;
elif Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
assert IsPrime(p) and p notin {2,3,5};

k:=GF(p); K<u>:=FunctionField(k); kk:=func<n|K!(k!n)>;
t:=kk(4)*(u^2+u-kk(6))/(u^2+kk(6));
ee:=-(kk(25)/kk(3))*t^2/
    (t^4-kk(25)*t^2+kk(1250)/kk(3));
R<M,U,v>:=PolynomialRing(K,3,"grevlex");
N:=(kk(3)*U+kk(6)*M)/kk(2);
Rc:=(kk(3)*U^2+kk(3)*v^2+(kk(2)/ee-kk(15))*M-N^2)/kk(2);
F3:=kk(2)*v^3+kk(2)*N*Rc-U^3-kk(6)*U*v^2-kk(22)*M;
F2:=Rc^2+kk(2)*N*v^3-kk(3)*U^2*v^2-kk(3)*v^4
    -(kk(1)/ee^2-kk(15))*M;
F1:=kk(2)*Rc*v^3-kk(3)*U*v^4-(kk(2)/ee+kk(6))*M;
I:=ideal<R|[R!Numerator(f):f in [F3,F2,F1]]>;

print "CONTACT6_M612_P8_CORE_MODULAR_QUOTIENT_ROOT";
print "PRIME",p,"RAW_DIMENSION",Dimension(I);
t0:=Cputime(); I:=Saturation(I,ideal<R|M>);
I:=Saturation(I,ideal<R|v>);
print "M_V_SATURATION_SECONDS",Cputime(t0),"DIMENSION",Dimension(I);
disc:=U^2-kk(4)*v^2;
Apre,qpre:=quo<R|I>;
print "PRE_DISCRIMINANT_LENGTH",Dimension(Apre),
      "DISCRIMINANT_UNIT",IsUnit(qpre(disc));
J:=I+ideal<R|disc>;
if Dimension(J) eq 0 then
    AJ,qJ:=quo<R|J>;
    print "DISCRIMINANT_BOUNDARY_LENGTH",Dimension(AJ);
end if;

// The open quotient is expected to have length 39.  Full saturation is
// cheap in the M-quotient model.
t1:=Cputime(); Iopen:=Saturation(I,ideal<R|disc>);
print "DISCRIMINANT_SATURATION_SECONDS",Cputime(t1),
      "OPEN_DIMENSION",Dimension(Iopen);
A,mp:=quo<R|Iopen>; d:=Dimension(A);
print "OPEN_QUOTIENT_LENGTH",d;

for c in [1..6] do
    zpoly:=M+kk(c)*U+kk(c*c+1)*v;
    f:=MinimalPolynomial(mp(zpoly)); fac:=Factorization(f);
    print "LINEAR_FORM",c,"MINPOLY_DEGREE",Degree(f),
          "SQUAREFREE",GCD(f,Derivative(f)) eq 1,
          "FACTOR_DEGREES_MULTIPLICITIES",
          [<Degree(fe[1]),fe[2]>:fe in fac];
    if Degree(f) eq d and GCD(f,Derivative(f)) eq 1 then break; end if;
end for;
print "CONTACT6_M612_P8_CORE_MODULAR_QUOTIENT_ROOT_DONE";
quit;
