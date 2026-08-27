// Lift the four rational points of Phi38 to full contact data:
// common c1 root of N0,N1; then a, c0, f; check disc(f) and the torsion classes.
Qdec := RationalField();
R<d,e,c1> := PolynomialRing(Qdec, 3);
F<dd,ee,cc1> := FieldOfFractions(R);

b  := (ee-dd)/2;
c2 := (ee+dd)/2;
r  := (dd*ee+7)/5;

Ra<A> := PolynomialRing(F);
Pa<xa> := PolynomialRing(FieldOfFractions(Ra));
h := 1 - (7/2)*xa + (FieldOfFractions(Ra)!A)*xa^2 + (FieldOfFractions(Ra)!b)*xa^3;
num := h^2 + (xa-1)^7;
f := &+[Coefficient(num,i+2)*xa^i : i in [0..5]];
e3 := Coefficient(f,3) - 2*(FieldOfFractions(Ra)!(cc1*c2)) - (FieldOfFractions(Ra)!(10*r^2));
e3n := Numerator(e3);
aval := -Coefficient(e3n,0)/Coefficient(e3n,1);

PF<xx> := PolynomialRing(F);
hF := 1 - (7/2)*xx + aval*xx^2 + F!b*xx^3;
numF := hF^2 + (xx-1)^7;
fF := &+[Coefficient(numF,i+2)*xx^i : i in [0..5]];
Rc<C0> := PolynomialRing(F);
e2 := Coefficient(fF,2) - (F!(cc1^2)) - 2*C0*(F!c2) - F!(-10*r^3);
c0val := -Coefficient(e2,0)/Coefficient(e2,1);
E1 := Coefficient(fF,1) - 2*c0val*(F!cc1) - F!(5*r^4);
E0 := Coefficient(fF,0) - c0val^2 - F!(-r^5);
N0 := Numerator(E1); N1 := Numerator(E0);

pts := [<-1,-3>, <1,-2>, <2,-1>, <3,1>];
QT<t> := PolynomialRing(Qdec);
for pt in pts do
  d0 := Qdec!pt[1]; e0 := Qdec!pt[2];
  printf "==== point (d,e) = (%o,%o), r = %o, b = %o, c2 = %o\n",
     d0, e0, (d0*e0+7)/5, (e0-d0)/2, (e0+d0)/2;
  n0 := UnivariatePolynomial(Evaluate(N0, [d0,e0,R.3]));
  n1 := UnivariatePolynomial(Evaluate(N1, [d0,e0,R.3]));
  printf "  deg N0=%o  deg N1=%o\n", Degree(n0), Degree(n1);
  if n0 eq 0 or n1 eq 0 then printf "  one N vanishes identically!\n"; end if;
  g := GCD(n0, n1);
  printf "  gcd degree %o : %o\n", Degree(g), g;
  rts := [rt[1] : rt in Roots(g)];
  printf "  rational c1 roots: %o\n", rts;
  for c1v in rts do
    av := Evaluate(Numerator(aval),[d0,e0,c1v])/Evaluate(Denominator(aval),[d0,e0,c1v]);
    c0v := F!c0val;
    c0n := Evaluate(Numerator(c0val),[d0,e0,c1v])/Evaluate(Denominator(c0val),[d0,e0,c1v]);
    bv := (e0-d0)/2; c2v := (e0+d0)/2; rv := (d0*e0+7)/5;
    hh := 1 - (7/2)*t + av*t^2 + bv*t^3;
    nn := hh^2 + (t-1)^7;
    fq := &+[Coefficient(nn,i+2)*t^i : i in [0..5]];
    qq := c0n + c1v*t + c2v*t^2;
    printf "  c1=%o: a=%o c0=%o\n", c1v, av, c0n;
    printf "    f = %o\n", fq;
    printf "    check f-q^2-(t-r)^5 = %o\n", fq - qq^2 - (t-rv)^5;
    printf "    disc(f) = %o\n", Discriminant(fq);
    printf "    h(1) = %o  (zero would kill the 7-class)\n", Evaluate(hh,1);
    printf "    q(r)^2 - f(r) = %o; f(r) = %o\n", Evaluate(qq,rv)^2-Evaluate(fq,rv), Evaluate(fq,rv);
  end for;
end for;
quit;
