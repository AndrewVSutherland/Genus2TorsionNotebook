// Reconstruct the residual plane curve Phi38(d,e) of the Z/35 nondegenerate
// contact-7 + contact-5 lane (notes/agent_Z35_next_route.md), then analyze it.
// h = 1-(7/2)x+a*x^2+b*x^3, f = (h^2+(x-1)^7)/x^2, impose f - q^2 = (x-r)^5,
// q = c0+c1*x+c2*x^2; d=c2-b, e=c2+b, r=(d*e+7)/5.
Qdec := RationalField();
R<d,e,c1> := PolynomialRing(Qdec, 3);
F<dd,ee,cc1> := FieldOfFractions(R);
P<x> := PolynomialRing(F);

b  := (ee-dd)/2;
c2 := (ee+dd)/2;
r  := (dd*ee+7)/5;

// f depends on a; keep a symbolic via a second layer: solve E3 for a first.
// Work in F(a) -- use another poly ring.
Ra<A> := PolynomialRing(F);
Pa<xa> := PolynomialRing(FieldOfFractions(Ra));

h := 1 - (7/2)*xa + (FieldOfFractions(Ra)!A)*xa^2 + (FieldOfFractions(Ra)!b)*xa^3;
num := h^2 + (xa-1)^7;
// divide by x^2
assert Coefficient(num,0) eq 0 and Coefficient(num,1) eq 0;
f := &+[Coefficient(num,i+2)*xa^i : i in [0..5]];

// E3: coeff of x^3 of f - q^2 - (x-r)^5, q = c0+c1 x+c2 x^2.
// coeff3(q^2) = 2*c1*c2 (no c0? q^2 coeff3 = 2 c1 c2). coeff3((x-r)^5)=10 r^2.
e3 := Coefficient(f,3) - 2*(FieldOfFractions(Ra)!(cc1*c2)) - (FieldOfFractions(Ra)!(10*r^2));
// e3 is linear in A: e3 = e31*A + e30
e3n := Numerator(e3);  // in Ra = F[A]
assert Degree(e3n) eq 1;
aval := -Coefficient(e3n,0)/Coefficient(e3n,1);   // element of F

// now specialize f at A = aval
PF<xx> := PolynomialRing(F);
hF := 1 - (7/2)*xx + aval*xx^2 + F!b*xx^3;
numF := hF^2 + (xx-1)^7;
assert Coefficient(numF,0) eq 0 and Coefficient(numF,1) eq 0;
fF := &+[Coefficient(numF,i+2)*xx^i : i in [0..5]];

// E2: coeff2( f - q^2 - (x-r)^5 ) = 0, q^2 coeff2 = c1^2 + 2 c0 c2 -> solve c0
Rc<C0> := PolynomialRing(F);
e2 := Coefficient(fF,2) - (F!(cc1^2)) - 2*C0*(F!c2) - F!(-10*r^3);
assert Degree(e2) eq 1;
c0val := -Coefficient(e2,0)/Coefficient(e2,1);

// E1: coeff1: q^2 coeff1 = 2 c0 c1 ; (x-r)^5 coeff1 = 5 r^4
E1 := Coefficient(fF,1) - 2*c0val*(F!cc1) - F!(5*r^4);
// E0: coeff0: q^2 coeff0 = c0^2 ; (x-r)^5 coeff0 = -r^5
E0 := Coefficient(fF,0) - c0val^2 - F!(-r^5);

N0 := Numerator(E1); N1 := Numerator(E0);
// clear rational content
N0 := N0 * LCM([Denominator(c) : c in Coefficients(N0)]);
N0 := N0 / GCD([Numerator(c) : c in Coefficients(N0)]);
N1 := N1 * LCM([Denominator(c) : c in Coefficients(N1)]);
N1 := N1 / GCD([Numerator(c) : c in Coefficients(N1)]);
N0 := R!N0; N1 := R!N1;
printf "N0 deg/terms %o %o\n", TotalDegree(N0), #Terms(N0);
printf "N1 deg/terms %o %o\n", TotalDegree(N1), #Terms(N1);

res := Resultant(N0, N1, c1);
fac := Factorization(res);
printf "resultant factors (deg, mult):\n";
Phi := R!1;
for t in fac do
  printf "  deg %o  mult %o  terms %o\n", TotalDegree(t[1]), t[2], #Terms(t[1]);
  if TotalDegree(t[1]) gt 4 then Phi := t[1]; end if;
end for;
printf "Phi total degree %o, terms %o\n", TotalDegree(Phi), #Terms(Phi);
S<D,E> := PolynomialRing(Qdec,2);
PhiDE := Evaluate(Phi, [S.1, S.2, 0]);
assert TotalDegree(PhiDE) eq TotalDegree(Phi); // c1 eliminated
// irreducibility over Q
facQ := Factorization(PhiDE);
printf "Phi38 over Q: %o factors, degrees %o\n", #facQ, [TotalDegree(t[1]) : t in facQ];

// genus over a large prime first (fast certificate of the Q-genus generically)
for p in [10007, 32003] do
  Sp<Dp,Ep> := PolynomialRing(GF(p),2);
  Phip := Sp!0;
  for t in Terms(PhiDE) do Phip +:= Sp!t; end for;
  facp := Factorization(Phip);
  printf "mod %o: %o factors, degrees %o\n", p, #facp, [<TotalDegree(t[1]),t[2]> : t in facp];
  for t in facp do
    Cp := Curve(AffineSpace(Sp), t[1]);
    gp := Genus(ProjectiveClosure(Cp));
    printf "  mod %o component deg %o: genus %o\n", p, TotalDegree(t[1]), gp;
  end for;
end for;

// save Phi38 for reuse (repo-relative; run from the repo root, or override
// with magma -b outfile:=/path/to/file)
if not assigned outfile then outfile := "data/claude_prod_04_35_phi38_poly.txt"; end if;
PrintFile(outfile, Sprint(PhiDE) : Overwrite := true);
quit;
