// The q0=c^2, q1=1, q2=-1/4 family: verify [15] structure, simplicity, and set up the [30]/[2,30] cut.
SetColumns(0);
Qx<x> := PolynomialRing(Rationals());
function buildF(q0,q1)
  q2 := -1/4;
  Q := q2*x^2 + q1*x + q0;
  Qp := Q - x;
  h := Qp - x*Q;
  return h^2 + 4*Q^2*Qp;
end function;

// (i) direct order check of D = [(2c,4c^2)+(-2c,4c^2)-K] for several c
print "=== (i) order of D and full torsion for c = 1,2,3,7,1/2 ===";
for c in [1,2,3,7,1/2] do
  F := buildF(c^2, 1);
  den := LCM([Denominator(co) : co in Coefficients(F)]);
  Fi := F*den^2;  // integral-ish
  C := HyperellipticCurve(Fi);
  J := Jacobian(C);
  // points (2c, 4c^2) scale: y^2 = F -> (den*y)^2 = Fi: y-coord scales by den
  D := elt<J | x^2 - 4*c^2, Qx!(4*c^2*den)>;
  printf "c=%o: Order(D)=%o, torsion=%o\n", c, Order(D), Invariants(TorsionSubgroup(J));
end for;

// (ii) symbolic contact identity over Q(c): h3^2 - F = A3^2 (x^2-4c^2)^3
print "=== (ii) symbolic contact over Q(c) ===";
Kc<c> := FunctionField(Rationals());
P4<A0,A1,A2,A3> := PolynomialRing(Kc, 4);
P4x<X> := PolynomialRing(P4);
q2 := -1/4;
Qq := q2*X^2 + X + c^2;
Qp := Qq - X;
hh := Qp - X*Qq;
FF := hh^2 + 4*Qq^2*Qp;
H3 := A3*X^3 + A2*X^2 + A1*X + A0;
EE := H3^2 - FF - A3^2*(X^2 - 4*c^2)^3;
eqs := [Coefficient(EE, i) : i in [0..6]];
I := ideal<P4 | eqs>;
GB := GroebnerBasis(I);
print "Groebner basis:", GB;
print "dim:", Dimension(I);
if Dimension(I) eq 0 then
  V := Variety(I);
  print "rational solutions (A0,A1,A2,A3):", V;
end if;

// (iii) simplicity certificates: L-poly at good primes, irreducible deg 4 AND 12th-power transform irreducible
print "=== (iii) simplicity certs ===";
function simplecert(F, pmax)
  den := LCM([Denominator(co) : co in Coefficients(F)]);
  Fi := F*den^2;
  co := Coefficients(Fi);
  g := GCD([Integers()!cc : cc in co]);
  sq := 1;
  for pf in Factorization(g) do sq *:= pf[1]^(pf[2] div 2); end for;
  Fi := Fi div sq^2;
  C := HyperellipticCurve(Fi);
  disc := Integers()!Discriminant(C);
  certs := [];
  for p in PrimesInInterval(5, pmax) do
    if disc mod p eq 0 then continue; end if;
    Cp := ChangeRing(C, GF(p));
    chi := Numerator(ZetaFunction(Cp));  // reciprocal L-poly
    chi := Qx!chi;
    if not IsIrreducible(chi) then continue; end if;
    // 12th power transform: charpoly of Frob^12
    K := NumberField(chi);
    al := K.1;
    chi12 := MinimalPolynomial(al^12);
    if Degree(chi12) eq 4 and IsIrreducible(chi12) then
      Append(~certs, p);
      if #certs ge 3 then break; end if;
    end if;
  end for;
  return certs;
end function;

for pr in [<1,1>, <4,1>, <9,1>, <2,1>, <-36, 11/2>, <-2, 16/7>, <49/4, 1>] do
  F := buildF(pr[1], pr[2]);
  certs := simplecert(F, 200);
  printf "q0=%o q1=%o: simplicity cert primes %o %o\n", pr[1], pr[2], certs,
    (#certs ge 2 select "SIMPLE (multi-prime)" else "NO CERT");
end for;

// (iv) the [30]/[2,30] cut: rational-root curve F(r; c)=0
print "=== (iv) rational-root curve F(r,c)=0 ===";
A2<r,cc> := AffineSpace(Rationals(), 2);
q2 := -1/4;
Qr := q2*r^2 + r + cc^2;
Qpr := Qr - r;
hr := Qpr - r*Qr;
Fr := hr^2 + 4*Qr^2*Qpr;
Croot := Curve(A2, Numerator(Fr));
print "is irreducible:", IsIrreducible(Croot);
print "genus:", Genus(Croot);
print "singular points:", SingularPoints(Croot);
pts := PointSearch(Croot, 1000);
print "point search to height 1000:", pts;
quit;
