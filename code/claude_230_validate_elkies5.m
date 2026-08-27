// validate230.m — calibration for target (2,30)
// (a) reproduce repo fact: contact5xcontact6 member (u,s)=(125,5415) has J(Q)_tors = [30], simple
// (b) reproduce literature fact: Elkies (2,2,2,10) curve torsion
// (c) validate Elkies universal 5-torsion threefold: random (q0,q1,q2) give order-5 class;
//     involution (q0,q1,q2)->(q2,1-q1,q0) preserves the curve (G2 invariants)
QQ := Rationals();
P<x> := PolynomialRing(QQ);

print "=== (a) repo [30] member (u,s)=(125,5415) ===";
u := QQ!125; s := QQ!5415;
CC := (u^2+1)/(2*u); cc := (u^2-1)/(2*u);
q := (15*u^5 + 90*u^4 + 20*u^3*s - 6*u^2*s^2 + 231*u^3 + 2*u^2*s - 15*u*s^2
      + 90*u^2 - 20*u*s + 15*u - 2*s)
     / (u^6 + 6*u^4*s - 2*u^4 + 15*u^3*s - u*s^3 + u^2);
A := (s+q)/2; e := (s-q)/2; B := (15-s*q)/2; d := (B*CC+3)/cc;
K := -2*A - 6;
h6 := x^3 + A*x^2 + B*x + CC;
h5 := e*x^2 + d*x + cc;
f := h6^2 - (x-1)^6;
print "identity f = h5^2 - K*x^5:", f eq h5^2 - K*x^5;
den := LCM([Denominator(co) : co in Coefficients(f)]);
fZ := P!(den^2*f);   // y -> den*y : isomorphic Jacobian, integral model
C30 := HyperellipticCurve(fZ);
J30 := Jacobian(C30);
D6 := J30![x-1, P!(den*Evaluate(h6,1))];
D5 := J30![x, P!(den*cc)];
printf "Order(D5)=%o Order(D6)=%o Order(D5+D6)=%o\n", Order(D5), Order(D6), Order(D5+D6);
T30 := TorsionSubgroup(J30);
print "TorsionSubgroup invariants:", Invariants(T30);
for pr in [7,11,13,17,19,23,29,31,37] do
  if Valuation(Integers()!Discriminant(fZ), pr) eq 0 and Valuation(den, pr) eq 0 then
    Pp := PolynomialRing(GF(pr));
    Cp := HyperellipticCurve(Pp!fZ);
    Lp := LPolynomial(Cp);
    chi := P!Reverse(Coefficients(Lp));
    if Degree(chi) eq 4 and IsIrreducible(chi) then
      Pxy<X,T> := PolynomialRing(QQ, 2);
      chi12 := UnivariatePolynomial(Resultant(Evaluate(chi, X), T - X^12, X));
      if Degree(chi12) eq 4 and IsIrreducible(chi12) then
        printf "simplicity cert: p=%o Lp=%o IRREDUCIBLE + 12th-transform irreducible\n", pr, Lp;
        break;
      end if;
    end if;
  end if;
end for;

print "=== (b) Elkies (2,2,2,10) curve ===";
f80 := x*(x+1)*(x-1)*(3*x-7)*(8*x-13)*(24*x+25);
J80 := Jacobian(HyperellipticCurve(f80));
print "torsion invariants:", Invariants(TorsionSubgroup(J80));

print "=== (c) Elkies universal 5-torsion curve, random triples ===";
for tr in [[QQ|1,2,3],[QQ|2,-1,1],[QQ|-1,3,2],[QQ|1/2,3,-2],[QQ|5,7,-3]] do
  q0 := tr[1]; q1 := tr[2]; q2 := tr[3];
  Q := q2*x^2 + q1*x + q0;
  Qp := Q - x;
  h := Qp - x*Q;      // L'Q' - LQ with L=x, L'=1
  fE := Q^2*Qp;
  F := h^2 + 4*fE;    // completed square
  C0 := HyperellipticCurve(fE, h);
  Csq := HyperellipticCurve(F);
  JE := Jacobian(Csq);
  b := h mod Q;
  D := JE![Q, b];
  // involution triple
  Q2 := q0*x^2 + (1-q1)*x + q2;
  h2 := (Q2 - x) - x*Q2; F2 := h2^2 + 4*Q2^2*(Q2-x);
  printf "triple %o: genus=%o Order(D)=%o involutionG2match=%o\n",
    tr, Genus(C0), Order(D),
    G2Invariants(Csq) eq G2Invariants(HyperellipticCurve(F2));
end for;
quit;
