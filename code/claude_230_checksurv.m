// checksurv.m — exact torsion + simplicity certificate for sieve survivors
QQ := Rationals();
P<x> := PolynomialRing(QQ);
trips := [
 [QQ|-16/1, 4/1, -1/4],
 [QQ|-12/1, 11/3, -1/4],
 [QQ|-9/1, 3/1, -1/4],
 [QQ|-9/1, 27/4, -1/4],
 [QQ|-4/1, 3/1, -1/4],
 [QQ|-1/1, -1/1, -1/4],
 [QQ|13/1, -39/4, -1/4],
 [QQ|60/1, -45/1, -1/4],
 [QQ|1/2, 31/4, -1/4],
 [QQ|-23/3, 23/4, -1/4],
 [QQ|-9/4, -1/2, -1/4],
 [QQ|-1/4, 1/2, -1/4],
 [QQ|-1/4, -24/7, -1/4],
 [QQ|-1/4, -6/7, -1/4],
 [QQ|-1/4, 13/7, -1/4],
 [QQ|-1/4, 31/7, -1/4],
 [QQ|-1/4, -31/8, -1/4],
 [QQ|-1/4, 39/8, -1/4],
 [QQ|-38/5, 25/8, -1/4],
 [QQ|32/5, -24/5, -1/4],
 [QQ|34/5, 45/4, -1/4],
 [QQ|36/5, -27/5, -1/4],
 [QQ|-49/6, 49/8, -1/4],
 [QQ|-5/6, -1/5, -1/4],
 [QQ|19/6, -19/8, -1/4],
 [QQ|12/7, -56/1, -1/4],
 [QQ|12/7, 19/3, -1/4],
 [QQ|39/7, 37/8, -1/4] ];
for tr in trips do
  q0 := tr[1]; q1 := tr[2]; q2 := tr[3];
  Q := q2*x^2 + q1*x + q0;
  Qp := Q - x; h := Qp - x*Q; F := h^2 + 4*Q^2*Qp;
  if Degree(F) lt 5 or not IsSeparable(F) then
    printf "q=%o : DEGENERATE\n", tr; continue;
  end if;
  den := LCM([Denominator(co) : co in Coefficients(F)]);
  fZ := P!(den^2*F);
  C := HyperellipticCurve(fZ);
  J := Jacobian(C);
  T := TorsionSubgroup(J);
  inv := Invariants(T);
  // simplicity certificate
  cert := 0;
  DD := Integers()!Discriminant(fZ);
  for pr in PrimesInInterval(7,97) do
    if DD mod pr ne 0 then
      Pp := PolynomialRing(GF(pr));
      Lp := LPolynomial(HyperellipticCurve(Pp!fZ));
      chi := P!Reverse(Coefficients(Lp));
      if Degree(chi) eq 4 and IsIrreducible(chi) then
        Pxy<X,T2> := PolynomialRing(QQ, 2);
        chi12 := UnivariatePolynomial(Resultant(Evaluate(chi, X), T2 - X^12, X));
        if Degree(chi12) eq 4 and IsIrreducible(chi12) then cert := pr; break; end if;
      end if;
    end if;
  end for;
  fac := [<Degree(g[1]), g[2]> : g in Factorization(fZ)];
  printf "q=%o : tors=%o cert_p=%o facF=%o\n", tr, inv, cert, fac;
end for;
quit;
