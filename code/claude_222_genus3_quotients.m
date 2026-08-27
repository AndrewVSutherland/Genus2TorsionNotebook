// claude_222_genus3_quotients.m — quotient each of the two genus-3 blockers
// by its unique involution; identify the quotient (genus 0/1/2) and, when
// elliptic or genus 2, set up the rank computation that closes X(Q).
SetColumns(0);
SetMemoryLimit(24*10^9);
Q := Rationals();
A2<s,x> := AffineSpace(Q, 2);

c1 := (s-1)^2*x^3 + (2*s^4-4*s^3+5*s^2-4*s+1)*x^2
    + (s^6-2*s^5+3*s^4-2*s^3+2*s^2-2*s+1)*x - s^2*(s^2-s+1)^2;
X1 := Curve(A2, c1);

FF<T> := FunctionField(Q);
uu := (-11/100*T^3 + 919/50*T^2 - 8017/25*T + 40154/25)/(T^3 - 23*T^2 + 168*T - 396);
vv := (2401/1600*T^4 - 12593/200*T^3 + 205647/200*T^2 - 385757/50*T + 2253001/100)/(T^4 - 34*T^3 + 421*T^2 - 2244*T + 4356);
tt := (-21609/40000*T^6 + 435561/10000*T^5 - 2897123/2000*T^4 + 6390763/250*T^3 - 126932427/500*T^2 + 844854361/625*T - 1894773841/625)/(T^6 - 46*T^5 + 865*T^4 - 8520*T^3 + 46440*T^2 - 133056*T + 156816);
PxF<X> := PolynomialRing(FF);
FlT := X^6+2*X^5+(2*tt+3)*X^4+2*X^3+(tt^2+1)*X^2+2*tt*(1-tt)*X+tt^2;
Q4 := FlT div (X^2 + FF!uu*X + FF!vv);
cfs := Coefficients(Q4);
den := LCM([Denominator(cq) : cq in cfs]);
R2<t2,r2> := PolynomialRing(Q, 2);
num := &+[ Evaluate(Numerator(cfs[i]*den), t2) * r2^(i-1) : i in [1..#cfs] ];
A2b<tb,rb> := AffineSpace(Q,2);
comp9 := 0;
for fa in Factorization(Evaluate(num, [tb, rb])) do
  if TotalDegree(fa[1]) ge 8 then comp9 := fa[1]; end if;
end for;
X2 := Curve(A2b, comp9);

for pair in [* <X1,"X1">, <X2,"X2"> *] do
  C := pair[1]; nm := pair[2];
  printf "==== %o ====\n", nm;
  Cp := ProjectiveClosure(C);
  G := AutomorphismGroup(Cp);
  printf "%o: #Aut = %o\n", nm, #G;
  if #G lt 2 then continue; end if;
  CQ, qm := CurveQuotient(G);
  gq := Genus(CQ);
  printf "%o: quotient genus = %o\n", nm, gq;
  if gq eq 1 then
    pts := PointSearch(CQ, 10^5);
    printf "  quotient points found: %o\n", #pts;
    if #pts gt 0 then
      E, _ := EllipticCurve(CQ, pts[1]);
      Em, _ := MinimalModel(E);
      printf "  ELLIPTIC QUOTIENT (minimal): %o\n", aInvariants(Em);
      printf "  conductor = %o\n", Conductor(Em);
      rl, ru := RankBounds(Em);
      printf "  rank bounds = [%o, %o]\n", rl, ru;
      printf "  torsion = %o\n", Invariants(TorsionSubgroup(Em));
      if ru eq 0 then
        printf "  RANK 0 => %o(Q) maps to a FINITE set: enumerate torsion fibers to close %o\n", nm, nm;
      end if;
    end if;
  elif gq eq 2 then
    ok, H2 := IsHyperelliptic(CQ);
    if ok then
      printf "  GENUS-2 QUOTIENT: %o\n", H2;
      J := Jacobian(H2);
      printf "  rank bound = %o\n", RankBound(J);
    else
      print "  genus-2 quotient not hyperelliptic-recognized";
    end if;
  elif gq eq 0 then
    printf "  quotient rational => %o is hyperelliptic; computing model\n", nm;
    ok, HC := IsHyperelliptic(Cp);
    if ok then printf "  hyperelliptic model: %o\n", HC; end if;
  end if;
end for;
print "QUOTIENTS_DONE";
quit;
