// Lane 2: Cab (the c=2 slice in symmetric coordinates) has GENUS 0.  Parametrize it,
// then the two square conditions become an explicit (2,2)-cover of P^1:
//    X :  u^2 = d1(t) = a(t)^2 - 4 b(t),   v^2 = d2(t) = (4-a(t))^2 - 4/b(t).
//    u^2 = d1  <->  three rational z  ([2,2,14]);   BOTH  <->  five ([2,2,2,14], order 112).
SetColumns(0);
Q := Rationals();
A2<a,b> := AffineSpace(Q,2);
F := -3*a*b^2 + 10*b^2 + 2*a^2*b - 8*a*b + 3*a - 2;
C := Curve(A2,F);  PC := ProjectiveClosure(C);
printf "genus(Cab) = %o\n", Genus(PC);
pt := PC ! [-22/9, 7/27, 1];
prm := Parametrization(PC, Place(pt));
printf "parametrization: %o\n", DefiningEquations(prm);
P1<t,s> := Domain(prm);
eqs := DefiningEquations(prm);
Rt<T> := FunctionField(Q);
at := Evaluate(eqs[1],[T,1])/Evaluate(eqs[3],[T,1]);
bt := Evaluate(eqs[2],[T,1])/Evaluate(eqs[3],[T,1]);
printf "a(T) = %o\nb(T) = %o\n", at, bt;
printf "check on Cab: %o\n", -3*at*bt^2 + 10*bt^2 + 2*at^2*bt - 8*at*bt + 3*at - 2 eq 0;
d1 := at^2 - 4*bt;
d2 := (4-at)^2 - 4/bt;
printf "d1(T) = %o\nd2(T) = %o\n", d1, d2;
// squarefree part of the square class of a rational function
function SqPart(f)
  n := Numerator(f); d := Denominator(f);
  g := n*d;
  P := Parent(g);
  out := P!1;
  for tp in Factorisation(g) do
    if IsOdd(tp[2]) then out *:= tp[1]; end if;
  end for;
  return out * (LeadingCoefficient(g) eq 0 select 1 else 1);
end function;
lc1 := LeadingCoefficient(Numerator(d1)*Denominator(d1));
lc2 := LeadingCoefficient(Numerator(d2)*Denominator(d2));
g1 := SqPart(d1); g2 := SqPart(d2);
// re-attach the leading constant (square class matters!)
c1 := (Numerator(d1)*Denominator(d1)) div (g1 * ((Numerator(d1)*Denominator(d1)) div g1));
g1 := g1 * LeadingCoefficient(Numerator(d1)*Denominator(d1)) / LeadingCoefficient(g1)^1;
g2 := g2 * LeadingCoefficient(Numerator(d2)*Denominator(d2)) / LeadingCoefficient(g2)^1;
// (the quotient (num*den)/squarefreepart is a square times a constant; fix the constant below)
q1 := (Numerator(d1)*Denominator(d1)) / g1;  q2 := (Numerator(d2)*Denominator(d2)) / g2;
printf "check q1 is a square in Q(T): %o ; q2: %o\n", IsSquare(q1), IsSquare(q2);
printf "g1 = %o   (deg %o)\n", g1, Degree(g1);
printf "g2 = %o   (deg %o)\n", g2, Degree(g2);
printf "factor g1: %o\n", Factorisation(g1);
printf "factor g2: %o\n", Factorisation(g2);
printf "gcd(g1,g2) = %o\n", GCD(g1,g2);
P<X> := PolynomialRing(Q);
G1 := P ! [Coefficient(g1,i) : i in [0..Degree(g1)]];
G2 := P ! [Coefficient(g2,i) : i in [0..Degree(g2)]];
for pr in [<G1,"u^2=g1 (three rational z, [2,2,14])">, <G2,"v^2=g2">] do
  GG := pr[1];
  if Degree(GG) in {3,4} then
    EE := MinimalModel(EllipticCurve(HyperellipticCurve(GG)));
    printf "%o -> %o  conductor %o  MW %o\n", pr[2], EE, Conductor(EE), Invariants(MordellWeilGroup(EE));
  elif Degree(GG) in {5,6} then
    HH := HyperellipticCurve(GG);
    printf "%o -> hyperelliptic genus %o\n", pr[2], Genus(HH);
  end if;
end for;
// the (2,2)-cover X : u^2 = g1(T), v^2 = g2(T)
KT := FunctionField(Q); TT := KT.1;
PK<Y> := PolynomialRing(KT);
F1 := FunctionField(Y^2 - Evaluate(g1,TT));
PK1<Y1> := PolynomialRing(F1);
F2 := FunctionField(Y1^2 - F1!Evaluate(g2,TT));
printf "GENUS of the order-112 curve X (u^2=g1, v^2=g2) = %o\n", Genus(F2);
printf "PARAM_DONE\n";
quit;
