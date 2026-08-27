// claude_222_genus3_round2.m — closing machinery for the two genus-3 blockers.
// X1: bielliptic over E1 = 92-class rank-1 curve => J(X1) ~ E1 x Prym.
//     Compute L_p(X1)/L_p(E1) for good p: the Prym's quadratic Euler factors.
//     Census the Prym's real-subfield discs (RM/split/generic signature) and
//     print (p, a_p) pairs for an LMFDB isogeny-class lookup.
// X2: its genus-1 quotient Y2 had NO points to 1e5. Print Y2's model and test
//     local solvability; if Y2(Q_v) = {} for some v, X2(Q) = {} — theorem.
SetColumns(0);
SetMemoryLimit(24*10^9);
Q := Rationals();
A2<s,x> := AffineSpace(Q, 2);
c1 := (s-1)^2*x^3 + (2*s^4-4*s^3+5*s^2-4*s+1)*x^2
    + (s^6-2*s^5+3*s^4-2*s^3+2*s^2-2*s+1)*x - s^2*(s^2-s+1)^2;
X1 := Curve(A2, c1);
E1 := EllipticCurve([0,0,0,-1,1]);
printf "E1 conductor %o\n", Conductor(E1);

// ---- X1 point counts over F_p via the function field ----
PQ<u> := PolynomialRing(Q);
for p in PrimesInInterval(3, 80) do
  try
    k := GF(p);
    A2k<sk,xk> := AffineSpace(k, 2);
    c1k := (sk-1)^2*xk^3 + (2*sk^4-4*sk^3+5*sk^2-4*sk+1)*xk^2
        + (sk^6-2*sk^5+3*sk^4-2*sk^3+2*sk^2-2*sk+1)*xk - sk^2*(sk^2-sk+1)^2;
    X1k := Curve(A2k, c1k);
    FFk := AlgorithmicFunctionField(FunctionField(X1k));
    Z := ZetaFunction(FFk);
    L := Numerator(Z);
    // L in Z[T], degree 6; divide by Euler factor of E1
    PT<T> := Parent(L);
    aE := p + 1 - #ChangeRing(E1, GF(p));
    LE := 1 - aE*T + p*T^2;
    ok, LA := IsDivisibleBy(L, PT!LE);
    if ok then
      printf "p=%o L_Prym=%o\n", p, LA;
    else
      printf "p=%o E1-factor does NOT divide (bad p for the fibration?) L=%o\n", p, L;
    end if;
  catch e
    printf "p=%o failed\n", p;
  end try;
end for;

// ---- X2 quotient model + local solvability ----
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
G2 := AutomorphismGroup(ProjectiveClosure(X2));
Y2, q2 := CurveQuotient(G2);
printf "Y2 ambient: %o\n", AmbientSpace(Y2);
printf "Y2 equations:\n%o\n", DefiningEquations(Y2);
printf "Y2 genus %o\n", Genus(Y2);
// try to recognize a hyperelliptic/genus-one model for local solvability
try
  ok2, HY, mp2 := IsHyperelliptic(Y2);
  if ok2 then
    printf "Y2 hyperelliptic model: %o\n", HY;
    f2, h2 := HyperellipticPolynomials(HY);
    F2 := 4*f2 + h2^2;
    printf "Y2 as y^2 = %o\n", F2;
    for p in [2,3,5,7,11,13,17,19,23] cat [q : q in PrimeDivisors(Integers()!Discriminant(HY))] do
      printf "  locally solvable at %o: %o\n", p, IsLocallySolvable(HY, p);
    end for;
    // Jacobian = elliptic curve (genus 1 with degree-2 model): rank of Jac
    E2j := Jacobian(GenusOneModel(HY));
    E2m, _ := MinimalModel(E2j);
    printf "Jac(Y2) minimal: %o conductor %o\n", aInvariants(E2m), Conductor(E2m);
    rl, ru := RankBounds(E2m);
    printf "Jac(Y2) rank bounds [%o,%o] torsion %o\n", rl, ru, Invariants(TorsionSubgroup(E2m));
  else
    print "Y2 not recognized hyperelliptic; raw model printed above";
  end if;
catch e
  printf "Y2 model handling failed: %o\n", e`Object;
end try;
print "ROUND2_DONE";
quit;
