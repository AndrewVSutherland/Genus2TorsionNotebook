// claude_222_genus3_analysis.m — structure of the two genus-3 covers that are
// the last in-Flynn hopes for [2,22]:
//   X1: branch-point cubic-root cover  c_-(s,x) = 0  (notes/order222_from_order11.md;
//       c_+ = c_-(-s,x) is isomorphic via s -> -s)
//   X2: the degree-9 genus-3 component of the Flynn43 extra-root cover
// For each: confirm genus; Frobenius charpolys of the Jacobian at many good p
// (degree-6 L-polynomial factorization patterns reveal isogeny factors);
// automorphism group; quotients by involutions (elliptic/genus-2 quotients =>
// rank machinery / Chabauty routes).
SetColumns(0);
SetMemoryLimit(24*10^9);
Q := Rationals();
A2<s,x> := AffineSpace(Q, 2);

c1 := (s-1)^2*x^3 + (2*s^4-4*s^3+5*s^2-4*s+1)*x^2
    + (s^6-2*s^5+3*s^4-2*s^3+2*s^2-2*s+1)*x - s^2*(s^2-s+1)^2;
X1 := Curve(A2, c1);
printf "X1 irred %o absirred %o genus %o\n", IsIrreducible(X1), IsAbsolutelyIrreducible(X1), Genus(X1);

// X2: rebuild the extra-root cover of the Flynn43 family and take the deg-9 factor
FF<T> := FunctionField(Q);
uu := (-11/100*T^3 + 919/50*T^2 - 8017/25*T + 40154/25)/(T^3 - 23*T^2 + 168*T - 396);
vv := (2401/1600*T^4 - 12593/200*T^3 + 205647/200*T^2 - 385757/50*T + 2253001/100)/(T^4 - 34*T^3 + 421*T^2 - 2244*T + 4356);
tt := (-21609/40000*T^6 + 435561/10000*T^5 - 2897123/2000*T^4 + 6390763/250*T^3 - 126932427/500*T^2 + 844854361/625*T - 1894773841/625)/(T^6 - 46*T^5 + 865*T^4 - 8520*T^3 + 46440*T^2 - 133056*T + 156816);
PxF<X> := PolynomialRing(FF);
FlT := X^6+2*X^5+(2*tt+3)*X^4+2*X^3+(tt^2+1)*X^2+2*tt*(1-tt)*X+tt^2;
Q4 := FlT div (X^2 + FF!uu*X + FF!vv);
cfs := Coefficients(Q4);
den := LCM([Denominator(c) : c in cfs]);
R2<t2,r2> := PolynomialRing(Q, 2);
num := &+[ Evaluate(Numerator(cfs[i]*den), t2) * r2^(i-1) : i in [1..#cfs] ];
A2b<tb,rb> := AffineSpace(Q,2);
comp9 := 0;
for fa in Factorization(Evaluate(num, [tb, rb])) do
  if TotalDegree(fa[1]) ge 8 then comp9 := fa[1]; end if;
end for;
assert comp9 ne 0;
X2 := Curve(A2b, comp9);
printf "X2 irred %o absirred %o genus %o\n", IsIrreducible(X2), IsAbsolutelyIrreducible(X2), Genus(X2);

// Frobenius factorization patterns for both
for pair in [* <X1, "X1">, <X2, "X2"> *] do
  C := pair[1]; nm := pair[2];
  printf "==== %o L-polynomial factorization degrees ====\n", nm;
  for p in PrimesInInterval(11, 60) do
    try
      Cp := ChangeRing(C, GF(p));
      if not IsIrreducible(Cp) then printf "p=%o reducible reduction\n", p; continue; end if;
      Z := ZetaFunction(Cp);
      L := Numerator(Z);
      degs := Sort([Degree(fe[1])^^fe[2] : fe in Factorization(L)]);
      printf "p=%o Lfact=%o\n", p, degs;
    catch e
      printf "p=%o failed\n", p;
    end try;
  end for;
  // automorphisms of the function field
  try
    Fn := AlgorithmicFunctionField(FunctionField(C));
    A := Automorphisms(Fn);
    printf "%o: #automorphisms(FF) = %o\n", nm, #A;
  catch e
    printf "%o: automorphism computation failed: %o\n", nm, e`Object;
  end try;
end for;
print "GENUS3_ANALYSIS_DONE";
quit;
