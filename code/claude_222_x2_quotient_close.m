// claude_222_x2_quotient_close.m — decide Y2(Q) (genus-1 quotient of the
// Flynn43 extra-root cover X2).  Y2 = intersection of two quadrics in P^3
// (degree-4 genus-one model).  Minimise+Reduce the model, then:
//   - local solvability at all bad primes (and small primes);
//   - Jacobian elliptic curve: conductor, rank bounds, torsion;
//   - point search on the reduced model.
// Outcomes: a point => X2 has rational points upstairs (pull back!);
// no local point at some v => X2(Q) empty (THEOREM);
// ELS + rank 0 + no points => nontrivial Sha[4]-torsor (near-theorem, then
// decide by comparing against E's finitely many rational points via descent).
SetColumns(0);
SetMemoryLimit(24*10^9);
Q := Rationals();

// rebuild X2 and its quotient (as in round 2)
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

md := GenusOneModel(Y2);
printf "raw model degree %o\n", Degree(md);
mdm, tr1 := Minimise(md);
mdr, tr2 := Reduce(mdm);
printf "REDUCED MODEL:\n%o\n", mdr;
E := Jacobian(mdr);
Em, _ := MinimalModel(E);
printf "Jac(Y2) = %o, conductor %o\n", aInvariants(Em), Conductor(Em);
rl, ru := RankBounds(Em);
printf "rank bounds [%o, %o], torsion %o\n", rl, ru, Invariants(TorsionSubgroup(Em));
bad := PrimeDivisors(Integers()!Discriminant(Em)) cat [2,3,5,7];
bad := Sort(SetToSequence(SequenceToSet(bad)));
els := true;
for p in bad do
  ok := IsLocallySolvable(mdr, p);
  printf "locally solvable at %o: %o\n", p, ok;
  if not ok then els := false; end if;
end for;
printf "everywhere locally solvable (tested set): %o\n", els;
C4 := Curve(mdr);
pts := PointSearch(C4, 10^7);
printf "points on reduced model to 1e7: %o\n", #pts;
if #pts gt 0 then
  printf "POINTS FOUND: %o\n", [pts[i] : i in [1..Minimum(#pts,5)]];
  print "=> X2 may have rational points: pull back through the quotient!";
end if;
print "X2CLOSE_DONE";
quit;
