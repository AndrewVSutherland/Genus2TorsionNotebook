SetColumns(0);
Q := Rationals();
R2<v1,v2> := PolynomialRing(Q,2);
g := func<t | t^5 + 4*t^4 + 10*t^3 + 215/16*t^2 + 35/4*t + 35/16>;
// psi(v) = g(v)/(v+1)^4 ; X2 : psi(v1)=psi(v2), v1 <> v2
N := g(v1)*(v2+1)^4 - g(v2)*(v1+1)^4;
F2 := ExactQuotient(N, v1-v2);
print "F2 =", F2;
A2 := AffineSpace(R2);
C2 := Curve(A2, F2);
print "X2: irreducible?", IsIrreducible(F2), " dims", [Degree(F2,v1),Degree(F2,v2)];
PC2 := ProjectiveClosure(C2);
print "X2 geometric genus =", Genus(C2);
// factor check
print "factorization of F2:", Factorization(F2);
quit;
