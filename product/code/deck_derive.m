// deck_derive.m — derive the deck correspondence of the E26 family:
// factor numerator(j(t) - j(u)) in Q[t,u]; components other than t-u are the
// deck curves (j-equal parameter pairs).  Pure symbolic, ~seconds.
// Usage: cd product/code && magma -b deck_derive.m > ../logs/deck_derive.log
SetColumns(0);
SetMemoryLimit(2*10^9);

QQ := Rationals();
R2<t,u> := PolynomialRing(QQ, 2);
F2 := FieldOfFractions(R2);

// j-invariant of E26(t) from the Weierstrass roots (2-division values)
function Jof(x)
    r1 := (-2*x+10)/((x+3)*(x-3));
    r2 := (-x^3+7*x^2-11*x+5)/(4*(x+3)*(x-3)^2);
    r3 := (-2*x^2+4*x-2)/((x+3)^2*(x-3));
    // lambda for the Legendre-type form via cross-ratio of roots
    lam := (r3-r1)/(r2-r1);
    return 256*(lam^2-lam+1)^3/(lam^2*(lam-1)^2);
end function;

jt := Jof(F2!t);
ju := Jof(F2!u);
num := Numerator(jt - ju);
fac := Factorization(R2!num);
printf "== deck components of j(t)=j(u) for E26 ==\n";
for f in fac do
    printf "COMPONENT deg=(%o,%o) mult=%o: %o\n", Degree(f[1], t), Degree(f[1], u), f[2], f[1];
end for;
printf "DECK_DERIVE_DONE\n";
quit;
