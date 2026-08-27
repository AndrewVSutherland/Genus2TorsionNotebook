SetColumns(0); SetMemoryLimit(3*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);
fD := x^6 - 8*x^3 + 25;
K<z> := CyclotomicField(3);
PKt<T> := PolynomialRing(K);
PT := PolynomialRing(Integers());
print "== GL2-type check: quartic L-factors over K = Q(zeta3) ==";
for p in [7,19,37,43,61,73,79] do
    fp := PolynomialRing(GF(p))!fD;
    if not IsSquarefree(fp) then continue; end if;
    Lp := Numerator(ZetaFunction(HyperellipticCurve(fp)));
    facQ := Factorization(PT!Lp);
    facK := Factorization(PKt![K!c : c in Coefficients(Lp)]);
    degsQ := [<Degree(g[1]),g[2]> : g in facQ];
    degsK := [<Degree(g[1]),g[2]> : g in facK];
    printf "p=%o : over Z %o ; over K %o", p, degsQ, degsK;
    if #facK eq 2 and Degree(facK[1][1]) eq 2 then
        // a_p in Z[zeta]: extract traces, check conjugacy
        a1 := -Coefficient(facK[1][1],1)/Coefficient(facK[1][1],2);
        a2 := -Coefficient(facK[2][1],1)/Coefficient(facK[2][1],2);
        conj := hom<K -> K | z^2>;
        printf "  traces %o, %o conjugate? %o", a1/p^0, a2, conj(K!a1) eq K!a2;
    end if;
    printf "\n";
end for;
print "== CM probe: inert primes p = 2 mod 3 ==";
for p in [5,11,17,23,29,41,47,53] do
    fp := PolynomialRing(GF(p))!fD;
    if not IsSquarefree(fp) then continue; end if;
    Lp := Numerator(ZetaFunction(HyperellipticCurve(fp)));
    printf "p=%o : trace(=-c1) %o ; L factors over Z: %o\n", p,
        -Coefficient(Lp,1), [<Degree(g[1]),g[2]> : g in Factorization(PT!Lp)];
end for;
print "GL2CHECK_DONE";
quit;
