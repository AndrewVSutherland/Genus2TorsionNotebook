// Lane 8 (2026-07-25, third session): try to certify f70's NON-SIMPLICITY
// unconditionally, by exhibiting a degree-3 elliptic subfield of Q(C).
//
// If J(C) is the (3,3)-glue of E1 x E2 then C carries degree-3 maps C -> E1 and
// C -> E2, i.e. Q(C) has index-3 subfields of genus 1.  Exhibiting one is an
// UNCONDITIONAL proof that J(C) is not simple (a fortiori not geometrically
// simple, since the map is defined over Q).  Contrast with the L-polynomial
// evidence, which is Faltings-conditional (no effective Sturm bound was used).
SetColumns(0);
SetMemoryLimit(4*10^9);
QQ := Rationals();
R<x> := PolynomialRing(QQ);

g := 3168*x^5 + 697*x^4 - 23220*x^3 + 37620*x^2 - 23328*x + 5184;   // 144*f70
C := HyperellipticCurve(g);
printf "C : y^2 = %o\n", g;
printf "genus %o   TorsionSubgroup %o\n", Genus(C), Invariants(TorsionSubgroup(Jacobian(C)));

F := FunctionField(C);
printf "function field: %o\n", F;

t0 := Cputime();
ok := true;
S := [];
try
    S := Subfields(F);
catch e
    ok := false;
    printf "Subfields FAILED: %o\n", e`Object;
end try;
printf "Subfields time %o  success %o\n", Cputime(t0), ok;
if ok then
    printf "number of subfields found: %o\n", #S;
    for s in S do
        K := s[1];
        d := Degree(F, K);
        gg := -1;
        try gg := Genus(K); catch e gg := -1; end try;
        printf "  subfield of index %o, genus %o\n", d, gg;
    end for;
end if;

// independent: the two elliptic curves and their mod-3 representations
E1 := EllipticCurve([1,0,0,-45,81]);
E2 := EllipticCurve([1,0,0,-5774401,5346023177]);
printf "\nE1 cond %o torsion %o ; E2 cond %o torsion %o\n",
       Conductor(E1), Invariants(TorsionSubgroup(E1)), Conductor(E2), Invariants(TorsionSubgroup(E2));
d1 := DivisionPolynomial(E1, 3);
d2 := DivisionPolynomial(E2, 3);
printf "3-division polynomials: deg %o, %o\n", Degree(d1), Degree(d2);
printf "  factor degrees E1: %o\n", [<Degree(t[1]),t[2]> : t in Factorization(d1)];
printf "  factor degrees E2: %o\n", [<Degree(t[1]),t[2]> : t in Factorization(d2)];
K1 := SplittingField(d1);  K2 := SplittingField(d2);
printf "  splitting field of psi_3: deg %o and %o ; ISOMORPHIC: %o\n",
       Degree(K1), Degree(K2), IsIsomorphic(K1, K2);
printf "  disc(E1) sqfree core %o ; disc(E2) sqfree core %o\n",
       Squarefree(Integers()!Discriminant(E1)), Squarefree(Integers()!Discriminant(E2));

printf "SEARCH_DONE f70subfield\n";
quit;
