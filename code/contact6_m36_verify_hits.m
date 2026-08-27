//////////////////////////////////////////////////////////////////////
//  Verify [6,6] hits found in the contact-6 extra-root search.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
P<x> := PolynomialRing(Q);

curves := [
    <"eps=1, r=21 and r=-1/3, a=187/21, b=-23/7",
     -111132*x^5 + 2646000*x^4 - 7101864*x^3 + 11226096*x^2 + 4630500*x>,
    <"eps=1, r=4/3, a=-1/42, b=-13/7",
     7112448*x^5 - 36091440*x^4 + 68732496*x^3 - 58231404*x^2 + 18522000*x>
];

for item in curves do
    label := item[1];
    f := item[2];
    print "HIT", label;
    print "f =", f;
    print "factorization =", Factorization(f);

    C := HyperellipticCurve(f);
    J := Jacobian(C);
    G, phi := TorsionSubgroup(J);
    print "torsion_invariants =", Invariants(G);

    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67] do
        fp := ChangeRing(f, GF(p));
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
            continue;
        end if;
        Lp := LPolynomial(ChangeRing(C, GF(p)));
        print "p", p, "L_p", Lp, "factorization", Factorization(Lp);
    end for;
    print "";
end for;

quit;
