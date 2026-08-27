//////////////////////////////////////////////////////////////////////
//  Frobenius simplicity certificate for Elkies' printed order-32 curve.
//
//  The reconstructed genus-0 family contains this member.  An irreducible
//  degree-4 L-polynomial at one good prime certifies Q-simplicity for the
//  printed specialization, hence the reconstructed component is not contained
//  in the decomposable locus.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
P<x> := PolynomialRing(Q);

f := (15*x - 1)*(1056*x^4 + 156183*x^3 + 26297*x^2 + 649*x - 121);
C := HyperellipticCurve(f);

print "Elkies printed N=32 curve";
print "f =", f;

found := false;
for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71] do
    try
        Cp := ChangeRing(C, GF(p));
        Lp := LPolynomial(Cp);
        fac := Factorization(Lp);
        print "p", p, "L", Lp, "factorization", fac;
        if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
            print "simple certificate at p =", p;
            found := true;
            break;
        end if;
    catch e
        print "p", p, "skipped";
    end try;
end for;

if not found then
    print "no certificate found in tested primes";
end if;
