//////////////////////////////////////////////////////////////////////
// Local sieve for quadratic twists of the order-96 record curve.
//
// Every quadratic twist retains the six rational Weierstrass points and
// geometric simplicity.  Thus a twist with torsion [2,2,2,24] would be a
// target.  At a good odd prime p there are only two reductions: the curve
// and its nonsquare twist.  If neither Jacobian order is divisible by 192,
// the entire quadratic-twist route is ruled out at that prime.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetLogFile("results/target_22224_record_twist_sieve.log" : Overwrite := true);
Q := Rationals();
P<x> := PolynomialRing(Q);

F := 3027600*x^6 + 2950382280*x^5 + 602288814361*x^4
     - 63417934304484*x^3 - 2122595910966478*x^2
     + 128056619498204124*x + 3322970988364151397;

print "TARGET_22224_RECORD_TWIST_SIEVE_START";
decisive := [];
for p in PrimesInInterval(5,199) do
    k := GF(p);
    R<X> := PolynomialRing(k);
    fp := R!F;
    if Degree(fp) ne 6 or Discriminant(fp) eq 0 then
        continue;
    end if;
    eta := k!2;
    while IsSquare(eta) do
        eta +:= 1;
    end while;
    np := #Jacobian(HyperellipticCurve(fp));
    nm := #Jacobian(HyperellipticCurve(eta*fp));
    plus := np mod 192 eq 0;
    minus := nm mod 192 eq 0;
    print "TWIST_LOCAL",p,"plus_order",np,"minus_order",nm,
          "plus_192",plus,"minus_192",minus,
          "v2",Valuation(np,2),Valuation(nm,2),
          "v3",Valuation(np,3),Valuation(nm,3);
    if not plus and not minus then
        Append(~decisive,p);
        print "TWIST_ROUTE_OBSTRUCTION",p,np,nm;
        break;
    end if;
end for;
print "TARGET_22224_RECORD_TWIST_SIEVE_DONE","decisive",decisive;
quit;
