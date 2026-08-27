// [2,10]-family from Elkies (ConM 796) Thm 2.3 / Cor 2.4 + 3 rational
// Weierstrass points.  Chart: (L,L') = (x,1), Q quadratic;
// g = 4Q^3 + (x^2-6x+1)Q^2 + 2x(x-1)Q + x^2;  T = [D-K], D = {Q=0,y=0}
// has order 5.  Weierstrass point at x=lam:  g(lam)=0  <=>
// lam^2 q^2 + 2 lam q(2q-1)(q-1) + (q-1)^2 = 0 with q = Q(lam)/lam,
// solvable over Q iff q(q-1) = m^2:  q = t^2/(t^2-1), m = t/(t^2-1),
// lam = (q-1)(-(2q-1)+2m)/q.  Three t's -> three (lam_i, q_i) ->
// interpolate Q (3 conditions, unique) -> g has 3 rational Weierstrass
// roots -> J(Q) >= Z/5 x (Z/2)^2 = [2,10].
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);

function WPdata(t)
    q := t^2/(t^2-1);
    m := t/(t^2-1);
    lam := (q-1)*(-(2*q-1)+2*m)/q;
    return lam, q;
end function;

function BuildG(ts)
    lams := []; qs := [];
    for t in ts do
        lam, q := WPdata(t);
        Append(~lams, lam); Append(~qs, q);
    end for;
    if #Set(lams) ne 3 then return false, P!0, lams; end if;
    // interpolate quadratic Q with Qq(lam_i) = q_i * lam_i
    Qq := Interpolation(lams, [ qs[i]*lams[i] : i in [1..3] ]);
    g := 4*Qq^3 + (x^2-6*x+1)*Qq^2 + 2*x*(x-1)*Qq + x^2;
    return true, g, lams;
end function;

function TorsOf(g)
    den := LCM([Denominator(c) : c in Coefficients(g)]);
    gi := den^2 * g;
    C := HyperellipticCurve(gi);
    J := Jacobian(C);
    return Invariants(TorsionSubgroup(J));
end function;

print "PART 1: validation on random triples";
triples := [ [2,3,5], [2,3,7], [3,5,7], [2,5,9], [2,3,-5], [4,7,9] ];
for ts in triples do
    ok, g, lams := BuildG([Q!t : t in ts]);
    if not ok then printf "  ts=%o DEGENERATE\n", ts; continue; end if;
    // assert the three Weierstrass roots
    assert &and[ Evaluate(g, l) eq 0 : l in lams ];
    if Discriminant(g) eq 0 then printf "  ts=%o singular\n", ts; continue; end if;
    tor := TorsOf(g);
    printf "  ts=%o roots ok, torsion %o\n", ts, tor;
end for;
printf "ROUND6_PART1_DONE\n";
quit;
