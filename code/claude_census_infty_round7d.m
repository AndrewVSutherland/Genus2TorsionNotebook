// Certification of the [2,2,10]-family: Elkies chart, W-points from
// t1 (free), t2 = 3, and the +-s pair with common q; the antidiagonal
// cubic a^3+5a^2-(5/9)ab^2+(1/3)b^2 = 0 is parametrized by
// t1 = -(45+3*tau^2)/(9-5*tau^2), s = tau*t1.
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);

function WPq(t) return t^2/(t^2-1), t/(t^2-1); end function;

function BuildG(tau)
    t1 := -(45+3*tau^2)/(9-5*tau^2);
    s := tau*t1;
    if t1 in {Q!0,1,-1,3,-3} or s in {Q!0,1,-1} then return false,P!0; end if;
    q1, m1 := WPq(t1); lam1 := (q1-1)*(-(2*q1-1)+2*m1)/q1;
    q2, m2 := WPq(Q!3); lam2 := (q2-1)*(-(2*q2-1)+2*m2)/q2;
    q, m := WPq(s);
    Cq := P![ (q-1)^2, 2*q*(2*q-1)*(q-1), q^2 ];
    if Evaluate(Cq, lam1) eq 0 or Evaluate(Cq, lam2) eq 0 then return false,P!0; end if;
    kap1 := (q1-q)*lam1/Evaluate(Cq, lam1);
    kap2 := (q2-q)*lam2/Evaluate(Cq, lam2);
    if kap1 ne kap2 then return false, P!0; end if;   // must vanish on the cubic
    Qq := q*x + kap1*Cq;
    g := 4*Qq^3 + (x^2-6*x+1)*Qq^2 + 2*x*(x-1)*Qq + x^2;
    if Degree(g) lt 5 or Discriminant(g) eq 0 then return false, P!0; end if;
    return true, g;
end function;

function StrictPrime(gi)
    for p in [pp : pp in PrimesUpTo(600) | pp ge 11] do
        if Valuation(Q!LCM([Denominator(c) : c in Coefficients(gi)]), p) ne 0 then continue; end if;
        gp := PolynomialRing(GF(p))!gi;
        if Degree(gp) lt 5 or not IsSquarefree(gp) then continue; end if;
        lp := LPolynomial(HyperellipticCurve(gp));
        cp := Parent(lp)!Reverse(Coefficients(lp));
        if not IsIrreducible(cp) then continue; end if;
        K<pi> := NumberField(cp);
        ok := true;
        for n in [1..12] do
            if Degree(MinimalPolynomial(pi^n)) ne 4 then ok := false; break; end if;
        end for;
        if ok then return p; end if;
    end for;
    return 0;
end function;

print "fibers tau = 2, 5, 7, 1/2, 4:";
invs := []; nex := 0;
for tau in [Q!2, 5, 7, 1/2, 4] do
    ok, g := BuildG(tau);
    if not ok then printf "  tau=%o degenerate\n", tau; continue; end if;
    // count rational Weierstrass roots as a sanity check
    nroots := &+[ r[2] : r in Roots(g) ];
    den := LCM([Denominator(c) : c in Coefficients(g)]);
    gi := den^2*g;
    C := HyperellipticCurve(gi);
    tor := Invariants(TorsionSubgroup(Jacobian(C)));
    sp := StrictPrime(gi);
    gg := G2Invariants(C);
    Append(~invs, gg);
    if tor eq [2,2,10] and sp gt 0 then nex +:= 1; end if;
    printf "  tau=%o: rational W-roots %o, torsion %o, strict p=%o\n", tau, nroots, tor, sp;
end for;
printf "distinct G2: %o of %o; exact+strict: %o\n", #Set(invs), #invs, nex;
printf "ROUND7D_DONE\n";
quit;
