// Certification of the [2,10]-family (Elkies 3-Weierstrass chart),
// slice (t1,t2,t3) = (t, 2, 3), t free: exactness, strictness, moduli.
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);

function WPdata(t)
    q := t^2/(t^2-1); m := t/(t^2-1);
    return (q-1)*(-(2*q-1)+2*m)/q, q;
end function;

function BuildG(ts)
    lams := []; qs := [];
    for t in ts do lam, q := WPdata(t); Append(~lams, lam); Append(~qs, q); end for;
    if #Set(lams) ne 3 then return false, P!0; end if;
    Qq := Interpolation(lams, [ qs[i]*lams[i] : i in [1..3] ]);
    g := 4*Qq^3 + (x^2-6*x+1)*Qq^2 + 2*x*(x-1)*Qq + x^2;
    if Discriminant(g) eq 0 then return false, P!0; end if;
    return true, g;
end function;

// strict prime: irreducible quartic charpoly with [Q(pi^n):Q]=4, n<=12
function StrictPrime(gi)
    for p in [pp : pp in PrimesUpTo(600) | pp ge 11] do
        if Valuation(Q!LCM([Denominator(c) : c in Coefficients(gi)]), p) ne 0 then continue; end if;
        gp := PolynomialRing(GF(p))!gi;
        if Degree(gp) lt 5 or not IsSquarefree(gp) then continue; end if;
        Cp := HyperellipticCurve(gp);
        lp := LPolynomial(Cp);
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

print "slice (t,2,3), fibers t = 5,7,9,11,13:";
invs := [];
nex := 0;
for t in [5,7,9,11,13] do
    ok, g := BuildG([Q!t, Q!2, Q!3]);
    if not ok then printf "  t=%o degenerate\n", t; continue; end if;
    den := LCM([Denominator(c) : c in Coefficients(g)]);
    gi := den^2 * g;
    C := HyperellipticCurve(gi);
    J := Jacobian(C);
    tor := Invariants(TorsionSubgroup(J));
    sp := StrictPrime(gi);
    gg := G2Invariants(C);
    Append(~invs, gg);
    if tor eq [2,10] and sp gt 0 then nex +:= 1; end if;
    printf "  t=%o torsion %o strict p=%o\n", t, tor, sp;
end for;
printf "distinct G2-invariants among fibers: %o of %o\n", #Set(invs), #invs;
printf "exact+strict fibers: %o\n", nex;
printf "ROUND6B_DONE\n";
quit;
