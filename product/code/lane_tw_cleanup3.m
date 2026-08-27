// lane_tw_cleanup3.m — final step of the twisted-diagonal cleanup: decide
// J(Q)_tors exactly by HALVING analysis.  From lane_tw_cleanup2.m:
// J(Q)[2] = (Z/2)^3 (tw12) resp. (Z/2)^4 (tw10), and J(Q) = J(Q)_tors is a
// finite abelian 2-group (rank 0 unconditional; order divides 128 resp. 32).
// A finite abelian 2-group G exceeds G[2] iff some nonzero xi in G[2] is
// halvable in G (any order-4 element eta has 2*eta in G[2] nonzero).  So:
// construct all 2-torsion classes as factor-subset Mumford divisors and run
// DivisionPoints(xi, 2) on each.  If none halve: J(Q) = J(Q)[2]; the
// reduced representative of a factor-subset class has degree = deg of the
// factor product (<= 3), and the classes [P - infinity] with P in C(Q) are
// exactly the degree<=1 ones => C(Q) = {infinity} + rational Weierstrass
// points, i.e. the known degenerate points EXHAUST C(Q), unconditionally.
// Usage: cd product/code && magma -b lane_tw_cleanup3.m > ../logs/lane_tw_cleanup3.log
SetColumns(0);
SetMemoryLimit(4*10^9);

QQ := Rationals();
Px<x> := PolynomialRing(QQ);
Pz<T> := PolynomialRing(QQ);

tw12 := -1296*T^8 + 5184*T^7 - 9072*T^6 + 9072*T^5 - 5580*T^4 + 2088*T^3 - 432*T^2 + 36*T;
tw10 := 32*T^7 - 160*T^6 + 256*T^5 - 156*T^4 + 16*T^3 + 16*T^2 - 4*T;

function MonicOdd(f)
    c := LeadingCoefficient(f);
    g := Evaluate(f, Parent(f).1/c) * c^6;
    return Px! [ QQ! Coefficient(g, i) : i in [0..7] ];
end function;

f12flip := Px! [ Coefficient(tw12, 8-i) : i in [0..8] ];
CURVES := [* <"tw12", MonicOdd(f12flip)>, <"tw10", MonicOdd(Px! [ Coefficient(tw10, i) : i in [0..7] ])> *];

for ent in CURVES do
    name := ent[1]; f := ent[2];
    printf "\n== %o ==\n", name;
    fac := [ g[1] : g in Factorization(f) ];
    k := #fac;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    // subsets of factors, nonempty, deg of product <= 3 (each class has a
    // unique such representative among {S, complement-with-infinity})
    idx := [1..k];
    reps := [];
    for S in Subsets(Set(idx)) do
        if #S eq 0 then continue; end if;
        g := &*[ fac[i] : i in S ];
        if Degree(g) le 3 then Append(~reps, g); end if;
    end for;
    printf "%o: %o nontrivial 2-torsion representatives (expect %o)\n", name, #reps, 2^(k-1) - 1;
    nhalf := 0; nfail := 0;
    for g in reps do
        okp := true; pt := J!0;
        try pt := elt<J | g, Px!0, Degree(g)>; catch e okp := false; end try;
        if not okp then
            try pt := J! [g, Px!0]; catch e okp := false; end try;
        end if;
        if not okp then printf "BUILDFAIL %o rep deg %o: %o\n", name, Degree(g), g; nfail +:= 1; continue; end if;
        error if not (2*pt eq J!0), "representative not 2-torsion";
        halves := [];
        okd := true;
        try halves := DivisionPoints(pt, 2); catch e okd := false; end try;
        if not okd then printf "DIVFAIL %o rep %o\n", name, g; nfail +:= 1; continue; end if;
        if #halves gt 0 then
            nhalf +:= 1;
            printf "HALVABLE %o rep %o: %o rational halves e.g. %o\n", name, g, #halves, halves[1];
        end if;
    end for;
    if nfail eq 0 and nhalf eq 0 then
        lin := [ g : g in fac | Degree(g) eq 1 ];
        printf "TORSION_FINAL %o: no 2-torsion class halvable => J(Q) = J(Q)[2] = (Z/2)^%o\n", name, k-1;
        printf "CQ_FINAL %o: C(Q) = {infinity} + %o rational Weierstrass points (x = %o) - the known degenerate points EXHAUST C(Q) unconditionally\n",
            name, #lin, [ -Coefficient(l,0)/Coefficient(l,1) : l in lin ];
    else
        printf "CLEANUP3_INCONCLUSIVE %o: halvable=%o buildfail/divfail=%o\n", name, nhalf, nfail;
    end if;
end for;
printf "TWCLEANUP3_DONE\n";
quit;
