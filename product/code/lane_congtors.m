// lane_congtors.m — treasure hunt over Frengley's LMFDB (N,r)-congruent pair
// lists: for N-congruent NON-isogenous pairs (E,F) with an ANTI-isometry
// available (some s with r*s^2 = -1 mod N) and gcd(N, #T1 * #T2) = 1, the
// (N,N)-gluing A = (ExF)/graph is a principally polarized abelian surface
// whose rational torsion CONTAINS T1 x T2 (the product injects: the graph
// meets the rational torsion trivially).  For non-geometrically-isogenous
// pairs A is a Jacobian.  So any pair with T1 x T2 not in KNOWN is a
// candidate NEW GROUP (construction of the genus-2 curve to follow).
// Usage: cd product/code && magma -b lane_congtors.m > ../logs/lane_congtors.log
SetColumns(0);
SetMemoryLimit(4*10^9);
load "split_lab.m";  // for KNOWN

// files with an anti-isometric scaling available: r*s^2 = -1 mod N solvable
USABLE := [ <5,1>, <7,3>, <8,7>, <10,1>, <11,2>, <12,11>, <13,1>, <14,3>, <16,3> ];

function TorsOf(lbl)
    E := 0; ok := true;
    try E := EllipticCurve(lbl); catch e ok := false; end try;
    if not ok then return false, 0, [Integers()|]; end if;
    return true, E, Invariants(TorsionSubgroup(E));
end function;

for nr in USABLE do
    N := nr[1]; r := nr[2];
    fn := Sprintf("../N-congruences/examples/lmfdb/%o-%o.m", N, r);
    ok := true; s := "";
    try s := Read(fn); catch e ok := false; end try;
    if not ok then printf "MISSING %o\n", fn; continue; end if;
    L := eval s;
    printf "== (N,r)=(%o,%o): %o pairs ==\n", N, r, #L;
    ncand := 0;
    for pr in L do
        ok1, E1, T1 := TorsOf(pr[1]);
        ok2, E2, T2 := TorsOf(pr[2]);
        if not (ok1 and ok2) then printf "LOOKUPFAIL %o %o\n", pr[1], pr[2]; continue; end if;
        n1 := IsEmpty(T1) select 1 else &*T1;
        n2 := IsEmpty(T2) select 1 else &*T2;
        if GCD(N, n1*n2) ne 1 then continue; end if;
        if n1*n2 lt 16 then continue; end if;   // report threshold
        prod := Invariants(AbelianGroup(T1 cat T2));
        isnew := not prod in KNOWN;
        // geometric isogeny (quadratic twist) check
        tw := IsQuadraticTwist(E1, E2);
        printf "PAIR N=%o %o %o T1=%o T2=%o prod=%o order=%o %o%o\n",
            N, pr[1], pr[2], T1, T2, prod, n1*n2,
            isnew select "NEW" else "known", tw select " TWISTPAIR" else "";
        ncand +:= 1;
    end for;
    printf "== (%o,%o): %o pairs with torsion product >= 16 and coprime to N ==\n", N, r, ncand;
end for;
printf "CONGTORS_DONE\n";
quit;
