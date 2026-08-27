// x7x10_glue.m — algebraic (3,3)-glue of the sieve's genuine non-isogenous
// 3-congruent pair x7(2/11) x x10(-1/2)/(1/3): product [70] (known) -- a
// PIPELINE CONTROL for congruence->gluing->torsion, and a [210]-shot if a
// 3-gain occurs.
SetColumns(0);
SetMemoryLimit(4*10^9);
load "split_lab.m";
RQ := Rationals();
b7 := (2/11)^3-(2/11)^2; c7 := (2/11)^2-(2/11);
E7 := EllipticCurve([1-c7, -b7, -b7, 0, 0]);
for u in [RQ|-1/2, 1/3] do
    b10 := u^3*(u-1)*(2*u-1)/(u^2-3*u+1)^2; c10 := -u*(u-1)*(2*u-1)/(u^2-3*u+1);
    E10 := EllipticCurve([1-c10, -b10, -b10, 0, 0]);
    printf "pair x7(2/11) x x10(%o): torsions %o %o\n", u,
        Invariants(TorsionSubgroup(E7)), Invariants(TorsionSubgroup(E10));
    L := [];
    try L := Genus2Elliptic3(E7, E10); catch e printf "GE3 error\n"; end try;
    printf "Genus2Elliptic3: %o curves\n", #L;
    for k in [1..#L] do
        C := L[k];
        ok, g := IntegralSextic(C);
        if not ok then continue; end if;
        einv := [Integers()|-1];
        try einv := ExactTorsion(g); catch e printf "torsfail\n"; continue; end try;
        printf "X7X10 RESULT %o invs=%o g=%o\n", k, einv, g;
    end for;
end for;
printf "X7X10_DONE\n";
quit;
