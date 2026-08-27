// Torsion of the 11 glued genus-2 curves in Frengley's N-congruences repo
// (optimal degree-N covers of (N,-1)-congruent pairs, N = 5,7,15,16,17).
SetColumns(0);
SetMemoryLimit(8*10^9);
P<x> := PolynomialRing(Rationals());
base := "~/.claude/jobs/a1db5dd4/tmp/N-congruences/genus-2-curves/";
for N in [5,7,15,16,17] do
    fs := eval Read(base cat Sprintf("%o.m", N));
    if Type(fs[1]) ne SeqEnum then fs := [fs]; end if;
    for i in [1..#fs] do
        f := P!fs[i];
        C := HyperellipticCurve(f);
        Cs := SimplifiedModel(ReducedMinimalWeierstrassModel(C));
        T := TorsionSubgroup(Jacobian(Cs));
        printf "N=%o ex=%o TORSION %o order %o disc %o\n", N, i, Invariants(T), #T,
            Factorization(Integers()!Discriminant(Cs));
    end for;
end for;
printf "FRENGLEY_TORSION_DONE\n";
quit;
