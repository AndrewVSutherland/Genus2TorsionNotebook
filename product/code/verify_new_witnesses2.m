// verify_new_witnesses2.m — second batch: the HLP §3.7 instantiation winners.
// Usage: magma -b verify_new_witnesses2.m > ../logs/verify_witnesses2.log
SetColumns(0);
SetMemoryLimit(4*10^9);
load "split_lab.m";  // run from product/code/

WITS := [*
  <"2.2.24_hlp37_t=241/81_u=1/3",
   [-1729978236,11223443268,-15881777387,10516435114,-574778279,-4134794160,581449680]>,
  <"2.2.24_hlp37_t=-95_u=3/7",
   [0,-32640753600,548428255729,365172573934,-179818919099,-75368404020,21208976100]>,
  <"2.2.4.8_hlp37_t=1/6_u=4/17",
   [949701503960100,812225828599020,-1568053362370059,-1006367755763986,860505249465645,289009358554092,-144295356865660]>
*];

for w in WITS do
    g := RQx!w[2];
    printf "== WITNESS %o ==\n", w[1];
    printf "g = %o\n", g;
    dz := Integers()!Discriminant(g);
    printf "disc(g) = %o digits\n", #IntegerToString(AbsoluteValue(dz));
    T := [Integers()|];
    try
        T := ExactTorsion(g);
    catch e
        printf "TORSION FAILED\n";
        continue;
    end try;
    printf "EXACT TORSION = %o\n", T;
    red, n := LpReducibleCount(g, 200);
    printf "L-poly reducible at %o/%o good primes < 200\n", red, n;
    C := HyperellipticCurve(g);
    okm := true;
    try C := ReducedMinimalWeierstrassModel(C); catch e okm := false; end try;
    L := [* *];
    try L := RichelotIsogenousSurfaces(C); catch e L := [* *]; end try;
    for X in L do
        if Type(X) eq SetCart then
            E1 := Component(X,1); E2 := Component(X,2);
            printf "FACTORS: E1 cond %o tors %o ainvs %o ; E2 cond %o tors %o ainvs %o\n",
                Conductor(E1), Invariants(TorsionSubgroup(E1)), aInvariants(MinimalModel(E1)),
                Conductor(E2), Invariants(TorsionSubgroup(E2)), aInvariants(MinimalModel(E2));
            break;
        end if;
    end for;
end for;
printf "SEARCH_DONE verify2\n";
quit;
