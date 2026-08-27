// verify_new_witnesses.m — independent rebuild + certification of the NEW
// exact split-torsion witnesses produced by the 2026-08-12 session lanes.
// For each curve: exact TorsionSubgroup on the integral y^2=g model, L-poly
// reducibility count (split evidence), and elliptic factor extraction via a
// rational Richelot to a product (SetCart) with factor conductors/torsions.
// Usage: magma -b verify_new_witnesses.m > ../logs/verify_witnesses.log
SetColumns(0);
SetMemoryLimit(4*10^9);
load "split_lab.m";  // run from product/code/

WITS := [*
  <"8.8_W1_x1fam_8(1/7)x8(-5/7)",
   [84285504,-4535664,-4045487,1800118,88597,88596,836]>,
  <"8.8_W2_glue2_210.e6x46410.ck6",
   [-3741868748800,0,246368752516,0,1077715393,0,814016]>,
  <"6.12_V1_glue2_30.a6x630.h2",
   [88825,81950,75207,-13354,-6347,396,132]>,
  <"6.12_V2_x1fam_12(1/7)x12(-2)",
   [-15993254144,0,413375730009,0,-178475623260,0,19299418112]>
*];

for w in WITS do
    g := RQx!w[2];
    printf "== WITNESS %o ==\n", w[1];
    printf "g = %o\n", g;
    dz := Integers()!Discriminant(g);
    printf "disc(g) = %o (%o digits)\n", dz, #IntegerToString(AbsoluteValue(dz));
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
    // elliptic factors via Richelot -> SetCart
    C := HyperellipticCurve(g);
    okmin := false;
    try C := ReducedMinimalWeierstrassModel(C); okmin := true; catch e okmin := false; end try;
    L := [* *];
    try L := RichelotIsogenousSurfaces(C); catch e L := [* *]; end try;
    found := false;
    for X in L do
        if Type(X) eq SetCart then
            E1 := Component(X,1); E2 := Component(X,2);
            printf "FACTORS: E1 cond %o tors %o ainvs %o ; E2 cond %o tors %o ainvs %o\n",
                Conductor(E1), Invariants(TorsionSubgroup(E1)), aInvariants(MinimalModel(E1)),
                Conductor(E2), Invariants(TorsionSubgroup(E2)), aInvariants(MinimalModel(E2));
            found := true;
            break;
        end if;
    end for;
    if not found then printf "FACTORS: no rational Richelot to a product from this model\n"; end if;
end for;
printf "SEARCH_DONE verify\n";
quit;
