// verify_pp_curves.m — independent verification of the Platonov-Petrunin
// order-36/48 curves as printed in Platonov, Russian Math. Surveys 69:1
// (2014), Section 6 (f36, f48,1, f48,2), and identification against the
// LMFDB alpha database split [36]/[48] rows.
SetColumns(0);
SetMemoryLimit(3*10^9);
load "split_lab.m";  // run from product/code/

CURVES := [*
  <"f36",  (x+2)*(3*x^2-6*x+8)*(3*x^3+6*x^2+3*x+4)> where x := RQx.1,
  <"f48_1",(x-2)*(x+2)*(x^4-10*x^2-3)> where x := RQx.1,
  <"f48_2",(x^2+x+2)*(x^4+3*x^3+21*x^2-27*x+18)> where x := RQx.1
*];

// DB split [36]/[48] rows: <label, fcoeffs, hcoeffs>
DB := [*
  <"4860.f1",  [9,-9,-1,-5,5,-1,1],   [0,1,1]>,
  <"24300.k1", [9,-9,-1,10,-7,1,1],   [0,1,1]>,
  <"170100.i1",[64,-152,98,-122,53,-27,9],[0,1,1]>,
  <"1764.a1",  [3,0,9,0,-4],          [0,1,0,1]>,
  <"2700.d2",  [5,0,7,0,10,0,2],      [0,1,0,1]>,
  <"5292.c2",  [1,1,5,10,10,6,2],     [0,1,1]>,
  <"5292.c4",  [0,16,8,-15,-5,3,1],   [0,1,1]>,
  <"31212.j1", [68,34,-34,-2,3,-3,1], [0,1,1]>
*];
dbinv := [* *];
for r in DB do
    Append(~dbinv, <r[1], G2Invariants(HyperellipticCurve(RQx!r[2], RQx!r[3]))>);
end for;

for cu in CURVES do
    g := cu[2];
    printf "== %o : g = %o ==\n", cu[1], g;
    T := ExactTorsion(g);
    printf "EXACT TORSION = %o\n", T;
    red, n := LpReducibleCount(g, 200);
    printf "L-poly reducible at %o/%o good primes < 200\n", red, n;
    gi := G2Invariants(HyperellipticCurve(g));
    match := "none (outside DB or different class)";
    for r in dbinv do
        if r[2] eq gi then match := r[1]; break; end if;
    end for;
    printf "G2-invariant match in DB split [36]/[48] rows: %o\n", match;
    C := HyperellipticCurve(g);
    okm := true;
    try C := ReducedMinimalWeierstrassModel(C); catch e okm := false; end try;
    L := [* *];
    try L := RichelotIsogenousSurfaces(C); catch e L := [* *]; end try;
    for X in L do
        if Type(X) eq SetCart then
            E1 := Component(X,1); E2 := Component(X,2);
            printf "FACTORS: cond %o tors %o  x  cond %o tors %o\n",
                Conductor(E1), Invariants(TorsionSubgroup(E1)),
                Conductor(E2), Invariants(TorsionSubgroup(E2));
            break;
        end if;
    end for;
end for;
printf "DONE\n";
quit;
