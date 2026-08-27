// claude_ov_b4_orbit6rank.m -- Lane 4 (route B4): settle the rational points of
// the GENUS-2 component of the Richelot-source locus of Flynn's order-11 family
// (the "D6"/transverse-partition orbit, degree 6 over the t-line).
//
// Model from claude_ov_b4_orbitcovers.m:
//   y^2 = 16x^6+3648x^5+346176x^4+17501312x^3+497180160x^2+7525173248x+47411310592
// equivalently Y^2 = g(x), Y = y/4,
//   g = x^6+228x^5+21636x^4+1093832x^3+31073760x^2+470323328x+2963206912.
// Points(bound 10^5) found 8 rational points; RankBounds returned [0,2] but the
// lower bound is only a search bound (differences of the 8 points are non-zero
// in J(Q) and J(Q)_tors is trivial, so the rank is >= 1).
//
// Goal: rank exactly, then Chabauty if rank = 1.
//
// Run: code/claude_magma_slot.sh -b MemGB:=24 code/claude_ov_b4_orbit6rank.m \
//        > results/claude_ov_b4_orbit6rank.log 2>&1 &

SetColumns(0);
if not assigned MemGB then MemGB := 8; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned PB then PB := 3000; elif Type(PB) eq MonStgElt then PB := StringToInteger(PB); end if;

Q := Rationals();
P<x> := PolynomialRing(Q);
g := x^6 + 228*x^5 + 21636*x^4 + 1093832*x^3 + 31073760*x^2 + 470323328*x + 2963206912;
H := HyperellipticCurve(g);
printf "H : y^2 = %o\n", g;
printf "genus %o  disc %o\n", Genus(H), Factorization(Integers()!Discriminant(H));
Hr, mr := ReducedModel(H);
printf "REDUCED MODEL: %o\n", Hr;
Hm, mm := ReducedMinimalWeierstrassModel(H);
printf "REDUCED MINIMAL WEIERSTRASS MODEL: %o\n", Hm;
printf "  disc %o\n", Factorization(Integers()!Discriminant(Hm));

pts := Points(Hm : Bound := 100000);
printf "Points(bound 10^5) on the reduced model: %o -> %o\n", #pts, pts;

J := Jacobian(Hm);
printf "TwoTorsionSubgroup invariants: %o\n", Invariants(TwoTorsionSubgroup(J));
printf "TorsionSubgroup: %o\n", Invariants(TorsionSubgroup(J));
rl, ru := RankBounds(J);
printf "RankBounds (default): [%o, %o]\n", rl, ru;
rb := RankBound(J);
printf "RankBound: %o\n", rb;

JP := Points(J : Bound := PB);
printf "Points(J : Bound := %o): %o\n", PB, #JP;
inf := [D : D in JP | Order(D) eq 0];
printf "points of infinite order: %o\n", #inf;
if #inf gt 0 then
    B := ReducedBasis(inf);
    printf "ReducedBasis size (= rank lower bound): %o\n", #B;
    for D in B do printf "  BASIS %o   height %o\n", D, Height(D); end for;
    if #B eq 1 and rb le 1 then
        printf "RANK = 1 -> CHABAUTY\n";
        ch := Chabauty(B[1]);
        printf "CHABAUTY SET (%o points): %o\n", #ch, ch;
        printf "index bound of <B[1]> in J(Q): rank bound %o\n", rb;
    elif #B eq 2 then
        printf "RANK >= 2 (= genus): classical Chabauty does NOT apply.\n";
    end if;
end if;

printf "ORBIT6RANK_DONE\n";
quit;
