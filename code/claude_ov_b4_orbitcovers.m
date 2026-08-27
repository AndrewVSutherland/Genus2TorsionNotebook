// claude_ov_b4_orbitcovers.m -- Lane 4 (route B4): the two components of the
// RICHELOT-SOURCE cover of Flynn's t-line.
//
// claude_ov_b4_galois15.m showed Gal(F_t/Q(t)) = S3 wr C2 (order 72), acting
// imprimitively on the six Weierstrass points with two blocks of size 3, and
// that its orbits on the 15 (2,2,2)-partitions have sizes 6 (the "transverse"
// partitions: every pair meets both blocks) and 9 (one pair inside each block
// plus the transverse leftover).  So the locus of Flynn members carrying a
// Galois-stable Richelot kernel is a cover of the t-line with exactly two
// components, of degrees 6 and 9.
//
// This script (a) computes the degree of t as a map from the (6,2) component
// E = 92.a1 to the t-line -- identifying which orbit the whole overnight
// campaign lived on -- and (b) tries to realise the two orbit covers as
// function fields and compute their genera.
//
// Run: code/claude_magma_slot.sh -b MemGB:=8 code/claude_ov_b4_orbitcovers.m \
//        > results/claude_ov_b4_orbitcovers.log 2>&1 &

SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();

// ---- (a) degree of t on the (6,2) component -------------------------------
printf "=== (a) the (6,2) component E = 92.a1 -> t-line ===\n";
Fu<u> := FunctionField(Q);
PY<Y> := PolynomialRing(Fu);
K<yy> := ext<Fu | Y^2 - (-u^3+3*u^2-2*u+1)>;
uu := K!u;
a2f := (uu-2)^4;
b1f := 2*(uu^5-5*uu^4+9*uu^3-8*uu^2+4*uu-2);
vv := (-b1f + 4*(uu-1)*yy) / (2*a2f);
RT<Tt> := PolynomialRing(K);
RTX<xx> := PolynomialRing(RT);
Ftt := xx^6 + 2*xx^5 + (2*Tt+3)*xx^4 + 2*xx^3 + (Tt^2+1)*xx^2 + 2*Tt*(1-Tt)*xx + Tt^2;
_, rm := Quotrem(Ftt, xx^2 + (RT!uu)*xx + (RT!vv));
gg := GCD(Coefficient(rm,1), Coefficient(rm,0));
tK := -Coefficient(gg,0)/Coefficient(gg,1);
printf "t on E92 has DEGREE %o (as a map E92 -> P^1)\n", Degree(tK);
printf "  genus of E92 function field: %o\n", Genus(K);

// ---- (a') degree of t on the (4,3) component ------------------------------
Fw<w> := FunctionField(Q);
tW := -((w-1)*(w^2-w+1)/w^2)^2;
printf "t on the (4,3) component (w-line) has DEGREE %o\n", Degree(tW);

// ---- (b) the two orbit covers ---------------------------------------------
printf "\n=== (b) the two orbit covers of the t-line ===\n";
Ft<t> := FunctionField(Q);
Px<x> := PolynomialRing(Ft);
F := x^6 + 2*x^5 + (2*t+3)*x^4 + 2*x^3 + (t^2+1)*x^2 + 2*t*(1-t)*x + t^2;
G, R, S := GaloisGroup(F);
printf "Gal = %o of order %o\n", GroupName(G), #G;

S6 := Sym(6);
parts := [];
for a in [2..6] do
    rest1 := [i : i in [2..6] | i ne a];
    b := rest1[1];
    for c in [i : i in rest1 | i ne b] do
        rest2 := [i : i in rest1 | i notin {b, c}];
        Append(~parts, { {1,a}, {b,c}, {rest2[1], rest2[2]} });
    end for;
end for;
parts := SetToSequence(SequenceToSet(parts));
function ActP(P, g) return { { i^g : i in pr } : pr in P }; end function;

seen := {}; orbs := [];
for i in [1..#parts] do
    if i in seen then continue; end if;
    orb := {};
    for g in G do Include(~orb, Index(parts, ActP(parts[i], S6!g))); end for;
    seen := seen join orb;
    Append(~orbs, orb);
end for;
printf "orbit sizes %o\n", Sort([#o : o in orbs]);

// also the orbits on the 15 PAIRS (= components of the quadratic-factor
// incidence variety, which is what the whole campaign so far worked inside)
prs := [ {i,j} : j in [i+1..6], i in [1..6] ];
seenp := {}; orbp := [];
for i in [1..#prs] do
    if i in seenp then continue; end if;
    o := {};
    for g in G do Include(~o, Index(prs, { k^(S6!g) : k in prs[i] })); end for;
    seenp := seenp join o;
    Append(~orbp, o);
end for;
printf "orbit sizes on the 15 PAIRS (= incidence-variety components): %o\n", Sort([#o : o in orbp]);

act := map< car<{1..#parts}, G> -> {1..#parts} |
            y :-> Index(parts, ActP(parts[y[1]], S6!y[2])) >;
GS := GSet(G, {1..#parts}, act);

for o in orbs do
    i0 := Min(o);
    Stab := Stabiliser(G, GS, i0);
    printf "\nORBIT SIZE %o : stabiliser order %o (index %o), %o\n",
        #o, #Stab, Index(G, Stab), GroupName(Stab);
    printf "  partition: %o\n", parts[i0];
    ok := true;
    try
        h := GaloisSubgroup(S, Stab);
        printf "  DEFINING POLYNOMIAL of the cover (degree %o):\n    %o\n", Degree(h), h;
        printf "  irreducible over Q(t): %o\n", IsIrreducible(h);
        LL := ext<Ft | h>;
        gg := Genus(LL);
        printf "  GENUS of the cover = %o\n", gg;
        printf "  constant field = %o\n", ConstantField(LL);
        // clear denominators and analyse the plane model
        dd := LCM([Denominator(c) : c in Coefficients(h)]);
        P2c<Tv,Yv> := PolynomialRing(Q, 2);
        hb := &+[ P2c!Evaluate(Numerator(dd*Coefficient(h,i)), Tv) * Yv^i
                  : i in [0..Degree(h)] ];
        printf "  PLANE MODEL (t,Y): %o\n", hb;
        Cv := Curve(AffineSpace(P2c), hb);
        PCv := ProjectiveClosure(Cv);
        printf "  plane model genus %o, absolutely irreducible %o\n",
            Genus(PCv), IsAbsolutelyIrreducible(Cv);
        ptsv := PointSearch(PCv, 5000);
        printf "  PointSearch(5000): %o points: %o\n", #ptsv, ptsv;
        if gg eq 2 then
            bhy, Hy, mhy := IsHyperelliptic(PCv);
            printf "  IsHyperelliptic: %o\n", bhy;
            if bhy then
                Hy := SimplifiedModel(Hy);
                printf "  GENUS-2 MODEL: %o\n", Hy;
                printf "  discriminant %o\n", Factorization(Integers()!Discriminant(Hy));
                ptsH := Points(Hy : Bound := 100000);
                printf "  Points(bound=10^5) on the genus-2 model: %o -> %o\n", #ptsH, ptsH;
                Jy := Jacobian(Hy);
                rl, ru := RankBounds(Jy);
                printf "  JACOBIAN RANK BOUNDS: [%o, %o]\n", rl, ru;
                printf "  TORSION of J: %o\n", Invariants(TorsionSubgroup(Jy));
                if ru eq 0 then
                    printf "  RANK 0 -> Chabauty0\n";
                    ch := Chabauty0(Jy);
                    printf "  CHABAUTY0 RATIONAL POINTS: %o -> %o\n", #ch, ch;
                elif rl eq 1 and ru eq 1 then
                    printf "  RANK 1 -> searching for a generator for Chabauty\n";
                    bas := Points(Jy : Bound := 5000);
                    inf := [P : P in bas | Order(P) eq 0];
                    printf "  points of infinite order found: %o\n", #inf;
                    if #inf gt 0 then
                        ch := Chabauty(inf[1]);
                        printf "  CHABAUTY (may be for a finite-index subgroup): %o -> %o\n", #ch, ch;
                    end if;
                end if;
            end if;
        end if;
    catch e
        ok := false;
        printf "  GaloisSubgroup/analysis failed: %o\n", e`Object;
    end try;
end for;

printf "ORBITCOVERS_DONE\n";
quit;
