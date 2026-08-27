// claude_ov_b4_galois15.m -- Lane 4 (route B4): SCOPE the Richelot-source
// problem over Flynn's order-11 family beyond the quadratic-factor incidence.
//
// A Richelot (2,2)-kernel over Q is a partition of the six Weierstrass points
// into three pairs that is stable under Gal(Qbar/Q) -- the three quadratics
// need NOT be individually rational.  The whole overnight Lane-4 campaign
// worked inside the "quadratic-factor incidence variety", i.e. inside the case
// where at least one of the three quadratics IS rational.  This script decides
// how much is left over: it computes G = Gal(F_t / Q(t)) and its orbits on the
// 15 partitions of the six roots into three pairs.  An orbit of size d is a
// degree-d component of the "Richelot-source" cover of the t-line; orbits of
// size 1 would mean every member is Richelot-able.
//
// Run: code/claude_magma_slot.sh -b MemGB:=8 code/claude_ov_b4_galois15.m \
//        > results/claude_ov_b4_galois15.log 2>&1 &

SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Ft<t> := FunctionField(Q);
Px<x> := PolynomialRing(Ft);
F := x^6 + 2*x^5 + (2*t+3)*x^4 + 2*x^3 + (t^2+1)*x^2 + 2*t*(1-t)*x + t^2;
printf "F irreducible over Q(t): %o\n", IsIrreducible(F);
printf "F factor degrees over Q(t): %o\n", [Degree(f[1]) : f in Factorization(F)];

G, R, S := GaloisGroup(F);
printf "GALOIS GROUP of F over Q(t): order %o, %o\n", #G, GroupName(G);
printf "  transitive: %o\n", IsTransitive(G);

// the 15 partitions of {1..6} into three pairs
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
printf "number of (2,2,2) partitions: %o\n", #parts;

function ActP(P, g)
    return { { i^g : i in pr } : pr in P };
end function;

orbs := [];
seen := {};
for i in [1..#parts] do
    if i in seen then continue; end if;
    orb := {};
    for g in G do
        P := ActP(parts[i], S6!g);
        j := Index(parts, P);
        Include(~orb, j);
    end for;
    seen := seen join orb;
    Append(~orbs, orb);
end for;
printf "ORBITS of Gal on the 15 partitions: sizes %o\n", Sort([#o : o in orbs]);
for o in orbs do
    printf "  orbit size %o : e.g. %o\n", #o, parts[Min(o)];
end for;
printf "STABLE PARTITIONS over Q(t) (orbit size 1): %o\n", #[o : o in orbs | #o eq 1];

// how many partitions have a RATIONAL quadratic in them, generically?
printf "\n(For reference: an orbit of size d is a degree-d component of the\n";
printf " Richelot-source cover of the t-line.  The quadratic-factor incidence\n";
printf " variety accounts only for partitions containing a rational pair.)\n";

printf "GALOIS15_DONE\n";
quit;
