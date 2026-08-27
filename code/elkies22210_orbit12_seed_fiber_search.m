/*
Special-fiber search for the first orbit-12 radicand quotient at

    s = 59/49.

The fiber contains the known split lift

    q = 8/7,  (r3,r4,r5) = (-9/7,-5/7,-1/7),

for which only the first of the four literal halving radicands is a
square.  This script computes the exact Mordell--Weil group of the
elliptic quotient, identifies the seed in those coordinates, and
searches a symmetric Mordell--Weil coefficient box (including every
torsion translate).  Every completely split lift is checked directly
against the CK equations, smoothness, and all four literal radicands.

Usage:

  magma -b coefficient_bound:=250 \
    code/elkies22210_orbit12_seed_fiber_search.m
*/

if not assigned coefficient_bound then
    coefficient_bound := 250;
elif Type(coefficient_bound) eq MonStgElt then
    coefficient_bound := StringToInteger(coefficient_bound);
end if;
require coefficient_bound ge 0 : "coefficient_bound must be nonnegative";

Q := Rationals();
P<z> := PolynomialRing(Q);
s := Q!59/49;

E := EllipticCurve([Q | 0, s^2-2*s-2, 0,
                    -4*s^2*(s+1), 4*s^2*(s+1)^2]);
Emin, minmap := MinimalModel(E);
lo, hi := RankBounds(Emin);
require lo eq hi : "rank is not proved";
MW, mwmap := MordellWeilGroup(Emin);
mw_invariants := Invariants(MW);
free_indices := [i : i in [1..#mw_invariants] |
                 mw_invariants[i] eq 0];
torsion_indices := [i : i in [1..#mw_invariants] |
                    mw_invariants[i] ne 0];
require #free_indices eq lo : "MW presentation and proved rank disagree";
returned_free_min := [mwmap(MW.i) : i in free_indices];
returned_free := [Inverse(minmap)(p) : p in returned_free_min];

torsion_vectors := [[]];
for i in torsion_indices do
    torsion_vectors := [v cat [a] : v in torsion_vectors,
                                    a in [0..mw_invariants[i]-1]];
end for;
returned_torsion := [];
for tv in torsion_vectors do
    coords := [0 : i in [1..#mw_invariants]];
    for j in [1..#torsion_indices] do
        coords[torsion_indices[j]] := tv[j];
    end for;
    Append(~returned_torsion,Inverse(minmap)(mwmap(MW!coords)));
end for;

// Pin a reproducible Mordell--Weil basis.  MordellWeilGroup may return
// different LLL-reduced bases in different processes, so certify these
// two fixed points against the full group it computes, then use them for
// the search box.
fixed_free := [
    E![Q!-288/343,Q!101160/16807,1],
    E![Q!-944/441,Q!367688/64827,1]
];
fixed_T := E![Q!0,Q!12744/2401,1];
assert Order(fixed_T) eq 4;

function ReturnedCoordinates(ep, returned_free, returned_torsion)
    for n1 in [-20..20] do
        for n2 in [-20..20] do
            for ti in [1..#returned_torsion] do
                if n1*returned_free[1]+n2*returned_free[2]+
                   returned_torsion[ti] eq ep then
                    return true,[n1,n2],ti;
                end if;
            end for;
        end for;
    end for;
    return false,[],0;
end function;

ok1, c1, ti1 := ReturnedCoordinates(fixed_free[1],returned_free,
                                     returned_torsion);
ok2, c2, ti2 := ReturnedCoordinates(fixed_free[2],returned_free,
                                     returned_torsion);
require ok1 and ok2 : "fixed basis not located in returned MW basis";
basis_change := Matrix(Integers(),2,2,c1 cat c2);
require Abs(Determinant(basis_change)) eq 1 :
        "fixed points do not project to a basis of the free quotient";

free := fixed_free;
free_min := [minmap(p) : p in free];
torsion := [i*fixed_T : i in [0..3]];

function IsRationalSquare(a)
    return IsSquare(a);
end function;

function LiteralRadicands(rr)
    aa := [x^2 : x in rr];
    return [
        -(aa[1]-aa[3])*(aa[1]-aa[4])*(aa[1]-aa[5]),
         (aa[3]-aa[2])*(aa[1]-aa[4])*(aa[1]-aa[5]),
         (aa[4]-aa[2])*(aa[1]-aa[3])*(aa[1]-aa[5]),
         (aa[5]-aa[2])*(aa[1]-aa[3])*(aa[1]-aa[4])
    ];
end function;

function SmoothCK(rr)
    if &or[x eq 0 : x in rr] then return false; end if;
    if &+[x : x in rr] ne 0 or &+[x^3 : x in rr] ne 0 then
        return false;
    end if;
    aa := [x^2 : x in rr];
    return &and[aa[i] ne aa[j] : i,j in [1..5] | i lt j];
end function;

function SplitRoots(q)
    cubic := z^3+(1+q)*z^2+s*z-(1+q)*(q-s);
    fac := Factorization(cubic);
    if #fac ne 3 or
       not &and[Degree(row[1]) eq 1 and row[2] eq 1 : row in fac] then
        return false, [];
    end if;
    roots := [-Coefficient(row[1],0)/Coefficient(row[1],1) : row in fac];
    return true, roots;
end function;

// Known split-cubic regression point.  The sign of v chooses one of
// the two points over q; both have the same splitting/radicand data.
qseed := Q!8/7;
vseed := Q!192/343;
Xseed := 2*s*(s+1)/qseed;
Yseed := 2*s*(s+1)*vseed/qseed^2;
seed := E![Xseed,Yseed,1];
seed_min := minmap(seed);
seed_roots := [Q!-9/7,Q!-5/7,Q!-1/7];
assert &and[Evaluate(z^3+(1+qseed)*z^2+s*z-
                     (1+qseed)*(qseed-s), a) eq 0 : a in seed_roots];
seed_G := LiteralRadicands([Q!1,qseed] cat seed_roots);
assert [IsRationalSquare(g) : g in seed_G] eq
       [true,false,false,false];

printf "ELKIES22210_ORBIT12_SEED_FIBER_SEARCH\n";
printf "s %o\n", s;
printf "coefficient_bound %o\n", coefficient_bound;
printf "E %o\n", E;
printf "minimal_model %o\n", Emin;
printf "rank_bounds %o %o\n", lo, hi;
printf "mw_invariants %o\n", mw_invariants;
printf "torsion_invariants %o\n",
       [mw_invariants[i] : i in torsion_indices];
printf "minimal_generators %o\n", free_min;
printf "generators_on_quotient_model %o\n", free;
printf "fixed_basis_in_returned_free_coordinates %o %o\n", c1,c2;
printf "fixed_basis_change_determinant %o\n", Determinant(basis_change);
printf "fixed_torsion_generator %o\n", fixed_T;
printf "seed_on_quotient_model %o\n", seed;
printf "seed_on_minimal_model %o\n", seed_min;
printf "seed_roots %o\n", seed_roots;
printf "seed_radicands %o\n", seed_G;

// Identify the seed in the basis returned by Generators.  Since the
// rank is small, a modest exact coefficient box is enough to find it.
seed_coordinates_found := false;
seed_coordinates := [];
identify_bound := 20;
if #free eq 1 then
    for n in [-identify_bound..identify_bound] do
        for ti in [1..#torsion] do
            if n*free[1]+torsion[ti] eq seed then
                seed_coordinates_found := true;
                seed_coordinates := [n,ti];
            end if;
        end for;
    end for;
elif #free eq 2 then
    for n1 in [-identify_bound..identify_bound] do
        for n2 in [-identify_bound..identify_bound] do
            for ti in [1..#torsion] do
                if n1*free[1]+n2*free[2]+torsion[ti] eq seed then
                    seed_coordinates_found := true;
                    seed_coordinates := [n1,n2,ti];
                end if;
            end for;
        end for;
    end for;
end if;
printf "seed_coordinates_found %o\n", seed_coordinates_found;
if seed_coordinates_found then
    printf "seed_coordinates_free_then_torsion_index %o\n", seed_coordinates;
    seed_fixed_coordinates :=
        [seed_coordinates[#free+1]-1] cat seed_coordinates[1..#free];
    printf "seed_coordinates_fixed_T_P1_P2 %o\n", seed_fixed_coordinates;
end if;

require #free in [1,2] : "search loop supports rank one or two";

function LexicographicallyAtMost(a,b)
    for i in [1..#a] do
        if a[i] lt b[i] then return true; end if;
        if a[i] gt b[i] then return false; end if;
    end for;
    return true;
end function;

// Encode the fixed-basis group element n1*P1+n2*P2+k*T by an integer.
// First retain one of P and -P.  Then apply only necessary local tests:
// if the complementary cubic splits over Q, at every good prime where
// q has finite reduction it must split completely over F_p.
B := coefficient_bound;
width := 2*B+1;
mw_combinations := width^2*4;
candidate_ids := {};
for n1 in [-B..B] do
    for n2 in [-B..B] do
        for k in [0..3] do
            if n1 eq 0 and n2 eq 0 and k eq 0 then continue; end if;
            coord := [n1,n2,k];
            negcoord := [-n1,-n2,(-k) mod 4];
            if not LexicographicallyAtMost(coord,negcoord) then continue; end if;
            id := ((n1+B)*width+(n2+B))*4+k;
            Include(~candidate_ids,id);
        end for;
    end for;
end for;
sign_representatives := #candidate_ids;

sieve_primes := [11,13,17,19,23,29,31,37,41,43,
                 47,53,61,67,71,73,79,83];
allowed_ids := candidate_ids;
printf "sieve_primes %o\n", sieve_primes;
printf "sign_representatives_before_sieve %o\n", sign_representatives;
for p in sieve_primes do
    F := GF(p);
    sp := F!59/F!49;
    Ep := EllipticCurve([F | 0,sp^2-2*sp-2,0,
                        -4*sp^2*(sp+1),4*sp^2*(sp+1)^2]);
    assert Discriminant(Ep) ne 0;
    P1p := Ep![F!free[1][1],F!free[1][2],1];
    P2p := Ep![F!free[2][1],F!free[2][2],1];
    Tp := Ep![F!fixed_T[1],F!fixed_T[2],1];
    P1multiples := [n*P1p : n in [-B..B]];
    P2multiples := [n*P2p : n in [-B..B]];
    Tmultiples := [k*Tp : k in [0..3]];
    Rp<zp> := PolynomialRing(F);
    next_ids := {};
    for id in allowed_ids do
        k := id mod 4;
        u := id div 4;
        n2 := (u mod width)-B;
        n1 := (u div width)-B;
        epp := P1multiples[n1+B+1]+P2multiples[n2+B+1]+
               Tmultiples[k+1];
        // At O or X=0 the rational function q has indeterminate or
        // infinite reduction, so this prime is inconclusive.
        if IsIdentity(epp) or epp[1] eq 0 then
            Include(~next_ids,id);
            continue;
        end if;
        qp := 2*sp*(sp+1)/epp[1];
        cubicp := zp^3+(1+qp)*zp^2+sp*zp-(1+qp)*(qp-sp);
        rootsp := Roots(cubicp);
        if #rootsp gt 0 and &+[row[2] : row in rootsp] eq 3 then
            Include(~next_ids,id);
        end if;
    end for;
    allowed_ids := next_ids;
    printf "sieve_after_p%o %o\n", p,#allowed_ids;
end for;

// Positive controls: the boundary q=s and the known smooth seed must
// survive every local split-cubic sieve.
boundary_id := ((0+B)*width+(0+B))*4+2;
assert boundary_id in allowed_ids;
if B ge 1 then
    // Lexicographic sign convention retains -seed.
    seed_id := ((-1+B)*width+(1+B))*4+2;
    assert seed_id in allowed_ids;
end if;

exact_points_constructed := 0;
finite_q_points := 0;
cubic_discriminant_squares := 0;
cubics_split := 0;
smooth_split_lifts := 0;
full_hits := 0;
max_square_mask_weight := 0;
allowed_sequence := Setseq(allowed_ids);
Sort(~allowed_sequence);
for id in allowed_sequence do
    k := id mod 4;
    u := id div 4;
    n2 := (u mod width)-B;
    n1 := (u div width)-B;
    coeffs := [n1,n2];
    ti := k+1;
    ep := n1*free[1]+n2*free[2]+torsion[ti];
    exact_points_constructed +:= 1;
    if IsIdentity(ep) or ep[1] eq 0 then continue; end if;
    finite_q_points +:= 1;
    q := Q!(2*s*(s+1)/ep[1]);

    cubic := z^3+(1+q)*z^2+s*z-(1+q)*(q-s);
    if not IsSquare(Discriminant(cubic)) then continue; end if;
    cubic_discriminant_squares +:= 1;
    split, roots := SplitRoots(q);
    if not split then continue; end if;
    cubics_split +:= 1;

    rr := [Q!1,q] cat roots;
    if not SmoothCK(rr) then
        printf "SPLIT_BOUNDARY coefficients=%o torsion_index=%o q=%o roots=%o\n",
               coeffs,ti,q,roots;
        continue;
    end if;
    smooth_split_lifts +:= 1;
    GG := LiteralRadicands(rr);
    sq := [IsRationalSquare(g) : g in GG];
    assert sq[1];
    weight := #[b : b in sq | b];
    max_square_mask_weight := Max(max_square_mask_weight,weight);
    printf "SPLIT_LIFT coefficients=%o torsion_index=%o q=%o roots=%o mask=%o G=%o\n",
           coeffs,ti,q,roots,sq,GG;
    if &and sq then
        full_hits +:= 1;
        printf "HIT r=%o\n",rr;
    end if;
end for;

printf "mw_combinations %o\n", mw_combinations;
printf "distinct_elliptic_points %o\n", mw_combinations;
printf "sign_representatives %o\n", sign_representatives;
printf "modular_sieve_survivors %o\n", #allowed_ids;
printf "exact_points_constructed %o\n", exact_points_constructed;
printf "finite_q_points %o\n", finite_q_points;
printf "cubic_discriminant_squares %o\n", cubic_discriminant_squares;
printf "cubics_split_completely %o\n", cubics_split;
printf "smooth_split_lifts %o\n", smooth_split_lifts;
printf "max_exact_square_mask_weight %o\n", max_square_mask_weight;
printf "full_hits %o\n", full_hits;
printf "DONE\n";

quit;
