/*
Mordell--Weil coefficient sieve on the known split first-cover fiber

    s = 59/49.

The quotient elliptic curve has rank two and torsion Z/4.  Its known split
point is -P1+P2.  Direct exact enumeration becomes expensive because the
rational coordinates grow quadratically in the coefficient bound.  This
script first requires the complementary cubic to split modulo several good
primes, then constructs exact points only for the surviving coefficient
tuples.

Usage:

  magma -b coefficient_bound:=200 \
    code/elkies22210_orbit12_seed_fiber_mw_sieve.m
*/

if not assigned coefficient_bound then
    coefficient_bound := 200;
elif Type(coefficient_bound) eq MonStgElt then
    coefficient_bound := StringToInteger(coefficient_bound);
end if;
require coefficient_bound ge 0 : "coefficient_bound must be nonnegative";

Q := Rationals();
s := Q!59/49;
E := EllipticCurve([Q | 0, s^2-2*s-2, 0,
                    -4*s^2*(s+1), 4*s^2*(s+1)^2]);

// Exact free generators returned by MordellWeilGroup/Generators.  Their
// independence and saturation are certified below by recomputing the full
// Mordell--Weil group and comparing its rank and torsion.
P1 := E![Q!864/1225,Q!182088/42875,1];
P2 := E![Q!-944/441,Q!-367688/64827,1];
T4 := E![Q!0,2*s*(s+1),1];
assert Order(T4) eq 4;

Emin, minmap := MinimalModel(E);
lo, hi := RankBounds(Emin);
assert lo eq 2 and hi eq 2;
MW, mwmap := MordellWeilGroup(Emin);
assert Invariants(MW) eq [4,0,0];

free_indices := [i : i in [1..#Invariants(MW)] |
                 Invariants(MW)[i] eq 0];
torsion_indices := [i : i in [1..#Invariants(MW)] |
                    Invariants(MW)[i] ne 0];
returned_free := [Inverse(minmap)(mwmap(MW.i)) : i in free_indices];
returned_torsion := [];
for k in [0..Invariants(MW)[torsion_indices[1]]-1] do
    coordinates := [0 : i in [1..#Invariants(MW)]];
    coordinates[torsion_indices[1]] := k;
    Append(~returned_torsion,Inverse(minmap)(mwmap(MW!coordinates)));
end for;

function FreeCoordinates(ep,free,torsion)
    for a,b in [-20..20] do
        for tt in torsion do
            if a*free[1]+b*free[2]+tt eq ep then
                return true,[a,b];
            end if;
        end for;
    end for;
    return false,[];
end function;

ok1,c1 := FreeCoordinates(P1,returned_free,returned_torsion);
ok2,c2 := FreeCoordinates(P2,returned_free,returned_torsion);
assert ok1 and ok2;
basis_change := Matrix(Integers(),2,2,c1 cat c2);
assert Abs(Determinant(basis_change)) eq 1;

// Positive control: -P1+P2 is the known q=8/7 point.
seed := -P1+P2;
assert Q!(2*s*(s+1)/seed[1]) eq Q!8/7;

function ReduceRational(a,F)
    return F!Numerator(a)/F!Denominator(a);
end function;

function ReducePoint(P,Ep,F)
    if IsIdentity(P) then return Ep!0; end if;
    return Ep![ReduceRational(P[1],F),ReduceRational(P[2],F),1];
end function;

function AllowedSplitQ(F)
    PF<z> := PolynomialRing(F);
    sf := ReduceRational(s,F);
    allowed := {};
    for q in F do
        f := z^3+(1+q)*z^2+sf*z-(1+q)*(q-sf);
        fac := Factorization(f);
        if &and[Degree(row[1]) eq 1 : row in fac] then
            Include(~allowed,q);
        end if;
    end for;
    return allowed;
end function;

function PassesPrime(ep,allowed,F)
    // If q=2s(s+1)/X is not integral at this prime, reduction of the monic
    // complementary cubic is not a safe obstruction, so retain the tuple.
    if IsIdentity(ep) or ep[1] eq 0 then return true; end if;
    sf := ReduceRational(s,F);
    qf := 2*sf*(sf+1)/ep[1];
    return qf in allowed;
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

// Avoid 2,3,5,7 (denominators of the model/generators), and automatically
// discard any accidental bad-reduction prime.
prime_candidates := [11,13,17,19,23,29,31,37,41,43,47,
                     53,59,61,67,71,73,79,83,89,97,
                     101,103,107,109,113,127,131,137,139,149,151,
                     157,163,167,173,179,181,191,193,197,199];
good_primes := [];
for p in prime_candidates do
    F := GF(p);
    discf := ReduceRational(256*s^2*(s+1)^4*(s^2+8*s+8),F);
    if discf eq 0 then continue; end if;
    Ep := EllipticCurve([F | 0,ReduceRational(s^2-2*s-2,F),0,
                         ReduceRational(-4*s^2*(s+1),F),
                         ReduceRational(4*s^2*(s+1)^2,F)]);
    assert Discriminant(Ep) ne 0;
    Append(~good_primes,p);
end for;

N := coefficient_bound;
candidates := [ <n1,n2,k> : n1,n2 in [-N..N], k in [0..3] ];
printf "ELKIES22210_ORBIT12_SEED_FIBER_MW_SIEVE\n";
printf "s %o\n",s;
printf "coefficient_bound %o\n",N;
printf "rank_bounds %o %o\n",lo,hi;
printf "mw_invariants %o\n",Invariants(MW);
printf "generators %o %o\n",P1,P2;
printf "basis_change_determinant %o\n",Determinant(basis_change);
printf "torsion_generator %o\n",T4;
printf "seed_coordinates [-1,1,0]\n";
printf "initial_coefficient_tuples %o\n",#candidates;

for p in good_primes do
    F := GF(p);
    Ep := EllipticCurve([F | 0,ReduceRational(s^2-2*s-2,F),0,
                         ReduceRational(-4*s^2*(s+1),F),
                         ReduceRational(4*s^2*(s+1)^2,F)]);
    p1 := ReducePoint(P1,Ep,F);
    p2 := ReducePoint(P2,Ep,F);
    tt := ReducePoint(T4,Ep,F);
    allowed := AllowedSplitQ(F);
    m1 := [n*p1 : n in [-N..N]];
    m2 := [n*p2 : n in [-N..N]];
    mt := [k*tt : k in [0..3]];
    next := [];
    for c in candidates do
        ep := m1[c[1]+N+1]+m2[c[2]+N+1]+mt[c[3]+1];
        if PassesPrime(ep,allowed,F) then Append(~next,c); end if;
    end for;
    candidates := next;
    printf "SIEVE p=%o allowed_q=%o survivors=%o\n",p,#allowed,#candidates;
end for;

if #candidates le 100 then
    printf "surviving_coefficient_tuples %o\n",candidates;
end if;

P<z> := PolynomialRing(Q);
seen_q := {};
finite_exact_q := 0;
exact_split := 0;
smooth_split := 0;
full_hits := 0;
max_mask_weight := 0;
for c in candidates do
    ep := c[1]*P1+c[2]*P2+c[3]*T4;
    if IsIdentity(ep) or ep[1] eq 0 then continue; end if;
    q := Q!(2*s*(s+1)/ep[1]);
    if q in seen_q then continue; end if;
    Include(~seen_q,q);
    finite_exact_q +:= 1;
    cubic := z^3+(1+q)*z^2+s*z-(1+q)*(q-s);
    fac := Factorization(cubic);
    if #fac ne 3 or
       not &and[Degree(row[1]) eq 1 and row[2] eq 1 : row in fac] then
        continue;
    end if;
    exact_split +:= 1;
    roots := [-Coefficient(row[1],0)/Coefficient(row[1],1) : row in fac];
    rr := [Q!1,q] cat roots;
    if not SmoothCK(rr) then continue; end if;
    smooth_split +:= 1;
    GG := LiteralRadicands(rr);
    mask := [IsSquare(g) : g in GG];
    max_mask_weight := Max(max_mask_weight,#[b : b in mask | b]);
    printf "SPLIT_LIFT coefficients=%o q=%o roots=%o mask=%o G=%o\n",
           c,q,roots,mask,GG;
    if &and mask then
        full_hits +:= 1;
        printf "HIT r=%o\n",rr;
    end if;
end for;

printf "modular_survivors %o\n",#candidates;
printf "distinct_finite_exact_q %o\n",finite_exact_q;
printf "cubics_split_completely %o\n",exact_split;
printf "smooth_split_lifts %o\n",smooth_split;
printf "max_exact_square_mask_weight %o\n",max_mask_weight;
printf "full_hits %o\n",full_hits;
printf "DONE\n";

quit;
