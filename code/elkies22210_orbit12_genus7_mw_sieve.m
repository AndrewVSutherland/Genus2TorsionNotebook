/*
Mordell--Weil coefficient-box sieve for the genus-7 discriminant cover

  v^2 = G(q),  z^2 = D(q)

on the seed fibre s=59/49.  The first equation is the rank-two elliptic
quotient.  The second is imposed modulo good primes before exact
reconstruction.  The pinned generators are a certified full basis of the
free part of E(Q); see elkies22210_orbit12_seed_fiber_mw_sieve.m.

Usage:

  magma -b coefficient_bound:=500 \
    code/elkies22210_orbit12_genus7_mw_sieve.m
*/

if not assigned coefficient_bound then
    coefficient_bound := 500;
elif Type(coefficient_bound) eq MonStgElt then
    coefficient_bound := StringToInteger(coefficient_bound);
end if;
require coefficient_bound ge 0 : "coefficient_bound must be nonnegative";
if not assigned n1_low then
    n1_low := -coefficient_bound;
elif Type(n1_low) eq MonStgElt then
    n1_low := StringToInteger(n1_low);
end if;
if not assigned n1_high then
    n1_high := coefficient_bound;
elif Type(n1_high) eq MonStgElt then
    n1_high := StringToInteger(n1_high);
end if;
require -coefficient_bound le n1_low and n1_low le n1_high and
        n1_high le coefficient_bound : "invalid n1 interval";

Q := Rationals();
s := Q!59/49;
E := EllipticCurve([Q | 0,s^2-2*s-2,0,
                    -4*s^2*(s+1),4*s^2*(s+1)^2]);
P1 := E![Q!864/1225,Q!182088/42875,1];
P2 := E![Q!-944/441,Q!-367688/64827,1];
T4 := E![Q!0,2*s*(s+1),1];
assert Order(T4) eq 4;

Emin,minmap := MinimalModel(E);
lo,hi := RankBounds(Emin);
MW,mwmap := MordellWeilGroup(Emin);
assert lo eq 2 and hi eq 2;
assert Invariants(MW) eq [4,0,0];

free_indices := [i:i in [1..3]|Invariants(MW)[i] eq 0];
torsion_index := [i:i in [1..3]|Invariants(MW)[i] ne 0][1];
returned_free := [Inverse(minmap)(mwmap(MW.i)):i in free_indices];
returned_torsion := [
    Inverse(minmap)(mwmap(MW![k*(i eq torsion_index select 1 else 0):i in [1..3]]))
    : k in [0..3]
];

function FreeCoordinates(ep,free,torsion)
    for a,b in [-20..20] do
        for tt in torsion do
            if a*free[1]+b*free[2]+tt eq ep then return true,[a,b]; end if;
        end for;
    end for;
    return false,[];
end function;

ok1,c1 := FreeCoordinates(P1,returned_free,returned_torsion);
ok2,c2 := FreeCoordinates(P2,returned_free,returned_torsion);
assert ok1 and ok2;
basis_change := Matrix(Integers(),2,2,c1 cat c2);
assert Abs(Determinant(basis_change)) eq 1;
assert Q!(2*s*(s+1)/(-P1+P2)[1]) eq Q!8/7;

function ReduceRational(a,F)
    return F!Numerator(a)/F!Denominator(a);
end function;

function ReducePoint(P,Ep,F)
    if IsIdentity(P) then return Ep!0; end if;
    return Ep![ReduceRational(P[1],F),ReduceRational(P[2],F),1];
end function;

function DiscPolynomial(q,sf)
    return (1+q)^2*sf^2-4*sf^3+4*(1+q)^4*(q-sf)
           -27*(1+q)^2*(q-sf)^2-18*(1+q)^2*sf*(q-sf);
end function;

function PassesDiscriminant(ep,allowed_zero,allowed_q,sf)
    if IsIdentity(ep) then
        // q=a/X reduces to zero when X has a pole.
        return allowed_zero;
    end if;
    if ep[1] eq 0 then
        // q=infinity.  D has square leading coefficient 4.
        return true;
    end if;
    qf := 2*sf*(sf+1)/ep[1];
    return qf in allowed_q;
end function;

function IsRationalSquare(a)
    return a ge 0 and IsSquare(Numerator(a)) and IsSquare(Denominator(a));
end function;

function LiteralRadicands(rr)
    aa := [x^2:x in rr];
    return [
        -(aa[1]-aa[3])*(aa[1]-aa[4])*(aa[1]-aa[5]),
         (aa[3]-aa[2])*(aa[1]-aa[4])*(aa[1]-aa[5]),
         (aa[4]-aa[2])*(aa[1]-aa[3])*(aa[1]-aa[5]),
         (aa[5]-aa[2])*(aa[1]-aa[3])*(aa[1]-aa[4])
    ];
end function;

function SmoothCK(rr)
    if &or[x eq 0:x in rr] then return false; end if;
    if &+[x:x in rr] ne 0 or &+[x^3:x in rr] ne 0 then return false; end if;
    aa := [x^2:x in rr];
    return &and[aa[i] ne aa[j]:i,j in [1..5]|i lt j];
end function;

prime_candidates := [11,13,17,19,23,29,31,37,41,43,47,
                     53,59,61,67,71,73,79,83,89,97,
                     101,103,107,109,113,127,131,137,139,149,151,
                     157,163,167,173,179,181,191,193,197,199,
                     211,223,227,229,233,239,241,251,257,263,269,271,
                     277,281,283,293,307,311,313,317,331,337,347,349];
good_primes := [];
for p in prime_candidates do
    F := GF(p);
    if p eq 7 then continue; end if;
    discf := ReduceRational(256*s^2*(s+1)^4*(s^2+8*s+8),F);
    if discf ne 0 then Append(~good_primes,p); end if;
end for;

N := coefficient_bound;
W := 2*N+1;
initial_count := 4*(n1_high-n1_low+1)*W;
candidates := [];
printf "ELKIES22210_ORBIT12_GENUS7_MW_SIEVE\n";
printf "s %o\n",s;
printf "coefficient_bound %o\n",N;
printf "n1_interval %o %o\n",n1_low,n1_high;
printf "rank_bounds %o %o\n",lo,hi;
printf "mw_invariants %o\n",Invariants(MW);
printf "basis_change_determinant %o\n",Determinant(basis_change);
printf "seed_coordinates [-1,1,0]\n";
printf "initial_coefficient_tuples %o\n",initial_count;

// Applying the first primes simultaneously avoids materialising the very
// large intermediate list (about half the box after just one prime).
batch_size := Min(12,#good_primes);
filters := [* *];
for prime_index in [1..batch_size] do
    p := good_primes[prime_index];
    F := GF(p);
    sf := ReduceRational(s,F);
    allowed_q := {q:q in F|IsSquare(DiscPolynomial(q,sf))};
    allowed_zero := F!0 in allowed_q;
    Ep := EllipticCurve([F|0,ReduceRational(s^2-2*s-2,F),0,
                        ReduceRational(-4*s^2*(s+1),F),
                        ReduceRational(4*s^2*(s+1)^2,F)]);
    p1 := ReducePoint(P1,Ep,F);
    p2 := ReducePoint(P2,Ep,F);
    tt := ReducePoint(T4,Ep,F);
    m1 := [n*p1:n in [-N..N]];
    m2 := [n*p2:n in [-N..N]];
    mt := [k*tt:k in [0..3]];
    Append(~filters,<sf,allowed_zero,allowed_q,m1,m2,mt>);
end for;

for n1 in [n1_low..n1_high] do
  for n2 in [-N..N] do
    for k in [0..3] do
        keep := true;
        for row in filters do
            ep := row[4][n1+N+1]+row[5][n2+N+1]+row[6][k+1];
            if not PassesDiscriminant(ep,row[2],row[3],row[1]) then
                keep := false; break;
            end if;
        end for;
        if keep then Append(~candidates,((n1+N)*W+(n2+N))*4+k); end if;
    end for;
  end for;
end for;
printf "SIEVE_BATCH primes=%o survivors=%o\n",
       good_primes[1..batch_size],#candidates;

for prime_index in [batch_size+1..#good_primes] do
    p := good_primes[prime_index];
    F := GF(p);
    sf := ReduceRational(s,F);
    allowed_q := {q:q in F|IsSquare(DiscPolynomial(q,sf))};
    allowed_zero := F!0 in allowed_q;
    Ep := EllipticCurve([F|0,ReduceRational(s^2-2*s-2,F),0,
                        ReduceRational(-4*s^2*(s+1),F),
                        ReduceRational(4*s^2*(s+1)^2,F)]);
    p1 := ReducePoint(P1,Ep,F);
    p2 := ReducePoint(P2,Ep,F);
    tt := ReducePoint(T4,Ep,F);
    m1 := [n*p1:n in [-N..N]];
    m2 := [n*p2:n in [-N..N]];
    mt := [k*tt:k in [0..3]];
    next := [Integers()|];
    for code in candidates do
        k := code mod 4;
        ij := code div 4;
        n2 := (ij mod W)-N;
        n1 := (ij div W)-N;
        ep := m1[n1+N+1]+m2[n2+N+1]+mt[k+1];
        if PassesDiscriminant(ep,allowed_zero,allowed_q,sf) then
            Append(~next,code);
        end if;
    end for;
    candidates := next;
    printf "SIEVE p=%o survivors=%o\n",p,#candidates;
end for;

if #candidates le 200 then
    decoded := [];
    for code in candidates do
        k := code mod 4; ij := code div 4;
        Append(~decoded,<(ij div W)-N,(ij mod W)-N,k>);
    end for;
    printf "surviving_coefficient_tuples %o\n",decoded;
end if;

P<z> := PolynomialRing(Q);
seen_q := {};
finite_exact_q := 0;
exact_D_squares := 0;
exact_split := 0;
smooth_split := 0;
full_hits := 0;
for code in candidates do
    k := code mod 4; ij := code div 4;
    c := <(ij div W)-N,(ij mod W)-N,k>;
    ep := c[1]*P1+c[2]*P2+c[3]*T4;
    if IsIdentity(ep) then
        q := Q!0;
    elif ep[1] eq 0 then
        continue; // q=infinity, a boundary point of the CK chart.
    else
        q := Q!(2*s*(s+1)/ep[1]);
    end if;
    if q in seen_q then continue; end if;
    Include(~seen_q,q);
    finite_exact_q +:= 1;
    D := DiscPolynomial(q,s);
    if not IsRationalSquare(D) then continue; end if;
    exact_D_squares +:= 1;
    cubic := z^3+(1+q)*z^2+s*z-(1+q)*(q-s);
    fac := Factorization(cubic);
    if #fac ne 3 or not &and[Degree(a[1]) eq 1 and a[2] eq 1:a in fac] then
        printf "DISCRIMINANT_ONLY coefficients=%o q=%o D=%o\n",c,q,D;
        continue;
    end if;
    exact_split +:= 1;
    roots := [-Coefficient(a[1],0)/Coefficient(a[1],1):a in fac];
    rr := [Q!1,q] cat roots;
    if not SmoothCK(rr) then
        printf "SPLIT_BOUNDARY coefficients=%o q=%o roots=%o\n",c,q,roots;
        continue;
    end if;
    smooth_split +:= 1;
    GG := LiteralRadicands(rr);
    mask := [IsRationalSquare(a):a in GG];
    printf "SPLIT_LIFT coefficients=%o q=%o roots=%o mask=%o G=%o\n",
           c,q,roots,mask,GG;
    if &and mask then full_hits +:= 1; printf "HIT r=%o\n",rr; end if;
end for;

printf "modular_survivors %o\n",#candidates;
printf "distinct_finite_exact_q %o\n",finite_exact_q;
printf "exact_discriminant_squares %o\n",exact_D_squares;
printf "exact_cubics_split_completely %o\n",exact_split;
printf "smooth_split_lifts %o\n",smooth_split;
printf "full_hits %o\n",full_hits;
printf "DONE\n";
quit;
