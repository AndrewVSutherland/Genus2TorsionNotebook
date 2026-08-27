//////////////////////////////////////////////////////////////////////
// Bounded search for either an independent rational 3-torsion class or
// a rational 9-torsion class on Lepr\'evost's direct Z/21 family.
//
// If J_t(Q) contains [3,21], then J_t(Q)[3] has rank two.  At every
// good prime p != 3 this injects into J_t(F_p)[3].  For each small p we
// therefore tabulate the residues t in F_p for which the finite
// Jacobian has 3-rank at least two.  Similarly, a cyclic order-63 class
// would force the finite 3-primary exponent to be at least 9.  A
// rational height search uses both necessary residue tests before any
// exact TorsionSubgroup call.
//
// Typical runs:
//   magma -b mode:="finite" prime_bound:=43 \
//       code/z21_extra3_residue_search.m
//   magma -b mode:="search" height:=1000 prime_bound:=43 max_exact:=20 \
//       code/z21_extra3_residue_search.m
//
// The script streams rational parameters and uses a 300 MB hard cap.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then mode := "search"; end if;
if not assigned height then height := 200; end if;
if not assigned prime_bound then prime_bound := 43; end if;
if not assigned max_exact then max_exact := 20; end if;
if not assigned max_print then max_print := 40; end if;
if not assigned progress_interval then progress_interval := 1000000; end if;
if not assigned MemMB then MemMB := 300; end if;

if Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if Type(max_print) eq MonStgElt then max_print := StringToInteger(max_print); end if;
if Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if Type(MemMB) eq MonStgElt then MemMB := StringToInteger(MemMB); end if;
SetMemoryLimit(MemMB*10^6);

Z := Integers();
Q := Rationals();
Qx<x> := PolynomialRing(Q);

function FamilyPolynomial(K, tt)
    P<X> := PolynomialRing(K);
    p2 := tt^14+4*tt^13+19*tt^12+32*tt^11+113*tt^10+188*tt^9
        +379*tt^8+448*tt^7+379*tt^6+188*tt^5+113*tt^4
        +32*tt^3+19*tt^2+4*tt+1;
    p1 := tt^10+4*tt^9+17*tt^8+24*tt^7+54*tt^6+56*tt^5
        +54*tt^4+24*tt^3+17*tt^2+4*tt+1;
    p0 := tt^6+4*tt^5+15*tt^4+16*tt^3+15*tt^2+4*tt+1;
    A := p2*X^2-2*(tt^2+1)^2*p1*X+(tt^2+1)^4*p0;
    k := 64*tt^4*(tt+1)^2*(tt^2+1)^3
       *(tt^4+2*tt^3+6*tt^2+2*tt+1)^3;
    return A^2-k*X^3*(X-1)^2;
end function;

function IsGoodModel(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function ThreeRank(invs)
    return #[n : n in invs | (Z!n) mod 3 eq 0];
end function;

function HasNineExponent(invs)
    return #invs gt 0 and (&or[(Z!n) mod 9 eq 0 : n in invs]);
end function;

function FiniteThreeRank(f)
    A, mp := AbelianGroup(Jacobian(HyperellipticCurve(f)));
    invs := Invariants(A);
    return ThreeRank(invs), invs;
end function;

function ResidueTable(p)
    K := GF(p);
    allowed := {};
    allowed9 := {};
    bad := {};
    counts := AssociativeArray();
    examples := AssociativeArray();
    for tt in K do
        t0 := Z!tt;
        f := FamilyPolynomial(K, tt);
        if not IsGoodModel(f) then
            Include(~bad, t0);
            continue;
        end if;
        r, invs := FiniteThreeRank(f);
        if IsDefined(counts, r) then
            counts[r] +:= 1;
        else
            counts[r] := 1;
            examples[r] := <t0, invs>;
        end if;
        if r ge 2 then Include(~allowed, t0); end if;
        if HasNineExponent(invs) then Include(~allowed9, t0); end if;
    end for;
    return allowed, allowed9, bad, counts, examples;
end function;

function ResidueOfFraction(a, b, p)
    if b mod p eq 0 then return false, 0; end if;
    K := GF(p);
    return true, Z!((K!a)/(K!b));
end function;

function IntegralSquareModel(f)
    L := Z!1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f,i)));
    end for;
    return Qx!(L^2*f), L;
end function;

primes := [p : p in PrimesUpTo(prime_bound) | p notin {2,3}];
allowed_by_p := AssociativeArray();
allowed9_by_p := AssociativeArray();
bad_by_p := AssociativeArray();

print "Z21_PLUS_INDEPENDENT_3_RESIDUE_SEARCH";
print "mode",mode,"height",height,"prime_bound",prime_bound,
      "primes",primes,"max_exact",max_exact,"memory_MB",MemMB;

for p in primes do
    allowed, allowed9, bad, counts, examples := ResidueTable(p);
    allowed_by_p[p] := allowed;
    allowed9_by_p[p] := allowed9;
    bad_by_p[p] := bad;
    print "PRIME",p,"allowed_rank_ge_2",#allowed,"bad",#bad,
          "allowed_exponent_9",#allowed9,
          "rank_counts",Sort([<r,counts[r]> : r in Keys(counts)]),
          "rank2_residues",Sort(Setseq(allowed)),
          "exponent9_residues",Sort(Setseq(allowed9));
end for;

if mode eq "finite" then
    print "DONE_FINITE_TABLES";
    quit;
end if;

checked := 0;
smooth := 0;
survivors := 0;
survivors9 := 0;
exact_tests := 0;
hits := [];
hits9 := [];
printed := 0;
kill_counts := AssociativeArray();
kill9_counts := AssociativeArray();
bad_skip_counts := AssociativeArray();

for b in [1..height] do
    for a in [-height..height] do
        if GCD(a,b) ne 1 or a in {0,-b} then continue; end if;
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "PROGRESS",checked,
                  "rank2_survivors",survivors,
                  "exponent9_survivors",survivors9,
                  "exact",exact_tests,"extra3_hits",#hits,
                  "order63_hits",#hits9;
        end if;

        killed := false;
        killed9 := false;
        usable := 0;
        skipped := [];
        for p in primes do
            ok, tt := ResidueOfFraction(a,b,p);
            if not ok or tt in bad_by_p[p] then
                Append(~skipped,p);
                if IsDefined(bad_skip_counts,p) then
                    bad_skip_counts[p] +:= 1;
                else
                    bad_skip_counts[p] := 1;
                end if;
                continue;
            end if;
            usable +:= 1;
            if not killed and tt notin allowed_by_p[p] then
                killed := true;
                if IsDefined(kill_counts,p) then
                    kill_counts[p] +:= 1;
                else
                    kill_counts[p] := 1;
                end if;
            end if;
            if not killed9 and tt notin allowed9_by_p[p] then
                killed9 := true;
                if IsDefined(kill9_counts,p) then
                    kill9_counts[p] +:= 1;
                else
                    kill9_counts[p] := 1;
                end if;
            end if;
            if killed and killed9 then break; end if;
        end for;
        if killed and killed9 then continue; end if;

        // A survivor has passed every available good-prime necessary test.
        ttQ := Q!a/Q!b;
        fQ := FamilyPolynomial(Q,ttQ);
        if not IsGoodModel(fQ) then continue; end if;
        smooth +:= 1;
        if not killed then survivors +:= 1; end if;
        if not killed9 then survivors9 +:= 1; end if;
        if printed lt max_print then
            print "SURVIVOR","t",ttQ,"usable_primes",usable,
                  "extra3",not killed,"nine",not killed9,
                  "skipped_primes",skipped;
            printed +:= 1;
        end if;

        if exact_tests ge max_exact then continue; end if;
        fI,L := IntegralSquareModel(fQ);
        G,mp := TorsionSubgroup(Jacobian(HyperellipticCurve(fI)));
        invs := Invariants(G);
        exact_tests +:= 1;
        print "EXACT","t",ttQ,"torsion",invs,"order",#G,
              "three_rank",ThreeRank(invs),
              "has_nine_exponent",HasNineExponent(invs),"y_scale",L;
        if ThreeRank(invs) ge 2 then
            Append(~hits,<ttQ,invs,fI>);
            print "INDEPENDENT_3_HIT","t",ttQ,"torsion",invs,
                  "integral_model",fI;
        end if;
        if HasNineExponent(invs) then
            Append(~hits9,<ttQ,invs,fI>);
            print "ORDER_63_HIT","t",ttQ,"torsion",invs,
                  "integral_model",fI;
        end if;
    end for;
end for;

print "DONE_SEARCH";
print "checked",checked,"smooth_survivors",smooth,
      "rank2_residue_survivors",survivors,
      "exponent9_residue_survivors",survivors9,
      "exact_tests",exact_tests,"extra3_hits",#hits,
      "order63_hits",#hits9;
print "KILL_COUNTS";
for p in Sort([p : p in Keys(kill_counts)]) do
    print p,kill_counts[p];
end for;
print "KILL9_COUNTS";
for p in Sort([p : p in Keys(kill9_counts)]) do
    print p,kill9_counts[p];
end for;
print "BAD_OR_DENOMINATOR_SKIP_COUNTS";
for p in Sort([p : p in Keys(bad_skip_counts)]) do
    print p,bad_skip_counts[p];
end for;
quit;
