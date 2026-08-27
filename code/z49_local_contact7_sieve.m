//////////////////////////////////////////////////////////////////////
//  Local exponent-49 sieve on the two-parameter contact-7 family.
//
//      h = 1 - (7/2)*x + a*x^2 + b*x^3,
//      f = (h^2 + (x - 1)^7)/x^2.
//
//  A rational point of order 49 injects into J(F_p) at every prime
//  p != 7 of good reduction.  Consequently the EXPONENT of J(F_p),
//  not merely its order, must be divisible by 49.  Bad, singular, and
//  nonintegral displayed reductions are retained as unresolved boundary
//  cases; only good reductions with exponent not divisible by 49 kill a
//  rational parameter pair.
//
//  Typical runs:
//      magma -b mode:=finite prime_bound:=23 \
//          code/z49_local_contact7_sieve.m
//      magma -b mode:=search height:=12 prime_bound:=23 \
//          code/z49_local_contact7_sieve.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "finite";
end if;
if not assigned height then
    height := 12;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 23;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_exact then
    max_exact := 200;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned progress_interval then
    progress_interval := 10000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Z := Integers();
Q := Rationals();
P<x> := PolynomialRing(Q);
filter_primes := [ p : p in PrimesUpTo(prime_bound) | p ge 3 and p ne 7 ];

function Increment(counts, key)
    if IsDefined(counts, key) then
        counts[key] +:= 1;
    else
        counts[key] := 1;
    end if;
    return counts;
end function;

function HasExponent49(invs)
    return exists{ n : n in invs | (Z!n) mod 49 eq 0 };
end function;

// Sanity check: order divisibility alone can be a false positive.
assert not HasExponent49([ 7, 7 ]);
assert HasExponent49([ 49 ]);
assert HasExponent49([ 2, 98 ]);

function Contact7PolynomialFinite(F, a, b)
    PF<X> := PolynomialRing(F);
    h := 1 - (F!7/F!2)*X + a*X^2 + b*X^3;
    num := h^2 + (X - 1)^7;
    if Coefficient(num, 0) ne 0 or Coefficient(num, 1) ne 0 then
        return false, PF!0, PF!0, "identity_failure";
    end if;
    f := ExactQuotient(num, X^2);
    if Degree(f) ne 5 then
        return false, f, h, "degree_drop";
    end if;
    if Discriminant(f) eq 0 then
        return false, f, h, "singular";
    end if;
    return true, f, h, "good";
end function;

function Contact7Polynomial(a, b)
    h := 1 - (Q!7/2)*x + a*x^2 + b*x^3;
    num := h^2 + (x - 1)^7;
    if Coefficient(num, 0) ne 0 or Coefficient(num, 1) ne 0 then
        return false, P!0, P!0;
    end if;
    return true, ExactQuotient(num, x^2), h;
end function;

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            q := Q!num/Q!den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function ResidueIndex(q, p)
    d := Denominator(q);
    if d mod p eq 0 then
        return false, 0;
    end if;
    F := GF(p);
    r := F!(Numerator(q) mod p) / F!(d mod p);
    return true, Z!r;
end function;

// The mask is true exactly for a good affine residue whose finite
// Jacobian exponent is NOT divisible by 49.  Such a residue rigorously
// kills a rational order-49 specialization.  All boundary residues are
// false and hence retained.
function BuildKillMask(p)
    F := GF(p);
    mask := [ false : i in [1..p^2] ];
    total := 0;
    good := 0;
    pass_order := 0;
    pass_exponent := 0;
    false_order_positives := 0;
    boundary := 0;
    invariant_counts := AssociativeArray();
    open_pass := [];

    for ai in [0..p-1] do
        for bi in [0..p-1] do
            total +:= 1;
            a := F!ai;
            b := F!bi;
            ok, f, h, reason := Contact7PolynomialFinite(F, a, b);
            idx := ai*p + bi + 1;
            if not ok then
                boundary +:= 1;
                continue;
            end if;

            good +:= 1;
            C := HyperellipticCurve(f);
            J := Jacobian(C);
            A, phi := AbelianGroup(J);
            invs := [ Z!n : n in Invariants(A) ];
            invariant_counts := Increment(invariant_counts, Sprint(invs));

            order_pass := (#A) mod 49 eq 0;
            exponent_pass := HasExponent49(invs);
            if order_pass then
                pass_order +:= 1;
            end if;
            if exponent_pass then
                pass_exponent +:= 1;
                Append(~open_pass, <ai, bi, invs>);
            else
                mask[idx] := true;
                if order_pass then
                    false_order_positives +:= 1;
                end if;
            end if;
        end for;
    end for;

    summary := <total, good, boundary, pass_order, pass_exponent,
                false_order_positives, invariant_counts, open_pass>;
    return mask, summary;
end function;

procedure PrintFiniteSummary(p, summary)
    print "p", p,
          "total", summary[1],
          "good", summary[2],
          "boundary", summary[3],
          "pass_order49", summary[4],
          "pass_exponent49", summary[5],
          "false_order_positives", summary[6];
    print " invariant_counts", summary[7];
    if #summary[8] le 30 then
        print " open_exponent49", summary[8];
    else
        print " open_exponent49_first30", summary[8][1..30];
    end if;
end procedure;

procedure RunFinite()
    print "Z49 LOCAL FINITE SCOUT ON CONTACT-7";
    print "criterion", "49 divides exponent of full AbelianGroup(J(F_p))";
    print "primes", filter_primes;
    for p in filter_primes do
        mask, summary := BuildKillMask(p);
        PrintFiniteSummary(p, summary);
    end for;
end procedure;

procedure RunSearch()
    print "Z49 LOCAL RATIONAL SIEVE ON CONTACT-7";
    print "criterion", "good reduction requires exponent divisible by 49";
    print "height", height, "prime_bound", prime_bound,
          "primes", filter_primes, "max_exact", max_exact;

    masks := [];
    for p in filter_primes do
        mask, summary := BuildKillMask(p);
        Append(~masks, mask);
        PrintFiniteSummary(p, summary);
    end for;

    params := RationalParametersOfHeight(height);
    checked := 0;
    smooth_marked7 := 0;
    survivors := [];
    exact_tests := 0;
    exact_hits := [];
    first_kill := AssociativeArray();
    unresolved_denominator_counts := AssociativeArray();

    for a in params do
        for b in params do
            checked +:= 1;
            if progress_interval gt 0 and checked mod progress_interval eq 0 then
                print "progress", checked, "survivors", #survivors,
                      "exact_tests", exact_tests, "exact_hits", #exact_hits;
            end if;

            killed := false;
            for j in [1..#filter_primes] do
                p := filter_primes[j];
                oka, ai := ResidueIndex(a, p);
                okb, bi := ResidueIndex(b, p);
                if not oka or not okb then
                    unresolved_denominator_counts :=
                        Increment(unresolved_denominator_counts, p);
                    continue;
                end if;
                idx := ai*p + bi + 1;
                if masks[j][idx] then
                    first_kill := Increment(first_kill, p);
                    killed := true;
                    break;
                end if;
            end for;
            if killed then
                continue;
            end if;

            Append(~survivors, <a,b>);
            ok, f, h := Contact7Polynomial(a, b);
            if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 then
                continue;
            end if;
            yP := Evaluate(h, Q!1);
            if yP eq 0 then
                continue;
            end if;
            C := HyperellipticCurve(f);
            J := Jacobian(C);
            D7 := J![x-1, yP];
            if Order(D7) ne 7 then
                continue;
            end if;
            smooth_marked7 +:= 1;

            if exact_tests ge max_exact then
                continue;
            end if;
            fI, L := IntegralModel(f);
            CI := HyperellipticCurve(fI);
            JI := Jacobian(CI);
            G, phi := TorsionSubgroup(JI);
            invs := [ Z!n : n in Invariants(G) ];
            exact_tests +:= 1;
            print "EXACT_SURVIVOR", "a", a, "b", b,
                  "torsion", invs, "integral_scale", L;
            if HasExponent49(invs) then
                Append(~exact_hits, <a,b,invs,fI>);
                print "HIT49", exact_hits[#exact_hits];
            end if;
        end for;
    end for;

    print "DONE Z49 LOCAL RATIONAL SIEVE";
    print "parameters", #params, "checked", checked,
          "mask_survivors", #survivors,
          "smooth_marked7_survivors", smooth_marked7,
          "exact_tests", exact_tests, "exact_hits", #exact_hits;
    print "FIRST_KILL", first_kill;
    print "UNRESOLVED_DENOMINATOR_VISITS", unresolved_denominator_counts;
    if #survivors le 200 then
        print "SURVIVORS", survivors;
    else
        print "SURVIVORS_FIRST200", survivors[1..200];
    end if;
end procedure;

if mode eq "finite" then
    RunFinite();
elif mode eq "search" then
    RunSearch();
else
    error "mode must be finite or search";
end if;

quit;
