//////////////////////////////////////////////////////////////////////
//  Finite-field subgroup scout for [3,18] in the rational-root
//  subfamily of the contact-9 construction.
//
//  This deliberately uses the full abelian invariants of J(F_p).
//  Divisibility of #J(F_p) by 54 is not enough to ensure that
//  Z/3 Z x Z/18 Z embeds.  At a good prime p, only the prime-to-p
//  part of rational torsion is required to inject.
//
//  Typical run:
//      magma -b prime_bound:=101 code/contact9_318_finite_scout.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned prime_bound then
    prime_bound := 101;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;

Z := Integers();
target := [ Z!3, Z!18 ];

function StripPPart(n, p)
    m := Z!n;
    while m mod p eq 0 do
        m div:= p;
    end while;
    return m;
end function;

function PrimeToPTarget(invs, p)
    result := [];
    for n in invs do
        m := StripPPart(n, p);
        if m gt 1 then
            Append(~result, m);
        end if;
    end for;
    return result;
end function;

// The rank of the l^k-torsion in a finite abelian group is the number
// of invariant factors divisible by l^k.  Comparing these ranks for
// every l and k is equivalent to subgroup containment.
function HasSubgroupEmbedding(ambient, required)
    if #required eq 0 then
        return true;
    end if;
    primes := PrimeDivisors(&*required);
    for ell in primes do
        max_exp := Max([ Valuation(n, ell) : n in required ]);
        for k in [1..max_exp] do
            need := #[ n : n in required | Valuation(n, ell) ge k ];
            have := #[ n : n in ambient | Valuation(Z!n, ell) ge k ];
            if have lt need then
                return false;
            end if;
        end for;
    end for;
    return true;
end function;

assert PrimeToPTarget(target, 2) eq [ 3, 9 ];
assert PrimeToPTarget(target, 3) eq [ 2 ];
assert PrimeToPTarget(target, 7) eq target;
assert HasSubgroupEmbedding([ 3, 18 ], target);
assert not HasSubgroupEmbedding([ 2, 36 ], target);
assert not HasSubgroupEmbedding([ 2, 27 ], target);

function Contact9RootPolynomialFinite(F, s, eps)
    PF<X> := PolynomialRing(F);
    if Characteristic(F) eq 2 then
        return false, PF!0, PF!0, F!0, F!0, "characteristic_2";
    end if;

    r := 1 - s^2;
    if r eq 0 then
        return false, PF!0, PF!0, F!0, r, "r_zero";
    end if;

    h0 := 1 - (F!9/F!2)*r + (F!63/F!8)*r^2
          - (F!105/F!16)*r^3;
    a := (eps*s^9 - h0)/r^4;
    h := 1 - (F!9/F!2)*X + (F!63/F!8)*X^2
         - (F!105/F!16)*X^3 + a*X^4;
    f := ExactQuotient(h^2 + (X - 1)^9, X^4);
    if Evaluate(f, r) ne 0 then
        return false, f, h, a, r, "root_identity";
    end if;
    if Degree(f) ne 5 then
        return false, f, h, a, r, "degree_drop";
    end if;
    if Discriminant(f) eq 0 then
        return false, f, h, a, r, "singular";
    end if;
    if Evaluate(h, F!1) eq 0 then
        return false, f, h, a, r, "marked_boundary";
    end if;
    return true, f, h, a, r, "good";
end function;

function Increment(counts, key)
    if IsDefined(counts, key) then
        counts[key] +:= 1;
    else
        counts[key] := 1;
    end if;
    return counts;
end function;

print "Contact-9 rational-root finite subgroup scout for [3,18]";
print "prime_bound", prime_bound, "target", target;
print "containment_test", "valuation ranks in full AbelianGroup invariants";

for p in [ q : q in PrimesUpTo(prime_bound) | q ge 3 ] do
    F := GF(p);
    required := PrimeToPTarget(target, p);
    total := 0;
    good := 0;
    pass_order := 0;
    pass_embedding := 0;
    false_order_positives := 0;
    failure_counts := AssociativeArray();
    invariant_counts := AssociativeArray();
    open_target := [];

    for s in F do
        for eps in [ F!-1, F!1 ] do
            total +:= 1;
            ok, f, h, a, r, reason :=
                Contact9RootPolynomialFinite(F, s, eps);
            if not ok then
                failure_counts := Increment(failure_counts, reason);
                continue;
            end if;

            good +:= 1;
            C := HyperellipticCurve(f);
            J := Jacobian(C);
            A, phi := AbelianGroup(J);
            invs := [ Z!n : n in Invariants(A) ];
            key := Sprint(invs);
            invariant_counts := Increment(invariant_counts, key);

            required_order := #required eq 0 select 1 else &*required;
            if (#A) mod required_order eq 0 then
                pass_order +:= 1;
            end if;
            embeds := HasSubgroupEmbedding(invs, required);
            if embeds then
                pass_embedding +:= 1;
                Append(~open_target,
                       <Z!s, Z!eps, Z!a, Z!r, invs>);
            elif (#A) mod required_order eq 0 then
                false_order_positives +:= 1;
            end if;
        end for;
    end for;

    print "p", p, "required", required, "total", total,
          "good", good, "pass_order", pass_order,
          "pass_embedding", pass_embedding,
          "false_order_positives", false_order_positives;
    if #Keys(failure_counts) gt 0 then
        print " failures", failure_counts;
    end if;
    print " invariant_counts", invariant_counts;
    if #open_target le 30 then
        print " open_target", open_target;
    else
        print " open_target_first30", open_target[1..30];
    end if;
end for;

quit;
