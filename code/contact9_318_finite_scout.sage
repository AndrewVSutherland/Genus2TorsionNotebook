#!/usr/bin/env sage
"""Independent Sage fallback for the contact-9 [3,18] finite scout.

Sage does not expose a finite genus-2 Jacobian's abelian invariants directly.
For these one-parameter, small-field diagnostics we enumerate the unique
reduced Mumford representatives, then recover the Sylow invariant factors
from the sizes of the kernels of multiplication by ell^k.
"""

import argparse
from collections import Counter


TARGET = [ZZ(3), ZZ(18)]


def strip_p_part(n, p):
    n = ZZ(n)
    while n % p == 0:
        n //= p
    return n


def prime_to_p_target(p):
    return [m for n in TARGET if (m := strip_p_part(n, p)) > 1]


def has_subgroup_embedding(ambient, required):
    """Test finite abelian subgroup containment prime by prime."""
    if not required:
        return True
    for ell in prime_divisors(prod(required)):
        max_exp = max(ZZ(n).valuation(ell) for n in required)
        for k in range(1, max_exp + 1):
            need = sum(ZZ(n).valuation(ell) >= k for n in required)
            have = sum(ZZ(n).valuation(ell) >= k for n in ambient)
            if have < need:
                return False
    return True


assert prime_to_p_target(2) == [3, 9]
assert prime_to_p_target(3) == [2]
assert prime_to_p_target(7) == TARGET
assert has_subgroup_embedding([3, 18], TARGET)
assert not has_subgroup_embedding([2, 36], TARGET)
assert not has_subgroup_embedding([2, 27], TARGET)


def contact9_root_polynomial(F, s, eps):
    R = PolynomialRing(F, "x")
    x = R.gen()
    r = 1 - s**2
    if r == 0:
        return None, None, None, r, "r_zero"

    h0 = 1 - F(9)/2*r + F(63)/8*r**2 - F(105)/16*r**3
    a = (eps*s**9 - h0)/r**4
    h = 1 - F(9)/2*x + F(63)/8*x**2 - F(105)/16*x**3 + a*x**4
    f, rem = (h**2 + (x - 1)**9).quo_rem(x**4)
    if rem:
        return f, h, a, r, "root_identity"
    if f(r) != 0:
        return f, h, a, r, "root_identity"
    if f.degree() != 5:
        return f, h, a, r, "degree_drop"
    if f.discriminant() == 0:
        return f, h, a, r, "singular"
    if h(1) == 0:
        return f, h, a, r, "marked_boundary"
    return f, h, a, r, "good"


def square_roots(a):
    if not a.is_square():
        return []
    if a == 0:
        return [a]
    root = a.sqrt()
    return [root, -root]


def reduced_jacobian_points(f):
    """Enumerate reduced Mumford pairs [u,v] for y^2=f, deg(f)=5."""
    F = f.base_ring()
    R = f.parent()
    x = R.gen()
    X = HyperellipticCurve(f).jacobian()(F)
    points_by_rep = {repr(X(0)): X(0)}

    def add_point(D):
        points_by_rep[repr(D)] = D

    # Degree-one reduced divisors.
    for a in F:
        u = x + a
        value = f(-a)
        for d in square_roots(value):
            add_point(X([u, R(d)]))

    # For u=x^2+a*x+b and v=c*x+d, reduction modulo u gives
    #   -a*c^2+2*c*d=A,  d^2-b*c^2=B,
    # where f mod u=A*x+B.  Eliminating d leaves a quadratic in c^2.
    T = PolynomialRing(F, "t")
    t = T.gen()
    for a in F:
        for b in F:
            u = x**2 + a*x + b
            rem = f % u
            A = rem[1]
            B = rem[0]
            candidates = set()

            if A == 0:
                for d in square_roots(B):
                    candidates.add((F(0), d))

            eq_t = (a**2 - 4*b)*t**2 + (2*a*A - 4*B)*t + A**2
            if eq_t:
                for tv, multiplicity in eq_t.roots():
                    if tv == 0:
                        continue
                    for c in square_roots(tv):
                        d = (A + a*c**2)/(2*c)
                        candidates.add((c, d))

            for c, d in candidates:
                v = c*x + d
                if (f - v**2) % u:
                    raise ArithmeticError("Mumford solution failed verification")
                try:
                    add_point(X([u, v]))
                except (ArithmeticError, ValueError, ZeroDivisionError):
                    # A non-coprime pair is not a reduced divisor.  The same
                    # class already occurs in lower degree.
                    pass
    return X, list(points_by_rep.values())


def jacobian_order_from_curve(f):
    C = HyperellipticCurve(f)
    n1, n2 = C.count_points(2)
    p = f.base_ring().cardinality()
    s1 = p + 1 - n1
    sum_alpha2 = p**2 + 1 - n2
    e2 = (s1**2 - sum_alpha2)//2
    return ZZ(1 - s1 + e2 - p*s1 + p**2)


def abelian_invariants(f):
    X, points = reduced_jacobian_points(f)
    expected = jacobian_order_from_curve(f)
    if len(points) != expected:
        raise ArithmeticError(
            f"Mumford enumeration has {len(points)} points, expected {expected}"
        )

    zero = X(0)
    sylow_exponents = {}
    for ell, total_exp in factor(expected):
        kernel_logs = []
        for k in range(1, total_exp + 1):
            scalar = ell**k
            kernel_size = sum(scalar*D == zero for D in points)
            log_size = ZZ(kernel_size).valuation(ell)
            if ell**log_size != kernel_size:
                raise ArithmeticError("multiplication kernel size is not an ell-power")
            kernel_logs.append(log_size)
            if log_size == total_exp:
                break

        ranks = []
        previous = 0
        for log_size in kernel_logs:
            ranks.append(log_size - previous)
            previous = log_size
        ranks.append(0)

        exponents = []
        for k in range(1, len(ranks)):
            exponents.extend([k]*(ranks[k - 1] - ranks[k]))
        sylow_exponents[ZZ(ell)] = sorted(exponents)

    rank = max((len(exps) for exps in sylow_exponents.values()), default=0)
    invs = [ZZ(1)]*rank
    for ell, exps in sylow_exponents.items():
        padded = [0]*(rank - len(exps)) + exps
        invs = [n*ell**e for n, e in zip(invs, padded)]
    return [n for n in invs if n > 1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-bound", type=int, default=43)
    parser.add_argument("--show-limit", type=int, default=30)
    args = parser.parse_args()

    print("Contact-9 rational-root finite subgroup scout for [3,18]")
    print(f"prime_bound {args.prime_bound} target {TARGET}")
    print("containment_test valuation ranks in full enumerated invariants")

    for p in prime_range(3, args.prime_bound + 1):
        F = GF(p)
        required = prime_to_p_target(p)
        failures = Counter()
        invariant_counts = Counter()
        open_target = []
        boundary_records = []
        good_nontarget = []
        total = good = pass_order = pass_embedding = false_order_positives = 0

        for s in F:
            for eps in (F(-1), F(1)):
                total += 1
                f, h, a, r, reason = contact9_root_polynomial(F, s, eps)
                if reason != "good":
                    failures[reason] += 1
                    boundary_records.append(
                        (ZZ(s), ZZ(eps), None if a is None else ZZ(a), ZZ(r), reason)
                    )
                    continue

                good += 1
                invs = abelian_invariants(f)
                invariant_counts[tuple(invs)] += 1
                required_order = prod(required, ZZ(1))
                group_order = prod(invs, ZZ(1))
                order_passes = group_order % required_order == 0
                embeds = has_subgroup_embedding(invs, required)
                pass_order += order_passes
                pass_embedding += embeds
                false_order_positives += order_passes and not embeds
                if embeds:
                    open_target.append((ZZ(s), ZZ(eps), ZZ(a), ZZ(r), invs))
                else:
                    good_nontarget.append((ZZ(s), ZZ(eps), ZZ(a), ZZ(r), invs))

        print(
            f"p {p} required {required} total {total} good {good} "
            f"pass_order {pass_order} pass_embedding {pass_embedding} "
            f"false_order_positives {false_order_positives}"
        )
        print(f" failures {dict(sorted(failures.items()))}")
        print(f" invariant_counts {dict(sorted(invariant_counts.items()))}")
        label = "open_target" if len(open_target) <= args.show_limit else "open_target_first"
        print(f" {label} {open_target[:args.show_limit]}")
        if pass_embedding == 0:
            print(f" boundary_records {boundary_records[:args.show_limit]}")
            print(f" good_nontarget {good_nontarget[:args.show_limit]}")


main()
