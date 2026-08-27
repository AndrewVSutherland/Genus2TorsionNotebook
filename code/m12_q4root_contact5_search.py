#!/usr/bin/env python3
"""Exact point-contact-5 search on the second rational M(12) chart.

The completed-square M(12) sextic is

    W = (T+1) Q4,
    T = a*x^2-x+r,
    Q4 = (x-r)^2*(T+1)+4*a*x^2*T.

This script uses the rational surface on which Q4 has a rational root v.
Put

    k = a*v^2/(v-r),  D = 1-k,  q = 2*k-1.

Then

    v-r = q^2/D,   a = k*q^2/(D*v^2).

After X=v/(x-v), multiplying the transformed sextic by the square D^4
gives the odd quintic F=P*Q, where

    P = 4*k*D^2*X^2 + (2*k*q^2-v*D)*X + k*q^2,
    R = -q^2*D*X^2 + (2*k*q^2-v*D)*X + k*q^2,
    Q = (q^2*X+v*D)^2*P + 4*k*q^2*D*(X+1)^2*R.

The X^4 coefficient of Q cancels, so Q is generically cubic.  The marked
order-12 point has X=-1.  A point (u,c) gives order five precisely when

    E3 = 8*A0^2*A3-A1*(4*A0*A2-A1^2) = 0,
    E4 = 64*A0^3*A4-(4*A0*A2-A1^2)^2 = 0,
    c^2=A0=F(u) != 0,

where Ai=F^(i)(u)/i!.  Thus a smooth chart point passing the exact tests is
an exact Z/60 candidate, to be certified independently in Magma.

This pure-Python implementation is intentionally dependency-free.  It is a
positive-controlled finite-field funnel followed by exact rational gcd and
square tests; it does not replace the final Magma torsion/simplicity proof.

With ``--parameter q --height H``, the searched box is exactly all reduced
rationals q=n/d and v=m/e with positive denominators and
``|n|,d,|m|,e <= H``.  This is a naive height box in (q,v), not in (k,a,r).
"""

from __future__ import annotations

import argparse
from collections import Counter
import math
import time
from fractions import Fraction


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def padd(a, b):
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return trim(out)


def psub(a, b):
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] -= x
    return trim(out)


def pscale(a, c):
    return trim([c * x for x in a])


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def pdivmod_q(a, b):
    """Polynomial division over Q (or any field represented exactly)."""
    a = trim(a)
    b = trim(b)
    if b == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    if len(a) < len(b):
        return [0], a
    q = [0] * (len(a) - len(b) + 1)
    while a != [0] and len(a) >= len(b):
        d = len(a) - len(b)
        c = Fraction(a[-1]) / Fraction(b[-1])
        q[d] = c
        for j in range(len(b)):
            a[d + j] -= c * b[j]
        a = trim(a)
    return trim(q), a


def pgcd_q(a, b):
    a, b = trim(a), trim(b)
    while b != [0]:
        _, r = pdivmod_q(a, b)
        a, b = b, r
    if a == [0]:
        return [0]
    return pscale(a, Fraction(1) / Fraction(a[-1]))


def positive_divisors(n):
    """Positive divisors of a nonzero integer, in increasing order."""
    n = abs(n)
    small, large = [], []
    for d in range(1, math.isqrt(n) + 1):
        if n % d == 0:
            small.append(d)
            if d * d != n:
                large.append(n // d)
    return small + list(reversed(large))


def rational_roots_q(f):
    """All distinct rational roots of an exact Q-polynomial.

    The usual linear-gcd case is constant time.  For the exceptional
    higher-degree gcd case, clear denominators and use the rational-root
    theorem.  This prevents a valid rational contact root from being lost
    merely because E3 and E4 have an additional common factor.
    """
    f = trim([Fraction(c) for c in f])
    if len(f) <= 1:
        return []
    if len(f) == 2:
        return [-f[0] / f[1]]

    roots = []
    while len(f) > 1 and f[0] == 0:
        roots.append(Fraction(0))
        f = trim(f[1:])
    if len(f) <= 1:
        return roots
    if len(f) == 2:
        root = -f[0] / f[1]
        if root not in roots:
            roots.append(root)
        return roots

    denominator_lcm = 1
    for c in f:
        denominator_lcm = math.lcm(denominator_lcm, c.denominator)
    integral = [int(c * denominator_lcm) for c in f]
    content = 0
    for c in integral:
        content = math.gcd(content, abs(c))
    integral = [c // content for c in integral]
    for num in positive_divisors(integral[0]):
        for den in positive_divisors(integral[-1]):
            for sign in (-1, 1):
                root = Fraction(sign * num, den)
                if root not in roots and peval(f, root) == 0:
                    roots.append(root)
    return roots


def pderiv(a):
    return trim([i * a[i] for i in range(1, len(a))]) if len(a) > 1 else [0]


def peval(a, x):
    z = 0
    for c in reversed(a):
        z = z * x + c
    return z


def chart_quintic(k, v):
    """Return low-to-high coefficients of the square-equivalent quintic."""
    D = 1 - k
    q = 2 * k - 1
    q2 = q * q
    P = [k * q2, 2 * k * q2 - v * D, 4 * k * D * D]
    R = [k * q2, 2 * k * q2 - v * D, -q2 * D]
    linear_sq = pmul([v * D, q2], [v * D, q2])
    one_plus_x_sq = [1, 2, 1]
    Q = padd(pmul(linear_sq, P),
             pscale(pmul(one_plus_x_sq, R), 4 * k * q2 * D))
    return trim(pmul(P, Q))


def taylor_polynomials(f):
    """Ai(u)=f^(i)(u)/i! as polynomials in u, for i=0,...,4."""
    ans = []
    for i in range(5):
        Ai = [0] * (len(f) - i)
        for j in range(i, len(f)):
            Ai[j - i] = f[j] * math.comb(j, i)
        ans.append(trim(Ai))
    return ans


def contact_polynomials(f):
    A0, A1, A2, A3, A4 = taylor_polynomials(f)
    q2 = psub(pscale(pmul(A0, A2), 4), pmul(A1, A1))
    E3 = psub(pscale(pmul(pmul(A0, A0), A3), 8), pmul(A1, q2))
    E4 = psub(pscale(pmul(pmul(pmul(A0, A0), A0), A4), 64), pmul(q2, q2))
    return trim(E3), trim(E4)


def is_square_fraction(x):
    if x <= 0:
        return False, None
    n, d = x.numerator, x.denominator
    sn, sd = math.isqrt(n), math.isqrt(d)
    if sn * sn == n and sd * sd == d:
        return True, Fraction(sn, sd)
    return False, None


def rational_values(height):
    vals = []
    for den in range(1, height + 1):
        for num in range(-height, height + 1):
            if math.gcd(num, den) == 1:
                vals.append(Fraction(num, den))
    # Denominator-positive reduced representation is already unique.
    return vals


def mod_inv(a, p):
    return pow(a % p, p - 2, p)


def residue(x, p):
    if x.denominator % p == 0:
        return None
    return (x.numerator % p) * mod_inv(x.denominator, p) % p


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for d in range(3, math.isqrt(n) + 1, 2):
        if n % d == 0:
            return False
    return True


def validate_primes(primes):
    if len(set(primes)) != len(primes):
        raise ValueError("finite-mask primes must be distinct")
    invalid = [p for p in primes if p <= 5 or not is_prime(p)]
    if invalid:
        raise ValueError(
            "finite-mask moduli must be primes > 5; invalid: "
            + ",".join(str(p) for p in invalid)
        )


def mod_trim(a, p):
    a = [x % p for x in a]
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def mod_divmod(a, b, p):
    a = mod_trim(a, p)
    b = mod_trim(b, p)
    if b == [0]:
        raise ZeroDivisionError
    if len(a) < len(b):
        return [0], a
    q = [0] * (len(a) - len(b) + 1)
    ib = mod_inv(b[-1], p)
    while a != [0] and len(a) >= len(b):
        d = len(a) - len(b)
        c = a[-1] * ib % p
        q[d] = c
        for j in range(len(b)):
            a[d + j] = (a[d + j] - c * b[j]) % p
        a = mod_trim(a, p)
    return mod_trim(q, p), a


def mod_gcd(a, b, p):
    a, b = mod_trim(a, p), mod_trim(b, p)
    while b != [0]:
        _, r = mod_divmod(a, b, p)
        a, b = b, r
    if a == [0]:
        return [0]
    return mod_trim([(x * mod_inv(a[-1], p)) % p for x in a], p)


def finite_contact(k, v, p):
    """Return (state, roots): state is bad, killed, or allowed.

    This affine enumeration is rigorous for the validated primes p>5.
    If f=c5*X^5+... is a good quintic, the degree-12 leading coefficient
    of E3 is 5*c5^3, hence is nonzero modulo p.  Thus E3 and E4 cannot have
    a common projective root at u=infinity.  A rational contact abscissa
    whose denominator is divisible by p is therefore impossible at such a
    good reduction.  Moreover, A0 cannot reduce to zero: E3=E4=0 and A0=0
    force A1=0, contradicting squarefreeness.  It is consequently safe to
    enumerate only finite u with nonzero square A0.
    """
    D, q = (1 - k) % p, (2 * k - 1) % p
    if k % p == 0 or v % p == 0 or D == 0 or q == 0:
        return "bad", []
    f = mod_trim(chart_quintic(k % p, v % p), p)
    if len(f) != 6 or mod_gcd(f, pderiv(f), p) != [1]:
        return "bad", []
    roots = []
    for u in range(p):
        A = []
        for i in range(5):
            Ai = 0
            for j in range(i, 6):
                Ai += f[j] * math.comb(j, i) * pow(u, j - i, p)
            A.append(Ai % p)
        A0, A1, A2, A3, A4 = A
        q2c = (4 * A0 * A2 - A1 * A1) % p
        e3 = (8 * A0 * A0 * A3 - A1 * q2c) % p
        e4 = (64 * A0 * A0 * A0 * A4 - q2c * q2c) % p
        if e3 == 0 and e4 == 0 and A0 != 0 and pow(A0, (p - 1) // 2, p) == 1:
            roots.append(u)
    return ("allowed", roots) if roots else ("killed", [])


def finite_masks(primes):
    masks = {}
    for p in primes:
        allowed, bad = set(), set()
        good = killed = root_count = 0
        for k in range(p):
            for v in range(p):
                state, roots = finite_contact(k, v, p)
                key = k * p + v
                if state == "bad":
                    bad.add(key)
                elif state == "allowed":
                    allowed.add(key)
                    good += 1
                    root_count += len(roots)
                else:
                    good += 1
                    killed += 1
        masks[p] = (allowed, bad)
        print(f"MASK p={p} allowed_pairs={len(allowed)} good_pairs={good} "
              f"killed_pairs={killed} contact_roots={root_count} bad={len(bad)}",
              flush=True)
    return masks


def passes_masks(k_res, v_res, primes, masks):
    good_checked = 0
    for p in primes:
        kr, vr = k_res[p], v_res[p]
        if kr is None or vr is None:
            continue
        key = kr * p + vr
        allowed, bad = masks[p]
        if key in bad:
            continue
        good_checked += 1
        if key not in allowed:
            return False, good_checked
    return True, good_checked


def exact_candidates(k, v):
    if k == 0 or v == 0 or k == 1 or 2 * k == 1:
        return [], "chart_boundary"
    f = chart_quintic(k, v)
    # On the open chart, degree < 5 means that the chosen Q4 root ceased to
    # be a simple branch point.  The transformed model is not a smooth
    # genus-2 quintic, so it cannot be a target specialization.
    if len(f) != 6:
        return [], "degree_boundary"
    if len(pgcd_q(f, pderiv(f))) != 1:
        return [], "singular"
    E3, E4 = contact_polynomials(f)
    g = pgcd_q(E3, E4)
    if len(g) == 1:
        return [], "no_common_root"
    roots = rational_roots_q(g)
    if not roots:
        return [], f"gcd_degree_{len(g)-1}_no_rational_root"

    candidates = []
    for u in roots:
        if peval(E3, u) != 0 or peval(E4, u) != 0:
            raise AssertionError("reported gcd root is not a common root")
        f0 = peval(f, u)
        square, c = is_square_fraction(f0)
        if not square or f0 == 0:
            continue

        # Exact reconstruction of h=c+d*(X-u)+e*(X-u)^2.  The other
        # square root -c gives -h and the inverse order-5 class, so one
        # sign is enough for existence and subsequent certification.
        A0, A1, A2, _, _ = [peval(A, u) for A in taylor_polynomials(f)]
        dd = A1 / (2 * c)
        ee = (A2 - dd * dd) / (2 * c)
        h = [c - dd * u + ee * u * u, dd - 2 * ee * u, ee]
        contact = psub(f, pmul(h, h))
        # f-h^2 = lc(f)*(X-u)^5.
        rhs = pscale(
            [(-u) ** (5 - i) * math.comb(5, i) for i in range(6)],
            f[-1],
        )
        if trim(contact) != trim(rhs):
            raise AssertionError("contact reconstruction failed")
        candidates.append((u, c, f, h))
    if candidates:
        return candidates, "hit"
    return [], "nonsquare"


def positive_control():
    # h=2+3X+5X^2 and f=h^2+7X^5 has exact contact at u=0.
    h = [Fraction(2), Fraction(3), Fraction(5)]
    f = padd(pmul(h, h), [0, 0, 0, 0, 0, Fraction(7)])
    E3, E4 = contact_polynomials(f)
    g = pgcd_q(E3, E4)
    if peval(E3, 0) != 0 or peval(E4, 0) != 0 or peval(f, 0) != 4:
        raise AssertionError("contact covariant positive control failed")
    if len(g) < 2:
        raise AssertionError("positive-control gcd has no root")
    if Fraction(0) not in rational_roots_q(g):
        raise AssertionError("positive-control rational root was not recovered")

    # Exercise the exceptional higher-degree-gcd rational-root path.
    test_gcd = pmul([Fraction(0), Fraction(1)],
                    [Fraction(1), Fraction(0), Fraction(1)])
    if rational_roots_q(test_gcd) != [Fraction(0)]:
        raise AssertionError("higher-degree rational-root control failed")

    # Compact-chart algebraic control at (k,v)=(2,1).
    chart_f = chart_quintic(Fraction(2), Fraction(1))
    if len(chart_f) != 6 or len(pgcd_q(chart_f, pderiv(chart_f))) != 1:
        raise AssertionError("compact-chart smooth quintic control failed")
    if peval(chart_f, Fraction(-1)) != 12100:
        raise AssertionError("marked D12 square control failed")

    # Corrected p=19 affine mask has this explicit contact residue.
    state, roots = finite_contact(15, 2, 19)
    if state != "allowed" or 0 not in roots:
        raise AssertionError("finite p=19 contact control failed")
    print("CONTACT_SELF_TEST_PASS", "Q4_CHART_PASS", "P19_MASK_PASS",
          flush=True)


def parse_primes(s):
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--height", type=int, default=30)
    ap.add_argument("--primes", default="7,11,13,17,19,23,29")
    ap.add_argument("--parameter", choices=("q", "k"), default="q",
                    help="height coordinate; q=2*k-1 is the branch-free normalization")
    ap.add_argument("--progress", type=int, default=500000)
    ap.add_argument("--max-hits", type=int, default=20)
    args = ap.parse_args()

    if args.height < 1:
        ap.error("--height must be positive")
    if args.progress < 0:
        ap.error("--progress must be nonnegative")
    if args.max_hits < 1:
        ap.error("--max-hits must be positive")

    primes = parse_primes(args.primes)
    try:
        validate_primes(primes)
    except ValueError as exc:
        ap.error(str(exc))
    positive_control()
    masks = finite_masks(primes)
    # Most selective first reduces pair-scan work without changing rigor.
    primes.sort(key=lambda p: len(masks[p][0]) / max(1, p * p - len(masks[p][1])))
    print("MASK_ORDER", primes, flush=True)

    vals = rational_values(args.height)
    base_vals = ([(x + 1) / 2 for x in vals]
                 if args.parameter == "q" else vals)
    residues = {p: [residue(x, p) for x in base_vals] for p in primes}
    print(f"SEARCH parameter={args.parameter} height={args.height} "
          f"rational_values={len(vals)} "
          f"pairs={len(vals)**2}", flush=True)
    print("HEIGHT_BOX numerator_abs<=height denominator<=height "
          "reduced_denominator_positive; v uses this box and "
          + ("q uses this box with k=(q+1)/2"
             if args.parameter == "q" else "k uses this box"),
          flush=True)

    pairs = survivors = smooth_tests = hits = 0
    state_counts = Counter()
    stopped_early = False
    t0 = time.time()
    for ik, k in enumerate(base_vals):
        k_res = {p: residues[p][ik] for p in primes}
        for iv, v in enumerate(vals):
            pairs += 1
            v_res = {p: residues[p][iv] for p in primes}
            ok, good_checked = passes_masks(k_res, v_res, primes, masks)
            if ok:
                survivors += 1
                candidates, state = exact_candidates(k, v)
                state_counts[state] += 1
                if state not in {"chart_boundary", "degree_boundary", "singular"}:
                    smooth_tests += 1
                for u, c, f, h in candidates:
                    hits += 1
                    D, q = 1 - k, 2 * k - 1
                    d = q * q / D
                    a = k * q * q / (D * v * v)
                    r = v - d
                    print("CONTACT5_CANDIDATE", "k", k, "v", v, "u", u,
                          "c", c, "a", a, "r", r, "f", f, "h", h,
                          flush=True)
                    if hits >= args.max_hits:
                        stopped_early = True
                        break
            if args.progress and pairs % args.progress == 0:
                print(f"PROGRESS pairs={pairs} survivors={survivors} "
                      f"smooth_tests={smooth_tests} "
                      f"hits={hits} seconds={time.time()-t0:.1f}", flush=True)
            if stopped_early:
                break
        if stopped_early:
            break

    chart_boundary = state_counts["chart_boundary"]
    affine_survivors = survivors - chart_boundary
    complete = pairs == len(vals) ** 2 and not stopped_early
    print("EXACT_STATE_COUNTS",
          " ".join(f"{name}={state_counts[name]}"
                   for name in sorted(state_counts)), flush=True)
    print(f"DONE pairs={pairs} survivors={survivors} "
          f"chart_boundary={chart_boundary} "
          f"affine_survivors={affine_survivors} "
          f"smooth_tests={smooth_tests} hits={hits} "
          f"complete={complete} seconds={time.time()-t0:.1f}",
          flush=True)


if __name__ == "__main__":
    main()
