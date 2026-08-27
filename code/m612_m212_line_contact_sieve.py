#!/usr/bin/env python3
"""Local cubic-contact masks for the simple M(2,12) line.

The line is

    a = (1-r)/4,              equivalently z = +/- r,

in the M(2,12) chart.  Its completed-square sextic is normalized to an
odd quintic with a rational branch point at X=0.  For a smooth odd model
f(X), every rational order-3 direction has a cubic-contact presentation

    h(X)^2 - f(X) = kappa * q(X)^3,

with q=X^2+U*X+V monic.  We enumerate these identities over F_p by a
triangular coefficient solve.  The marked order-12 class supplies one
q-support.  Hence a second distinct q-support is a necessary condition for
an independent rational 3-torsion direction.

This is deliberately independent of Magma; it is intended to make the
one-variable local masks while long Magma jobs occupy the available slots.

Examples:

    python3 code/m612_m212_line_contact_sieve.py --primes 5,7,11,13,17,19,23
    python3 code/m612_m212_line_contact_sieve.py --prime-bound 59 \
        --write-masks data/m612_m212_line_masks_p59.txt
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def trim(a: list[int], p: int) -> list[int]:
    a = [x % p for x in a]
    while a and a[-1] == 0:
        a.pop()
    return a


def padd(a: list[int], b: list[int], p: int) -> list[int]:
    n = max(len(a), len(b))
    return trim(
        [
            (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(n)
        ],
        p,
    )


def pmul(a: list[int], b: list[int], p: int) -> list[int]:
    if not a or not b:
        return []
    c = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            c[i + j] = (c[i + j] + ai * bj) % p
    return trim(c, p)


def pscale(a: list[int], c: int, p: int) -> list[int]:
    return trim([c * x for x in a], p)


def ppow(a: list[int], n: int, p: int) -> list[int]:
    out = [1]
    base = a
    while n:
        if n & 1:
            out = pmul(out, base, p)
        base = pmul(base, base, p)
        n >>= 1
    return out


def pdivrem(a: list[int], b: list[int], p: int) -> tuple[list[int], list[int]]:
    a = trim(a, p)
    b = trim(b, p)
    q = [0] * max(1, len(a) - len(b) + 1)
    ib = pow(b[-1], -1, p)
    while a and len(a) >= len(b):
        k = len(a) - len(b)
        c = a[-1] * ib % p
        q[k] = c
        for j, bj in enumerate(b):
            a[k + j] = (a[k + j] - c * bj) % p
        a = trim(a, p)
    return trim(q, p), a


def pgcd(a: list[int], b: list[int], p: int) -> list[int]:
    while trim(b, p):
        a, b = b, pdivrem(a, b, p)[1]
    a = trim(a, p)
    return pscale(a, pow(a[-1], -1, p), p) if a else []


def derivative(a: list[int], p: int) -> list[int]:
    return trim([i * a[i] for i in range(1, len(a))], p)


def original_sextic(r: int, p: int) -> list[int]:
    """W=(x-r)^2(T+1)^2+4*a*x^2*T*(T+1), low coefficients first."""
    a = (1 - r) * pow(4, -1, p) % p
    T = [r, -1, a]
    Tp1 = [(r + 1) % p, -1, a]
    h = pmul([-r, 1], Tp1, p)
    W = padd(
        pmul(h, h, p),
        pscale(pmul(pmul([0, 0, 1], T, p), Tp1, p), 4 * a, p),
        p,
    )
    return W


def is_good_genus2(W: list[int], p: int) -> bool:
    W = trim(W, p)
    return len(W) - 1 in (5, 6) and len(pgcd(W, derivative(W, p), p)) == 1


def odd_model(r: int, p: int) -> list[int]:
    """Return a degree-5 model f with f(0)=0, isomorphic over F_p.

    The rational roots 2 and u=2r/(2-r) are sent respectively to infinity
    and zero by x=(2X-u)/(X-1).  At r=2 the original model already has odd
    degree, and translating the root x=2 to zero is enough.
    """
    W = original_sextic(r, p)
    if r % p == 2 % p:
        # f(X)=W(X+2)
        out: list[int] = []
        for i, ci in enumerate(W):
            out = padd(out, pscale(pmul([1], ppow([2, 1], i, p), p), ci, p), p)
        return trim(out, p)

    u = 2 * r * pow(2 - r, -1, p) % p
    out: list[int] = []
    for i, ci in enumerate(W):
        # (X-1)^6 * ci*((2X-u)/(X-1))^i
        term = pmul(ppow([-u, 2], i, p), ppow([-1, 1], 6 - i, p), p)
        out = padd(out, pscale(term, ci, p), p)
    return trim(out, p)


def contact_supports(f: list[int], p: int) -> set[tuple[int, int]]:
    """Enumerate q=(X^2+U X+V) supporting cubic-contact identities."""
    ff = trim(f, p) + [0] * (7 - len(trim(f, p)))
    f0, f1, f2, f3, f4, f5, f6 = ff[:7]
    supports: set[tuple[int, int]] = set()

    for U in range(p):
        for V in range(p):
            qpoly = [V, U, 1]
            if len(pgcd(qpoly, f, p)) > 1:
                continue

            # Low-to-high coefficients of (X^2+U X+V)^3.
            q3 = [
                V**3,
                3 * U * V**2,
                3 * (U**2 * V + V**2),
                U**3 + 6 * U * V,
                3 * (U**2 + V),
                3 * U,
                1,
            ]
            q3 = [x % p for x in q3]
            found = False

            for m in range(p):
                kappa = (m * m - f6) % p
                if m:
                    inv2m = pow(2 * m, -1, p)
                    n = (f5 + kappa * q3[5]) * inv2m % p
                    s = (f4 + kappa * q3[4] - n * n) * inv2m % p
                    t = (f3 + kappa * q3[3] - 2 * n * s) * inv2m % p
                    if (s * s + 2 * n * t - f2 - kappa * q3[2]) % p:
                        continue
                    if (2 * s * t - f1 - kappa * q3[1]) % p:
                        continue
                    if (t * t - f0 - kappa * q3[0]) % p:
                        continue
                    found = True
                    break

                # The triangular solve has a small exceptional m=0 branch.
                if (-f5 - kappa * q3[5]) % p:
                    continue
                for n in range(p):
                    if (n * n - f4 - kappa * q3[4]) % p:
                        continue
                    if n:
                        s_values = [
                            (f3 + kappa * q3[3]) * pow(2 * n, -1, p) % p
                        ]
                    else:
                        if (f3 + kappa * q3[3]) % p:
                            continue
                        s_values = range(p)
                    for s in s_values:
                        if n:
                            t_values = [
                                (f2 + kappa * q3[2] - s * s)
                                * pow(2 * n, -1, p)
                                % p
                            ]
                        elif (s * s - f2 - kappa * q3[2]) % p == 0:
                            t_values = range(p)
                        else:
                            continue
                        for t in t_values:
                            if (2 * s * t - f1 - kappa * q3[1]) % p:
                                continue
                            if (t * t - f0 - kappa * q3[0]) % p:
                                continue
                            found = True
                            break
                        if found:
                            break
                    if found:
                        break
                if found:
                    break

            if found:
                supports.add((U, V))
    return supports


def prime_mask(p: int) -> tuple[list[int], list[int], dict[int, int]]:
    """Return affine allowed residues, bad residues, and contact counts."""
    allowed: list[int] = []
    bad: list[int] = []
    counts: dict[int, int] = {}
    for r in range(p):
        W = original_sextic(r, p)
        if not is_good_genus2(W, p):
            bad.append(r)
            allowed.append(r)
            continue
        f = odd_model(r, p)
        if len(f) - 1 != 5 or not f or f[0] % p != 0:
            raise RuntimeError(f"odd normalization failed at p={p}, r={r}: {f}")
        n = len(contact_supports(f, p))
        counts[r] = n
        # One support is the marked 3-direction.  A rational independent
        # direction necessarily supplies a second distinct support.
        if n >= 2:
            allowed.append(r)
    return allowed, bad, counts


def primes_up_to(n: int) -> list[int]:
    out = []
    for q in range(2, n + 1):
        if all(q % d for d in range(2, int(math.isqrt(q)) + 1)):
            out.append(q)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--primes", help="comma-separated primes; overrides --prime-bound")
    ap.add_argument("--prime-bound", type=int, default=43)
    ap.add_argument("--write-masks", type=Path)
    args = ap.parse_args()

    if args.primes:
        primes = [int(x) for x in args.primes.split(",") if x]
    else:
        primes = [p for p in primes_up_to(args.prime_bound) if p not in (2, 3)]

    lines = [
        "# p : affine_allowed ; affine_bad ; good_contact_counts(r=count)",
        "# The projective infinity residue is always conservatively allowed.",
    ]
    for p in primes:
        allowed, bad, counts = prime_mask(p)
        good_extra = [r for r, n in sorted(counts.items()) if n >= 2]
        count_text = ",".join(f"{r}={counts[r]}" for r in sorted(counts))
        line = (
            f"{p} : {','.join(map(str, allowed))} ; {','.join(map(str, bad))} ; "
            f"{count_text}"
        )
        lines.append(line)
        print(
            f"p={p} allowed={allowed} bad={bad} "
            f"good_extra_contact={good_extra} "
            f"projective_density={(len(allowed)+1)/(p+1):.6f}"
        )

    if args.write_masks:
        args.write_masks.parent.mkdir(parents=True, exist_ok=True)
        args.write_masks.write_text("\n".join(lines) + "\n")
        print(f"wrote {args.write_masks}")


if __name__ == "__main__":
    main()
