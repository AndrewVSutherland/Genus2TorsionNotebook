#!/usr/bin/env python3
"""Raw p=13 boundary intersection for A(2,2,2,8) plus 3-torsion.

The full A(2,2,2,8) cover in signed branch coordinates (a,b,c,d)
requires the four radicands

    abcd,
    a(a+b)(a+c)(a+d),
    b(b+a)(b+c)(b+d),
    c(c+a)(c+b)(c+d)

to be squares.  At p=13 the smooth cover has no rational 3-contact
point, so this script allows zero radicands and singular branch reduction,
then intersects the full cover with the raw cubic-contact equations.

It is deliberately finite and dependency-free.  Run from torsion_jac:

    python3 code/target_22224_full_boundary13.py
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import product


P = 13
SQUARES = {(x * x) % P for x in range(P)}


def inv(x: int) -> int:
    return pow(x % P, P - 2, P)


def trim(f: list[int]) -> list[int]:
    while len(f) > 1 and f[-1] % P == 0:
        f.pop()
    return [x % P for x in f]


def poly_mod(a: list[int], b: list[int]) -> list[int]:
    a = trim(a[:])
    b = trim(b[:])
    ib = inv(b[-1])
    while len(a) >= len(b) and not (len(a) == 1 and a[0] == 0):
        coeff = a[-1] * ib % P
        shift = len(a) - len(b)
        for i, bi in enumerate(b):
            a[i + shift] = (a[i + shift] - coeff * bi) % P
        a = trim(a)
    return a


def poly_gcd(a: list[int], b: list[int]) -> list[int]:
    a, b = trim(a), trim(b)
    while not (len(b) == 1 and b[0] == 0):
        a, b = b, poly_mod(a, b)
    ia = inv(a[-1])
    return [(x * ia) % P for x in a]


def contact_witnesses() -> dict[tuple[int, int, int, int], list[tuple[int, int, int]]]:
    """Return all raw (L,U,v) witnesses indexed by (e1,e2,e3,e4)."""
    out: dict[tuple[int, int, int, int], list[tuple[int, int, int]]] = defaultdict(list)
    i2, i4, i16 = inv(2), inv(4), inv(16)
    for L in range(1, P):
        M = L * L % P
        iM = inv(M)
        for U, v, e1 in product(range(P), repeat=3):
            A = (2 * M * e1 + 6 * (U * U + v * v) - (M + 3 * U) ** 2) % P
            e2 = (
                ((M + 3 * U) * A + 8 * v**3 - 4 * U**3 - 24 * U * v * v)
                * i4
                * iM
            ) % P
            e3 = (
                (A * A + 16 * (M + 3 * U) * v**3 - 48 * (U * U * v * v + v**4))
                * i16
                * iM
            ) % P
            e4 = ((A * v**3 - 6 * U * v**4) * i2 * iM) % P
            out[(e1, e2, e3, e4)].append((L, U, v))
    return out


def radicands(vals: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    a, b, c, d = vals
    return (
        a * b * c * d % P,
        a * (a + b) * (a + c) * (a + d) % P,
        b * (b + a) * (b + c) * (b + d) % P,
        c * (c + a) * (c + b) * (c + d) % P,
    )


def curve_key(vals: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    A, B, C, D = (x * x % P for x in vals)
    return (
        (A + B + C + D) % P,
        (A * B + A * C + A * D + B * C + B * D + C * D) % P,
        (A * B * C + A * B * D + A * C * D + B * C * D) % P,
        A * B * C * D % P,
    )


def labels(vals: tuple[int, int, int, int]) -> tuple[str, ...]:
    ans: list[str] = []
    for i, x in enumerate(vals, 1):
        if x == 0:
            ans.append(f"Z{i}")
    for i in range(4):
        for j in range(i + 1, 4):
            if (vals[i] - vals[j]) % P == 0:
                ans.append(f"E{i+1}{j+1}+")
            if (vals[i] + vals[j]) % P == 0:
                ans.append(f"E{i+1}{j+1}-")
    return tuple(sorted(ans)) if ans else ("open",)


def projective_representatives():
    for vals in product(range(P), repeat=4):
        if vals == (0, 0, 0, 0):
            continue
        first = next(x for x in vals if x)
        if first == 1:
            yield vals


def cover_flags(key: tuple[int, int, int, int], witness: tuple[int, int, int]) -> tuple[str, ...]:
    _, U, v = witness
    flags: list[str] = []
    if v == 0:
        flags.append("v=0")
    if (U * U - 4 * v * v) % P == 0:
        flags.append("discq=0")
    e1, e2, e3, e4 = key
    # Constant-term-first coefficients.
    q = [v * v % P, U, 1]
    f = [0, e4, e3, e2, e1, 1]
    if len(poly_gcd(q, f)) > 1:
        flags.append("gcd(q,f)")
    return tuple(flags) if flags else ("contact_open",)


def main() -> None:
    witnesses = contact_witnesses()
    base_counts = Counter()
    contact_counts = Counter()
    contact_open_counts = Counter()
    flag_counts = Counter()
    component_bases: dict[str, set[tuple[int, int, int, int]]] = defaultdict(set)
    component_open_bases: dict[str, set[tuple[int, int, int, int]]] = defaultdict(set)
    total_projective = total_cover = total_open_cover = 0
    bases_with_contact: set[tuple[int, int, int, int]] = set()
    bases_with_open_contact: set[tuple[int, int, int, int]] = set()
    samples: list[tuple] = []

    for vals in projective_representatives():
        total_projective += 1
        rads = radicands(vals)
        if any(r not in SQUARES for r in rads):
            continue
        total_cover += 1
        labs = labels(vals)
        sig = "+".join(labs)
        base_counts[sig] += 1
        if labs == ("open",):
            total_open_cover += 1
        key = curve_key(vals)
        ws = witnesses.get(key, ())
        if not ws:
            continue
        bases_with_contact.add(vals)
        for lab in labs:
            component_bases[lab].add(vals)
        for w in ws:
            fl = cover_flags(key, w)
            fsig = "+".join(fl)
            contact_counts[sig] += 1
            flag_counts[fsig] += 1
            is_open = fl == ("contact_open",)
            if is_open:
                bases_with_open_contact.add(vals)
                contact_open_counts[sig] += 1
                for lab in labs:
                    component_open_bases[lab].add(vals)
                if len(samples) < 40:
                    samples.append((vals, labs, rads, key, w))

    print("TARGET_22224_FULL_BOUNDARY13")
    print("contact_keys", len(witnesses))
    print("projective_bases", total_projective)
    print("full_cover_bases_including_boundary", total_cover)
    print("smooth_full_cover_bases", total_open_cover)
    print("bases_with_raw_contact", len(bases_with_contact))
    print("bases_with_contact_open", len(bases_with_open_contact))
    print("raw_contact_presentations", sum(contact_counts.values()))
    print("contact_open_presentations", sum(contact_open_counts.values()))
    print("flag_counts")
    for key, count in sorted(flag_counts.items()):
        print(" ", key, count)
    print("component_base_counts")
    all_components = sorted(set(component_bases) | set(component_open_bases))
    for lab in all_components:
        print(" ", lab, len(component_bases[lab]), len(component_open_bases[lab]))
    print("live_signature_counts")
    for sig, count in sorted(contact_open_counts.items()):
        print(" ", sig, "bases", base_counts[sig], "open_contacts", count)
    print("samples")
    for sample in samples:
        print(" ", sample)


if __name__ == "__main__":
    main()
