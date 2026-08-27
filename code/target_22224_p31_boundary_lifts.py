#!/usr/bin/env python3
"""Exhaustive normalized p=31 boundary lifts for the q-square cover.

Every finite double-square boundary base is lifted to 31^2.  Internal
collision classes that remain deep modulo 31^2 are then exhaustively lifted
to 31^3.  External collisions are checked against

    S - x R = -(x-1)^3/x,

which explains their automatic cubic thickness.  This is a base-cover and
stable-depth classifier; it does not pretend that a singular special fibre
is a finite Jacobian.
"""

from collections import Counter, defaultdict
from pathlib import Path

p = 31


def inv(a, m):
    return pow(a, -1, m)


def val(a, m, cap):
    a %= m
    if a == 0:
        return cap
    v = 0
    while a % p == 0:
        a //= p
        v += 1
    return min(v, cap)


def square_mod_pk(a, k):
    m = p**k
    a %= m
    if a == 0:
        return True
    v = val(a, m, k)
    if v % 2:
        return False
    return pow((a // (p**v)) % p, (p - 1) // 2, p) == 1


def data(A, B, k):
    m = p**k
    A %= m
    B %= m
    C = inv(A * B % m, m)
    x, y, z = A * A % m, B * B % m, C * C % m
    R = (x + y + z - 3) % m
    S = (inv(x, m) + inv(y, m) + inv(z, m) - 3) % m
    w = S * inv(R, m) % m if R % p else None
    return x, y, z, w, R, S


labels = ("x=y", "x=z", "x=w", "y=z", "y=w", "z=w")
pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))


def signature(A, B):
    x, y, z, w, R, S = data(A, B, 1)
    if R == 0 or S == 0:
        return "sheet"
    q = (x, y, z, w)
    hit = [labels[i] for i, (a, b) in enumerate(pairs) if q[a] == q[b]]
    return "+".join(hit) if hit else "smooth"


def main():
    root = Path("results")
    log_path = root / "target_22224_p31_boundary_lifts.log"
    tsv_path = root / "target_22224_p31_boundary_lifts_summary.tsv"
    bases = []
    for A in range(1, p):
        for B in range(1, p):
            x, y, z, w, R, S = data(A, B, 1)
            if not square_mod_pk(R, 1) or not square_mod_pk(S, 1):
                continue
            sig = signature(A, B)
            if sig != "smooth":
                bases.append((A, B, sig))

    counts = Counter()
    per_base = defaultdict(Counter)
    deep_internal = []
    external_identity_fail = 0
    for A0, B0, sig in bases:
        for a1 in range(p):
            A = A0 + p * a1
            for b1 in range(p):
                B = B0 + p * b1
                x, y, z, w, R, S = data(A, B, 2)
                if not square_mod_pk(R, 2) or not square_mod_pk(S, 2):
                    counts[(sig, "p2_nonsquare")] += 1
                    per_base[(A0, B0)]["p2_nonsquare"] += 1
                    continue
                counts[(sig, "p2_compatible")] += 1
                per_base[(A0, B0)]["p2_compatible"] += 1
                if sig == "sheet":
                    rv, sv = val(R, p**2, 2), val(S, p**2, 2)
                    counts[(sig, f"vR{rv}_vS{sv}")] += 1
                    continue
                idx = labels.index(sig)
                q = (x, y, z, w)
                depth = val(q[pairs[idx][0]] - q[pairs[idx][1]], p**2, 2)
                counts[(sig, f"depth{depth}")] += 1
                per_base[(A0, B0)][f"depth{depth}"] += 1
                if sig in ("x=y", "x=z", "y=z") and depth == 2:
                    deep_internal.append((A, B, sig))
                if sig in ("x=w", "y=w", "z=w"):
                    # Select the appropriate square coordinate and verify the
                    # cubic identity exactly modulo p^2.
                    xx = {"x=w": x, "y=w": y, "z=w": z}[sig]
                    lhs = (S - xx * R) % (p**2)
                    rhs = (-(xx - 1) ** 3 * inv(xx, p**2)) % (p**2)
                    if lhs != rhs:
                        external_identity_fail += 1

    # Exhaustive p^3 children only above internal p^2 classes that remained
    # on the collision divisor.  The other internal classes already have
    # exact thickness one.
    p3 = Counter()
    for A2, B2, sig in deep_internal:
        idx = labels.index(sig)
        for a2 in range(p):
            A = A2 + p**2 * a2
            for b2 in range(p):
                B = B2 + p**2 * b2
                x, y, z, w, R, S = data(A, B, 3)
                assert square_mod_pk(R, 3) and square_mod_pk(S, 3)
                q = (x, y, z, w)
                depth = val(q[pairs[idx][0]] - q[pairs[idx][1]], p**3, 3)
                p3[(sig, depth)] += 1

    # For every external p^2 lift, the cubic identity forces depth at least
    # three at the canonical p^3 child.  This verifies the implementation as
    # well as the symbolic argument without expanding 168 million children.
    external_p3_checks = 0
    external_p3_fail = 0
    for A0, B0, sig in bases:
        if sig not in ("x=w", "y=w", "z=w"):
            continue
        for a1 in range(p):
            A = A0 + p * a1
            for b1 in range(p):
                B = B0 + p * b1
                x, y, z, w, R, S = data(A, B, 3)
                q = (x, y, z, w)
                idx = labels.index(sig)
                external_p3_checks += 1
                if val(q[pairs[idx][0]] - q[pairs[idx][1]], p**3, 3) < 3:
                    external_p3_fail += 1

    with tsv_path.open("w") as f:
        f.write("A0\tB0\tsignature\tp2_compatible\tp2_nonsquare\tdepth1\tdepth2plus\n")
        for A0, B0, sig in bases:
            c = per_base[(A0, B0)]
            f.write(f"{A0}\t{B0}\t{sig}\t{c['p2_compatible']}\t{c['p2_nonsquare']}\t{c['depth1']}\t{c['depth2']}\n")

    lines = []
    lines.append(f"P31_BOUNDARY_BASES total {len(bases)} signatures {dict(Counter(s for _,_,s in bases))}")
    for key in sorted(counts):
        lines.append(f"P2 {key[0]} {key[1]} {counts[key]}")
    lines.append(f"P2_EXTERNAL_IDENTITY failures {external_identity_fail}")
    lines.append(f"P3_INTERNAL_DEEP parents {len(deep_internal)} children {sum(p3.values())}")
    for key in sorted(p3):
        lines.append(f"P3 {key[0]} depth{key[1]} {p3[key]}")
    lines.append(f"P3_EXTERNAL_CANONICAL checks {external_p3_checks} failures {external_p3_fail} depth_at_least_3_by_cubic_identity 1")
    lines.append("SCOPE p2 is exhaustive for every normalized finite p31 boundary base; p3 is exhaustive above internal depth>=2 classes and symbolic+canonical-check for external cubic-thick classes")
    lines.append(f"OUTPUT {tsv_path}")
    text = "\n".join(lines) + "\n"
    print(text, end="")
    log_path.write_text(text)


if __name__ == "__main__":
    main()
