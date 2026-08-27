#!/usr/bin/env python3
"""Mine the primitive off-rectangle tor2228 bank for low-dimensional structure.

The input bank is projective, so every row is first made primitive.  We retain
only the genuinely off-rectangle rows and compute:

* shared coordinate triples and pairs;
* the inverse-Cremona/Adam first-leg coordinates for every ordering;
* symmetric pair-partition invariants of the reciprocal coordinates;
* exact low-degree interpolation tests in several natural charts.

This is an exact rational computation; no floating-point equality is used.
"""

from __future__ import annotations

from ast import literal_eval
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import gcd, isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BANK = ROOT / "data/tor2228_bank.txt"
OUT = ROOT / "results/target_22224_offrectangle_new_curves_analysis.txt"


def primitive(v: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for z in v:
        g = gcd(g, abs(z))
    if g == 0:
        return v
    w = tuple(z // g for z in v)
    for z in w:
        if z:
            return w if z > 0 else tuple(-x for x in w)
    return w


def rectangle(v: tuple[int, int, int, int]) -> bool:
    a, b, c, d = v
    return a * b == c * d or a * c == b * d or a * d == b * c


def is_square_int(n: int) -> bool:
    if n < 0:
        return False
    r = isqrt(n)
    return r * r == n


def is_square_q(x: Q) -> bool:
    return x >= 0 and is_square_int(x.numerator) and is_square_int(x.denominator)


def sqrt_q(x: Q) -> Q | None:
    if not is_square_q(x):
        return None
    return Q(isqrt(x.numerator), isqrt(x.denominator))


def load_rows() -> list[tuple[int, int, int, int]]:
    seen: set[tuple[int, int, int, int]] = set()
    rows: list[tuple[int, int, int, int]] = []
    for raw in BANK.read_text().splitlines():
        if not raw.startswith("["):
            continue
        z = literal_eval(raw)
        if not isinstance(z, list) or len(z) != 4:
            continue
        v = primitive(tuple(int(x) for x in z))
        if v in seen:
            continue
        seen.add(v)
        if not rectangle(v):
            rows.append(v)
    return rows


def square_roots(v: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    a, b, c, d = v
    rad = (
        a * b * c * d,
        a * (a + b) * (a + c) * (a + d),
        b * (a + b) * (b + c) * (b + d),
        c * (a + c) * (b + c) * (c + d),
    )
    assert all(is_square_int(z) for z in rad), (v, rad)
    return tuple(isqrt(z) for z in rad)  # type: ignore[return-value]


def first_leg_inverse(v: tuple[int, int, int, int]):
    """Inverse Adam coordinates for this *ordered* projective tuple.

    Returns (u,v,w2,w) when the symmetric ten-factor character is square,
    otherwise (u,v,w2,None).  Here w is the nonnegative rational square root.
    """
    x = [Q(1, z) for z in v]
    s = sum(x)
    x = [z / s for z in x]
    A, B, C = x[0] + x[1], x[0] + x[2], x[0] + x[3]
    u = A / (1 - A)
    vv = B / (1 - B)
    w2 = C * u * vv / (1 - C)
    return u, vv, w2, sqrt_q(w2)


def qheight(x: Q) -> int:
    return max(abs(x.numerator), x.denominator)


def canonical_adam(v: tuple[int, int, int, int]):
    candidates = []
    for p in permutations(v):
        u, vv, w2, w = first_leg_inverse(p)
        if w is None:
            continue
        for ws in (w, -w):
            score = (max(qheight(u), qheight(vv), qheight(ws)),
                     qheight(u) + qheight(vv) + qheight(ws),
                     u, vv, ws, p)
            candidates.append((score, (u, vv, w2, ws, p)))
    return min(candidates)[1] if candidates else None


def partition_invariants(v: tuple[int, int, int, int]):
    """Canonical ratios of complementary pair sums of reciprocal coords."""
    x = [Q(1, z) for z in v]
    parts = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    vals = []
    for ij, kl in parts:
        r = (x[ij[0]] + x[ij[1]]) / (x[kl[0]] + x[kl[1]])
        vals.append(min(r, 1 / r))
    return tuple(sorted(vals))


def monomials(nvars: int, degree: int):
    return [e for e in product(range(degree + 1), repeat=nvars) if sum(e) <= degree]


def nullspace(matrix: list[list[Q]]) -> list[list[Q]]:
    if not matrix:
        return []
    A = [row[:] for row in matrix]
    m, n = len(A), len(A[0])
    pivots: list[int] = []
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, m) if A[i][c]), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        z = A[r][c]
        A[r] = [q / z for q in A[r]]
        for i in range(m):
            if i != r and A[i][c]:
                z = A[i][c]
                A[i] = [A[i][j] - z * A[r][j] for j in range(n)]
        pivots.append(c)
        r += 1
        if r == m:
            break
    free = [c for c in range(n) if c not in pivots]
    out = []
    for f in free:
        vec = [Q(0) for _ in range(n)]
        vec[f] = Q(1)
        for i in range(len(pivots) - 1, -1, -1):
            c = pivots[i]
            vec[c] = -sum(A[i][j] * vec[j] for j in free)
        out.append(vec)
    return out


def interpolate(points: list[tuple[Q, ...]], degree: int):
    mons = monomials(len(points[0]), degree)
    M = []
    for pt in points:
        M.append([product_value(pt, e) for e in mons])
    return mons, nullspace(M)


def product_value(pt: tuple[Q, ...], e: tuple[int, ...]) -> Q:
    z = Q(1)
    for x, n in zip(pt, e):
        z *= x**n
    return z


def format_poly(mons: list[tuple[int, ...]], coeffs: list[Q], names: tuple[str, ...]) -> str:
    terms = []
    for c, e in zip(coeffs, mons):
        if not c:
            continue
        factors = []
        for name, n in zip(names, e):
            if n == 1:
                factors.append(name)
            elif n > 1:
                factors.append(f"{name}^{n}")
        mon = "*".join(factors) or "1"
        terms.append(f"({c})*{mon}")
    return " + ".join(terms) if terms else "0"


def repeated_features(label: str, mapping: dict[object, set[int]], rows, lines, minimum=2):
    hits = [(k, sorted(ix)) for k, ix in mapping.items() if len(ix) >= minimum]
    hits.sort(key=lambda z: (-len(z[1]), str(z[0])))
    lines.append(f"\n[{label}] repeated={len(hits)}")
    for k, ix in hits[:100]:
        lines.append(f"  value={k} rows={ix} tuples={[rows[i] for i in ix]}")


def main() -> None:
    rows = load_rows()
    lines = [f"OFFRECTANGLE_ROWS {len(rows)}"]

    triple_map: dict[object, set[int]] = defaultdict(set)
    pair_map: dict[object, set[int]] = defaultdict(set)
    coord_map: dict[object, set[int]] = defaultdict(set)
    for i, row in enumerate(rows):
        for z in row:
            coord_map[z].add(i)
        for p in combinations(row, 2):
            pair_map[p].add(i)
        for t in combinations(row, 3):
            triple_map[t].add(i)
    repeated_features("SHARED_TRIPLES", triple_map, rows, lines)
    repeated_features("SHARED_PAIRS", pair_map, rows, lines)

    # Projective diagonal-torus coincidences.  A 3+1 ratio pattern is a
    # shared-coordinate fiber after rescaling; a 2+2 pattern is a pair-scaling
    # fiber.  We search all coordinate matchings, not only the sorted one.
    ratio31 = []
    ratio22 = []
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            found31 = None
            found22 = None
            for w in permutations(rows[j]):
                ratios = tuple(Q(w[k], rows[i][k]) for k in range(4))
                counts = sorted(Counter(ratios).values(), reverse=True)
                if counts == [3, 1] and found31 is None:
                    found31 = (w, ratios)
                if counts == [2, 2] and found22 is None:
                    found22 = (w, ratios)
            if found31 is not None:
                ratio31.append((i, j, found31))
            if found22 is not None:
                ratio22.append((i, j, found22))
    lines.append(f"\n[PROJECTIVE_RATIO_3PLUS1] count={len(ratio31)}")
    for item in ratio31:
        lines.append(f"  {item} tuples=({rows[item[0]]},{rows[item[1]]})")
    lines.append(f"\n[PROJECTIVE_RATIO_2PLUS2] count={len(ratio22)}")
    for item in ratio22:
        lines.append(f"  {item} tuples=({rows[item[0]]},{rows[item[1]]})")

    firstleg = []
    standard = []
    all_adam_features: dict[str, dict[object, set[int]]] = {
        name: defaultdict(set) for name in
        ("u", "v", "w2", "uv", "u_plus_v", "u_minus_v", "u_over_v",
         "pair_uv", "pair_u_w2", "pair_v_w2")
    }
    partition_map: dict[object, set[int]] = defaultdict(set)
    for i, row in enumerate(rows):
        roots = square_roots(row)
        pinv = partition_invariants(row)
        partition_map[pinv].add(i)
        u0, v0, w20, w0 = first_leg_inverse(row)
        standard.append((u0, v0, w20, w0))
        can = canonical_adam(row)
        lines.append(
            f"ROW {i:02d} tuple={row} roots={roots} partitions={pinv} "
            f"standard=(u={u0},v={v0},w2={w20},w={w0}) canonical={can}"
        )
        if can is None:
            continue
        firstleg.append(i)
        u, vv, w2, w, p = can
        feats = {
            "u": u, "v": vv, "w2": w2, "uv": u * vv,
            "u_plus_v": u + vv, "u_minus_v": u - vv,
            "u_over_v": u / vv, "pair_uv": (u, vv),
            "pair_u_w2": (u, w2), "pair_v_w2": (vv, w2),
        }
        for name, value in feats.items():
            all_adam_features[name][value].add(i)

    lines.append(f"\nFIRST_LEG_SQUARE {len(firstleg)} indices={firstleg}")
    repeated_features("PARTITION_INVARIANTS", partition_map, rows, lines)
    for name, mapping in all_adam_features.items():
        repeated_features(f"ADAM_{name}", mapping, rows, lines)

    # Exact interpolation in the natural normalized branch chart.
    branch_points = [(Q(b, a), Q(c, a), Q(d, a)) for a, b, c, d in rows]
    for deg in (1, 2, 3):
        mons, ns = interpolate(branch_points, deg)
        lines.append(f"\nBRANCH_INTERPOLATION degree={deg} monomials={len(mons)} nullity={len(ns)}")
        for vec in ns[:3]:
            lines.append("  " + format_poly(mons, vec, ("r", "s", "t")))

    # Standard ordered Adam chart, only where its w^2 is a rational square.
    adam_points = [(u, v, w2) for u, v, w2, w in standard if w is not None]
    if adam_points:
        for deg in (1, 2, 3):
            mons, ns = interpolate(adam_points, deg)
            lines.append(f"\nADAM_INTERPOLATION degree={deg} points={len(adam_points)} monomials={len(mons)} nullity={len(ns)}")
            for vec in ns[:3]:
                lines.append("  " + format_poly(mons, vec, ("u", "v", "W")))

    # Fiber equations attached to every exactly shared coordinate triple.
    lines.append("\nSHARED_TRIPLE_FIBERS")
    for triple, ix in sorted(triple_map.items(), key=lambda z: (-len(z[1]), z[0])):
        if len(ix) < 2:
            continue
        variable = []
        for i in sorted(ix):
            row = rows[i]
            rest = list(row)
            for z in triple:
                rest.remove(z)
            assert len(rest) == 1
            variable.append(rest[0])
        base = variable[0]
        ratios = [Q(z, base) for z in variable]
        lines.append(f"  fixed={triple} variable={variable} ratios={ratios} square_ratios={[sqrt_q(z) for z in ratios]}")
        lines.append(
            "    genus5 fiber: d=d0*T^2 and "
            + ", ".join(f"Y_{z}^2=({z}+d0*T^2)/({z}+d0)" for z in triple)
        )
        for x, y in combinations(triple, 2):
            lines.append(
                f"    elliptic quotient E_{x}_{y}: Z^2=({x}+d0*T^2)*({y}+d0*T^2)"
            )

    OUT.write_text("\n".join(lines) + "\n")
    print(f"WROTE {OUT}")
    print(f"OFFRECTANGLE {len(rows)} FIRST_LEG_SQUARE {len(firstleg)}")
    print("SHARED_TRIPLES", [(k, sorted(v)) for k, v in triple_map.items() if len(v) > 1])


if __name__ == "__main__":
    main()
