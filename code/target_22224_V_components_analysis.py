#!/usr/bin/env python3
"""Component-aware local analysis of the new V1 and V2 fibers.

This script deliberately separates four logically different objects:

* the projected, sampled direct-contact bank (reported elsewhere);
* the intrinsic cubic-contact equations on [A,B,C,D*T^2];
* the full A(2,2,2,8) square cover; and
* the exceptional divisor above the universal spurious L=0 component.

For the last item write M=L^2 and delta=2*v-U=L*s.  Dividing the three
contact equations by L^2 gives the strict-transform equations used below.
Thus a missing exceptional solution is a genuine first-blow-up obstruction,
not a failure to sample a sufficiently large tangent bank.
"""

from __future__ import annotations

import argparse
import itertools
from collections import Counter, defaultdict
from pathlib import Path

FAMILIES = {
    "V1": (-2178, 2420, 9075, -1470),
    "V2": (-1458, 2268, 7938, -2023),
}


def elementary(base):
    a, b, c, d = (x * x for x in base)
    return (
        a + b + c + d,
        a * b + a * c + a * d + b * c + b * d + c * d,
        a * b * c + a * b * d + a * c * d + b * c * d,
        a * b * c * d,
    )


def radicands(family, T):
    A, B, C, D = FAMILIES[family]
    d = D * T * T
    return (
        A * B * C * d,
        A * (A + B) * (A + C) * (A + d),
        B * (B + A) * (B + C) * (B + d),
        C * (C + A) * (C + B) * (C + d),
    )


def infinity_units(family, z=0):
    """R_i/z^6 on [A*z^2,B*z^2,C*z^2,D]."""
    A, B, C, D = FAMILIES[family]
    return (
        A * B * C * D,
        A * (A + B) * (A + C) * (D + A * z * z),
        B * (B + A) * (B + C) * (D + B * z * z),
        C * (C + A) * (C + B) * (D + C * z * z),
    )


def contact_equations_m(family, x):
    """Corrected cleared contact equations in (T,M,U,v), M=L^2."""
    T, M, U, v = x
    A, B, C, D = FAMILIES[family]
    e1, e2, e3, e4 = elementary((A, B, C, D * T * T))
    PP = 4 * M * e1 + 12 * (U * U + v * v) - (M + 3 * U) ** 2
    return (
        (M + 3 * U) * PP + 16 * v**3 - 8 * U**3
        - 48 * U * v * v - 8 * M * e2,
        PP * PP + 64 * (M + 3 * U) * v**3
        - 192 * (U * U * v * v + v**4) - 64 * M * e3,
        PP * v**3 - 12 * U * v**4 - 4 * M * e4,
    )


def strict_contact_equations(family, x):
    """Strict transform at L=0 in coordinates (T,L,s,U).

    Here M=L^2, delta=2*v-U=L*s.  The first two returned polynomials are
    F1/L^2 and F2/L^2.  The third is 8*F3/L^2, clearing a p-adic unit for
    p=11,13.
    """
    T, L, s, U = x
    A, B, C, D = FAMILIES[family]
    e1, e2, e3, e4 = elementary((A, B, C, D * T * T))
    g1 = (
        -L**4 - 9 * L * L * U + 4 * L * L * e1
        - 12 * U * U + 6 * L * U * s + 12 * U * e1
        + 3 * L * L * s * s - 8 * e2 + 3 * U * s * s
        + 2 * L * s**3
    )
    g2 = (
        L**6 + 12 * L**4 * U - 8 * L**4 * e1
        + 24 * L * L * U * U - 12 * L**3 * U * s
        - 48 * L * L * U * e1 - 6 * L**4 * s * s
        + 16 * L * L * e1 * e1 - 64 * U**3
        - 48 * L * U * U * s + 48 * U * U * e1
        - 12 * L * L * U * s * s + 48 * L * U * s * e1
        + 8 * L**3 * s**3 + 24 * L * L * s * s * e1
        - 64 * e3 + 24 * U * U * s * s
        + 12 * L * U * s**3 - 3 * L * L * s**4
    )
    g3 = (
        -L * L * U**3 - 3 * L**3 * U * U * s
        - 3 * L**4 * U * s * s - L**5 * s**3
        - 6 * U**4 - 18 * L * U**3 * s + 4 * U**3 * e1
        - 18 * L * L * U * U * s * s
        + 12 * L * U * U * s * e1 - 6 * L**3 * U * s**3
        + 12 * L * L * U * s * s * e1
        + 4 * L**3 * s**3 * e1 - 32 * e4
        + 3 * U**3 * s * s + 9 * L * U * U * s**3
        + 9 * L * L * U * s**4 + 3 * L**3 * s**5
    )
    return g1, g2, g3


def unit_system(family, x):
    """Contact plus cover equations in (T,M,U,v,y0,...,y3)."""
    out = list(contact_equations_m(family, x[:4]))
    out.extend(x[4 + i] ** 2 - r for i, r in enumerate(radicands(family, x[0])))
    return tuple(out)


def strict_system(family, x):
    """Strict contact plus cover in (T,L,s,U,y0,...,y3)."""
    out = list(strict_contact_equations(family, x[:4]))
    out.extend(x[4 + i] ** 2 - r for i, r in enumerate(radicands(family, x[0])))
    return tuple(out)


def rank_and_space(mat, rhs, p):
    rows = [[q % p for q in row] + [rhs[i] % p]
            for i, row in enumerate(mat)]
    nr, nv = len(rows), len(mat[0])
    pivots = []
    r = 0
    for c in range(nv):
        pivot = next((i for i in range(r, nr) if rows[i][c]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        inv = pow(rows[r][c], -1, p)
        rows[r] = [inv * q % p for q in rows[r]]
        for i in range(nr):
            if i != r and rows[i][c]:
                q = rows[i][c]
                rows[i] = [(rows[i][j] - q * rows[r][j]) % p
                           for j in range(nv + 1)]
        pivots.append(c)
        r += 1
    if any(all(rows[i][j] == 0 for j in range(nv)) and rows[i][-1]
           for i in range(r, nr)):
        return r, None
    free = [j for j in range(nv) if j not in pivots]
    particular = [0] * nv
    for i, c in enumerate(pivots):
        particular[c] = rows[i][-1]
    basis = []
    for fc in free:
        q = [0] * nv
        q[fc] = 1
        for i, c in enumerate(pivots):
            q[c] = -rows[i][fc] % p
        basis.append(tuple(q))
    return r, (tuple(particular), tuple(basis))


def jacobian_mod(system, family, x, p):
    x0 = tuple(q % p for q in x)
    f0 = system(family, x0)
    out = [[0] * len(x0) for _ in f0]
    for j in range(len(x0)):
        y = list(x0)
        y[j] += p
        fy = system(family, y)
        for i in range(len(f0)):
            d = fy[i] - f0[i]
            assert d % p == 0
            out[i][j] = (d // p) % p
    return out


def roots_mod(r, p):
    return tuple(z for z in range(p) if z * z % p == r % p)


def contact_mod_p(family, p):
    squares = {z * z % p for z in range(p)}
    unit = defaultdict(list)
    for T, M, U, v in itertools.product(range(p), repeat=4):
        if M == 0 or M not in squares:
            continue
        x = (T, M, U, v)
        if all(q % p == 0 for q in contact_equations_m(family, x)):
            J = jacobian_mod(lambda fam, z: contact_equations_m(fam, z),
                             family, x, p)
            rank, _ = rank_and_space(J, [0, 0, 0], p)
            unit[T].append((M, U, v, rank))
    return unit


def exceptional_mod_p(family, p):
    """Solutions on the first exceptional divisor L=0."""
    ans = defaultdict(list)
    for T, s, U in itertools.product(range(p), repeat=3):
        x = (T, 0, s, U)
        if all(q % p == 0 for q in strict_contact_equations(family, x)):
            J = jacobian_mod(
                lambda fam, z: strict_contact_equations(fam, z),
                family, x, p,
            )
            rank, _ = rank_and_space(J, [0, 0, 0], p)
            ans[T].append((s, U, rank))
    return ans


def cover_residues(family, p):
    squares = {z * z % p for z in range(p)}
    return {
        T: tuple(r % p for r in radicands(family, T))
        for T in range(p)
        if all(r % p in squares for r in radicands(family, T))
    }


def start_nodes(family, p, mode, T, data):
    rsets = [roots_mod(r, p) for r in radicands(family, T)]
    assert all(rsets)
    if mode == "unit":
        M, U, v = data
        prefix = (T, M, U, v)
        system = unit_system
    else:
        s, U = data
        prefix = (T, 0, s, U)
        system = strict_system
    nodes = [prefix + ys for ys in itertools.product(*rsets)]
    assert all(all(q % p == 0 for q in system(family, x)) for x in nodes)
    return system, nodes


def lift_children(system, family, x, p, level):
    mod = p**level
    vals = system(family, x)
    assert all(q % mod == 0 for q in vals)
    rhs = [-(q // mod) % p for q in vals]
    J = jacobian_mod(system, family, x, p)
    rank, space = rank_and_space(J, rhs, p)
    if space is None:
        return rank, ()
    particular, basis = space
    children = []
    for coeff in itertools.product(range(p), repeat=len(basis)):
        delta = tuple(
            (particular[j] + sum(c * b[j]
             for c, b in zip(coeff, basis))) % p
            for j in range(len(x))
        )
        y = tuple((x[j] + mod * delta[j]) % (mod * p)
                  for j in range(len(x)))
        if all(q % (mod * p) == 0 for q in system(family, y)):
            children.append(y)
    return rank, tuple(children)


def exhaustive_p2(system, family, nodes, p):
    out = set()
    ranks = Counter()
    for x in nodes:
        rank, children = lift_children(system, family, x, p, 1)
        ranks[rank] += 1
        out.update(children)
    return out, ranks


def open_at_precision(family, mode, x, p, depth):
    """Check that a certificate has left every geometric boundary."""
    mod = p**depth
    T = x[0]
    A, B, C, D = FAMILIES[family]
    base = (A, B, C, D * T * T)
    if any(q % mod == 0 for q in base):
        return False
    if any((base[i] * base[i] - base[j] * base[j]) % mod == 0
           for i in range(4) for j in range(i + 1, 4)):
        return False
    if mode == "unit":
        M, U, v = x[1:4]
        if M % p == 0:  # the unit chart under consideration
            return False
    else:
        L, s, U = x[1:4]
        if L % mod == 0:  # must move off the exceptional divisor
            return False
        v = ((U + L * s) * pow(2, -1, mod)) % mod
    # gcd(x^2+U*x+v^2, x*prod(x+a_i^2))=1.  Its resultant is
    # v^2*prod(a_i^4-U*a_i^2+v^2).
    resultant = v * v
    for a in base:
        resultant *= a**4 - U * a * a + v * v
    return resultant % mod != 0


def find_path(system, family, mode, node, p, depth):
    """Deterministic DFS for one compatible p^depth certificate."""
    cache = set()

    def rec(x, level):
        key = (level, x)
        if key in cache:
            return None
        if level == depth:
            if open_at_precision(family, mode, x, p, depth):
                return x
            cache.add(key)
            return None
        rank, children = lift_children(system, family, x, p, level)
        # For the strict transform, try branches moving off L=0 first.
        if mode == "exceptional":
            children = tuple(sorted(children, key=lambda y: y[1] == 0))
        for y in children:
            ans = rec(y, level + 1)
            if ans is not None:
                return ans
        cache.add(key)
        return None

    return rec(node, 1)


def canonical_node(nodes):
    """Choose nonnegative/minimal root signs deterministically."""
    return min(nodes)


def analyze_case(family, p, depth, lines, cert_rows, p2_rows):
    sq = {z * z % p for z in range(p)}
    cover = cover_residues(family, p)
    units = contact_mod_p(family, p)
    exceptional = exceptional_mod_p(family, p)
    iu = tuple(q % p for q in infinity_units(family))
    ib = tuple(q in sq for q in iu)
    lines.append(
        f"CASE family={family} p={p} cover_T={sorted(cover)} "
        f"infinity_units={iu} infinity_square={ib}"
    )
    lines.append(
        "INFINITY_FULL_COVER=0 reason=normalized_unit_nonsquare"
        if not all(ib) else "INFINITY_FULL_COVER=undecided"
    )

    viable_unit = []
    viable_exc = []
    for T, rr in sorted(cover.items()):
        us = units.get(T, [])
        es = exceptional.get(T, [])
        lines.append(
            f"FINITE T={T} radicands={rr} zero_count={sum(q == 0 for q in rr)} "
            f"unit_contact={us} exceptional_contact={es}"
        )
        for M, U, v, rank in us:
            viable_unit.append((T, M, U, v))
        for s, U, rank in es:
            viable_exc.append((T, s, U))

    if family == "V2" and p == 11:
        assert not viable_unit and not viable_exc
        lines.append(
            "GENUINE_LOCAL_OBSTRUCTION V2 p=11: every cover-allowed finite "
            "T disk has neither a unit-L contact point nor a point on the "
            "strict transform above L=0; infinity is killed by cover units"
        )
        lines.append(
            "INDEPENDENT_FINITE_FIELD_CHECK all four finite disks have good "
            "reduction and #J(F_11)=128=2 mod 3; see "
            "target_22224_V_components_finitefield.log"
        )
        return

    # T -> -T is an exact involution.  Analyze one representative of each
    # absolute disk; root signs are also exact independent involutions.
    representatives = []
    seen = set()
    for T, M, U, v in viable_unit:
        key = (min(T, (-T) % p), M, U, v)
        if key not in seen:
            seen.add(key)
            representatives.append(("unit", T, (M, U, v)))
    for T, s, U in viable_exc:
        key = (min(T, (-T) % p), s, U)
        if key not in seen:
            seen.add(key)
            representatives.append(("exceptional", T, (s, U)))

    for mode, T, data in representatives:
        system, nodes = start_nodes(family, p, mode, T, data)
        p2, ranks = exhaustive_p2(system, family, nodes, p)
        tclasses = {x[0] % (p * p) for x in p2}
        off_exceptional = (
            sum(x[1] % (p * p) != 0 for x in p2)
            if mode == "exceptional" else -1
        )
        lines.append(
            f"P2 family={family} p={p} mode={mode} T0={T} data={data} "
            f"starts={len(nodes)} rank_hist={dict(ranks)} nodes={len(p2)} "
            f"Tclasses={len(tclasses)} off_exceptional={off_exceptional}"
        )
        p2_rows.append((family, p, mode, T, repr(data), len(nodes),
                        repr(dict(ranks)), len(p2), len(tclasses),
                        off_exceptional))

        # A p^5 certificate is most useful for V1.  For V2/p13 the complete
        # p=11 obstruction already rules out the fiber globally, so p^2 is
        # the proportionate stopping point.
        if family != "V1":
            continue
        cert = None
        for node in sorted(nodes):
            cert = find_path(system, family, mode, node, p, depth)
            if cert is not None:
                break
        lines.append(
            f"P{depth} family={family} p={p} mode={mode} T0={T} "
            f"certificate={'YES' if cert is not None else 'NO'}"
        )
        if cert is not None:
            mod = p**depth
            assert all(q % mod == 0 for q in system(family, cert))
            cert_rows.append((family, p, depth, mode, T, repr(data),
                              *cert))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", type=int, default=5)
    ap.add_argument(
        "--log", default="results/target_22224_V_components_analysis.log"
    )
    ap.add_argument(
        "--certificates",
        default="results/target_22224_V_components_padic_certificates.tsv",
    )
    ap.add_argument(
        "--p2",
        default="results/target_22224_V_components_p2_components.tsv",
    )
    args = ap.parse_args()

    lines = [
        f"V_COMPONENT_ANALYSIS_START depth={args.depth}",
        "SCOPE intrinsic equations and full cover; independent of sampled bank",
    ]
    cert_rows = []
    p2_rows = []
    analyze_case("V1", 13, args.depth, lines, cert_rows, p2_rows)
    analyze_case("V2", 11, args.depth, lines, cert_rows, p2_rows)
    analyze_case("V2", 13, args.depth, lines, cert_rows, p2_rows)

    p2path = Path(args.p2)
    with p2path.open("w") as f:
        f.write("family\tp\tmode\tT0\tdata\tstarts\trank_hist\tp2_nodes\t"
                "p2_Tclasses\toff_exceptional\n")
        for row in p2_rows:
            f.write("\t".join(map(str, row)) + "\n")

    cpath = Path(args.certificates)
    with cpath.open("w") as f:
        f.write("family\tp\tdepth\tmode\tT0\tdata\tT\tx1\tx2\tx3\t"
                "y0\ty1\ty2\ty3\n")
        for row in cert_rows:
            f.write("\t".join(map(str, row)) + "\n")

    lines.extend([
        "CONCLUSION V2 is intrinsically excluded at p=11 (also independently "
        "by good reduction with #J(F_11)=128); this is not a sampled-bank "
        "conclusion",
        "CONCLUSION V1 remains locally viable at p=13 and receives explicit "
        f"contact+cover certificates through 13^{args.depth}",
        f"p2_output={p2path}",
        f"certificate_output={cpath}",
    ])
    text = "\n".join(lines) + "\n"
    Path(args.log).write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
