#!/usr/bin/env python3
"""First-order 13-adic boundary analysis for full A(2,2,2,8)+3.

This is deliberately the full three-dimensional order-8 cover, not the
point-on-the-curve K3 slice.  In signed square-branch coordinates its four
cover equations are

    r0^2 = abcd,
    r1^2 = a(a+b)(a+c)(a+d),
    r2^2 = b(b+a)(b+c)(b+d),
    r3^2 = c(c+a)(c+b)(c+d).

They are coupled to the three eliminated cubic-contact equations for a
rational 3-torsion class.  We enumerate the normalized mod-13 incidence
solutions (including contact and branch boundary), compute their affine
linear lift spaces to 13^2, and classify base boundary signatures and
Jacobian ranks.  This is a Hensel-seed diagnostic, not a proof that every
p-adic chart is represented: the contact normalization L=1/m is assumed to
be a unit.
"""

from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path

P = 13


def inv(x: int, p: int = P) -> int:
    return pow(x % p, -1, p)


def radicands(base: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    a, b, c, d = base
    return (
        a * b * c * d,
        a * (a + b) * (a + c) * (a + d),
        b * (b + a) * (b + c) * (b + d),
        c * (c + a) * (c + b) * (c + d),
    )


def coeffs(base: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    A, B, C, D = (x * x for x in base)
    return (
        A + B + C + D,
        A * B + A * C + A * D + B * C + B * D + C * D,
        A * B * C + A * B * D + A * C * D + B * C * D,
        A * B * C * D,
    )


def contact_equations(
    base: tuple[int, int, int, int], L: int, U: int, v: int
) -> tuple[int, int, int]:
    e1, e2, e3, e4 = coeffs(base)
    M = L * L
    aux = 2 * M * e1 + 6 * (U * U + v * v) - (M + 3 * U) ** 2
    return (
        (M + 3 * U) * aux + 8 * v**3 - 4 * e2 * M - 4 * U**3 - 24 * U * v * v,
        aux * aux
        + 16 * (M + 3 * U) * v**3
        - 16 * e3 * M
        - 48 * (U * U * v * v + v**4),
        aux * v**3 - 2 * e4 * M - 6 * U * v**4,
    )


def equations(z: tuple[int, ...], norm_index: int) -> list[int]:
    base = tuple(z[:4])
    roots = z[4:8]
    L, U, v = z[8:11]
    return [
        *(roots[i] * roots[i] - radicands(base)[i] for i in range(4)),
        *contact_equations(base, L, U, v),
        base[norm_index] - 1,
    ]


def roots_mod(value: int, p: int = P) -> tuple[int, ...]:
    return tuple(x for x in range(p) if (x * x - value) % p == 0)


SQRTS = {x: roots_mod(x) for x in range(P)}


def normalized_bases(p: int = P):
    for base in product(range(p), repeat=4):
        first = next((i for i, x in enumerate(base) if x), None)
        if first is None or base[first] != 1:
            continue
        rr = tuple(x % p for x in radicands(base))
        if all(SQRTS[x] for x in rr):
            yield base, first, rr


def contact_is_open(base: tuple[int, int, int, int], L: int, U: int, v: int, p: int = P) -> bool:
    if L % p == 0 or v % p == 0 or (U * U - 4 * v * v) % p == 0:
        return False
    # q=x^2+Ux+v^2 must be coprime to the completely split f.
    if v * v % p == 0:
        return False
    for x in base:
        A = x * x % p
        if (A * A - U * A + v * v) % p == 0:
            return False
    return True


def mod_solutions(p: int = P):
    sols: list[tuple[tuple[int, ...], int]] = []
    base_count = 0
    for base, norm_index, rr in normalized_bases(p):
        base_count += 1
        contacts = []
        for L in range(1, p):
            for U in range(p):
                for v in range(p):
                    if all(e % p == 0 for e in contact_equations(base, L, U, v)) and contact_is_open(base, L, U, v, p):
                        contacts.append((L, U, v))
        if not contacts:
            continue
        for roots in product(*(SQRTS[x] for x in rr)):
            for contact in contacts:
                sols.append((tuple(base + roots + contact), norm_index))
    return base_count, sols


def rref_lift_space(z: tuple[int, ...], norm_index: int, p: int = P):
    f0 = equations(z, norm_index)
    nrows, nvars = len(f0), len(z)
    rhs = [(-(f // p)) % p for f in f0]
    mat = [[0] * nvars for _ in range(nrows)]
    for j in range(nvars):
        zz = list(z)
        zz[j] += p
        fj = equations(tuple(zz), norm_index)
        for i in range(nrows):
            mat[i][j] = ((fj[i] - f0[i]) // p) % p
    aug = [mat[i] + [rhs[i]] for i in range(nrows)]
    row = 0
    pivots: list[int] = []
    for col in range(nvars):
        pivot = next((r for r in range(row, nrows) if aug[r][col] % p), None)
        if pivot is None:
            continue
        aug[row], aug[pivot] = aug[pivot], aug[row]
        scale = inv(aug[row][col], p)
        aug[row] = [(x * scale) % p for x in aug[row]]
        for r in range(nrows):
            if r != row and aug[r][col] % p:
                q = aug[r][col] % p
                aug[r] = [(aug[r][c] - q * aug[row][c]) % p for c in range(nvars + 1)]
        pivots.append(col)
        row += 1
    for r in range(row, nrows):
        if all(aug[r][c] % p == 0 for c in range(nvars)) and aug[r][-1] % p:
            return row, None
    free = [c for c in range(nvars) if c not in pivots]
    y0 = [0] * nvars
    for r, col in enumerate(pivots):
        y0[col] = aug[r][-1] % p
    basis = []
    for fc in free:
        vec = [0] * nvars
        vec[fc] = 1
        for r, col in enumerate(pivots):
            vec[col] = (-aug[r][fc]) % p
        basis.append(tuple(vec))
    return row, (tuple(y0), basis)


def signature(z: tuple[int, ...], p: int = P) -> tuple[str, ...]:
    base = z[:4]
    sq = [(x * x) % p for x in base]
    out = []
    for i, x in enumerate(base):
        if x % p == 0:
            out.append(f"Z{i+1}")
    for i in range(4):
        for j in range(i + 1, 4):
            if sq[i] == sq[j]:
                sign = "+" if (base[i] - base[j]) % p == 0 else "-"
                out.append(f"E{i+1}{j+1}{sign}")
    for i, r in enumerate(z[4:8]):
        if r % p == 0:
            out.append(f"R{i}=0")
    L, U, v = z[8:11]
    if v % p == 0:
        out.append("v=0")
    if (U * U - 4 * v * v) % p == 0:
        out.append("discq=0")
    # f=x*prod(x+a_i^2) is completely split, so gcd(q,f)>1 is
    # equivalent to one of these five evaluations vanishing.
    if v % p == 0:
        out.append("G0")
    for i, x in enumerate(base):
        A = x * x
        if (A * A - U * A + v * v) % p == 0:
            out.append(f"G{i+1}")
    return tuple(out)


def form_zero_on_space(form, y0, basis, p: int = P) -> bool:
    dot = lambda y: sum(a * b for a, b in zip(form, y)) % p
    return dot(y0) == 0 and all(dot(v) == 0 for v in basis)


def analyze():
    base_count, sols = mod_solutions(P)
    ranks = Counter()
    lift_summary = Counter()
    sig_counts = Counter()
    sig_lifts: dict[tuple[str, ...], Counter] = defaultdict(Counter)
    examples = {}
    for z, ni in sols:
        sig = signature(z)
        sig_counts[sig] += 1
        rank, space = rref_lift_space(z, ni, P)
        ranks[rank] += 1
        if space is None:
            lift_summary["no_p2_lift"] += 1
            sig_lifts[sig]["no_p2_lift"] += 1
            continue
        lift_summary["has_p2_lift"] += 1
        sig_lifts[sig]["has_p2_lift"] += 1
        y0, basis = space
        # Can every base zero, every simple +/- collision, and every
        # q--f intersection acquire a nonzero first coefficient?  There
        # are fewer than 13 such forms, so non-identical affine
        # hyperplanes cannot cover the lift space.
        forms = []
        for i, x in enumerate(z[:4]):
            if x % P == 0:
                f = [0] * 11
                f[i] = 1
                forms.append(tuple(f))
        for i in range(4):
            for j in range(i + 1, 4):
                if (z[i] * z[i] - z[j] * z[j]) % P == 0 and (z[i] or z[j]):
                    f = [0] * 11
                    f[i] = 2 * z[i]
                    f[j] = -2 * z[j]
                    forms.append(tuple(x % P for x in f))
        L, U, v = z[8:11]
        if v % P == 0:
            # q(0)=v^2 has zero differential here.  It cannot be made a
            # p-adic unit at first order in this integral chart.
            forms.append(tuple([0] * 11))
        for i, x in enumerate(z[:4]):
            A = x * x
            if (A * A - U * A + v * v) % P == 0:
                f = [0] * 11
                f[i] = 2 * x * (2 * A - U)
                f[9] = -A
                f[10] = 2 * v
                forms.append(tuple(y % P for y in f))
        if forms and all(not form_zero_on_space(f, y0, basis) for f in forms):
            lift_summary["all_base_and_gcd_boundaries_can_move_first_order"] += 1
            sig_lifts[sig]["all_base_and_gcd_boundaries_can_move_first_order"] += 1
            examples.setdefault(sig, z)
    return base_count, sols, ranks, lift_summary, sig_counts, sig_lifts, examples


def primitive_tuple(row: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    g = math.gcd(math.gcd(abs(row[0]), abs(row[1])), math.gcd(abs(row[2]), abs(row[3])))
    return tuple(x // g for x in row)


def normalize_mod(row: tuple[int, int, int, int], modulus: int) -> tuple[int, int, int, int]:
    vals = tuple(x % modulus for x in row)
    first = next(x for x in vals if math.gcd(x, modulus) == 1)
    ii = pow(first, -1, modulus)
    return tuple(x * ii % modulus for x in vals)


def bank_deep13_counts(bank_path: Path, sols) -> tuple[dict, list, list]:
    live_bases = {tuple(z[:4]) for z, _ni in sols}
    rows = []
    for line in bank_path.read_text().splitlines():
        line = line.strip()
        if not (line.startswith("[") and line.endswith("]")):
            continue
        vals = tuple(int(x.strip()) for x in line[1:-1].split(","))
        if len(vals) == 4:
            rows.append(vals)
    primitive = [primitive_tuple(row) for row in rows]
    live = []
    deep = []
    for idx, row in enumerate(primitive, 1):
        norm13 = normalize_mod(row, P)
        if norm13 not in live_bases:
            continue
        zeros = [i for i, x in enumerate(row) if x % P == 0]
        minus = [
            (i, j)
            for i in range(4)
            for j in range(i + 1, 4)
            if (row[i] + row[j]) % P == 0
        ]
        live.append((idx, row, norm13, tuple(zeros), tuple(minus)))
        ok = False
        if len(zeros) == 1:
            z = zeros[0]
            # The p^2 lift calculation forces both the zero and the live
            # opposite-pair collision to persist through 13^2.
            ok = row[z] % (P * P) == 0 and any(
                i != z and j != z and (row[i] + row[j]) % (P * P) == 0
                for i, j in minus
            )
        elif len(zeros) == 3:
            # abcd is a rational square on the full cover.  Three exact
            # valuation-one coordinates would give odd total valuation,
            # so at least one must lie one level deeper.
            ok = any(row[z] % (P * P) == 0 for z in zeros)
        if ok:
            deep.append((idx, row, norm13, tuple(zeros), tuple(minus)))
    return {
        "rows": len(rows),
        "unique_primitive": len(set(primitive)),
        "mod13_contact_open_live": len(live),
        "deep_mod169_survivors": len(deep),
    }, live, deep


def fixed_base_contact_lifts(row, sols, depth: int = 5):
    """Hensel-lift contact coordinates while holding one rational base fixed."""
    base1 = normalize_mod(row, P)
    seeds = sorted({tuple(z[8:11]) for z, _ni in sols if tuple(z[:4]) == base1})
    history = [(1, seeds)]
    for k in range(2, depth + 1):
        pk = P**k
        basek = normalize_mod(row, pk)
        step = P ** (k - 1)
        lifted = []
        for seed in seeds:
            for delta in product(range(P), repeat=3):
                W = tuple(seed[i] + step * delta[i] for i in range(3))
                if all(e % pk == 0 for e in contact_equations(basek, *W)):
                    lifted.append(W)
        seeds = sorted(set(lifted))
        history.append((k, seeds))
        if not seeds:
            break
    return history


def write_report(path: Path, bank_path: Path):
    base_count, sols, ranks, summary, sigs, bysig, examples = analyze()
    lines = [
        "Full A(2,2,2,8)+3 normalized 13-adic boundary analysis",
        f"normalized_full_cover_bases_mod13 {base_count}",
        f"contact_open_incidence_solutions_mod13 {len(sols)}",
        f"jacobian_ranks {dict(ranks)}",
        f"p2_lift_summary {dict(summary)}",
        "",
        "Boundary/contact signatures:",
    ]
    for sig, count in sigs.most_common():
        lines.append(f"  {sig}: mod13={count} lifts={dict(bysig[sig])} example={examples.get(sig)}")
    lines += [
        "",
        "Scope:",
        "  Full four-radicand A(2,2,2,8) cover; not the K3 sublocus.",
        "  Projective base normalized by its first nonzero coordinate.",
        "  Contact chart assumes L=1/m is a 13-adic unit.",
        "  Linear lift spaces solve all four cover equations, all three",
        "  eliminated contact equations, and the normalization through 13^2.",
    ]
    if bank_path.exists():
        counts, live, deep = bank_deep13_counts(bank_path, sols)
        lines += ["", f"supplied_bank_deep13_counts {counts}", "bank_mod13_live_rows:"]
        for item in live:
            lines.append(f"  {item}")
        lines.append("fixed_base_contact_hensel_lifts:")
        seen_deep = set()
        for _idx, row, _norm, _zeros, _minus in deep:
            if row in seen_deep:
                continue
            seen_deep.add(row)
            lines.append(f"  base={row} history={fixed_base_contact_lifts(row, sols)}")
    path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--output",
        default="results/target_22224_full_family_halving_padic13.txt",
    )
    ap.add_argument(
        "--bank",
        default="data/tor2228_bank.txt",
    )
    args = ap.parse_args()
    out = Path(args.output)
    write_report(out, Path(args.bank))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
