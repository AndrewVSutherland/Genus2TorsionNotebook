#!/usr/bin/env python3
"""Projective p-adic boundary lifts for the Orbit-01 unit-circle chart.

Write t_i=X_i/Y_i, i=2,...,5.  The two affine equations

    E4=E1+5,  E3=E2-10

are homogenized on (P1)^4.  A finite chart has Y_i=1 and variable X_i;
an infinity chart has X_i=1 and variable Y_i, whose residue is zero.
Starting from every F_p point, this script lifts the exact two equations by
their affine tangent system.  It also detects a genuinely open p-adic branch:
if a mod-p^k lift has every smooth-open factor nonzero mod p^k and the base
Jacobian has rank 2, multivariable Hensel gives an open Q_p point.
"""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations, product


def xy(mask, v):
    xs, ys = [], []
    for i, z in enumerate(v):
        if mask >> i & 1:
            xs.append(1)
            ys.append(z)
        else:
            xs.append(z)
            ys.append(1)
    return xs, ys


def equations(mask, v):
    x, y = xy(mask, v)
    x2 = [z * z for z in x]
    y2 = [z * z for z in y]
    py = 1
    px = 1
    for z in y2:
        py *= z
    for z in x2:
        px *= z

    # E4-E1-5=0 after multiplication by product Y_i^2.
    f1 = px - 5 * py
    for i in range(4):
        term = x2[i]
        for j in range(4):
            if j != i:
                term *= y2[j]
        f1 -= term

    # E3-E2+10=0 after the same clearing denominator.
    f2 = 10 * py
    for i in range(4):
        term = y2[i]
        for j in range(4):
            if j != i:
                term *= x2[j]
        f2 += term
    for i, j in combinations(range(4), 2):
        term = x2[i] * x2[j]
        for ell in range(4):
            if ell != i and ell != j:
                term *= y2[ell]
        f2 -= term
    return f1, f2


def jacobian(mask, v, p):
    inv2 = pow(2, -1, p)
    rows = [[], []]
    for j in range(4):
        vp = list(v)
        vm = list(v)
        vp[j] += 1
        vm[j] -= 1
        fp = equations(mask, vp)
        fm = equations(mask, vm)
        for i in range(2):
            rows[i].append((fp[i] - fm[i]) * inv2 % p)
    return rows


def jacobian_exact(mask, v):
    """Exact derivative matrix (each variable occurs only quadratically)."""
    rows = [[], []]
    for j in range(4):
        vp = list(v)
        vm = list(v)
        vp[j] += 1
        vm[j] -= 1
        fp = equations(mask, vp)
        fm = equations(mask, vm)
        for i in range(2):
            assert (fp[i] - fm[i]) % 2 == 0
            rows[i].append((fp[i] - fm[i]) // 2)
    return rows


def vp(z, p):
    if z == 0:
        return 10**9
    z = abs(z)
    e = 0
    while z % p == 0:
        z //= p
        e += 1
    return e


def quantitative_hensel_witness(mask, v, p):
    """Return data for the square-Jacobian quantitative Hensel lemma.

    We fix two of the four variables.  If d=v_p(det J_I),
    min_i v_p(F_i)>2d, Newton-Hensel gives a root whose correction has
    valuation at least min_i v_p(F_i)-d.  Requiring every open factor to
    have smaller valuation proves that the limiting root remains open.
    """
    f = equations(mask, v)
    vf = min(vp(f[0], p), vp(f[1], p))
    j = jacobian_exact(mask, v)
    best = None
    for a, b in combinations(range(4), 2):
        det = j[0][a] * j[1][b] - j[0][b] * j[1][a]
        d = vp(det, p)
        if d >= 10**9 or vf <= 2 * d:
            continue
        correction_v = vf - d
        factor_v = max(vp(z, p) for z in open_factors(mask, v))
        if factor_v < correction_v:
            cand = (d, vf, correction_v, factor_v, (a, b), det)
            if best is None or cand[:2] < best[:2]:
                best = cand
    return best


def rank_mod_p(a, p):
    a = [[z % p for z in row] for row in a]
    row = 0
    for col in range(len(a[0])):
        piv = next((i for i in range(row, len(a)) if a[i][col]), None)
        if piv is None:
            continue
        a[row], a[piv] = a[piv], a[row]
        inv = pow(a[row][col], -1, p)
        a[row] = [z * inv % p for z in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                c = a[i][col]
                a[i] = [(a[i][j] - c * a[row][j]) % p for j in range(len(a[0]))]
        row += 1
        if row == len(a):
            break
    return row


def affine_solutions(matrix, rhs, p):
    m, n = len(matrix), len(matrix[0])
    a = [[z % p for z in matrix[i]] + [rhs[i] % p] for i in range(m)]
    pivots = []
    row = 0
    for col in range(n):
        piv = next((i for i in range(row, m) if a[i][col]), None)
        if piv is None:
            continue
        a[row], a[piv] = a[piv], a[row]
        inv = pow(a[row][col], -1, p)
        a[row] = [z * inv % p for z in a[row]]
        for i in range(m):
            if i != row and a[i][col]:
                c = a[i][col]
                a[i] = [(a[i][j] - c * a[row][j]) % p for j in range(n + 1)]
        pivots.append(col)
        row += 1
        if row == m:
            break
    for i in range(row, m):
        if all(a[i][j] == 0 for j in range(n)) and a[i][n]:
            return []
    free = [j for j in range(n) if j not in pivots]
    out = []
    for vals in product(range(p), repeat=len(free)):
        z = [0] * n
        for j, val in zip(free, vals):
            z[j] = val
        for i, col in enumerate(pivots):
            z[col] = (a[i][n] - sum(a[i][j] * z[j] for j in free)) % p
        out.append(tuple(z))
    return out


def base_points(p):
    out = []
    for mask in range(16):
        inf = [i for i in range(4) if mask >> i & 1]
        fin = [i for i in range(4) if not (mask >> i & 1)]
        for fv in product(range(p), repeat=len(fin)):
            v = [0] * 4
            for i, z in zip(fin, fv):
                v[i] = z
            f = equations(mask, v)
            if f[0] % p == 0 and f[1] % p == 0:
                out.append((mask, tuple(v)))
    return out


def lift(states, q, p):
    out = []
    for mask, v, base_rank in states:
        f = equations(mask, v)
        assert f[0] % q == 0 and f[1] % q == 0
        j = jacobian(mask, v, p)
        rhs = [-(f[0] // q), -(f[1] // q)]
        for d in affine_solutions(j, rhs, p):
            vv = tuple(v[i] + q * d[i] for i in range(4))
            ff = equations(mask, vv)
            assert ff[0] % (q * p) == 0 and ff[1] % (q * p) == 0
            out.append((mask, vv, base_rank))
    return out


def open_factors(mask, v):
    """Numerators whose nonvanishing is exactly the CK smooth open."""
    x, y = xy(mask, v)
    vals = []
    # Collision with r1^2=1, and r_i=0.
    for i in range(4):
        vals.extend((x[i], y[i], y[i] * y[i] - x[i] * x[i]))
    # r_i-r_j and r_i+r_j numerators.
    for i, j in combinations(range(4), 2):
        vals.append(y[i] * y[i] * x[j] * x[j] - x[i] * x[i] * y[j] * y[j])
        vals.append(y[i] * y[i] * y[j] * y[j] - x[i] * x[i] * x[j] * x[j])
    return vals


def signature(mask, v, p, k):
    sig = []
    for z in open_factors(mask, v):
        z %= p**k
        if z == 0:
            sig.append(k)
        else:
            e = 0
            while z % p == 0:
                z //= p
                e += 1
            sig.append(e)
    return tuple(sorted(e for e in sig if e))


def main():
    ap = ArgumentParser()
    ap.add_argument("--p", type=int, default=3, choices=(3, 7))
    ap.add_argument("--max-k", type=int, default=3)
    ap.add_argument("--top", type=int, default=12)
    args = ap.parse_args()
    p = args.p
    base = base_points(p)
    states = [
        (mask, v, rank_mod_p(jacobian(mask, v, p), p))
        for mask, v in base
    ]
    print("ELKIES22210_ORBIT01_BOUNDARY_LIFTS")
    print("p", p, "max_k", args.max_k)
    print("level 1 states", len(states))
    print("infinity_count", dict(sorted(Counter(bin(mask).count("1") for mask, _, _ in states).items())))
    print("jacobian_rank", dict(sorted(Counter(rank for _, _, rank in states).items())))

    for k in range(2, args.max_k + 1):
        states = lift(states, p ** (k - 1), p)
        resolved = [s for s in states if all(z % (p**k) for z in open_factors(s[0], s[1]))]
        hensel_open = [s for s in resolved if s[2] == 2]
        quantitative = []
        for s in resolved:
            proof = quantitative_hensel_witness(s[0], s[1], p)
            if proof is not None:
                quantitative.append((s, proof))
        print("level", k, "modulus", p**k, "states", len(states))
        print("resolved_open", len(resolved), "rank2_resolved_open", len(hensel_open))
        print("quantitative_hensel_open", len(quantitative))
        print("signatures", len(Counter(signature(m, v, p, k) for m, v, _ in states)))
        for sig, count in Counter(signature(m, v, p, k) for m, v, _ in states).most_common(args.top):
            print("SIGNATURE", count, sig)
        if hensel_open:
            mask, v, rank = hensel_open[0]
            print("HENSEL_OPEN_WITNESS", "mask", mask, "variables", v,
                  "rank", rank, "equations", equations(mask, v),
                  "factors_mod", tuple(z % (p**k) for z in open_factors(mask, v)))
            # One witness is enough to prove local nonemptiness; deeper
            # enumeration is no longer useful for an obstruction search.
            break
        if quantitative:
            (mask, v, rank), proof = quantitative[0]
            print("QUANTITATIVE_HENSEL_OPEN_WITNESS", "mask", mask,
                  "variables", v, "base_rank", rank,
                  "proof_d_vf_correction_factor_cols_det", proof,
                  "equations", equations(mask, v),
                  "factors", tuple(open_factors(mask, v)))
            break
    print("DONE")


if __name__ == "__main__":
    main()
