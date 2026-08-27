#!/usr/bin/env python3
"""Normalized p-adic boundary classifier for the exact orbit-12 cover.

The marked class is {r1^2,r2^2}.  In every projective residue class we
normalize the first unit among r1,...,r5 to 1; hence the five charts cover
the cases in which either or both marked coordinates are nonunits.

The four exact Stoll--Zarhin radicands are

  G0 = -(a1-a3)(a1-a4)(a1-a5),
  G3 =  (a3-a2)(a1-a4)(a1-a5),
  G4 =  (a4-a2)(a1-a3)(a1-a5),
  G5 =  (a5-a2)(a1-a3)(a1-a4).

We retain normalized solutions of sum(r_i)=sum(r_i^3)=0 modulo p^k
for which all four G's are squares modulo p^k.  Lifting the two CK
equations is done by their exact affine Hensel system, so only p^2
(normally) tangent digits are enumerated at each level.

This square-residue classifier deliberately records G=0 mod p^k as
"deep": such a residue can overcount genuine Q_p square roots until its
valuation resolves.  The output separates these from resolved valuation-2
radicands.

Run:
  python3 code/elkies22210_orbit12_boundary_lifts.py --p 11 --max-k 3
"""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from itertools import product


def radicands(r, modulus=None):
    a = [x * x for x in r]
    g = [
        -(a[0] - a[2]) * (a[0] - a[3]) * (a[0] - a[4]),
        (a[2] - a[1]) * (a[0] - a[3]) * (a[0] - a[4]),
        (a[3] - a[1]) * (a[0] - a[2]) * (a[0] - a[4]),
        (a[4] - a[1]) * (a[0] - a[2]) * (a[0] - a[3]),
    ]
    return tuple(x % modulus for x in g) if modulus else tuple(g)


def square_residues(q):
    return {x * x % q for x in range(q)}


def cover_ok(r, q, squares):
    return all(g in squares for g in radicands(r, q))


def affine_solutions(matrix, rhs, p):
    """Return all solutions to matrix*x=rhs over F_p."""
    m = len(matrix)
    n = len(matrix[0])
    aug = [[x % p for x in matrix[i]] + [rhs[i] % p] for i in range(m)]
    pivots = []
    row = 0
    for col in range(n):
        pivot = next((i for i in range(row, m) if aug[i][col] % p), None)
        if pivot is None:
            continue
        aug[row], aug[pivot] = aug[pivot], aug[row]
        inv = pow(aug[row][col], -1, p)
        aug[row] = [(v * inv) % p for v in aug[row]]
        for i in range(m):
            if i == row or aug[i][col] == 0:
                continue
            c = aug[i][col]
            aug[i] = [(aug[i][j] - c * aug[row][j]) % p for j in range(n + 1)]
        pivots.append(col)
        row += 1
        if row == m:
            break
    for i in range(row, m):
        if all(aug[i][j] == 0 for j in range(n)) and aug[i][n] != 0:
            return [], len(pivots)
    free = [j for j in range(n) if j not in pivots]
    out = []
    for vals in product(range(p), repeat=len(free)):
        x = [0] * n
        for j, v in zip(free, vals):
            x[j] = v
        for i, col in enumerate(pivots):
            x[col] = (aug[i][n] - sum(aug[i][j] * x[j] for j in free)) % p
        out.append(tuple(x))
    return out, len(pivots)


def base_points(p):
    sq = square_residues(p)
    states = []
    for chart in range(5):
        # Earlier entries are nonunits, hence zero modulo p; the first
        # unit is normalized to 1.  Later entries are arbitrary.
        for tail in product(range(p), repeat=4 - chart):
            r = [0] * chart + [1] + list(tail)
            if sum(r) % p or sum(x ** 3 for x in r) % p:
                continue
            if cover_ok(r, p, sq):
                states.append((tuple(r), chart))
    return states


def lift_deltas(r, chart, q, p):
    free = [i for i in range(5) if i != chart]
    s1 = sum(r)
    s3 = sum(x ** 3 for x in r)
    assert s1 % q == 0 and s3 % q == 0
    matrix = [
        [1 for _ in free],
        [(3 * (r[i] % p) ** 2) % p for i in free],
    ]
    rhs = [-(s1 // q), -(s3 // q)]
    return affine_solutions(matrix, rhs, p)[0]


def ck_rank_mod_p(r, chart, p):
    free = [i for i in range(5) if i != chart]
    matrix = [
        [1 for _ in free],
        [(3 * (r[i] % p) ** 2) % p for i in free],
    ]
    _, rank = affine_solutions(matrix, [0, 0], p)
    return rank


def lift_level(states, p, k_from, base_ids):
    q = p ** k_from
    qnext = q * p
    sq = square_residues(qnext)
    lifted = []
    child_ids = []
    per_parent = Counter()
    for idx, ((r, chart), base_id) in enumerate(zip(states, base_ids)):
        count = 0
        free = [i for i in range(5) if i != chart]
        for delta in lift_deltas(r, chart, q, p):
            rr = list(r)
            for i, d in zip(free, delta):
                rr[i] += q * d
            rr = tuple(rr)
            assert rr[chart] == 1
            assert sum(rr) % qnext == 0
            assert sum(x ** 3 for x in rr) % qnext == 0
            if cover_ok(rr, qnext, sq):
                lifted.append((rr, chart))
                child_ids.append(base_id)
                count += 1
        per_parent[count] += 1
    return lifted, child_ids, per_parent


def base_stratum(r, p):
    zero = tuple(i + 1 for i, x in enumerate(r) if x % p == 0)
    a = [(x * x) % p for x in r]
    coll = tuple(
        f"{i+1}{j+1}" for i in range(5) for j in range(i + 1, 5) if a[i] == a[j]
    )
    rz = radicands(r, p)
    gzero = tuple(i for i, x in enumerate(rz) if x == 0)
    marked = ("U" if r[0] % p else "N") + ("U" if r[1] % p else "N")
    return f"chart={next(i+1 for i,x in enumerate(r) if x%p)} marked={marked} Z={zero} C={coll} G0={gzero}"


def valuation_signature(r, p, k):
    q = p ** k
    sig = []
    for g in radicands(r, q):
        if g == 0:
            sig.append(f"D{k}")
            continue
        v = 0
        while g % p == 0:
            v += 1
            g //= p
        sig.append(f"V{v}")
    return tuple(sig)


def curve_open_mod_q(r, q):
    """Certify that 0,r_1^2,...,r_5^2 are distinct modulo q."""
    if any(x % q == 0 for x in r):
        return False
    a = [(x * x) % q for x in r]
    return all((a[i] - a[j]) % q != 0
               for i in range(5) for j in range(i + 1, 5))


def main():
    ap = ArgumentParser()
    ap.add_argument("--p", type=int, default=11)
    ap.add_argument("--max-k", type=int, default=3, choices=(1, 2, 3))
    ap.add_argument("--top", type=int, default=30)
    args = ap.parse_args()
    p = args.p

    states = base_points(p)
    base_ids = list(range(len(states)))
    strata = [base_stratum(r, p) for r, _ in states]
    marked_for_base = [s.split("marked=")[1][:2] for s in strata]
    print("ELKIES22210_ORBIT12_BOUNDARY_LIFTS")
    print("p", p, "max_k", args.max_k)
    print("level 1 modulus", p, "states", len(states))
    print("marked_status", dict(Counter(marked_for_base)))
    tangent_counts = Counter(
        len(lift_deltas(r, chart, p, p)) for r, chart in states
    )
    print("ck_tangent_lift_counts", dict(sorted(tangent_counts.items())))
    print("base_strata", len(set(strata)))
    for key, count in Counter(strata).most_common(args.top):
        print("BASE_STRATUM", count, key)

    for k in range(1, args.max_k):
        parents = states
        parent_ids = base_ids
        states, base_ids, per_parent = lift_level(parents, p, k, parent_ids)
        print("level", k + 1, "modulus", p ** (k + 1), "states", len(states))
        print("per_parent_lift_counts", dict(sorted(per_parent.items())))
        print("marked_status", dict(Counter(marked_for_base[i] for i in base_ids)))
        by_base_stratum = Counter(strata[i] for i in base_ids)
        print("surviving_base_strata", len(by_base_stratum))
        for key, count in by_base_stratum.most_common(args.top):
            print("LIFTED_STRATUM", k + 1, count, key)
        vsig = Counter(valuation_signature(r, p, k + 1) for r, _ in states)
        print("valuation_signatures", len(vsig))
        for sig, count in vsig.most_common(args.top):
            print("VALUATION", k + 1, count, sig)
        deep = sum(count for sig, count in vsig.items() if any(x.startswith("D") for x in sig))
        resolved = len(states) - deep
        print("resolved", resolved, "deep", deep)
        resolution_marked = Counter()
        resolved_strata = Counter()
        deep_strata = Counter()
        resolved_sample = None
        qnow = p ** (k + 1)
        for (r, chart), base_id in zip(states, base_ids):
            sig = valuation_signature(r, p, k + 1)
            status = "deep" if any(x.startswith("D") for x in sig) else "resolved"
            resolution_marked[(marked_for_base[base_id], status)] += 1
            if status == "resolved":
                resolved_strata[strata[base_id]] += 1
                if resolved_sample is None and curve_open_mod_q(r, qnow):
                    resolved_sample = (r, chart, base_id, sig)
            else:
                deep_strata[strata[base_id]] += 1
        print("resolution_by_marked", dict(sorted(resolution_marked.items())))
        print("resolved_base_strata", len(resolved_strata),
              "deep_base_strata", len(deep_strata))
        for key, count in resolved_strata.most_common(args.top):
            print("RESOLVED_STRATUM", k + 1, count, key)
        if resolved_sample is not None:
            rr0, chart0, bid0, sig0 = resolved_sample
            gs0 = radicands(rr0, qnow)
            root_map = {}
            for z in range(qnow):
                root_map.setdefault(z * z % qnow, z)
            roots0 = tuple(root_map[g] for g in gs0)
            rank0 = ck_rank_mod_p(rr0, chart0, p)
            even_finite0 = all(s.startswith("V") and int(s[1:]) % 2 == 0
                               for s in sig0)
            open0 = curve_open_mod_q(rr0, qnow)
            certified0 = rank0 == 2 and even_finite0 and open0
            aa0 = tuple((x * x) % qnow for x in rr0)
            diffs0 = tuple((aa0[i] - aa0[j]) % qnow
                           for i in range(5) for j in range(i + 1, 5))
            print("RESOLVED_SAMPLE", "r", rr0, "chart", chart0 + 1,
                  "base_stratum", strata[bid0], "radicands", gs0,
                  "roots", roots0, "valuation", sig0,
                  "ck_rank_mod_p", rank0)
            print("OPENNESS_SAMPLE", "r_squares", aa0,
                  "pair_differences", diffs0)
            print("GENUINE_QP_BRANCH", str(certified0).lower(),
                  "ck_hensel_rank", rank0,
                  "all_radicands_even_finite", even_finite0,
                  "curve_open_modulus", open0)
        if states:
            print("sample", states[0], "base_stratum", strata[base_ids[0]],
                  "valuation", valuation_signature(states[0][0], p, k + 1))
    print("DONE")


if __name__ == "__main__":
    main()
