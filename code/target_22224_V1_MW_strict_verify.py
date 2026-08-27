#!/usr/bin/env python3
"""Certify the strict-transform contact branches over V1 MW p=13 hits.

For fixed T, the exceptional equations have an invertible 3x3 Jacobian in
(L,s,U).  Ordinary multivariate Hensel therefore gives a unique infinite
Q_13 solution above each mod-13 point.  A finite nonzero valuation of L and
of the q/f resultant proves that this solution leaves the universal L=0 and
gcd boundaries, respectively.
"""

from __future__ import annotations

import csv
from pathlib import Path

import target_22224_V_components_analysis as core

P = 13
DEPTH = 12
MOD = P**DEPTH
FAMILY = "V1"


def valuation_mod(q):
    q %= MOD
    if q == 0:
        return DEPTH
    v = 0
    while q % P == 0:
        q //= P
        v += 1
    return v


def fixed_t_lift(T, s0):
    z = (0, s0, 2)  # (L,s,U)
    x0 = (T % P, *z)
    Jfull = core.jacobian_mod(
        lambda fam, x: core.strict_contact_equations(fam, x),
        FAMILY,
        x0,
        P,
    )
    J = [row[1:] for row in Jfull]
    rank, _ = core.rank_and_space(J, [0, 0, 0], P)
    assert rank == 3
    for level in range(1, DEPTH):
        mod = P**level
        next_mod = mod * P
        vals = core.strict_contact_equations(FAMILY, (T % next_mod, *z))
        assert all(q % mod == 0 for q in vals)
        rhs = [-(q // mod) % P for q in vals]
        rank, space = core.rank_and_space(J, rhs, P)
        assert rank == 3 and space is not None and not space[1]
        delta = space[0]
        z = tuple((z[i] + mod * delta[i]) % next_mod for i in range(3))
    assert all(q % MOD == 0 for q in
               core.strict_contact_equations(FAMILY, (T % MOD, *z)))
    return z


def invariants(T, z):
    L, s, U = z
    inv2 = pow(2, -1, MOD)
    vv = ((U + L * s) * inv2) % MOD
    M = L * L % MOD
    original = core.contact_equations_m(FAMILY, (T, M, U, vv))
    assert all(q % MOD == 0 for q in original)
    A, B, C, D = core.FAMILIES[FAMILY]
    base = (A, B, C, D * T * T)
    resultant = vv * vv
    for a in base:
        resultant *= a**4 - U * a * a + vv * vv
    diffs = [base[i] * base[i] - base[j] * base[j]
             for i in range(4) for j in range(i)]
    return valuation_mod(L), valuation_mod(resultant), max(
        valuation_mod(q) for q in diffs
    )


def main():
    inp = Path("results/target_22224_V1_MW_p13_components.tsv")
    out = Path("results/target_22224_V1_MW_strict_verify.tsv")
    log = Path("results/target_22224_V1_MW_strict_verify.log")
    rows = [r for r in csv.DictReader(inp.open(), delimiter="\t")
            if r["component"] != "none"]
    output = []
    lines = [
        f"V1_MW_STRICT_VERIFY_START rows={len(rows)} depth={DEPTH}",
        "FIXED_T_JACOBIAN_RANK=3 (formal Hensel certificate)",
    ]
    for row in rows:
        T = int(row["T_mod13pow12"])
        for s0 in (6, 7):
            z = fixed_t_lift(T, s0)
            vL, vres, vcollision = invariants(T, z)
            genuine = vL < DEPTH and vres < DEPTH and vcollision < DEPTH
            assert genuine
            output.append((row["m"], row["n"], row["torsion_coset"],
                           row["component"], s0, T, *z, vL, vres,
                           vcollision, int(genuine)))
            lines.append(
                f"ROW m={row['m']} n={row['n']} component={row['component']} "
                f"s0={s0} v13(L)={vL} v13(resultant)={vres} "
                f"max_collision_v={vcollision} genuine={genuine}"
            )
    with out.open("w") as f:
        f.write("m\tn\ttorsion_coset\tcomponent\ts0\tT_mod13pow12\t"
                "L\ts\tU\tv13_L\tv13_resultant\tmax_collision_v\tgenuine\n")
        for row in output:
            f.write("\t".join(map(str, row)) + "\n")
    lines.append(f"V1_MW_STRICT_VERIFY_DONE genuine_rows={len(output)} output={out}")
    text = "\n".join(lines) + "\n"
    log.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
