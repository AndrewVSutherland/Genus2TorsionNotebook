#!/usr/bin/env python3
"""Focused CRT-lattice search in one compatible orbit-12 p-adic disk.

The selected local cover branches are represented by primitive CK tuples
modulo 11^3, 19^2, and 23^2.  Their pullbacks to the complete rational
(t,m)-chart are combined by CRT.  At p=19 the chart has a base point:
all five polynomial chart coordinates have a common factor 19.  Therefore
the primitive p^2 branch requires t,m modulo 19^3.  This script performs
that one-digit lift before forming the lattice.

For a residue X modulo M, reduced fractions n/d in the selected disk obey

    n - X*d = 0 (mod M),

so they are lattice points in <(M,0),(X,1)>.  The bounded search enumerates
these fractions separately for t and m, crosses only the resulting short
lists, and tests all four exact Stoll--Zarhin radicands.
"""

from argparse import ArgumentParser
from functools import reduce
from math import gcd, isqrt, prod


LOCAL_MODULI = (11**3, 19**2, 23**2)
LOCAL_SEEDS = (
    (1, 242, 959, 9, 120),
    (1, 248, 98, 3, 11),
    (1, 392, 118, 8, 10),
)
LOCAL_G0_ROOTS = (22, 97, 49)


def crt(residues, moduli):
    modulus = prod(moduli)
    return sum(
        r * (modulus // q) * pow(modulus // q, -1, q)
        for r, q in zip(residues, moduli)
    ) % modulus


def chart_coordinates(t, m, modulus=None):
    r = (
        1 + t * (t + 2) * m,
        t * m * (m - t - 2),
        -1 + m + t * (t + 1) * m * m,
        1 + t - m - t * m * m,
        -(1 + t) * (1 + t * m * m),
    )
    return tuple(x % modulus for x in r) if modulus else r


def inverse_chart(r, modulus):
    t = (r[2] + r[3]) * pow((r[0] + r[1]) % modulus, -1, modulus)
    m = (r[2] + r[0] + r[1]) * pow(r[0] % modulus, -1, modulus)
    return t % modulus, m % modulus


def ck_ok(r, modulus):
    return sum(r) % modulus == 0 and sum(x**3 for x in r) % modulus == 0


def radicands(r):
    a = [x * x for x in r]
    return (
        -(a[0] - a[2]) * (a[0] - a[3]) * (a[0] - a[4]),
        (a[2] - a[1]) * (a[0] - a[3]) * (a[0] - a[4]),
        (a[3] - a[1]) * (a[0] - a[2]) * (a[0] - a[4]),
        (a[4] - a[1]) * (a[0] - a[2]) * (a[0] - a[3]),
    )


def lift_p19_seed():
    """Lift the printed p^2 seed to p^3, fixing its first three entries."""
    p = 19
    q = p * p
    qq = q * p
    r = list(LOCAL_SEEDS[1])
    assert ck_ok(r, q)
    rhs1 = -(sum(r) // q) % p
    rhs3 = -(sum(x**3 for x in r) // q) % p
    # d4+d5=rhs1 and 3*r4^2*d4+3*r5^2*d5=rhs3 mod p.
    coefficient = 3 * (r[3] * r[3] - r[4] * r[4]) % p
    right = (rhs3 - 3 * r[4] * r[4] * rhs1) % p
    d4 = right * pow(coefficient, -1, p) % p
    d5 = (rhs1 - d4) % p
    r[3] += q * d4
    r[4] += q * d5
    assert ck_ok(r, qq)
    assert tuple(x % q for x in r) == LOCAL_SEEDS[1]
    return tuple(r), (d4, d5)


def selected_chart_disk():
    lifted19, deltas19 = lift_p19_seed()
    chart_moduli = (11**3, 19**3, 23**2)
    chart_seeds = (LOCAL_SEEDS[0], lifted19, LOCAL_SEEDS[2])
    tm = tuple(inverse_chart(r, q) for r, q in zip(chart_seeds, chart_moduli))

    # At 11 and 23 the chart scale is a unit.  At 19 all polynomial chart
    # coordinates have valuation one; divide by 19 to recover primitive p^2
    # data.  These checks guard the otherwise easy-to-miss lost digit.
    scales = []
    for index, (r, q, (t, m)) in enumerate(zip(chart_seeds, chart_moduli, tm)):
        rr = chart_coordinates(t, m, q)
        if index != 1:
            lam = rr[0] * pow(r[0], -1, q) % q
            assert gcd(lam, (11, 19, 23)[index]) == 1
            assert rr == tuple(lam * x % q for x in r)
            scales.append((lam, 0))
        else:
            assert all(x % 19 == 0 for x in rr)
            primitive = tuple((x // 19) % (19**2) for x in rr)
            lam = primitive[0]
            assert primitive == tuple(lam * x % (19**2) for x in r)
            scales.append((lam, 1))

    modulus = prod(chart_moduli)
    t_residue = crt(tuple(z[0] for z in tm), chart_moduli)
    m_residue = crt(tuple(z[1] for z in tm), chart_moduli)
    q_local = tuple(r[1] * pow(r[0], -1, q) % q
                    for r, q in zip(chart_seeds, chart_moduli))
    q_residue = crt(q_local, chart_moduli)

    # Pull one sign of the G0 square root to the chart.  At 19, lift the
    # primitive root from p^2 to p^3, then divide the chart coordinates by
    # their common p; a degree-six radicand root scales cubically.
    g019 = radicands(lifted19)[0] % (19**3)
    lifted_roots19 = [
        (LOCAL_G0_ROOTS[1] + 19**2 * e) % (19**3) for e in range(19)
        if (LOCAL_G0_ROOTS[1] + 19**2 * e) ** 2 % (19**3) == g019
    ]
    assert len(lifted_roots19) == 1
    chart_y = (
        LOCAL_G0_ROOTS[0] * scales[0][0] ** 3 % (11**3),
        (lifted_roots19[0] % (19**2)) * scales[1][0] ** 3 % (19**2),
        LOCAL_G0_ROOTS[2] * scales[2][0] ** 3 % (23**2),
    )
    # With Z=Y/19^3, the first-cover equation is
    # (19^3 Z)^2=G0(R(t,m)).  At 19, chart_y already means Y/19^3.
    root_moduli = LOCAL_MODULI
    z_roots = (
        chart_y[0] * pow(19**3, -1, 11**3) % (11**3),
        chart_y[1],
        chart_y[2] * pow(19**3, -1, 23**2) % (23**2),
    )
    signed_z_crt = []
    for s11 in (1, -1):
        for s19 in (1, -1):
            for s23 in (1, -1):
                signs = (s11, s19, s23)
                residues = tuple(s * z % q for s, z, q in zip(signs, z_roots, root_moduli))
                signed_z_crt.append((signs, crt(residues, root_moduli)))
    return {
        "lifted19": lifted19,
        "deltas19": deltas19,
        "moduli": chart_moduli,
        "tm": tm,
        "scales": tuple(scales),
        "modulus": modulus,
        "t": t_residue,
        "m": m_residue,
        "q_local": q_local,
        "q": q_residue,
        "p19_g0_lifted_root": lifted_roots19[0],
        "chart_y": chart_y,
        "z_roots": z_roots,
        "z_modulus": prod(root_moduli),
        "signed_z_crt": tuple(signed_z_crt),
    }


def nearest_integer(num, den):
    if num >= 0:
        return (2 * num + den) // (2 * den)
    return -nearest_integer(-num, den)


def gauss_reduce(v1, v2):
    dot = lambda a, b: a[0] * b[0] + a[1] * b[1]
    while True:
        if dot(v2, v2) < dot(v1, v1):
            v1, v2 = v2, v1
        mu = nearest_integer(dot(v1, v2), dot(v1, v1))
        if mu == 0:
            return v1, v2
        v2 = (v2[0] - mu * v1[0], v2[1] - mu * v1[1])


def fractions_in_disk(residue, modulus, height):
    if 2 * height >= modulus:
        raise ValueError("implementation requires 2*height < CRT modulus")
    out = []
    for denominator in range(1, height + 1):
        if gcd(denominator, modulus) != 1:
            continue
        numerator = residue * denominator % modulus
        if numerator > modulus // 2:
            numerator -= modulus
        if abs(numerator) > height:
            continue
        if gcd(abs(numerator), denominator) == 1:
            out.append((numerator, denominator))
    return out


def cleared_chart(t_fraction, m_fraction):
    a, b = t_fraction
    c, d = m_fraction
    bb, dd = b * b, d * d
    r = [
        bb * dd + a * (a + 2 * b) * c * d,
        a * c * (b * c - d * (a + 2 * b)),
        -bb * dd + bb * c * d + a * (a + b) * c * c,
        b * (a + b) * dd - bb * c * d - a * b * c * c,
        -(a + b) * (b * dd + a * c * c),
    ]
    common = reduce(gcd, (abs(x) for x in r))
    return tuple(x // common for x in r)


SMALL_MODULI = (64, 63, 65, 11, 13, 17, 19, 23, 29, 31)
SMALL_SQUARES = tuple({x * x % q for x in range(q)} for q in SMALL_MODULI)


def possible_square(n):
    return n > 0 and all(
        n % q in squares for q, squares in zip(SMALL_MODULI, SMALL_SQUARES)
    )


def exact_positive_square(n):
    return possible_square(n) and isqrt(n) ** 2 == n


def run_search(height):
    disk = selected_chart_disk()
    modulus = disk["modulus"]
    t_values = fractions_in_disk(disk["t"], modulus, height)
    m_values = fractions_in_disk(disk["m"], modulus, height)
    q_values = fractions_in_disk(disk["q"], modulus, height)

    counts = {
        "pairs": 0,
        "nonzero": 0,
        "smooth": 0,
        "g0_positive": 0,
        "all_positive": 0,
        "all_small_mod_squares": 0,
        "g0_exact_square": 0,
    }
    masks = {}
    hits = []
    for tv in t_values:
        for mv in m_values:
            counts["pairs"] += 1
            r = cleared_chart(tv, mv)
            if any(x == 0 for x in r):
                continue
            counts["nonzero"] += 1
            assert sum(r) == 0 and sum(x**3 for x in r) == 0
            a = [x * x for x in r]
            if len(set(a)) != 5:
                continue
            counts["smooth"] += 1
            gs = radicands(r)
            if gs[0] > 0:
                counts["g0_positive"] += 1
            if exact_positive_square(gs[0]):
                counts["g0_exact_square"] += 1
            if any(g <= 0 for g in gs):
                continue
            counts["all_positive"] += 1
            if all(possible_square(g) for g in gs):
                counts["all_small_mod_squares"] += 1
            mask = sum(
                1 << i for i, g in enumerate(gs) if exact_positive_square(g)
            )
            masks[mask] = masks.get(mask, 0) + 1
            if mask == 15:
                hits.append((tv, mv, r))
    return disk, t_values, m_values, q_values, counts, masks, hits


def main():
    ap = ArgumentParser()
    ap.add_argument("--height", type=int, default=2_000_000)
    args = ap.parse_args()
    disk, t_values, m_values, q_values, counts, masks, hits = run_search(args.height)

    print("ELKIES22210_ORBIT12_CHART_CRT_LATTICE")
    print("local_moduli", LOCAL_MODULI)
    print("local_seeds", LOCAL_SEEDS)
    print("p19_lift_deltas", disk["deltas19"])
    print("p19_lifted_seed_mod_6859", disk["lifted19"])
    print("chart_moduli", disk["moduli"])
    print("chart_tm_residues", disk["tm"])
    print("chart_scales_and_valuations", disk["scales"])
    print("crt_modulus", disk["modulus"])
    print("crt_t_residue", disk["t"])
    print("crt_m_residue", disk["m"])
    print("marked_ratio_q_local", disk["q_local"])
    print("crt_q_residue", disk["q"])
    print("p19_primitive_g0_root_mod_6859", disk["p19_g0_lifted_root"])
    print("chart_first_cover_y_data", disk["chart_y"])
    print("z_equals_y_over_19_cubed_residues", disk["z_roots"])
    print("z_crt_modulus", disk["z_modulus"])
    print("signed_z_crt", disk["signed_z_crt"])
    for label in ("t", "m", "q"):
        basis = gauss_reduce((disk["modulus"], 0), (disk[label], 1))
        print(label + "_lattice_basis", basis)
    print("height", args.height)
    print("t_fractions", len(t_values))
    print("m_fractions", len(m_values))
    print("q_fractions", len(q_values))
    print("short_q_fractions", sorted(
        q_values, key=lambda z: (max(abs(z[0]), z[1]), z[1], z[0])
    )[:10])
    print("parameter_pairs", counts["pairs"])
    for key in (
        "nonzero", "smooth", "g0_positive", "all_positive",
        "all_small_mod_squares", "g0_exact_square",
    ):
        print(key, counts[key])
    print("exact_masks", " ".join(f"{k}:{masks[k]}" for k in sorted(masks)))
    print("hits", len(hits))
    for hit in hits:
        print("HIT", hit)
    print("DONE")


if __name__ == "__main__":
    main()
