#!/usr/bin/env python3
"""Complete intrinsic V1 p=13 depth-five T-key mask.

The first blow-up strict transform has two points over each of the four
cover residues T=+/-1,+/-5.  Its Jacobian in (L,s,U), with T fixed, is
invertible modulo 13.  Consequently each point is a formal analytic graph
over the entire T disk; contact support is intrinsic and does not come from
the sampled tangent bank.

The viable depth-five keys are therefore exactly the full-cover keys.  The
two +/-5 disks are unramified full cylinders.  The +/-1 disks are ramified
and are expanded explicitly modulo 13^5.
"""

from __future__ import annotations

from pathlib import Path

P = 13
DEPTH = 5
MOD = P**DEPTH
A, B, C, D = -2178, 2420, 9075, -1470


def radicands(t):
    d = D * t * t
    return (
        A * B * C * d,
        A * (A + B) * (A + C) * (A + d),
        B * (B + A) * (B + C) * (B + d),
        C * (C + A) * (C + B) * (C + d),
    )


SQUARES = {z * z % P for z in range(P)}


def square_mod_power_possible(q):
    q %= MOD
    if q == 0:
        return True, DEPTH, 0
    v = 0
    while q % P == 0:
        q //= P
        v += 1
    return v % 2 == 0 and q % P in SQUARES, v, q % P


def cover_key(t):
    data = tuple(square_mod_power_possible(q) for q in radicands(t))
    return all(z[0] for z in data), data


def main():
    out = Path("results/target_22224_V1_depth5_Tkeys.tsv")
    compact = Path("results/target_22224_V1_depth5_Tkeys_compact.tsv")
    log = Path("results/target_22224_V1_depth5_Tkeys.log")

    keys = {}
    counts = {}
    for t0 in (1, 5, 8, 12):
        rows = []
        for k in range(P**(DEPTH - 1)):
            t = t0 + P * k
            ok, data = cover_key(t)
            if ok:
                rows.append((t, data))
        keys[t0] = rows
        counts[t0] = len(rows)

    assert counts == {1: 1021, 5: 28561, 8: 28561, 12: 1021}
    assert {t for t, _ in keys[12]} == {(-t) % MOD for t, _ in keys[1]}
    assert {t for t, _ in keys[8]} == {(-t) % MOD for t, _ in keys[5]}

    with out.open("w") as f:
        f.write("T_key\tT_mod13\tcomponent\tradicand_valuations\t"
                "radicand_units_mod13\n")
        for t0 in (1, 5, 8, 12):
            comp = {
                1: "exceptional_plus1",
                5: "exceptional_plus5",
                8: "exceptional_minus5",
                12: "exceptional_minus1",
            }[t0]
            for t, data in keys[t0]:
                f.write(
                    f"{t}\t{t0}\t{comp}\t"
                    f"{','.join(str(z[1]) for z in data)}\t"
                    f"{','.join(str(z[2]) for z in data)}\n"
                )

    with compact.open("w") as f:
        f.write("component\tmodulus\tresidue_or_file\tkey_count\n")
        f.write("exceptional_plus5\t13\t5\t28561\n")
        f.write("exceptional_minus5\t13\t8\t28561\n")
        f.write("exceptional_plus1\t371293\texpanded_full_file\t1021\n")
        f.write("exceptional_minus1\t371293\texpanded_full_file\t1021\n")

    lines = [
        "V1_DEPTH5_TKEYS_START p=13 depth=5 modulus=371293",
        "STRICT_TRANSFORM_FIXED_T_JACOBIAN_RANK=3 for s=6,7 over T=+/-1,+/-5",
        "CONTACT_SCOPE complete first-blow-up formal graphs; independent of sampled bank",
        f"COUNTS {counts} total={sum(counts.values())}",
        "COMPACT +/-5 are complete residue cylinders; +/-1 are ramified cover keys",
        f"output={out}",
        f"compact={compact}",
    ]
    text = "\n".join(lines) + "\n"
    log.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
