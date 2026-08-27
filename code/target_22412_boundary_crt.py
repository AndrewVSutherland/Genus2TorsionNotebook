#!/usr/bin/env python3
"""Boundary and CRT attack for the target J(Q)_tors = [2,2,4,12].

There are two independent computations in this file.

1.  For the standard second-J[2]-halving cover

        y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2),

    with the fixed partition [a,b | c,d], enumerate every projective
    mod-p base point on the cover and every normalized lift to p^2.  The
    four cross differences must have a common Q_p squareclass.  Dividing
    each known difference by its forced power of p gives an exact p^2
    classification into killed, resolved, and deep residue classes.

2.  Search the genuine q=(x+t)^2 slice of A(2,2,2,12).  In projective
    branch coordinates this slice is

        [A, B, C, D] = [A, B, 1/(AB), sigma/rho],

        rho^2   = A^2+B^2+C^2-3,
        sigma^2 = A^-2+B^-2+C^-2-3.

    Exact finite masks at many primes impose both square equations and at
    least one of the three second-half partitions.  A bounded rational
    search then tests the surviving pairs exactly.  This search is on a
    two-dimensional subfamily of the global A(2,2,2,12) cover; it is not a
    rescan of tor2244.txt and is not exhaustive for the global family.

No third-party Python modules are required.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import gcd, isqrt
from pathlib import Path
from time import perf_counter


PARTITIONS = (
    ("12|34", (0, 1, 2, 3)),
    ("13|24", (0, 2, 1, 3)),
    ("14|23", (0, 3, 1, 2)),
)


def legendre_or_zero(z: int, p: int) -> int:
    z %= p
    if z == 0:
        return 0
    return 1 if pow(z, (p - 1) // 2, p) == 1 else -1


def is_square_mod_p_allow_zero(z: int, p: int) -> bool:
    return legendre_or_zero(z, p) >= 0


def branch_squares(vals: tuple[int, int, int, int] | list[int]) -> list[int]:
    return [z * z for z in vals]


def cross_differences(
    vals: tuple[int, int, int, int] | list[int],
) -> tuple[int, int, int, int]:
    """Cross differences for the fixed partition [1,2 | 3,4]."""
    aa, bb, cc, dd = branch_squares(vals)
    return cc - aa, dd - aa, cc - bb, dd - bb


def second_half_radicands_from_squares(
    squares: list[int | Fraction], order: tuple[int, int, int, int]
) -> list[int | Fraction]:
    aa, bb, cc, dd = [squares[i] for i in order]
    return [
        (cc - aa) * (cc - bb),
        (cc - aa) * (dd - aa),
        (cc - bb) * (dd - bb),
        (dd - aa) * (dd - bb),
    ]


def fixed_cover_mod_p(vals: tuple[int, int, int, int], p: int) -> bool:
    diffs = cross_differences(vals)
    radicands = (
        diffs[0] * diffs[2],
        diffs[0] * diffs[1],
        diffs[2] * diffs[3],
        diffs[1] * diffs[3],
    )
    return all(is_square_mod_p_allow_zero(r, p) for r in radicands)


def projective_representatives(p: int):
    """Unique reps whose first nonzero coordinate is 1."""
    for vals in product(range(p), repeat=4):
        pivot = next((i for i, z in enumerate(vals) if z), None)
        if pivot is not None and vals[pivot] == 1:
            yield vals, pivot


def signed_boundary_labels(vals: tuple[int, int, int, int], p: int) -> tuple[str, ...]:
    labels: list[str] = []
    for i, z in enumerate(vals):
        if z % p == 0:
            labels.append(f"Z{i + 1}")
    for i in range(4):
        for j in range(i + 1, 4):
            if (vals[i] - vals[j]) % p == 0:
                labels.append(f"E{i + 1}{j + 1}+")
            if (vals[i] + vals[j]) % p == 0:
                labels.append(f"E{i + 1}{j + 1}-")
    return tuple(labels)


def unsigned_boundary_labels_from_squares(squares: list[int], p: int) -> tuple[str, ...]:
    labels: list[str] = []
    for i, z in enumerate(squares):
        if z % p == 0:
            labels.append(f"Z{i + 1}")
    for i in range(4):
        for j in range(i + 1, 4):
            if (squares[i] - squares[j]) % p == 0:
                labels.append(f"E{i + 1}{j + 1}")
    return tuple(labels)


def squareclass_mod_pk(z: int, p: int, modulus: int) -> tuple[int, int] | None:
    """Known Q_p squareclass from a residue modulo p^k; None means deep."""
    z %= modulus
    if z == 0:
        return None
    valuation = 0
    while z % p == 0:
        valuation += 1
        z //= p
    return valuation % 2, legendre_or_zero(z, p)


def classify_cross_differences_p2(
    vals: tuple[int, int, int, int] | list[int], p: int
) -> tuple[str, tuple[tuple[int, int] | None, ...]]:
    modulus = p * p
    squareclasses = tuple(
        squareclass_mod_pk(delta, p, modulus) for delta in cross_differences(vals)
    )
    known = {sc for sc in squareclasses if sc is not None}
    if len(known) > 1:
        return "killed", squareclasses
    if None in squareclasses:
        return "deep", squareclasses
    return "resolved", squareclasses


def format_counter(counter: Counter, keys: tuple[str, ...] | None = None) -> str:
    if keys is None:
        keys = tuple(sorted(counter))
    return " ".join(f"{key}={counter[key]}" for key in keys)


def run_p2_boundary_classifier(p: int, emit) -> dict:
    t0 = perf_counter()
    bases: list[tuple[tuple[int, int, int, int], int, tuple[str, ...]]] = []
    projective = 0
    open_bases = 0
    component_bases: Counter = Counter()
    signature_bases: Counter = Counter()

    for vals, pivot in projective_representatives(p):
        projective += 1
        labels = signed_boundary_labels(vals, p)
        if not labels:
            open_bases += 1
        if not fixed_cover_mod_p(vals, p):
            continue
        bases.append((vals, pivot, labels))
        signature_bases["+".join(labels) if labels else "open"] += 1
        for label in labels:
            component_bases[label] += 1

    lift_status: Counter = Counter()
    base_outcomes: Counter = Counter()
    component_lifts: dict[str, Counter] = defaultdict(Counter)
    component_outcomes: dict[str, Counter] = defaultdict(Counter)
    signature_outcomes: dict[str, Counter] = defaultdict(Counter)
    squareclass_patterns: Counter = Counter()

    for vals, pivot, labels in bases:
        free = [i for i in range(4) if i != pivot]
        local: Counter = Counter()
        for digits in product(range(p), repeat=3):
            lift = list(vals)
            for i, digit in zip(free, digits):
                lift[i] += p * digit
            status, pattern = classify_cross_differences_p2(lift, p)
            lift_status[status] += 1
            local[status] += 1
            squareclass_patterns[pattern] += 1
            for label in labels:
                component_lifts[label][status] += 1

        if local["resolved"] and local["deep"]:
            outcome = "resolved_and_deep"
        elif local["resolved"]:
            outcome = "resolved_only"
        elif local["deep"]:
            outcome = "deep_only"
        else:
            outcome = "all_killed"
        base_outcomes[outcome] += 1
        signature = "+".join(labels) if labels else "open"
        signature_outcomes[signature][outcome] += 1
        for label in labels:
            component_outcomes[label][outcome] += 1

    emit("P2_BOUNDARY_CLASSIFIER")
    emit(f"prime {p}")
    emit(f"projective_bases {projective}")
    emit(f"ambient_open_bases {open_bases}")
    emit(f"modp_cover_bases {len(bases)}")
    emit(f"modp_cover_open_bases {sum(1 for _, _, labels in bases if not labels)}")
    emit(f"modp_cover_boundary_bases {sum(1 for _, _, labels in bases if labels)}")
    emit(f"normalized_p2_lifts {len(bases) * p**3}")
    emit(
        "p2_lift_status "
        + format_counter(lift_status, ("killed", "deep", "resolved"))
    )
    emit(
        "p2_compatible_lifts "
        + str(lift_status["deep"] + lift_status["resolved"])
    )
    emit(
        "modp_base_outcomes "
        + format_counter(
            base_outcomes,
            ("all_killed", "deep_only", "resolved_only", "resolved_and_deep"),
        )
    )
    emit("P2_COMPONENT_MATRIX")
    for label in sorted(component_bases):
        emit(
            f"{label} bases={component_bases[label]} "
            f"lifts({format_counter(component_lifts[label], ('killed','deep','resolved'))}) "
            f"base_outcomes({format_counter(component_outcomes[label], ('all_killed','deep_only','resolved_only','resolved_and_deep'))})"
        )
    emit("P2_SIGNATURE_OUTCOMES")
    for signature in sorted(signature_bases):
        emit(
            f"{signature} bases={signature_bases[signature]} "
            + format_counter(
                signature_outcomes[signature],
                ("all_killed", "deep_only", "resolved_only", "resolved_and_deep"),
            )
        )
    emit("P2_TOP_SQUARECLASS_PATTERNS")
    for pattern, count in squareclass_patterns.most_common(20):
        emit(f"{count} {pattern}")
    emit(f"P2_BOUNDARY_SECONDS {perf_counter() - t0:.3f}")
    return {
        "bases": len(bases),
        "lift_status": lift_status,
        "base_outcomes": base_outcomes,
    }


def qslice_data_mod_p(a: int, b: int, p: int):
    if a % p == 0 or b % p == 0:
        return None
    c = pow((a * b) % p, -1, p)
    r = (a * a + b * b + c * c - 3) % p
    s = (pow(a, -2, p) + pow(b, -2, p) + pow(c, -2, p) - 3) % p
    return c, r, s


def any_partition_cover_mod_p(squares: list[int], p: int) -> tuple[str, ...]:
    hits: list[str] = []
    for name, order in PARTITIONS:
        radicands = second_half_radicands_from_squares(squares, order)
        if all(is_square_mod_p_allow_zero(int(z), p) for z in radicands):
            hits.append(name)
    return tuple(hits)


def qslice_mask(p: int) -> tuple[set[tuple[int, int]], Counter]:
    """Necessary affine mask; R=0 is retained as unresolved chart boundary."""
    allowed: set[tuple[int, int]] = set()
    stats: Counter = Counter()
    for a in range(1, p):
        for b in range(1, p):
            c, r, s = qslice_data_mod_p(a, b, p)
            if legendre_or_zero(r, p) < 0 or legendre_or_zero(s, p) < 0:
                stats["square_equation_kill"] += 1
                continue
            if r == 0:
                allowed.add((a, b))
                stats["rho_boundary_retained"] += 1
                continue
            d2 = s * pow(r, -1, p) % p
            squares = [a * a % p, b * b % p, c * c % p, d2]
            hits = any_partition_cover_mod_p(squares, p)
            if hits:
                allowed.add((a, b))
                stats["cover_retained"] += 1
                if s == 0:
                    stats["sigma_boundary_retained"] += 1
            else:
                stats["cover_kill"] += 1
    stats["allowed"] = len(allowed)
    return allowed, stats


def run_qslice_finite(primes: list[int], emit) -> dict[int, dict]:
    results = {}
    emit("QSQUARE_FINITE_INTERSECTION")
    for p in primes:
        stats: Counter = Counter()
        cover_signatures: Counter = Counter()
        partition_hits: Counter = Counter()
        samples: list[tuple] = []
        for a in range(1, p):
            for b in range(1, p):
                c, r, s = qslice_data_mod_p(a, b, p)
                lr = legendre_or_zero(r, p)
                ls = legendre_or_zero(s, p)
                if lr < 0 or ls < 0:
                    stats["square_kill"] += 1
                    continue
                if r == 0 or s == 0:
                    if r == 0:
                        stats["rho_boundary"] += 1
                    if s == 0:
                        stats["sigma_boundary"] += 1
                    continue
                stats["slice_open"] += 1
                d2 = s * pow(r, -1, p) % p
                squares = [a * a % p, b * b % p, c * c % p, d2]
                labels = unsigned_boundary_labels_from_squares(squares, p)
                if labels:
                    stats["branch_boundary"] += 1
                else:
                    stats["branch_smooth"] += 1
                hits = any_partition_cover_mod_p(squares, p)
                if not hits:
                    continue
                stats["cover_any"] += 1
                stats["cover_partition_hits"] += len(hits)
                for hit in hits:
                    partition_hits[hit] += 1
                signature = "+".join(labels) if labels else "smooth"
                cover_signatures[signature] += 1
                if labels:
                    stats["cover_boundary"] += 1
                else:
                    stats["cover_smooth"] += 1
                    if len(samples) < 5:
                        samples.append((a, b, c, d2, hits))
        emit(
            f"QSQUARE_FINITE p={p} "
            + format_counter(
                stats,
                (
                    "slice_open",
                    "rho_boundary",
                    "sigma_boundary",
                    "branch_smooth",
                    "branch_boundary",
                    "cover_any",
                    "cover_smooth",
                    "cover_boundary",
                    "cover_partition_hits",
                ),
            )
        )
        emit(f"QSQUARE_FINITE_PARTITIONS p={p} {dict(sorted(partition_hits.items()))}")
        emit(f"QSQUARE_FINITE_COVER_SIGNATURES p={p} {dict(cover_signatures.most_common())}")
        emit(f"QSQUARE_FINITE_SMOOTH_SAMPLES p={p} {samples}")
        results[p] = {
            "stats": stats,
            "signatures": cover_signatures,
            "partitions": partition_hits,
        }
    return results


def positive_height_rationals(height: int) -> list[tuple[Fraction, int, int]]:
    vals: list[tuple[Fraction, int, int]] = []
    for den in range(1, height + 1):
        for num in range(1, height + 1):
            if gcd(num, den) == 1:
                vals.append((Fraction(num, den), num, den))
    vals.sort(key=lambda z: z[0])
    return vals


def rational_square_root(q: Fraction) -> tuple[bool, Fraction]:
    if q < 0:
        return False, Fraction(0)
    sn = isqrt(q.numerator)
    sd = isqrt(q.denominator)
    if sn * sn != q.numerator or sd * sd != q.denominator:
        return False, Fraction(0)
    return True, Fraction(sn, sd)


def fraction_text(q: Fraction) -> str:
    return f"{q.numerator}/{q.denominator}"


def qslice_exact_data(a: Fraction, b: Fraction):
    c = 1 / (a * b)
    r = a * a + b * b + c * c - 3
    s = 1 / (a * a) + 1 / (b * b) + 1 / (c * c) - 3
    okr, rho = rational_square_root(r)
    oks, sigma = rational_square_root(s)
    if not okr or not oks or rho == 0:
        return None
    d2 = s / r
    return c, r, s, rho, sigma, d2


def exact_cover_partitions(squares: list[Fraction]) -> tuple[str, ...]:
    hits: list[str] = []
    for name, order in PARTITIONS:
        radicands = second_half_radicands_from_squares(squares, order)
        if all(rational_square_root(Fraction(z))[0] and z != 0 for z in radicands):
            hits.append(name)
    return tuple(hits)


def run_qslice_search(
    height: int,
    primes: list[int],
    candidate_path: Path,
    emit,
) -> dict:
    t0 = perf_counter()
    vals = positive_height_rationals(height)
    emit(f"QSQUARE_CRT_SEARCH_START height={height} positive_values={len(vals)}")

    masks: dict[int, set[tuple[int, int]]] = {}
    mask_stats: dict[int, Counter] = {}
    for p in primes:
        masks[p], mask_stats[p] = qslice_mask(p)
    # Low density first is a substantial speedup.  The ordering is reported
    # so the first-kill counts remain reproducible.
    ordered_primes = sorted(primes, key=lambda p: len(masks[p]) / (p - 1) ** 2)
    emit(f"QSQUARE_MASK_ORDER {ordered_primes}")
    for p in ordered_primes:
        emit(
            f"QSQUARE_MASK p={p} size={len(masks[p])} ambient={(p-1)**2} "
            + format_counter(
                mask_stats[p],
                (
                    "square_equation_kill",
                    "cover_kill",
                    "rho_boundary_retained",
                    "sigma_boundary_retained",
                    "cover_retained",
                ),
            )
        )

    residues: dict[int, list[int | None]] = {}
    for p in ordered_primes:
        row: list[int | None] = []
        for _, num, den in vals:
            if den % p == 0:
                row.append(None)
            else:
                row.append((num % p) * pow(den % p, -1, p) % p)
        residues[p] = row

    counts: Counter = Counter()
    first_kill: Counter = Counter()
    boundary_bypass: Counter = Counter()
    candidates: list[dict] = []
    local_survivors: list[tuple[Fraction, Fraction]] = []
    exact_double_samples: list[tuple[Fraction, Fraction, Fraction, tuple[str, ...]]] = []
    exact_collision_signatures: Counter = Counter()

    for i, (a, an, ad) in enumerate(vals):
        for j in range(i, len(vals)):
            b, bn, bd = vals[j]
            # A <= B <= C=1/(AB), without constructing C for every pair.
            if an * bn * bn > ad * bd * bd:
                break
            counts["ordered_pairs"] += 1

            survived = True
            for p in ordered_primes:
                ra = residues[p][i]
                rb = residues[p][j]
                if ra is None or rb is None or ra == 0 or rb == 0:
                    boundary_bypass[p] += 1
                    continue
                if (ra, rb) not in masks[p]:
                    first_kill[p] += 1
                    survived = False
                    break
            if not survived:
                continue
            counts["mask_survivors"] += 1
            if len(local_survivors) < 20:
                local_survivors.append((a, b))

            exact = qslice_exact_data(a, b)
            if exact is None:
                counts["exact_double_square_fail"] += 1
                continue
            c, r, s, rho, sigma, d2 = exact
            counts["exact_double_square"] += 1
            squares = [a * a, b * b, c * c, d2]
            smooth = len(set(squares)) == 4 and all(z != 0 for z in squares)
            if smooth:
                counts["smooth_double_square"] += 1
            else:
                counts["collision_double_square"] += 1
                collision_labels = []
                for ii in range(4):
                    for jj in range(ii + 1,4):
                        if squares[ii] == squares[jj]:
                            collision_labels.append(f"E{ii+1}{jj+1}")
                signature = tuple(collision_labels)
                exact_collision_signatures["+".join(signature)] += 1
            if len(exact_double_samples) < 40:
                exact_double_samples.append((a,b,c,tuple(collision_labels) if not smooth else ()))
            hits = exact_cover_partitions(squares)
            if not hits:
                counts["exact_second_half_fail"] += 1
                continue
            counts["exact_second_half"] += 1
            if not smooth:
                counts["singular_second_half"] += 1
                continue

            # Recover the direct A(2,2,2,12) parameters.
            scale = 1 / (2 * rho)
            t = scale * scale
            u = 2 * t
            aa = scale * a
            bb = scale * b
            cc = scale * c
            dd = 2 * scale * scale * sigma
            v = 3 * t * t + dd * dd / 2
            candidate = {
                "A": a,
                "B": b,
                "C": c,
                "rho": rho,
                "sigma": sigma,
                "a": aa,
                "b": bb,
                "c": cc,
                "d": dd,
                "u": u,
                "t": t,
                "v": v,
                "partitions": hits,
            }
            candidates.append(candidate)
            counts["smooth_candidates"] += 1

    with candidate_path.open("w") as out:
        out.write("# A B C rho sigma a b c d u t v partitions\n")
        for cand in candidates:
            fields = [
                cand[key]
                for key in ("A", "B", "C", "rho", "sigma", "a", "b", "c", "d", "u", "t", "v")
            ]
            out.write(" ".join(fraction_text(q) for q in fields))
            out.write(" " + ",".join(cand["partitions"]) + "\n")

    emit(
        "QSQUARE_CRT_COUNTS "
        + format_counter(
            counts,
            (
                "ordered_pairs",
                "mask_survivors",
                "exact_double_square_fail",
                "exact_double_square",
                "smooth_double_square",
                "collision_double_square",
                "exact_second_half_fail",
                "exact_second_half",
                "singular_second_half",
                "smooth_candidates",
            ),
        )
    )
    emit(f"QSQUARE_FIRST_KILL {dict((p, first_kill[p]) for p in ordered_primes)}")
    emit(f"QSQUARE_BOUNDARY_BYPASS {dict((p, boundary_bypass[p]) for p in ordered_primes)}")
    emit(
        "QSQUARE_MASK_SURVIVOR_SAMPLES "
        + str([(fraction_text(a), fraction_text(b)) for a, b in local_survivors])
    )
    emit(f"QSQUARE_EXACT_COLLISION_SIGNATURES {dict(exact_collision_signatures)}")
    emit(
        "QSQUARE_EXACT_DOUBLE_SAMPLES "
        + str([
            (fraction_text(a),fraction_text(b),fraction_text(c),signature)
            for a,b,c,signature in exact_double_samples
        ])
    )
    emit(f"QSQUARE_CANDIDATE_FILE {candidate_path}")
    emit(f"QSQUARE_CRT_SECONDS {perf_counter() - t0:.3f}")
    return {
        "counts": counts,
        "first_kill": first_kill,
        "candidates": candidates,
    }


def parse_primes(text: str) -> list[int]:
    return [int(z) for z in text.split(",") if z]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=160)
    parser.add_argument("--boundary-prime", type=int, default=11)
    parser.add_argument(
        "--finite-primes",
        default="5,7,11,13,17,19,23,29,31,37,41,43,47,53",
    )
    parser.add_argument(
        "--search-primes",
        default="11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97",
    )
    parser.add_argument(
        "--log",
        type=Path,
        default=Path("results/target_22412_boundary_crt.log"),
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        default=Path("results/target_22412_boundary_crt_candidates.txt"),
    )
    parser.add_argument("--skip-boundary", action="store_true")
    parser.add_argument("--skip-finite", action="store_true")
    parser.add_argument("--skip-search", action="store_true")
    args = parser.parse_args()

    args.log.parent.mkdir(parents=True, exist_ok=True)
    args.candidates.parent.mkdir(parents=True, exist_ok=True)
    with args.log.open("w") as log:
        def emit(line: str) -> None:
            print(line, flush=True)
            log.write(line + "\n")
            log.flush()

        emit("TARGET_22412_BOUNDARY_CRT_START")
        emit(
            f"CONFIG height={args.height} boundary_prime={args.boundary_prime} "
            f"finite_primes={args.finite_primes} search_primes={args.search_primes}"
        )
        if not args.skip_boundary:
            run_p2_boundary_classifier(args.boundary_prime, emit)
        if not args.skip_finite:
            run_qslice_finite(parse_primes(args.finite_primes), emit)
        if not args.skip_search:
            run_qslice_search(
                args.height,
                parse_primes(args.search_primes),
                args.candidates,
                emit,
            )
        emit("TARGET_22412_BOUNDARY_CRT_DONE")


if __name__ == "__main__":
    main()
