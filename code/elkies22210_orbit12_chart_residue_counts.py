#!/usr/bin/env python3
"""Project the certified orbit-12 boundary lifts to the (t,m) chart.

This is a diagnostic, not a sufficient global sieve.  It distinguishes all
truncated square-residue states (including deep zero radicands) from states
whose four valuations have already resolved to finite even values.
"""

import importlib.util
from pathlib import Path


def load_boundary_module():
    path = Path(__file__).with_name("elkies22210_orbit12_boundary_lifts.py")
    spec = importlib.util.spec_from_file_location("orbit12_boundary", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def resolved(module, r, p, k):
    signature = module.valuation_signature(r, p, k)
    return all(s.startswith("V") and int(s[1:]) % 2 == 0 for s in signature)


def main():
    module = load_boundary_module()
    print("ELKIES22210_ORBIT12_CHART_RESIDUE_COUNTS")
    for p, k in ((11, 3), (19, 2), (23, 2)):
        states = module.base_points(p)
        base_ids = list(range(len(states)))
        for level in range(1, k):
            states, base_ids, _ = module.lift_level(
                states, p, level, base_ids
            )
        q = p**k
        chart_states = []
        chart_resolved = []
        for r, chart in states:
            # r1=1 is the normalization chart and r1+r2 a unit is the
            # integral inverse-chart locus used here.
            if chart != 0 or (r[0] + r[1]) % p == 0:
                continue
            t = (r[2] + r[3]) * pow((r[0] + r[1]) % q, -1, q) % q
            m = (r[2] + r[0] + r[1]) * pow(r[0] % q, -1, q) % q
            chart_states.append((t, m))
            if resolved(module, r, p, k):
                chart_resolved.append((t, m))
        all_pairs = set(chart_states)
        resolved_pairs = set(chart_resolved)
        print(
            "p", p, "k", k, "modulus", q,
            "all_cover_states", len(states),
            "inverse_chart_states", len(chart_states),
            "distinct_tm", len(all_pairs),
            "resolved_states", len(chart_resolved),
            "distinct_resolved_tm", len(resolved_pairs),
            "ambient_tm", q * q,
        )
    print("DONE")


if __name__ == "__main__":
    main()

