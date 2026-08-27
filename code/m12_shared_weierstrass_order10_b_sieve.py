#!/usr/bin/env python3
"""Sieve compact-M(12) parameter b, with no height bound on z.

The unique non-boundary resultant factor from
``m12_shared_weierstrass_order10_param.m`` has

    lc_z = -64/25*(b-1/2)^6,
    const_z = b^42*(b-1)^22.

Thus away from 2 and 5 a point with good b cannot hide at z=0 or z=infinity,
and the projection to b of the affine finite-field masks is rigorous.
"""

from argparse import ArgumentParser
from math import gcd

from m12_shared_weierstrass_order10_sieve import projection_mask


def search(height, primes):
    bmasks = {}
    for p in primes:
        mask, _ = projection_mask(p, 0)
        bmasks[p] = {b for b, _ in mask}
        print(f"BMASK p={p} allowed_open_b={len(bmasks[p])}/{p} "
              f"values={sorted(bmasks[p])}")

    pass_counts = [0] * (len(primes) + 1)
    survivors = []
    tested = 0
    for den in range(1, height + 1):
        for num in range(-height, height + 1):
            if gcd(abs(num), den) != 1:
                continue
            if num == 0 or num == den or 2*num == den:
                continue
            tested += 1
            pass_counts[0] += 1
            ok = True
            for j, p in enumerate(primes):
                if den % p == 0:
                    pass_counts[j + 1] += 1
                    continue
                br = num * pow(den, -1, p) % p
                if br in (0, 1, pow(2, -1, p)):
                    pass_counts[j + 1] += 1
                    continue
                if br not in bmasks[p]:
                    ok = False
                    break
                pass_counts[j + 1] += 1
            if ok:
                survivors.append((num, den))

    print(f"BSEARCH H={height} tested={tested}")
    print("PASS_COUNTS", list(zip(["start"] + primes, pass_counts)))
    print(f"B_SURVIVORS {len(survivors)}")
    for num, den in survivors:
        print(f"  b={num}/{den}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("height", type=int)
    parser.add_argument("--primes", type=int, nargs="+", required=True)
    args = parser.parse_args()
    search(args.height, args.primes)
