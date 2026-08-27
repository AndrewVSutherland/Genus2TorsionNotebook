#!/usr/bin/env python3
"""Low-memory modular sieve for the exact shared-Weierstrass [60] lane.

The square cover in ``m12_shared_weierstrass_order10_param.m`` is
parametrized by

    R = -z^2/b,  a = 2*b^4*(b-1)/z^5.

For fixed (b,z), the two remaining equations have degrees 2 and 4 in w.
This script builds exact finite-field projection masks by exhausting w, then
scans bounded-height rational (b,z) pairs.  It never stores a three-variable
table and uses only O(sum p^2) memory.

Examples:
  python3 code/m12_shared_weierstrass_order10_sieve.py local 7 11 13
  python3 code/m12_shared_weierstrass_order10_sieve.py search 50 \
      --primes 7 11 13 17 19 23 29 31
"""

from argparse import ArgumentParser
from math import gcd


def equations(b, w, z, p):
    """Return the two residual equations modulo p on the open chart."""
    inv = lambda x: pow(x % p, -1, p)
    rp = -z * z * inv(b) % p
    a = 2 * pow(b, 4, p) * (b - 1) * inv(pow(z, 5, p)) % p
    s = a * a % p

    g1 = (4*b**6*w**2 - 20*b**5*w**2 + 24*b**5*w
          + 41*b**4*w**2 - 72*b**4*w - 44*b**3*w**2
          + 4*b**4 + 78*b**3*w + 26*b**2*w**2 - 12*b**3
          - 36*b**2*w - 8*b*w**2 + 9*b**2 + 6*b*w + w**2) % p
    g2 = (-8*b**5*w**2 + 24*b**4*w**2 + 16*b**4*w
          - 26*b**3*w**2 - 32*b**3*w + 12*b**2*w**2
          + 12*b**2*w - 2*b*w**2 + 4*b*w - 6*b - 2*w) % p
    g3 = (4*b**4*w**2 - 4*b**3*w**2 + 24*b**3*w
          + b**2*w**2 - 24*b**2*w + 4*b**2 + 6*b*w
          + 4*b + 1) % p
    g4 = -4*b % p

    ccap = (g4 - 5*s*rp) * inv(2) % p
    dcap = (s*g3 - ccap*ccap + 10*s*s*rp*rp) * inv(2*s) % p
    e2 = (s*g2 - 2*ccap*dcap - 10*s*s*rp**3) % p
    e1 = (s*g1 - dcap*dcap + 5*s*s*rp**4) % p
    return e2, e1


def projection_mask(p, samples=8):
    if p in (2, 5):
        raise ValueError("use odd primes different from 5")
    allowed = set()
    triples = []
    for b in range(p):
        if b == 0 or b == 1 or (2*b-1) % p == 0:
            continue
        for z in range(1, p):
            roots = []
            for w in range(p):
                if equations(b, w, z, p) == (0, 0):
                    roots.append(w)
                    if len(triples) < samples:
                        triples.append((b, w, z, (-z*z*pow(b, -1, p)) % p))
            if roots:
                allowed.add((b, z))
    return allowed, triples


def rational_values(height):
    vals = []
    for den in range(1, height + 1):
        for num in range(-height, height + 1):
            if gcd(abs(num), den) == 1:
                vals.append((num, den))
    return vals


def reduced_residue(value, p, kind):
    num, den = value
    if den % p == 0:
        return None
    x = num * pow(den, -1, p) % p
    # Parametrization denominators, and the degree-five model, can lose
    # information on these reductions.  Passing them is conservative.
    if kind == "b" and (x == 0 or x == 1 or (2*x-1) % p == 0):
        return None
    if kind == "z" and x == 0:
        return None
    return x


def search(height, primes):
    masks = {}
    for p in primes:
        mask, triples = projection_mask(p)
        masks[p] = mask
        print(f"MASK p={p} allowed={len(mask)}/{p*p} samples={triples}")

    vals = rational_values(height)
    bvals = [v for v in vals if v not in ((0, 1), (1, 1), (1, 2))]
    zvals = [v for v in vals if v != (0, 1)]
    residues_b = {p: [reduced_residue(v, p, "b") for v in bvals]
                  for p in primes}
    residues_z = {p: [reduced_residue(v, p, "z") for v in zvals]
                  for p in primes}

    survivors = []
    tested = 0
    pass_counts = [0] * (len(primes) + 1)
    for ib, bv in enumerate(bvals):
        for iz, zv in enumerate(zvals):
            tested += 1
            pass_counts[0] += 1
            ok = True
            for j, p in enumerate(primes):
                br = residues_b[p][ib]
                zr = residues_z[p][iz]
                if br is not None and zr is not None and (br, zr) not in masks[p]:
                    ok = False
                    break
                pass_counts[j + 1] += 1
            if ok:
                survivors.append((bv, zv))
    print(f"SEARCH H={height} values={len(vals)} bvalues={len(bvals)} "
          f"zvalues={len(zvals)} tested={tested}")
    print("PASS_COUNTS", list(zip(["start"] + primes, pass_counts)))
    print(f"SURVIVORS {len(survivors)}")
    for bv, zv in survivors:
        print(f"  b={bv[0]}/{bv[1]} z={zv[0]}/{zv[1]}")


def main():
    parser = ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)
    loc = sub.add_parser("local")
    loc.add_argument("primes", type=int, nargs="+")
    loc.add_argument("--samples", type=int, default=8)
    sea = sub.add_parser("search")
    sea.add_argument("height", type=int)
    sea.add_argument("--primes", type=int, nargs="+", required=True)
    args = parser.parse_args()
    if args.mode == "local":
        for p in args.primes:
            mask, triples = projection_mask(p, args.samples)
            print(f"p={p} allowed_base={len(mask)}/{p*p} samples={triples}")
    else:
        search(args.height, args.primes)


if __name__ == "__main__":
    main()
