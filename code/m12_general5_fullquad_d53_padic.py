#!/usr/bin/env python3
"""7-adic graph search on the unique smooth d=5/3 branch modulo 7.

At p=7 the fixed-d slice has the single smooth/open point

    (b,w,c,d,s,t,a)=(3,5,6,4,1,4,5).

The residual Jacobian in (b,w,s,t,a) is invertible, so c is an etale local
parameter.  For every small rational c congruent to 6 modulo 7 this script
Newton-lifts the other five coordinates modulo 7^N, rationally reconstructs
them, and applies the complete exact identity check.  The computation stores
only a few integers regardless of the height or precision.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent


def import_file(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


D53 = import_file("fullquad_d53_padic_search", "m12_general5_fullquad_d53_search.py")
FQ = D53.FQ
P = 7
BASE_ALL = (3, 5, 6, 4, 1, 4, 5)  # b,w,c,d,s,t,a modulo 7
VAR_COLS = (0, 1, 4, 5, 6)
BASE_X = (3, 5, 1, 4, 5)


def matrix_inverse(A, p: int):
    n = len(A)
    M = [[x % p for x in row]+[1 if i == j else 0 for j in range(n)]
         for i, row in enumerate(A)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if M[i][j]), None)
        if pivot is None:
            raise ArithmeticError("singular Hensel Jacobian")
        M[j], M[pivot] = M[pivot], M[j]
        z = pow(M[j][j], -1, p)
        M[j] = [(z*x) % p for x in M[j]]
        for i in range(n):
            if i != j and M[i][j]:
                z = M[i][j]
                M[i] = [(x-z*y) % p for x, y in zip(M[i], M[j])]
    return [row[n:] for row in M]


def matvec(A, x, p: int):
    return [sum(a*b for a, b in zip(row, x)) % p for row in A]


def hensel_inverse():
    rr = FQ.dual_residuals(BASE_ALL, P)
    if any(x.v for x in rr):
        raise ArithmeticError("base point does not satisfy residual equations")
    J = [[x.g[j] for j in VAR_COLS] for x in rr]
    expected = [
        [2, 3, 2, 5, 6],
        [2, 3, 3, 4, 5],
        [6, 2, 2, 1, 4],
        [5, 3, 2, 6, 1],
        [1, 3, 3, 3, 1],
    ]
    if J != expected:
        raise ArithmeticError(f"unexpected Jacobian {J}")
    return J, matrix_inverse(J, P)


JACOBIAN, JINV = hensel_inverse()


def rat_residue(x: Fraction, modulus: int):
    return x.numerator*pow(x.denominator, -1, modulus) % modulus


def lift_c(c: Fraction, precision: int):
    if c.denominator % P == 0 or rat_residue(c, P) != 6:
        return None
    x = list(BASE_X)
    pn = P
    for _n in range(1, precision):
        q = pn*P
        b, w, s, t, a = x
        vals = (b, w, rat_residue(c, q), rat_residue(D53.D_FIXED, q), s, t, a)
        residuals = FQ.residual_values(vals, q)
        if any(z % pn for z in residuals):
            raise ArithmeticError("lost the Hensel congruence")
        rhs = [(-(z//pn)) % P for z in residuals]
        digit = matvec(JINV, rhs, P)
        x = [(z+pn*e) % q for z, e in zip(x, digit)]
        pn = q
    vals = (x[0], x[1], rat_residue(c, pn), rat_residue(D53.D_FIXED, pn),
            x[2], x[3], x[4])
    if any(FQ.residual_values(vals, pn)):
        raise ArithmeticError("final Hensel control failed")
    return tuple(x), pn


def rational_cs(height: int):
    for den in range(1, height+1):
        if den % P == 0:
            continue
        for num in range(-height, height+1):
            if math.gcd(abs(num), den) != 1:
                continue
            c = Fraction(num, den)
            if rat_residue(c, P) == 6:
                yield c


def scan(height: int, precision: int):
    total = reconstructed = 0
    exact_rows, hits = set(), []
    for c in rational_cs(height):
        total += 1
        lifted, modulus = lift_c(c, precision)
        coords = tuple(D53.rational_reconstruction(z, modulus) for z in lifted)
        if any(z is None for z in coords):
            continue
        reconstructed += 1
        b, w, s, t, a = coords
        row = (b, w, c, D53.D_FIXED, s, t, a)
        if row in exact_rows:
            continue
        exact_rows.add(row)
        ans = D53.exact_match(row)
        if ans is not None:
            hit = row+(ans[0],)
            hits.append(hit)
            print("EXACT_7ADIC_HIT b,w,c,d,s,t,a,tau=", hit)
        if total % 10000 == 0:
            print(f"PADIC_PROGRESS c_values={total} reconstructed={reconstructed} hits={len(hits)}")
    print(f"PADIC_D53_DONE H={height} precision={precision} modulus={P**precision} "
          f"c_values={total} reconstructed={reconstructed} "
          f"unique_exact_rows={len(exact_rows)} exact_hits={len(hits)}")
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--height", type=int, default=100)
    ap.add_argument("--precision", type=int, default=20)
    args = ap.parse_args()
    if args.height < 1 or args.precision < 2:
        raise SystemExit("height >=1 and precision >=2 required")
    print("JACOBIAN_MOD7", JACOBIAN)
    scan(args.height, args.precision)


if __name__ == "__main__":
    main()
