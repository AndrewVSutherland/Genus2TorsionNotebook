#!/usr/bin/env python3
"""Item 1, Step 2: complete twisted 2-descent sweep.

For every primitive tuple in the corrected tor2244 list (torsion >= (2,2,4,4)):
test, for each of the three components of the (2,2,4,8) locus
   I:   4G = T_0      (class family D0 + J[2])
   II:  4G = T_AB     (class family H_AB + J[2])
   III: 4G = T_CD     (class family D0+H_AB + J[2])
and each of the 16 twists, whether delta(class+twist) is trivial.
Any hit = a curve with torsion containing (2,2,4,8).

For every tuple in tor2228 (torsion >= (2,2,2,8)): test divisibility of all
15 nonzero 2-torsion classes (any finite-pair hit = (2,2,4,4) part => (2,2,4,8)).

delta convention (Magma-validated, 0 mismatches on 600+ ground-truth checks):
component at root e = prod over finite support P of (x(P)-e), with f'(e)
replacing vanishing factors; for monic U: (-1)^deg U * U(e).
"""
from fractions import Fraction
import math, sys

def is_sq(q):
    if q < 0:
        return False
    if isinstance(q, Fraction):
        n, d = q.numerator, q.denominator
        r = math.isqrt(n * d)
        return r * r == n * d
    r = math.isqrt(q)
    return r * r == q

def isqrt_exact(n):
    r = math.isqrt(n)
    assert r * r == n, n
    return r

def delta_of_U(Ucoeffs, roots, fprime):
    """delta vector of monic U (list of coeffs low->high, monic), deg 1 or 2."""
    deg = len(Ucoeffs) - 1
    vec = []
    for i, e in enumerate(roots):
        v = sum(Fraction(c) * e**k for k, c in enumerate(Ucoeffs))
        if v != 0:
            vec.append((-1)**deg * v)
        else:
            # U1 = U/(x-e); for deg2 monic U with roots e, r2: U1 = x - r2, U1(e) = e - r2
            if deg == 2:
                r2 = -Ucoeffs[1] - e          # sum of roots = -Ucoeffs[1]
                v1 = e - r2
                if v1 != 0:
                    vec.append(fprime[i] * (-1) * v1)
                else:
                    vec.append(Fraction(1))
            else:  # deg 1, U = x - e
                vec.append(fprime[i])
    return vec

def trivial(vec):
    return all(is_sq(v) for v in vec)

def nsq(vec):
    return sum(1 for v in vec if is_sq(v))

def load(path):
    out = []
    for line in open(path):
        line = line.strip()
        if line.startswith('['):
            out.append(eval(line))
    return out

def analyze_tor2244(path):
    tuples = load(path)
    print(f"tor2244 corrected primitives: {len(tuples)}")
    famnames = ["I:D0(4G=T0)", "II:H_AB(4G=TAB)", "III:D0+H(4G=TCD)"]
    hists = [dict() for _ in range(3)]     # per family: best-twist #square-components histogram
    hits = [[] for _ in range(3)]
    twist_pass_counts = [dict() for _ in range(3)]  # which twists ever get >=4 sq comps
    for cnt, (a, b, c, d) in enumerate(tuples):
        A, B, C, D = a*a, b*b, c*c, d*d
        roots = [Fraction(0), Fraction(-A), Fraction(-B), Fraction(-C), Fraction(-D)]
        fprime = []
        for i, e in enumerate(roots):
            p = Fraction(1)
            for j, e2 in enumerate(roots):
                if j != i:
                    p *= (e - e2)
            fprime.append(p)
        # D0
        s1, s2 = a+b+c+d, a*b+a*c+a*d+b*c+b*d+c*d
        s3, s4 = a*b*c+a*b*d+a*c*d+b*c*d, a*b*c*d
        vD0 = delta_of_U([Fraction(s4), Fraction(-s2), Fraction(1)], roots, fprime)
        # H_AB (halving.m formulas)
        u0 = isqrt_exact((A-C)*(A-D)); v0 = isqrt_exact((B-C)*(B-D))
        w0 = isqrt_exact((A-C)*(B-C)); t0 = isqrt_exact((A-D)*(B-D))
        if u0*v0 != w0*t0:
            t0 = -t0
        assert u0*v0 == w0*t0
        rho = Fraction(a, b); sig = Fraction(A-C, w0); tau = Fraction(A-D, t0)
        s2g = rho+sig+tau+rho*sig+rho*tau+sig*tau
        s4g = rho*sig*tau
        dl = 1+s2g+s4g
        UH = [ (A*A + A*B*s2g + B*B*s4g)/dl, (2*A + s2g*(A+B) + 2*B*s4g)/dl, Fraction(1) ]
        vH = delta_of_U(UH, roots, fprime)
        vD0H = [x*y for x, y in zip(vD0, vH)]
        # twists: identity + 15 nonzero 2-torsion classes
        twists = [("1", [Fraction(1)]*5)]
        for i in range(5):
            twists.append((f"e{i}inf", delta_of_U([-roots[i], Fraction(1)], roots, fprime)))
        for i in range(5):
            for j in range(i+1, 5):
                U = [roots[i]*roots[j], -(roots[i]+roots[j]), Fraction(1)]
                twists.append((f"e{i}{j}", delta_of_U(U, roots, fprime)))
        for fi, base in enumerate([vD0, vH, vD0H]):
            best = 0
            for tname, tv in twists:
                prod = [x*y for x, y in zip(base, tv)]
                k = nsq(prod)
                if k > best:
                    best = k
                if k >= 4:
                    twist_pass_counts[fi][tname] = twist_pass_counts[fi].get(tname, 0) + 1
                if k == 5:
                    hits[fi].append(((a, b, c, d), tname))
                    print(f"*** HIT family {famnames[fi]} tuple {(a,b,c,d)} twist {tname}")
            hists[fi][best] = hists[fi].get(best, 0) + 1
        if (cnt+1) % 2000 == 0:
            print(f"  ...{cnt+1} tuples", file=sys.stderr)
    for fi in range(3):
        print(f"family {famnames[fi]}: best-twist #sq-components histogram {dict(sorted(hists[fi].items()))}")
        print(f"  twists reaching >=4 components: {twist_pass_counts[fi]}")
        print(f"  FULL HITS: {len(hits[fi])} {hits[fi][:10]}")

def analyze_tor2228(path):
    tuples = load(path)
    print(f"\ntor2228 primitives: {len(tuples)}")
    divisible_pairs = dict()
    hits = []
    for (a, b, c, d) in tuples:
        A, B, C, D = a*a, b*b, c*c, d*d
        roots = [Fraction(0), Fraction(-A), Fraction(-B), Fraction(-C), Fraction(-D)]
        fprime = []
        for i, e in enumerate(roots):
            p = Fraction(1)
            for j, e2 in enumerate(roots):
                if j != i:
                    p *= (e - e2)
            fprime.append(p)
        names = ["0", "A", "B", "C", "D"]
        for i in range(5):
            v = delta_of_U([-roots[i], Fraction(1)], roots, fprime)
            if trivial(v):
                key = f"e({names[i]},inf)"
                divisible_pairs[key] = divisible_pairs.get(key, 0) + 1
                if i != 0:
                    hits.append(((a, b, c, d), key))
        for i in range(5):
            for j in range(i+1, 5):
                U = [roots[i]*roots[j], -(roots[i]+roots[j]), Fraction(1)]
                v = delta_of_U(U, roots, fprime)
                if trivial(v):
                    key = f"e({names[i]},{names[j]})"
                    divisible_pairs[key] = divisible_pairs.get(key, 0) + 1
                    hits.append(((a, b, c, d), key))
                    print(f"*** 2244-part HIT on tor2228 tuple {(a,b,c,d)} class {key}")
    print(f"divisible 2-torsion classes across list: {divisible_pairs}")
    print(f"hits beyond T_0: {len(hits)}")

if __name__ == "__main__":
    import os, sys
    # resolve tuple banks relative to the repository checkout (script lives in code/);
    # optional overrides: claude_twist_sweep.py [tor2244_bank] [tor2228_bank]
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    bank2244 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(repo, "data", "tor2244_all_primitives.txt")
    bank2228 = sys.argv[2] if len(sys.argv) > 2 else os.path.join(repo, "data", "tor2228.txt")
    analyze_tor2244(bank2244)
    analyze_tor2228(bank2228)
