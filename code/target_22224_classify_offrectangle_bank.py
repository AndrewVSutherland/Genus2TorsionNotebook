#!/usr/bin/env python3
"""Classify primitive off-rectangle tor2228 rows against known P1 families."""

from ast import literal_eval
from fractions import Fraction as Q
from itertools import permutations, product
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BANK = ROOT / "data/tor2228_bank.txt"


def primitive(v):
    g = 0
    for x in v:
        g = gcd(g, abs(x))
    return tuple(x // g for x in v)


def rectangle(v):
    a, b, c, d = v
    return a * b == c * d or a * c == b * d or a * d == b * c


def family_tuple(fam, t):
    if fam == 1:
        return (-(t*t+t+1)**2, -4*t*(t+1)**2, 4*t*(t+1), 4*t*t*(t+1))
    if fam == 2:
        if t == 0:
            return None
        den=t**4+2*t**3-t*t-2*t+1
        if den == 0:
            return None
        return (-(t**4-2*t**3-t*t+2*t+1)/den, -1/t, Q(1), t)
    if fam == 4:
        den=(t*t-2*t-1)*(t*t+1)
        if t == 0 or den == 0:
            return None
        return (-t*(t+1)**2*(t-1)/den, -t*t, Q(1), t)
    raise ValueError(fam)


def proportional(x, y):
    if any(z == 0 for z in x) or any(z == 0 for z in y):
        return False
    lam = Q(x[0], 1) / y[0]
    return all(Q(a, 1) == lam*b for a, b in zip(x, y))


def memberships(v):
    out = set()
    for perm in permutations(v):
        for signs in product((-1, 1), repeat=4):
            z=tuple(signs[i]*perm[i] for i in range(4))
            # In all three formulae c=1 and d=t, so t=d/c.
            t=Q(z[3],z[2])
            for fam in (1,2,4):
                f=family_tuple(fam,t)
                if f is not None and proportional(z,f):
                    out.add((fam,t))
    return sorted(out)


def main():
    rows=[];seen=set()
    for raw in BANK.read_text().splitlines():
        if not raw.startswith("["):
            continue
        v=literal_eval(raw)
        if not isinstance(v,list) or len(v)!=4:
            continue
        v=primitive(v)
        if v in seen:
            continue
        seen.add(v)
        if not rectangle(v):
            rows.append(v)
    classified=[];unclassified=[]
    for v in rows:
        m=memberships(v)
        (classified if m else unclassified).append((v,m))
    print("OFFRECTANGLE",len(rows),"CLASSIFIED",len(classified),"UNCLASSIFIED",len(unclassified))
    for v,m in classified:
        print("KNOWN",v,[(f,str(t)) for f,t in m])
    for v,_ in unclassified:
        print("UNKNOWN",v)


if __name__ == "__main__":
    main()
