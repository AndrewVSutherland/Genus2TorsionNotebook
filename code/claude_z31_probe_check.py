#!/usr/bin/env python3
# Validate parsed delta (duplication) and FULL 4x4 BBB (biquadratics)
# individually against Magma Jacobian-arithmetic ground truth (probe31.log):
#   delta(kA) ∝ k2A ;  BBB(kA,kB) ∝ kApB⊗kAmB + kAmB⊗kApB ;  kummeqn(kA)=0
import json, sys, os

_repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_repo, 'data')
_cand = [f'{sdir}/flynn_polys.json', f'{sdir}/claude_z31_flynn_polys.json']
J = json.load(open(next(p for p in _cand if os.path.exists(p))))

def evalp(name, p, V):
    acc = 0
    for e, c in J['polys'][name]:
        t = c
        for i, ei in enumerate(e):
            for _ in range(ei):
                t = (t * V[i]) % p
        acc = (acc + t) % p
    return acc % p

def prop(a, b, p):
    """projective proportionality of vectors mod p (not both zero)"""
    if all(x == 0 for x in a) or all(x == 0 for x in b): return False
    lam = None
    for x, y in zip(a, b):
        if (x == 0) != (y == 0): return False
        if x != 0:
            l2 = (y * pow(x, p - 2, p)) % p
            if lam is None: lam = l2
            elif lam != l2: return False
    return True

_probe = sys.argv[2] if len(sys.argv) > 2 else os.path.join(_repo, 'results', 'claude_z31_kummer_probe.log')
nvec = dfail = bfail = efail = 0
for ln in open(_probe):
    c = ln.split()
    if not c or c[0] != 'VEC': continue
    p = int(c[1]); f = [int(x) for x in c[2:9]]
    kA   = [int(x) for x in c[9:13]]
    kB   = [int(x) for x in c[13:17]]
    kApB = [int(x) for x in c[17:21]]
    kAmB = [int(x) for x in c[21:25]]
    k2A  = [int(x) for x in c[25:29]]
    nvec += 1
    # kummeqn on all five points
    for pt in (kA, kB, kApB, kAmB, k2A):
        if evalp('kummeqn', p, f + pt + [0]*4) != 0:
            efail += 1; print(f'EQFAIL vec {nvec} p={p} pt={pt}')
    # duplication
    d = [evalp(f'delta{i}', p, f + kA + [0]*4) for i in range(1, 5)]
    if not prop(d, k2A, p):
        dfail += 1; print(f'DBLFAIL vec {nvec} p={p}: {d} vs {k2A}')
    # full 4x4 biquadratic identity, matrix proportionality
    M = [[evalp(f'BBB[{i},{j}]', p, f + kA + kB) for j in range(1, 5)]
         for i in range(1, 5)]
    Z, W = kApB, kAmB
    R = [[(Z[i]*W[j] + Z[j]*W[i]) % p for j in range(4)] for i in range(4)]
    flatM = [x for row in M for x in row]
    flatR = [x for row in R for x in row]
    if not prop(flatM, flatR, p):
        bfail += 1; print(f'BBFAIL vec {nvec} p={p}')
print(f'probe vectors: {nvec}  kummeqn fails: {efail}  '
      f'duplication fails: {dfail}  biquadratic fails: {bfail}')
sys.exit(0 if dfail + bfail + efail == 0 else 1)
