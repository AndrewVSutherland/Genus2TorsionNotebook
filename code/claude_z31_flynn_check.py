#!/usr/bin/env python3
# Validate the parsed Flynn polynomials against B1's Magma-generated TSVs:
#  - chartU.tsv : ladder(order)*D == 0, ladder(order/31)*D != 0 iff div31flag
#  - kummer.tsv : kappa(D) and kappa(31 D) match Magma's normalized vectors
#  - kummeqn(kappa) == 0 on every constructed point
import json, sys

import os
_repo_data = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
sdir = sys.argv[1] if len(sys.argv) > 1 else _repo_data
_cand = [f'{sdir}/flynn_polys.json', f'{sdir}/claude_z31_flynn_polys.json']
J = json.load(open(next(p for p in _cand if os.path.exists(p))))
VARS = J['vars']

def compile_poly(name):
    """poly -> python function(p, vals[15]) via exec, mod-p per term"""
    terms = J['polys'][name]
    lines = ['def _f(p, V):']
    # precompute needed powers
    maxe = {}
    for e, c in terms:
        for i, ei in enumerate(e):
            if ei > 0: maxe[i] = max(maxe.get(i, 0), ei)
    for i in sorted(maxe):
        lines.append(f'    v{i}_1 = V[{i}]')
        for ei in range(2, maxe[i] + 1):
            lines.append(f'    v{i}_{ei} = (v{i}_{ei-1} * V[{i}]) % p')
    lines.append('    acc = 0')
    for e, c in terms:
        facs = [f'v{i}_{ei}' for i, ei in enumerate(e) if ei > 0]
        expr = '*'.join([str(c)] + facs) if facs else str(c)
        lines.append(f'    acc = (acc + {expr}) % p')
    lines.append('    return acc % p')
    src = '\n'.join(lines)
    g = {}
    exec(src, g)
    return g['_f']

delta = [compile_poly(f'delta{i}') for i in range(1, 5)]
kummeqn = compile_poly('kummeqn')
bbrow1 = [compile_poly(f'BBB[{i},1]') for i in range(1, 5)]

def f_from_chart(p, u1, u0, v1, v0, w4, w3, w2, w1, w0):
    f6 = w4 % p
    f5 = (w3 + u1*w4) % p
    f4 = (w2 + u1*w3 + u0*w4) % p
    f3 = (w1 + u1*w2 + u0*w3) % p
    f2 = (w0 + u1*w1 + u0*w2 + v1*v1) % p
    f1 = (u1*w0 + u0*w1 + 2*v1*v0) % p
    f0 = (u0*w0 + v0*v0) % p
    return [f0, f1, f2, f3, f4, f5, f6]

def kappa_from_chart(p, u1, u0, v1, v0, f):
    s1, s2 = (-u1) % p, u0 % p
    f0, f1_, f2_, f3_, f4_, f5_, f6_ = f
    F0 = (2*f0 + f1_*s1 + 2*f2_*s2 + f3_*s1*s2 + 2*f4_*s2*s2
          + f5_*s1*s2*s2 + 2*f6_*s2*s2*s2) % p
    y1y2 = (v1*v1*u0 - v1*v0*u1 + v0*v0) % p
    den = (u1*u1 - 4*u0) % p
    # projective rep: (den, -u1*den, u0*den, F0 - 2*y1y2)
    return [den % p, (-u1*den) % p, (u0*den) % p, (F0 - 2*y1y2) % p], den % p

def is_zero_pt(k, p):
    return k[0] == 0 and k[1] == 0 and k[2] == 0 and k[3] != 0

def all_zero(k):
    return all(x == 0 for x in k)

def dbl(p, f, x):
    V = f + x + [0, 0, 0, 0]
    return [d(p, V) for d in delta]

def padd_col1(p, f, x, y, w):
    V = f + x + y
    B = [b(p, V) for b in bbrow1]
    w1 = w[0]
    return [(2*w1*B[i] - B[0]*w[i]) % p for i in range(4)]

def ladder(p, f, base, n):
    """kappa(n * D) from base = kappa(D); requires base[0] != 0 for pseudo-add"""
    if n == 1: return list(base)
    bits = bin(n)[2:]
    R0, R1 = list(base), dbl(p, f, base)
    for b in bits[1:]:
        if b == '0':
            R1 = padd_col1(p, f, R0, R1, base)
            R0 = dbl(p, f, R0)
        else:
            R0 = padd_col1(p, f, R0, R1, base)
            R1 = dbl(p, f, R1)
    return R0

def normalize(k, p):
    for x in k:
        if x != 0:
            inv = pow(x, p - 2, p)
            return [(inv * y) % p for y in k]
    return list(k)

def run_chartU(path, limit=None):
    nrows = eqfail = ordfail = d31fail = unk = 0
    sieve31_pass = 0
    for ln in open(path):
        c = ln.split()
        if not c: continue
        p = int(c[0]); u1, u0, v1, v0 = (int(x) for x in c[1:5])
        w4, w3, w2, w1, w0 = (int(x) for x in c[5:10])
        order, flag = int(c[10]), int(c[11])
        nrows += 1
        f = f_from_chart(p, u1, u0, v1, v0, w4, w3, w2, w1, w0)
        base, den = kappa_from_chart(p, u1, u0, v1, v0, f)
        if den == 0:
            unk += 1
            print(f'  UNKNOWN (den=0): p={p} row {nrows}')
            continue
        if kummeqn(p, f + base + [0]*4) != 0:
            eqfail += 1; print(f'  EQFAIL p={p} row {nrows}')
        kN = ladder(p, f, base, order)
        if not is_zero_pt(kN, p):
            ordfail += 1; print(f'  ORDFAIL p={p} order={order} row {nrows} kN={kN}')
        if flag:
            kM = ladder(p, f, base, order // 31)
            if is_zero_pt(kM, p) or all_zero(kM):
                d31fail += 1; print(f'  D31FAIL p={p} order={order} row {nrows}')
        k31 = ladder(p, f, base, 31)
        if is_zero_pt(k31, p):
            sieve31_pass += 1   # expect 0: no row has order in {1,31}
        if limit and nrows >= limit: break
    print(f'chartU: rows={nrows} unknown={unk} eqfail={eqfail} '
          f'ordfail={ordfail} d31fail={d31fail} sieve31_pass={sieve31_pass}')
    return eqfail + ordfail + d31fail

def run_kummer(path):
    nrows = kfail = k31fail = 0
    for ln in open(path):
        c = ln.split()
        if not c: continue
        p = int(c[0]); u1, u0, v1, v0 = (int(x) for x in c[1:5])
        w4, w3, w2, w1, w0 = (int(x) for x in c[5:10])
        kD_ref = [int(x) % p for x in c[12:16]]
        k31_ref = [int(x) % p for x in c[16:20]]
        nrows += 1
        f = f_from_chart(p, u1, u0, v1, v0, w4, w3, w2, w1, w0)
        base, den = kappa_from_chart(p, u1, u0, v1, v0, f)
        if den == 0:
            print(f'  row {nrows}: den=0, comparing k2,k3 only? skip'); continue
        if normalize(base, p) != kD_ref:
            kfail += 1; print(f'  KFAIL p={p} row {nrows}: {normalize(base,p)} vs {kD_ref}')
        k31 = ladder(p, f, base, 31)
        if normalize(k31, p) != k31_ref:
            k31fail += 1
            print(f'  K31FAIL p={p} row {nrows}: {normalize(k31,p)} vs {k31_ref}')
    print(f'kummer: rows={nrows} kappaD_fail={kfail} kappa31D_fail={k31fail}')
    return kfail + k31fail

if __name__ == '__main__':
    import os
    ddir = sys.argv[3] if len(sys.argv) > 3 else os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
    bad = run_kummer(os.path.join(ddir, 'claude_z31_vectors_kummer.tsv'))
    lim = int(sys.argv[2]) if len(sys.argv) > 2 else None
    bad += run_chartU(os.path.join(ddir, 'claude_z31_vectors_chartU.tsv'), lim)
    print('TOTAL FAILURES:', bad)
    sys.exit(0 if bad == 0 else 1)
