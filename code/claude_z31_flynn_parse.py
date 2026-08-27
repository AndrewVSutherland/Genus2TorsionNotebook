#!/usr/bin/env python3
# Mechanical parser for Flynn's genus-2 Kummer files (Maple format) -> exact
# integer multivariate polynomials, stored as {exponent-tuple: coeff} dicts.
# Variable order (fixed everywhere): f0..f6, k1..k4, l1..l4  (15 vars).
import json, sys, re

VARS = ['f0','f1','f2','f3','f4','f5','f6','k1','k2','k3','k4','l1','l2','l3','l4']
VIDX = {v: i for i, v in enumerate(VARS)}
NV = len(VARS)

# ---------------- polynomial engine: dict[tuple[int]*NV] -> int ----------------
def pnew(): return {}
def pconst(c):
    return {tuple([0]*NV): c} if c != 0 else {}
def pvar(name):
    e = [0]*NV; e[VIDX[name]] = 1
    return {tuple(e): 1}
def padd(a, b):
    r = dict(a)
    for e, c in b.items():
        s = r.get(e, 0) + c
        if s: r[e] = s
        elif e in r: del r[e]
    return r
def pneg(a): return {e: -c for e, c in a.items()}
def psub(a, b): return padd(a, pneg(b))
def pmul(a, b):
    r = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            s = r.get(e, 0) + ca * cb
            if s: r[e] = s
            elif e in r: del r[e]
    return r
def ppow(a, n):
    r = pconst(1)
    for _ in range(n): r = pmul(r, a)
    return r

# ---------------- tokenizer ----------------
TOKRE = re.compile(r'\s*(\*\*|\^|[0-9]+|[a-zA-Z][a-zA-Z0-9]*(\[[0-9],[0-9]\])?|[-+*()])')
def tokenize(s):
    toks, i = [], 0
    while i < len(s):
        m = TOKRE.match(s, i)
        if not m:
            raise ValueError(f'tokenize fail at {s[i:i+40]!r}')
        toks.append(m.group(1))
        i = m.end()
    return toks

# ---------------- recursive descent parser ----------------
class P:
    def __init__(self, toks, env):
        self.t, self.i, self.env = toks, 0, env
    def peek(self): return self.t[self.i] if self.i < len(self.t) else None
    def next(self): tok = self.t[self.i]; self.i += 1; return tok
    def expr(self):
        neg = False
        while self.peek() in ('+', '-'):
            if self.next() == '-': neg = not neg
        r = self.term()
        if neg: r = pneg(r)
        while self.peek() in ('+', '-'):
            op = self.next()
            neg2 = (op == '-')
            while self.peek() in ('+', '-'):
                if self.next() == '-': neg2 = not neg2
            t = self.term()
            r = psub(r, t) if neg2 else padd(r, t)
        return r
    def term(self):
        r = self.factor()
        while self.peek() == '*':
            self.next()
            r = pmul(r, self.factor())
        return r
    def factor(self):
        neg = False
        while self.peek() == '-':
            self.next(); neg = not neg
        b = self.base()
        if self.peek() in ('**', '^'):
            self.next()
            n = self.next()
            assert n.isdigit(), f'bad exponent {n}'
            b = ppow(b, int(n))
        return pneg(b) if neg else b
    def base(self):
        tok = self.next()
        if tok == '(':
            r = self.expr()
            assert self.next() == ')', 'missing )'
            return r
        if tok.isdigit():
            return pconst(int(tok))
        if tok in VIDX:
            return pvar(tok)
        if tok in self.env:
            return self.env[tok]
        raise ValueError(f'unknown symbol {tok}')

# ---------------- statement extraction ----------------
def statements(path):
    body = []
    for line in open(path):
        ls = line.strip()
        if ls.startswith('#'): continue
        if ls.startswith('interface('): continue
        if ls.startswith('with(linalg)'): continue
        body.append(line.rstrip('\n'))
    text = ' '.join(body)
    # split on ':' or ';' that are NOT part of ':='
    stmts, cur, i = [], '', 0
    while i < len(text):
        ch = text[i]
        if ch == ':' and i + 1 < len(text) and text[i+1] == '=':
            cur += ':='; i += 2; continue
        if ch in ':;':
            if cur.strip(): stmts.append(cur.strip())
            cur = ''; i += 1; continue
        cur += ch; i += 1
    if cur.strip(): stmts.append(cur.strip())
    return [s for s in stmts if ':=' in s]

def parse_file(path, env, skip_names=()):
    pending = []
    for s in statements(path):
        name, rhs = s.split(':=', 1)
        name, rhs = name.strip(), rhs.strip()
        if name in skip_names or name.startswith('Wg') or rhs.startswith('matrix('):
            continue
        pending.append((name, rhs))
    # multi-pass: resolve when all referenced symbols are known
    progress = True
    while pending and progress:
        progress, rest = False, []
        for name, rhs in pending:
            try:
                poly = P(tokenize(rhs), env).expr()
            except ValueError:
                rest.append((name, rhs)); continue
            env[name] = poly
            progress = True
        pending = rest
    if pending:
        raise ValueError(f'unresolved: {[n for n, _ in pending]}')

def kdeg(e):  return sum(e[7:11])
def ldeg(e):  return sum(e[11:15])
def fdeg(e):  return sum(e[0:7])

def main(sdir):
    env = {}
    parse_file(f'{sdir}/duplication', env)
    parse_file(f'{sdir}/defining.equations', env,
               skip_names=('Wg[%d,%d]' % (i, j) for i in range(1,5) for j in range(1,5)))
    parse_file(f'{sdir}/biquadratic.forms', env)

    # structural checks
    errs = 0
    for i in range(1, 5):
        d = env[f'delta{i}']
        assert all(kdeg(e) == 4 and ldeg(e) == 0 for e in d), f'delta{i} not k-quartic'
        print(f'delta{i}: {len(d)} terms, fdeg max {max(fdeg(e) for e in d)}')
    kq = env['kummeqn']
    assert all(kdeg(e) == 4 and ldeg(e) == 0 for e in kq), 'kummeqn not quartic'
    print(f'kummeqn: {len(kq)} terms')
    for i in range(1, 5):
        for j in range(1, 5):
            b = env[f'BBB[{i},{j}]']
            assert all(kdeg(e) == 2 and ldeg(e) == 2 for e in b), f'BBB[{i},{j}] not (2,2)'
    # symmetry BBB[i,j] == BBB[j,i]
    for i in range(1, 5):
        for j in range(i+1, 5):
            a, b = env[f'BBB[{i},{j}]'], env[f'BBB[{j},{i}]']
            if a != b:
                print(f'SYMMETRY-FAIL BBB[{i},{j}] != BBB[{j},{i}]'); errs += 1
    # k<->l swap symmetry of each entry
    def swap(poly):
        return {tuple(list(e[0:7]) + list(e[11:15]) + list(e[7:11])): c for e, c in poly.items()}
    for i in range(1, 5):
        for j in range(1, 5):
            a = env[f'BBB[{i},{j}]']
            if swap(a) != a:
                print(f'KL-SWAP-FAIL BBB[{i},{j}]'); errs += 1
    print('symmetry checks done, errs =', errs)

    out = {}
    for name in [f'delta{i}' for i in range(1,5)] + ['kummeqn'] + \
                [f'BBB[{i},{j}]' for i in range(1,5) for j in range(1,5)]:
        out[name] = [[list(e), c] for e, c in sorted(env[name].items())]
    with open(f'{sdir}/flynn_polys.json', 'w') as fh:
        json.dump({'vars': VARS, 'polys': out}, fh)
    print('wrote flynn_polys.json;', 'TOTAL terms:',
          sum(len(v) for v in out.values()))
    return 0 if errs == 0 else 1

if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
