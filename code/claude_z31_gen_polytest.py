#!/usr/bin/env python3
# Generate random polytest vectors: p f0..f6 x1..x4 y1..y4 d1..d4 b11,b21,b31,b41 eq
# Values evaluated straight from the parsed Flynn JSON (pure polynomial eval,
# arbitrary points — NOT required to lie on the Kummer surface).
import json, random, sys
sdir = sys.argv[1]
J = json.load(open(f'{sdir}/flynn_polys.json'))
def evalp(name, p, V):
    acc = 0
    for e, c in J['polys'][name]:
        t = c
        for i, ei in enumerate(e):
            for _ in range(ei):
                t = (t * V[i]) % p
        acc = (acc + t) % p
    return acc % p
random.seed(31007)
rows = []
for p in [101, 211, 1009, 2003, 4001, 4003, 10007]:
    for _ in range(60):
        f = [random.randrange(p) for _ in range(7)]
        x = [random.randrange(p) for _ in range(4)]
        y = [random.randrange(p) for _ in range(4)]
        V = f + x + y
        Vx = f + x + [0, 0, 0, 0]
        d = [evalp(f'delta{i}', p, Vx) for i in range(1, 5)]
        b = [evalp(f'BBB[{i},1]', p, V) for i in range(1, 5)]
        eq = evalp('kummeqn', p, Vx)
        rows.append(' '.join(map(str, [p] + f + x + y + d + b + [eq])))
open(f'{sdir}/polytest_vectors.txt', 'w').write('\n'.join(rows) + '\n')
print(f'wrote {len(rows)} vectors')
