#!/bin/bash
# claude_88_drive.sh — overnight (8,8) hunt driver: 8-way parallel Magma workers,
# one (m,n) base each, 900s timeout per base.  Logs in data/t88hunt/.
cd /home/claude/torsion_jac
PAR=${1:-8}; TMO=${2:-900}; HMAX=${3:-12}
OUT=data/t88hunt; mkdir -p $OUT
python3 - "$HMAX" > $OUT/bases.txt << 'EOF'
import sys
from math import gcd
H = int(sys.argv[1])
def rats(H, excl):
    out = []
    for b in range(1, H+1):
        for a in range(-H, H+1):
            if a == 0 or gcd(abs(a), b) != 1: continue
            for (p, q) in [(a, b), (b, a)] if b > 1 else [(a, b)]:
                if q == 0: continue
                v = (p, q)
                fv = p/q
                if fv in excl: continue
                if v not in out: out.append(v)
    return out
ms = rats(H, {0.0, 1.0, -1.0})
ns = rats(H, {0.0})
bases = []
for (a, b) in ms:
    for (c, d) in ns:
        h = max(abs(a), b) + max(abs(c), d)
        bases.append((h, f"{a}/{b}" if b > 1 else f"{a}", f"{c}/{d}" if d > 1 else f"{c}"))
bases.sort()
seen = set()
count = 0
for h, m, n in bases:
    if (m, n) in seen: continue
    seen.add((m, n))
    print(m, n)
    count += 1
    if count >= 6000: break
EOF
echo "$(wc -l < $OUT/bases.txt) bases queued" >&2
# skip the already-tested (3, 1/3)
grep -v "^3 1/3$" $OUT/bases.txt | \
xargs -P "$PAR" -L1 bash -c 'L="data/t88hunt/b_${0//\//o}_${1//\//o}.log"; timeout '"$TMO"' magma -b m:=$0 n:=$1 code/claude_88_worker.m > "$L" 2>&1'
echo "DRIVER_DONE" >&2
