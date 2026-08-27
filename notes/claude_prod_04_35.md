# Z/35 production campaign (top-10 #4): A_1(5) sweep to height 88 + Phi38 lane closure

Session: 2026-07-18 (resumed after API-limit abort; reused engine + completed H=16 run
from the aborted launch).  Strategy recap (3 lines):
(A) CRT sieve on the Elkies `A_1(5)` chart (`y^2+(Q'-xQ)y=Q^2 Q'`, universal 5-torsion),
    kill tables through `p<=61` in C, involution-quotiented, pushed from height 16 to 88;
(B) 3-adic filter sharpened via exact `genus2red` verdicts (measured, not deployed);
(C) the nondegenerate contact-7+contact-5 lane subsumed into one plane curve `Phi38` —
    irreducible, genus >= 7 (= 7 mod two primes), its 4 small rational points all degenerate.

## Task A: production sweep — clean to height 88, no candidate

Engine: `code/claude_prod35_sweep.c` (binary + all raw/derived outputs in scratchpad
`prod35/`).  Upgrades over the exploratory `c35_sweep.c`: tables `p<=43` (13 primes) or
`p<=61` (17, `nfp=17`), big-prime-first check order, degenerate planes `q0=0`, `q2=0`
excluded in-engine, involution `(q0,q1,q2)->(q2,1-q1,q0)` quotient (whole `q1`-slices
skipped, `i0<=i2` on the fixed slice `q1=1/2`), chunked by kept `q1`-slice.

Commands (scratchpad `prod35/`, ~40M triples/s/thread, nice -n 10):

```
./prod35 search 40 2 {0,1} 13   # 4.69e9  triples,  ~90 s/chunk
./prod35 search 64 2 {0,1} 17   # 8.00e10 triples,  ~17 min/chunk
./prod35 search 80 3 {0,1,2} 17 # 3.04e11 triples,  ~35 min/chunk
./prod35 search 88 3 {0,1,2} 17 # 5.31e11 triples,  ~15-55 min/chunk (load-dependent)
gp post{40,64,80,88}.gp         # exact 35 | #J(F_q) for good q <= 199 (prime-to-q part
                                # at q=5,7), degenerate triples (disc=0) dropped
```

Results (raw survivors -> degenerate + genuine -> exact kills -> final):

```
H=16 (prior session): 347  -> 344 degen + 3 genuine  -> kills 47:2 59:1          -> 0
H=40: 12633 -> 2101 degen + 10532 genuine -> 47:8591 53:1626 59:259 61:42 67:11 71:2 79:1 -> 0
H=64:  5877 -> 5442 degen + 435 genuine   -> 67:363 71:64 73:6 79:2             -> 0
H=80: 10002 -> 8516 degen + 1486 genuine  -> 67:1196 71:249 73:31 79:10         -> 0
H=88: 13327 -> 10279 degen + 3048 genuine -> 67:2494 71:464 73:69 79:19 83:2    -> 0
```

**No rational C35 point on the `A_1(5)` chart with coordinate height <= 88.**
(Box ~1.06e12 raw triples; quotient enumerated 9.2e11 across the height ladder.)
Kill decay is perfectly
geometric (~0.16-0.2/prime) at every height — no anomalous tail, i.e. no hint of a
near-miss locus.  The only positive-dimensional survivor families are genus-drop
degeneracies of the chart (e.g. the line `(q0,q1)=(-1/4,-1/2)`, `disc(f)` identically 0,
found by the sweep and discarded exactly).

## Task B: 3-adic sharpening — measured, verdict "not worth deploying in-engine"

Exact `genus2red(f,3)` verdicts on ALL 435 genuine H=64 survivors
(`code/claude_prod35_taskB_3adic.gp`):

```
class 1  C good reduction at 3, magic Weil poly (x^2+x+3)(x^2+3x+3): 0
class 2  C good reduction, wrong poly  (the only strictly killable class): 4   (0.9%)
class 3  J bad reduction at 3 (honest conservative pass):            422   (97%)
class 4  J good, C bad: compact type [I{0}-I{0}-n], J ~ E1 x E2:      13   (3%)
```

So the speced mask ("valuation patterns that force good reduction with wrong charpoly")
buys ~1%: the in-engine smooth-chart p=3 table already catches the bulk, and the
conservative pass is dominated by honest bad reduction.  Real (undeployed) levers:
elliptic-trace tests on class 4, Neron component-group divisibility on class 3.

Valuation-pattern classifier (125 patterns `v3(q0,q1,q2) in {-2..2}^3`, sampled):
only `v=(0,>=0,0)`-type patterns admit C-good reduction with the magic poly.
Structural find: **`v=(-1,*,-1)` forces compact-type J-good reduction (400/400 deep
samples, J ~ E1 x E2 at 3, never chart-good)** — this is the concrete mechanism behind
the dossier's "3-adically boundary but good" prediction; a C35 curve of this shape
needs trace pair (-1,-3) on the two elliptic pieces.  Survivor pattern census (435
genuine H=64): 53 at (0,0,0), 35 at (-1,0,0), 25 at (0,-1,0), long tail; 15 in the
`(-1,*,-1)` family.

## Task C: Phi38 — the nondegenerate contact lane closed (small height), genus >= 7

Scripts: `code/claude_prod35_phi38_derive.m`, `_analyze.m`, `_lift.m`.
Re-derived the residual system of `notes/agent_Z35_next_route.md` from scratch
(h = 1-(7/2)x+ax^2+bx^3, f=(h^2+(x-1)^7)/x^2, f - q^2 = (x-r)^5, d=c2-b, e=c2+b):
`Res_c1(N0,N1) = d^3 e^3 (d-e)^8 (d+e)^4 Phi38(d,e)` with Phi38 of degree 38 and
280 terms — exactly matching the notes.  New results:

1. **Phi38 is irreducible over Q**; its plane curve has **geometric genus >= 7
   over Q** (genus exactly 7 mod both 10007 and 32003; delta invariants only grow
   under specialization, so genus_Q >= 7 >= 2 rigorously — the exact value over Q
   is unconfirmed, see the correction addendum below) => finitely many rational
   points (Faltings).  The whole nondegenerate contact lane is this one curve.
2. `F_3`-points of the projective closure: exactly 6 — the four affine degenerate
   centers `(2,0),(0,1),(1,1),(2,2)` of the notes plus `(0:1:0),(1:0:0)` at infinity.
   Every rational point is 3-adically confined to those 6 disks.
3. Rational point search (all d with height <= 40, all rational e via `nfroots`):
   exactly 4 points, `(-1,-3), (1,-2), (2,-1), (3,1)` (pairs swapped by the q->-q
   symmetry `(d,e)->(-e,-d)`).  **All four lift to full rational contact data
   (rational common root c1 of N0,N1) but every lift has `disc(f)=0` AND `h(1)=0`**
   — singular quintic, degenerate 7-contact point: none gives a genus-2 curve.
   `(1,-2),(2,-1)` additionally have `r=1` (5-contact collides with 7-contact).
   Lift data recorded in `data/claude_prod_04_35_phi38.txt` with the polynomial.

**The nondegenerate contact-7+contact-5 lane is closed for all points with
height(d) <= 40, and is a curve of genus >= 7 (= 7 mod both test primes) with at
most finitely many rational points,
all confined to 6 mod-3 disks whose 4 affine centers are exactly the known degenerate
ones.**  Realistic continuations: two-cover/elliptic-quotient descent on Phi38, or
declare the lane exhausted and prioritize Route B (`A_1(7)` mirror chart).

## Verdict

No Z/35 realization found; frontier pushed from chart height 16 to **88** with clean
geometric kill statistics, and the independent contact lane reduced to a genus->=7
curve whose accessible rational points are all degenerate.  Z/35 on `A_1(5)` now
needs either much larger height (cost ~H^6: H=128 is ~17x H=80), Route B (`A_1(7)`
threefold, NDE 2003 — sieve for 5 | #J with ~0.23 pass rate, likely lower-height
exposure of `A_1(5,7)`), or descent machinery on `A_1(5,7)` itself.

## Resume state

- Engine + tables rebuild from source in ~1-2 min: `code/claude_prod35_sweep.c`
  (`gcc -O3 -march=native -o prod35 claude_prod35_sweep.c`).
- Scratchpad `prod35/`: `h{40,64,80}_c*.{txt,log}`, `h*_all.txt`, `post{40,64,80}.gp`
  (exact sieve), `taskB_3adic.gp`, `taskB2.gp`, `phi38*.{m,gp,out}`, `phi38_poly.txt`.
- To extend the sweep: `./prod35 search 96 4 {0,1,2,3} 17` (~7.0e11 quotient triples,
  ~1.2 h/chunk at 40M/s; MAXR=12000 supports H<=99; raise MAXR beyond that).
- Phi38 genus over Q: do NOT rerun `phi38_genusQ.m` (raw `Genus(ProjectiveClosure)`
  on the degree-38 plane model — ran 11h/70GB with no output; see the correction
  addendum below). If the exact Q-genus is ever needed, compute it via a lifted
  function-field model (function field mod a good prime, then lift the map data).
- No jackpot candidates arose; nothing pending exact verification.

## Correction (2026-07-18 evening): phi38_genusQ.m stray process
The optional exact-Q genus job reported "killed by PID after 35 min" was NOT dead: the kill
hit the magma wrapper, and the magma.exe child ran detached for 11h14m (to ~70 GB RSS, empty
output) before being terminated at 18:19 EDT. No result was produced and none is needed: the
genus facts stand as proven: genus >= 7 over Q rigorously, with equality mod both
test primes 10007 and 32003 (the exact Q-value remains UNCONFIRMED). If
the exact Q-genus is ever required for the Phi38 descent, compute it via a lifted function-
field model, not Genus(ProjectiveClosure) on the raw degree-38 plane model.
