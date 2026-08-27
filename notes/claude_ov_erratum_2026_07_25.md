# The `genus2red` conductor erratum: statement, validation, repository sweep

**Date:** 2026-07-25.  **Lane:** 7 (overnight campaign), three successive agents.
**Status:** closed.  Every affected number in the repository has been recomputed
and corrected; the sweep below lists what was checked and found *not* affected.

This note is the canonical reference for the erratum.  It is cited from
`code/claude_gen_torsion_table.py`, `notes/claude_generic_2214_found_2026_07_23.md`
and `notes/claude_generic_222_2214_plan_2026_07_23.md`.  The lane's other
results (certificates, the RM sixth curve, the two side rungs, the Richelot
audit) are in `notes/claude_ov_lane7_2026_07_25.md`; §3 there is the first
statement of the erratum and this note supersedes it.

---

## 1. The bug

PARI/GP's `genus2red` is documented as

> *"Determines the reduction at **p > 2** of the (proper, smooth) hyperelliptic
> curve `C/Q` of genus 2 ..."*

When it cannot determine the conductor exponent at 2 it emits the **sentinel
row `[2,-1]`** in the returned factorisation matrix and **omits the 2-part from
the returned conductor `N`**.  The returned `N` is then the ODD PART only, with
no warning on stdout.

Reproducers (`gp -q`):

```text
genus2red(x^5-x)     = [1,   Mat([2,-1]),      ...]   <- claims N = 1
genus2red(x^6-x)     = [625, [2,-1; 5,4],      ...]   <- N = 5^4, no 2-part
genus2red(x^5+x^2+x) = [229, [2,-1; 229,1],    ...]
```

**The omission is NOT unconditional**, which is why it went unnoticed: when
`genus2red` *can* resolve `p = 2` it returns the full conductor.  Control curve
`9450.b.2` has `v_2(N) = 1` and `genus2red` returns the full `9450`, with
reduction type `[I{2-0-0}] page 170` at 2 and **no sentinel row**.

**Operational test (use this one, not "is `N` odd?"):** look for a `[2,-1]` row
in the returned factorisation matrix.

`code/claude_ov_erratum_sentinel.gp` -> `results/claude_ov_erratum_sentinel.log`
applies the test to every `genus2red` conductor recorded in this repository:

```text
SENTSUMMARY  contact7 sentinel present in 11 of 11 ; blpC4 sentinel 1
```

so **all twelve recorded values were odd parts**.

## 2. The replacement method, and three validations of it

Replacement: `Conductor(ReducedMinimalWeierstrassModel(C))` in Magma.  Magma
prints `WARNING: Using Ogg's formula when v_2(D)>=12, no correctness guarantee`
whenever `v_2(disc_min) >= 12`, which is the regime of 5 of the 12 curves, so
the replacement itself had to be validated before anything was rewritten.

(Remark for the record: Ogg's formula `f = v(Delta) + 1 - m` is a **theorem** —
Ogg 1967 for elliptic curves, Saito 1988 in general, including residue
characteristic 2 — so Magma's disclaimer is about its own implementation, not
about the mathematics.  We did not audit that implementation; the guarantees
below are empirical.)

### 2.1 Control A — the paper's own rows, 14/14

`code/claude_ov_erratum_control.m` -> `results/claude_ov_erratum_control.log`.
Fourteen LMFDB genus-2 curves taken **verbatim from the equations displayed in
`paper/torsion_realizations.tex`**, so their published conductors are an
independent datum.

```text
CTRLSUMMARY  matched=14  mismatched=0  oddpart_matched=14  of 14
```

`v_2(N)` ranges over `0,1,2,3,4`.  Four of the fourteen (`10512.n.1`,
`19044.h.2`, `8136.c.1`, `28200.e.1`; `v_2(disc_min) = 12,12,12,16`) triggered
the Ogg warning and **all four were still exactly right**.
A gp-side companion (`code/claude_ov_erratum_control.gp` ->
`results/claude_ov_erratum_control_gp.log`) shows `genus2red` on the same
fourteen: it agrees on the odd part everywhere and drops the 2-part exactly
where the sentinel appears.

### 2.2 Control B — the big Ogg control, 1113/1113

`code/claude_ov_erratum_bigcontrol.m` -> `results/claude_ov_erratum_bigcontrol.log`
(12 forked children, ~1 minute).  Control set: **every** curve in the
production LMFDB `g2c_curves` table with `2^12 | abs_disc`, i.e. every published
curve that lands in the Ogg-warned regime — 1113 of them, pulled 2026-07-25 into
`data/claude_ov_erratum_g2c_v2disc12.csv` /
`data/claude_ov_erratum_g2c_control.m`.  Their conductors are
Booker–Sijsling–Sutherland–Voight data, independent of this project.

```text
BIGSUMMARY curves=1113 errors=0 matched=1113 mismatched=0
BIGSUMMARY Ogg-warned regime (v_2(disc_min)>=12): 1113 curves, 1113 matched, 0 mismatched
BIGSUMMARY v_2(published conductor) values seen: [2..19]
```

**1113/1113 in exactly the regime Magma disclaims**, with `v_2(N)` spanning
2 to 19.  This is the strongest single piece of evidence in the erratum.

### 2.3 Control C — the 2-part recomputed off the regular model, 12/12

`code/claude_ov_erratum_ss2.m` -> `results/claude_ov_erratum_ss2_all.log`,
`..._ss2_rerun04.log`, `..._ss2_ctrl.log`.  For each curve, build
`M := RegularModel(Cm, 2)`, force multiplicities with `IntersectionMatrix(M)`,
count the components `m` of the special fibre, and evaluate Ogg's formula
`f_2 = v_2(disc_min) + 1 - m` **independently of the `Conductor` call**.

```text
curve                            v_2(disc_min)   m    f_2 = v+1-m   Magma Conductor 2-part   fibre reduced   Neron comp. gp
control 12300.e.2 (pub v_2=2)          6          5        2          2 (published)              yes            [2,4]
blpC4 [2,22]                          28         27        2          2^2                        yes            [2,110]
contact7  1 (-10,-10/7,-1/2)          10          9        2          2^2                        yes            [2,14]
contact7  2 (-5,-15/8,-15/22)         14         13        2          2^2                        yes            [2,28]
contact7  3 (-3,-3/4,-3/5)  [RM]       7          6        2          2^2                        yes            [14]
contact7  4 (-15/8,-15/19,-1/2)       14         13        2          2^2                        yes            [2,28]
contact7  5 (-5/18,-10/49,4/17)       30          -        -          NOT DETERMINED (OOM)        -              -
contact7  6 (-4/9,-4/25,4/17)         30         26        5          2^5                        NO             [2,2,14]
contact7  7 (-511/61,-511/625,-1/2)    8          6        3          2^3                        yes            [14]
contact7  8 (-165/41,-33/16,-165/289) 21         20        2          2^2                        yes            [140]
contact7  9 (-164/297,-1/2,164/361)   26         25        2          2^2                        yes            [2,112]
contact7 10 (-17/50,-34/189,34/121)   14         11        4          2^4                        NO             [2,14]
contact7 11 (-1/2,-13/49,13/50)       14         13        2          2^2                        yes            [2,28]
```

**12 of 12 agree.**  Four of those twelve have `v_2(disc_min) < 12` (control,
curves 1, 3, 7) and therefore draw no Ogg warning, i.e. by Magma's own message
`Conductor` did *not* fall back on Ogg there — so for those four this compares
two different routes to the same exponent, and they agree 4/4.  (This is an
inference from the warning's wording, not from Magma's source, which was not
inspected.)  For the eight in the warned regime it re-derives the same
formula from the same regular model, which checks reproducibility of the
component count but is not logically independent.

The Néron component groups at 2 are a by-product and are recorded because
nothing else in the repo has them.

## 3. The corrected numbers

### 3.1 The eleven contact-7 `[2,2,14]` curves

`(s,t,u)` are the three-root parameters of the contact-7 chart; odd parts are
`genus2red`'s, 2-parts are Magma's.

```text
idx (s,t,u)                       N
 1  (-10,-10/7,-1/2)              2^2 * 9550095752925                        = 38200383011700   <- displayed in the paper
 2  (-5,-15/8,-15/22)             2^2 * 17934836205825                       = 71739344823300
 3  (-3,-3/4,-3/5)   [RM]         2^2 * 38025                                = 152100 = 390^2
 4  (-15/8,-15/19,-1/2)           2^2 * 14439206077220925                    = 57756824308883700
 5  (-5/18,-10/49,4/17)           2^? * 37247776656021702225                   2-PART NOT DETERMINED
 6  (-4/9,-4/25,4/17)             2^5 * 21161909161575                       = 677181093170400
 7  (-511/61,-511/625,-1/2)       2^3 * 1363558168459661985137828367867185325
 8  (-165/41,-33/16,-165/289)     2^2 * 11352674691759169446159049425
 9  (-164/297,-1/2,164/361)       2^2 * 10200972432773189299959691439475
10  (-17/50,-34/189,34/121)       2^4 * 3069529451653112094922275
11  (-1/2,-13/49,13/50)           2^2 * 69195300331841508825
```

Independent corroboration at curve 3: it is certified RM by `Q(sqrt 2)` and
isogenous over `Q` to `A_f` for the level-390 weight-2 newform with Hecke field
`Q(sqrt2)` (L-polynomial identity at 135/135 good primes `p < 800`,
`notes/claude_ov_lane7_2026_07_25.md` §4).  A GL2-type surface has
`cond = level^2`, so `N = 390^2 = 152100` — which forces `v_2(N) = 2` and agrees
with the table.  The newform has `a_2 = 1 != 0` (multiplicative reduction at 2),
also consistent.  The `genus2red` value `38025 = 195^2` would have been the
conductor of a level-195 form, and level 195 has **no** degree-2 Hecke field at
all, so the odd part was not merely incomplete but pointed at an impossible
level.

### 3.2 The BLP C4-corrected `[2,22]` curve

`code/claude_ov_erratum_blp22.m` -> `results/claude_ov_erratum_blp22.log`.
`y^2 = (x-9)(x+21)(x^2-80x+439)(x^2+50x+109)`, recorded in
`notes/claude_generic_222_2214_plan_2026_07_23.md` as "conductor odd part
`645^2`":

```text
N = 1664100 = 2^2 * 3^2 * 5^2 * 43^2 = 1290^2      (odd part 416025 = 645^2, confirmed)
torsion re-verified [2,22]
```

Independent corroboration: the curve is RM(`sqrt5`), hence of GL2-type, hence of
**square** conductor.  `2^k * 645^2` is a square only for even `k`, and
`v_2(disc_min) = 28` puts it in the Ogg-warned regime — yet the square-conductor
constraint plus Control C's `f_2 = 2` both give `2^2`.

### 3.3 The one number that is still missing

Curve 5, `(s,t,u) = (-5/18,-10/49,4/17)`, `v_2(disc_min) = 30`.

* `Conductor` OOMs: 45.8 GB (`results/claude_ov_erratum_cond5.log`), 40 GB
  (`results/claude_ov_lane7_condaudit_5.log`), 3 GB
  (`results/claude_ov_lane7_conductor2.log`).
* `RegularModel(Cm,2)` OOMs at a 24 GB cap
  (`results/claude_ov_erratum_ss2_all.log`: "Current total memory usage:
  21115.1MB, failed memory request: 3289.3MB").
* A retry at a 90 GB cap (`code/claude_magma_slot.sh -b IDX:=5 MemGB:=90
  code/claude_ov_erratum_ss2.m`) ran 9 min 22 s, reached **73.8 GB resident**
  with the RSS flat and the process CPU-bound (state `R`), and was **killed by
  PID** when claude-box hit load 42 with nine other lanes competing.  It did not
  OOM and it did not finish; `results/claude_ov_erratum_ss2_idx5_big.log` is
  empty because Magma's `-b` buffers stdout.  Killing it was a resource
  decision, not a mathematical one: the answer changes no claim (see below), and
  74 GB of 125 GB under load 42 was not defensible.
  For scale: curve 6 has the *same* `v_2(disc_min) = 30` and its
  `RegularModel(.,2)` finished in 0.05 s, so the cost is not explained by `v_2`
  alone — something about curve 5's 2-adic geometry specifically is expensive,
  and that is the same wall `Conductor` hits.

**Only the odd part `37247776656021702225` is known.**  This affects nothing:
the odd part alone is `3.7e19`, six orders of magnitude above the displayed
witness's `3.8e13`, so no ordering and no displayed value depends on it.

Resume command (needs a box with more RAM than claude-box, or a Liu-style
`p = 2` reduction implementation rather than a regular model):

```bash
cd /home/claude/torsion_jac
code/claude_magma_slot.sh -b IDX:=5 MemGB:=<N> code/claude_ov_erratum_ss2.m \
    > results/claude_ov_erratum_ss2_idx5_big.log 2>&1
code/claude_magma_slot.sh -b MemGB:=<N> code/claude_ov_erratum_cond5.m \
    > results/claude_ov_erratum_cond5.log 2>&1
```

## 4. The load-bearing claim, and why it survives

`paper/torsion_realizations.tex` says the displayed `[2,2,14]` curve is *"the
smallest of ten generic witnesses"*.  Conductors are compared, so the 2-parts
matter.  Only **one** comparison is at risk: the odd parts of curves 1 and 2 are
in ratio `17934836205825 / 9550095752925 = 1.878 < 2`, so

```text
f_2(curve 2) = 2  ->  N = 71739344823300 > 38200383011700   claim holds
f_2(curve 2) = 1  ->  N = 35869672411650 < 38200383011700   claim FLIPS
```

The full sorted list of the ten generic witnesses (curve 3 is the RM one and is
not among them):

```text
 1 (-10,-10/7,-1/2)          N = 2^2*9550095752925                        =                 38200383011700
 2 (-5,-15/8,-15/22)         N = 2^2*17934836205825                       =                 71739344823300
 6 (-4/9,-4/25,4/17)         N = 2^5*21161909161575                       =                677181093170400
 4 (-15/8,-15/19,-1/2)       N = 2^2*14439206077220925                    =              57756824308883700
 5 (-5/18,-10/49,4/17)       N >= 37247776656021702225                       (odd part only, 2-part unknown)
11 (-1/2,-13/49,13/50)       N = 2^2*69195300331841508825                 =          276781201327366035300
10 (-17/50,-34/189,34/121)   N = 2^4*3069529451653112094922275            =        49112471226449793518756400
 8 (-165/41,-33/16,-165/289) N = 2^2*11352674691759169446159049425        =     45410698767036677784636197700
 9 (-164/297,-1/2,164/361)   N = 2^2*10200972432773189299959691439475     =  40803889731092757199838765757900
 7 (-511/61,-511/625,-1/2)   N = 2^3*1363558168459661985137828367867185325
```

Every witness below curve 2 is already larger in its odd part alone, so no
2-power (bounded by `2^{v_2(disc_min)}`, and in fact `<= 2^5` here) can reorder
them; curve 5's unknown 2-part cannot move it either.

Curve 2 is measured `f_2 = 2` by `Conductor` and by Control C
(`v_2(disc_min) = 14`, `m = 13`).  Both are Ogg in this regime, so the guarantee
is the 1113/1113 + 14/14 empirical record, not a proof.  **The claim stands.**

## 5. A shortcut that does NOT work (recorded so nobody retries it)

The obvious way to get `f_2 >= 2` with no appeal to Ogg is: `f = 2u + t + delta`,
so `f = 1` forces `(u,t,delta) = (0,1,0)`, i.e. `J` semistable at 2, i.e. the
special fibre of the minimal regular model is **reduced**.  Find a multiplicity
`> 1` and you have proved `f_2 >= 2` outright.  That much is sound, and it does
fire at curves 6 and 10 (fibre NOT reduced, so `f_2 >= 2` is proved there
without Ogg) — but those two were never at risk.

The converse step — *reduced fibre => semistable => `f_2 = t = b_1(dual graph)`*
— is **WRONG**, and curve 7 is the counterexample:

```text
contact7 7 : all 6 multiplicities = 1 (reduced), dual graph v=6, e=7, b_1 = 2,
             but f_2 = 3  (v_2(disc_min) = 8 < 12, so Magma used its NON-Ogg path;
             the Ogg count off the same regular model also gives 8+1-6 = 3)
```

`f_2 = 3 > 2 = g` is impossible for a semistable abelian surface, so the fibre is
reduced but **not** semistable — multiplicity 1 in Magma's `RegularModel` does
not certify smooth, geometrically irreducible components with only nodal
crossings.  At the other nine reduced curves `b_1` happens to equal `f_2 = 2`,
which is why the shortcut looks convincing until curve 7.

**Conclusion: the `b_1` column is a cross-check, never a proof.**  Curve 2's
`f_2 = 2` rests on Ogg + the empirical controls, and that is the honest status.

## 6. Repository-wide sweep

Method: (i) `grep -rl genus2red` over `notes/ paper/ data/ code/ results/`;
(ii) grep for each of the twelve recorded odd-part values; (iii) grep
`-i conductor` over every `.md` in `notes/`, every `.tex` in `paper/`, and every
`.log` in `results/`, and classify each hit by provenance.

### 6.1 Affected — all corrected

| file | what was wrong | now |
|---|---|---|
| `paper/torsion_realizations.tex` (endomorphism-refinement paragraph) | `[2,2,14]` conductor printed as `3^2 5^2 7^2 13^2*19*23*37*317 = 9550095752925` | `2^2 3^2 5^2 7^2 13^2*19*23*37*317 = 38200383011700`; regenerated + PDF rebuilt (commit `7e08b09`) |
| `code/claude_gen_torsion_table.py` | hard-coded that odd value | corrected, plus a `CONDUCTOR PROVENANCE` comment naming the sentinel test (commit `65125ab`/`7e08b09`) |
| `notes/claude_generic_2214_found_2026_07_23.md` | six `[2,2,14]` conductors were odd parts; **and** the derived claim *"all conductors are odd non-squares (no GL2-type smell)"* | erratum box + corrected table inserted; the "odd non-squares" inference struck |
| `notes/claude_generic_222_2214_plan_2026_07_23.md` | BLP C4 `[2,22]` "conductor odd part `645^2`" | full `N = 1664100 = 1290^2` recorded inline |
| `data/claude_endz_certificates.txt` | `[2,2,14]` headers carried odd parts | all six carry `N` with its 2-part, each tagged with `v_2(disc_min)` and (this session) the `RegularModel` cross-check |
| `notes/claude_ov_lane7_2026_07_25.md` §3 | first statement of the erratum; §8 bullet 3 wrongly said curve 5's 2-part had been obtained | §9(a) of that note already corrects it; §3.3 of this note is the current status |

### 6.2 Checked and NOT affected

* **`paper/torsion_sources.tex`** — contains **no conductor values at all**.  It
  is the attribution table (year, claim, simplicity verdict, citation); the only
  `grep -i cond` hits are the words "second" and "Comp.".  Nothing to fix.
* **`notes/claude_torsion_refs_dossier.md`** — exactly one occurrence,
  *"all 35 minimal-conductor simple representatives already have
  `geom_end_alg='Q'`"*, a statement about the LMFDB production census.  No
  computed conductor anywhere in the file.  Nothing to fix.
* **`notes/claude_prod_04_35.md` + `code/claude_prod35_taskB_3adic.gp`** — these
  *do* call `genus2red`, but as `genus2red(fz, 3)`, purely for the **local
  reduction type at `p = 3`**, never for a conductor.  `p = 3 > 2` is inside the
  documented range.  Unaffected.
* **Every other `conductor` hit in `notes/`** (about 30 lines across
  `claude_top10_*`, `claude_prod_*`, `claude_ari_surface_22212.md`,
  `claude_ov_lane1_*`, `claude_table1_errata.md`, `order222_from_order11.md`,
  `claude_torsion_sources_audit.md`, `contact30_*`) is either
  (a) an **elliptic-curve** conductor from Magma/PARI `ell*` — 24, 30, 54, 66,
  92, 288, 4290, 18744222, 70702170, ... — for which `genus2red` is not used at
  all, or (b) an **LMFDB label component** (published data), or (c) a
  qualitative remark ("small conductor", "square conductor").  None came from
  `genus2red`.
* **Every conductor in the body of `paper/torsion_realizations.tex`** — all are
  the first component of an LMFDB label, i.e. published data, per the
  `CONDUCTOR PROVENANCE` comment in the generator.  The `[2,2,14]` preamble
  value was the *only* computed conductor in the paper.
* **`results/*.log`** — the only conductor prints are elliptic (92, 30) from
  `claude_222_*` and `claude_ari_surface_verify`.
* **`data/` numeric hits** on `38025` / `416025` are coordinate tuples in
  `tor2244*`, `claude_prod_06_224_*`, `claude_prod_07_2248_*`,
  `claude_sib_t5_nearmisses.txt`, `m2248_*` — not conductors.
* **`notes/contact30_c3root_genus6_frobenius_2026_07_11.md`** explicitly says no
  rigorous conductor was computed (genus 6 anyway).  Unaffected.

## 7. Which published claims were affected, and which were not

**Affected (one number and one inference):**

1. The `[2,2,14]` conductor displayed in the paper:
   `9550095752925 -> 38200383011700`.  Corrected, regenerated, PDF rebuilt.
   The *curve* is unchanged, the *row* is unchanged, and it is still the
   smallest of the ten generic witnesses (§4).
2. The internal inference in `notes/claude_generic_2214_found_2026_07_23.md`
   that the six `[2,2,14]` conductors were *"odd non-squares, no GL2-type
   smell"*.  This was **entirely an artefact of the missing 2-part**: all six
   are even, and the RM one is `152100 = 390^2`, a perfect square — precisely
   the GL2-type signature the note declared absent.  Symmetrically, the BLP
   `[2,22]` curve's `645^2` became `1290^2`.  This is the only place where the
   bug caused a wrong *scientific* statement rather than a wrong digit.

**Not affected:**

* **No torsion claim.**  Every `TorsionSubgroup` in the project is exact Magma
  on an integral model and never touched `genus2red`.
* **No simplicity / `End = Z` claim.**  Those are Frobenius root-power +
  linear-disjointness certificates at good odd primes; the conductor plays no
  role.  All ten `[2,2,14]` certificates and the RM certification stand as
  written.
* **No count of realized groups**, no "new to this paper" attribution, no table
  row: the erratum changes displayed conductor values only.
* **The paper's table body**, which is LMFDB data throughout.
* **`paper/torsion_sources.tex`** in its entirety.
* **The order-35 production line** (`genus2red` at `p = 3` only).

Nothing is retracted.  No curve, certificate or realization claim changes.

## 8. Reproduction

```bash
cd /home/claude/torsion_jac
gp -q code/claude_ov_erratum_sentinel.gp   > results/claude_ov_erratum_sentinel.log
gp -q code/claude_ov_erratum_control.gp    > results/claude_ov_erratum_control_gp.log
code/claude_magma_slot.sh -b MemGB:=8 code/claude_ov_erratum_control.m \
    > results/claude_ov_erratum_control.log 2>&1
code/claude_magma_slot.sh -b NCH:=12 MemGB:=4 code/claude_ov_erratum_bigcontrol.m \
    > results/claude_ov_erratum_bigcontrol.log 2>&1     # ~1 min, 1113 curves
code/claude_magma_slot.sh -b MemGB:=8 code/claude_ov_erratum_blp22.m \
    > results/claude_ov_erratum_blp22.log 2>&1
for k in 0 1 2 3 4 6 7 8 9 10 11; do \
  code/claude_magma_slot.sh -b IDX:=$k MemGB:=24 code/claude_ov_erratum_ss2.m; done \
    > results/claude_ov_erratum_ss2_all.log 2>&1        # IDX:=-1 is the control;
                                                        # Magma cannot parse a negative
                                                        # on the command line, run it as
                                                        # a separate invocation
```

Wall clock: everything except curve 5 finishes in under three minutes on
claude-box.  `RegularModel(.,2)` costs 0.05–0.1 s per curve at
`v_2(disc_min) <= 26` and blows past 90 GB at `v_2(disc_min) = 30`.

## 9. What was NOT checked

* **Magma's `Conductor` implementation** was not audited; the 2-parts rest on
  1113 + 14 published-conductor controls and the Ogg-vs-non-Ogg agreement of
  §2.3, not on a proof.
* **Curve 5's 2-part** is unknown (§3.3), and so is its exact conductor.  Four
  attempts (Conductor at 3 / 40 / 45.8 GB, RegularModel at 24 and 90 GB) all
  ended in an OOM or a deliberate kill.  Nothing in the paper, the certificates
  file or any ordering depends on it.
* **`genus2red`'s behaviour at odd primes** was assumed correct throughout: the
  odd parts were never independently recomputed (Magma's odd part agreed with
  gp's `N` at all eleven contact-7 curves, per
  `notes/claude_ov_lane7_2026_07_25.md` §3.1, which is the only check there is).
* **No other PARI intrinsic** was swept for analogous `p = 2` restrictions.
* **The `[2,-1]` sentinel test was applied only to conductors recorded in this
  repository.**  If `genus2red` was used interactively in a session transcript
  and a number copied into prose that no longer exists as a file, this sweep
  would not see it.
* **`results/` logs were classified, not recomputed** — a conductor buried in a
  log but never promoted to a note or the paper is out of scope by design.
