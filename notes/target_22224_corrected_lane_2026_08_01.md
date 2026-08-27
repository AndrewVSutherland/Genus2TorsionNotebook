# [2,2,2,24] corrected-mask lane: integer sweep to d=1e18, family lane to h4000

Date: 2026-08-01

## Objective and design

Resume the corrected `A(2,2,2,8)+3` program of
`notes/target_22224_full_family_halving_2026_07_18.md` §6 ("generate rational
points of the full four-radicand cover subject to the corrected local
conditions").  Post-erratum design (the old `Aaux` contact masks and the
mod-169 Type I/II derivation must not be used):

- **bad branch reduction at both 11 and 13** (corrected smooth target masks
  are empty at both primes: `results/..._corrected_complete_p11.tsv` and
  `_p13.tsv` are header-only);
- at `17,19,23,29,31`: boundary **or** the corrected complete smooth
  3-contact mask (`results/target_22224_full_family_halving_corrected_complete_p*.tsv`).

Tooling: the resumable generator committed (undocumented) in `c2a5348`:
`code/target_22224_a2228_deep13_generator.cpp` (integer/fiber lane) and
`.py` (rational-family lane, families filip1/filip2/adam).  The committed
ELF binary was built on another machine and does not run on the target host; rebuild
with

```sh
g++ -O3 -march=native -std=c++17 -o deep13_gen \
    code/target_22224_a2228_deep13_generator.cpp
```

## Pipeline validation (mandatory known-positive test)

The current generator (correctly, per
`notes/target_22224_product_surface_bielliptic_2026_07_18.md`) excludes
rectangle tuples (`ad=bc` up to permutation) and non-smooth tuples.  This is
why re-running the old "integer control" finds none of the 19 rows of
`results/target_22224_a2228_deep13_generator_integer_control.tsv` — those are
all rectangle/degenerate tuples from the pre-exclusion code version.  The
committed control file is therefore obsolete as a reference.

New validation (2026-08-01): with **permissive** masks at 17..31 (all
residues allowed; boundary-only at 11,13 unchanged) and `d in [1,10^4]`, the
generator finds 133 full-cover points, all `old=1` (bank rows), including
both planted deep near-misses at their predicted fiber coordinates:

```text
FULL_COVER_HIT fiber=123 k=7  tuple=[121,1919,3211,4949]   (d = 101*7^2)
FULL_COVER_HIT fiber=198 k=13 tuple=[1369,1711,2349,9971]
```

With the **real** corrected masks the same range yields 0 candidates — the
near-misses are correctly rejected (they are smooth non-3-target at 23 and
31, matching their recorded first 3-killing primes).  Funnel validated
end-to-end; saved as `results/target_22224_a2228_deep13_generator_validation_permissive.tsv`.

## Integer lane: complete negative through d = 10^18

All 221 bank fibers (three-coordinate fibers of the 619 primitive
`tor2228.txt` tuples), all 4 sign patterns, `d = sign * sf(|abc|) * k^2`:

| range | k tested | mask survivors | radicand-tested | new cover points |
|---|---:|---:|---:|---:|
| `d <= 10^15` (single job) | 9,046,413,020 | 238,255,274 | 238,255,166 | **0** |
| `d in (10^15, 10^18]` (4 shards) | 277,026,300,016 | 7,296,061,297 | 7,296,061,297 | **0** |

Logs: `results/target_22224_a2228_deep13_generator_integer_full_d1e15.tsv`
(header-only; run summary in this note),
`results/target_22224_a2228_deep13_generator_integer_full_d1e18_part{0..3}.{tsv,log}`.
Wall time: ~110 s for d<=1e15; 12–19 min per shard for the 1e18 stretch.

**Conclusion.**  No new rational point of the full four-radicand cover exists
with three coordinates from the bank fibers and fourth coordinate of absolute
value up to 10^18 satisfying the corrected local conditions.  Combined with
the earlier bank scan, the fiber-extension direction is now decisively
closed far beyond the old B=16384 bank bound; new cover points must move in
at least two coordinates simultaneously (off-fiber).

## Family lane: denominators 2001..4000

`code/target_22224_a2228_deep13_generator.py`, families filip1/filip2/adam,
`t=m/n` with `den in [2001,4000]`, `|m| <= 4000`: 48,006,000 pairs tested,
4,069 cover points, **1,747** survivors of the 18-prime local profile
(primes 5..71), max tuple height ~6.4e21.  Candidates:
`results/target_22224_a2228_deep13_generator_candidates_h2001_4000.tsv`;
log `results/..._family_h4000.log`.

Survivors passed to the exact verifier
(`code/target_22224_a2228_deep13_generator_verify.m`, reduction 3-primary
scan at primes 73..499, then exact `TorsionSubgroup` on any survivor):
log `results/target_22224_a2228_deep13_generator_verify_h4000.log`,
hits file `results/target_22224_a2228_deep13_generator_hits_h4000.tsv`.
RESULT (2026-08-01): `A2228_DEEP13_VERIFY_DONE rows 1747
reduction_survivors 0 exact_hits 0` — every candidate is killed by a good
prime with `3 nmid #J(F_p)`, elimination histogram
`73:1002, 79:401, 83:157, 89:70, 97:60, 101:18, 103:11, 107:4, 109:3,
113:6, 127:2, 131:6, 139:3, 149:2, 157:1, 163:1`.  Together with the h2000
block (493 candidates, same outcome) the three rational families are
3-primary-dead through parameter denominator 4000.

## Status / next steps

1. If the h4000 verify is negative (expected base rate), the family lane can
   be resumed cheaply (`--den-start 4001`, checkpoint
   `results/..._checkpoint_h4000.json`); each doubling of the height box
   quadruples cost (~5 min generator + Magma verify per block).
2. The integer lane should NOT be extended further in d alone (it is a
   union of thin quadratic sections of genus >= 2 curves; 10^18 is already
   deep in the tail).  The genuinely new directions, per the July notes:
   two-coordinate fiber moves (pair fibers P1..P7 of
   `notes/target_22224_offrectangle_new_curves_2026_07_18.md` with the
   Mordell-Weil residue sieve at each family's zero-target prime), and the
   `--surface` / `--rank8-charts` / `--deep-p5-charts` lanes of the
   generator, which tie into the direct-contact deep-13 tangent charts.
3. A corrected re-derivation of the 13-adic boundary strata (the analogue of
   the invalidated Type I/II classification, now with the fixed `Aaux`
   formula `P = 4*M*e1 + 12*(U^2+v^2) - (M+3*U)^2`) is still missing and
   would sharpen the masks from "boundary at 13" to specific 169-charts.

## Reproduction

```sh
cd ~/torsion_jac
g++ -O3 -march=native -std=c++17 -o /tmp/deep13_gen code/target_22224_a2228_deep13_generator.cpp
/tmp/deep13_gen --d-min 1 --d-max 1000000000000000 \
  --output results/target_22224_a2228_deep13_generator_integer_full_d1e15.tsv
# shards: --d-min 1000000000000001 --d-max 1000000000000000000 \
#   --fiber-start <lo> --fiber-stop <hi> --output ..._d1e18_part<k>.tsv
python3 code/target_22224_a2228_deep13_generator.py --height 4000 \
  --den-start 2001 --den-stop 4000 \
  --output results/target_22224_a2228_deep13_generator_candidates_h2001_4000.tsv \
  --checkpoint results/target_22224_a2228_deep13_generator_checkpoint_h4000.json
magma -b input_file:=results/target_22224_a2228_deep13_generator_candidates_h2001_4000.tsv \
  log_file:=results/target_22224_a2228_deep13_generator_verify_h4000.log \
  hits_file:=results/target_22224_a2228_deep13_generator_hits_h4000.tsv \
  code/target_22224_a2228_deep13_generator_verify.m
```
