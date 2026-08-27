---
name: validate-and-record-a-hit
description: The non-negotiable checklist between "the search printed HIT" and "we claim a result" - independent rebuild in a fresh session, exact torsion invariants, marked-class order, 2-rank, the simplicity certificate, the never-claim-a-family-from-its-seed rule (the spurious kappa=0 tale), plus the recording conventions (notes, logs, commits, ERRATUM). WHEN a search prints a hit, when writing up any result, when a claim needs retracting, or when tempted to announce a family. Trigger words - hit, validate, claim a result, reverify, record, ERRATUM, commit, family claim, false positive.
---

# Validate and record a hit

## When to use this

A search printed a `HIT` line. Between that moment and the words "we found a
geometrically simple `[...]` curve" sits this checklist. It exists because the
lab has real scars: nonsimple `[6,6]` curves that passed every cheap filter,
a "rank-1 family of `Z/24` curves" that evaporated on second look, and a
committed quotient-rank claim that needed a public ERRATUM. Run the checklist
mechanically; it is designed to be executable by a Sonnet-class model.

## The two-certificate rule (top line)

**A result = exact torsion + geometric simplicity certificate, both recorded.**
Anything less is a candidate. A `[6,6]` curve whose Jacobian splits is
worthless; a `[24]` curve is not a `[2,24]` curve. (Hub: `g2-torsion-lab`,
"two things must hold".)

## The checklist

1. **Rebuild independently, in a FRESH Magma session.** Reconstruct `f` from
   the logged **parameters** via the chart formulas
   (`named-charts-reference`), clear denominators (`IntModel`), and compare
   against the logged `f` string. Never validate the logged string alone —
   transcription and scaling bugs live exactly there (see the Mumford scaling
   rule, `magma-lab-conventions` §6). A fresh session also clears any state
   contamination from the long-running search process.

2. **Exact torsion.** `Invariants(TorsionSubgroup(J))` must equal the target
   invariant factors **exactly**. `[ 24 ]` is not `[ 2, 24 ]`: the extra
   `Z/2` must appear IN the invariant factors, not be inferred from 2-rank.
   The recorded convention is Magma invariant-factor lists, as in the real
   hit line (verbatim in `data/agent_a2_24_composite_h12_part1.log`; the
   curve is recorded in `notes/agent_a2_24_composite.md`):

   ```text
   HIT24 r=5 p=-5/2 t=-9/2 torsion=[ 24 ] 2rank=1 simple=true (q=17 chi=x^4 - 5*x^3 + 16*x^2 - 85*x + 289)
   ```

3. **Marked-class order, directly.** Re-verify the chart's marked divisor on
   the rebuilt model with the exact-order pattern (verbatim shape from
   `code/agent_a2_24_wsplit_3tors.m`):

   ```magma
   D8 := J![g8, v8];
   ok := (8*D8 eq O) and (4*D8 ne O);    // order exactly 8, not merely dividing 8
   ```

   Remember the `Lden` scaling of `v8` on the integral model
   (`magma-lab-conventions` §6).

4. **2-rank, if the target has a `[2,...]` prefix.** `TwoRank(fInt)` with the
   sextic/quintic accounting of `two-rank-and-factor-types`.

5. **Simplicity certificate.** Run `SimplicityCertificate(fInt)`
   (`simplicity-certificates`) and **record the witness `(pp, chi)`**. House
   standard for recorded results: confirm at MORE than one good prime — the
   two `Z/24` curves are recorded as, e.g.:

   ```text
   torsion [24]; SIMPLE (chi_17 = x^4-5x^3+16x^2-85x+289, irred; D4 at 8 primes)
   ```

   (`notes/agent_a2_24_composite.md` — "D4 at 8 primes" = the power-transform
   certificate succeeded at eight witness primes.)

6. **Sanity block.** `Discriminant(fInt) ne 0`; `Degree(fInt) in {5,6}`; f is
   NOT an even sextic (`Coefficient(f,i) eq 0` for `i in [1,3,5]` means
   bielliptic ⇒ split ⇒ reject regardless of everything above).

7. **Family claims need independent members — the kappa=0 tale.** If the hit
   suggests a FAMILY ("this component/curve parametrizes infinitely many
   ..."), generate **at least two additional members** and run steps 1–6 on
   each BEFORE claiming anything. The cautionary tale, in full
   (`notes/agent_a2_24_d0_derivation.md`): the `d=0` derivation produced an
   elliptic curve of **rank 1** passing through the known simple `Z/24` curve
   B — on paper, "infinitely many simple `Z/24`". Generating actual points
   killed it: **every member except B had torsion `[8]`**. The curve was the
   spurious `kappa = 0` branch (the perfect-square degeneration of the
   contact identity), and B lay on it by coincidence; the genuine cover,
   extracted properly, was genus-1 **rank-0**. A seed point lying on a
   component proves nothing about the component. No independent members, no
   family claim — ever.

8. **Record it** (next section), and if a **target fell**, update the
   state-of-play in the hub skill (`skills/g2-torsion-lab/SKILL.md`)
   and check `paper/` tracking.

## Recording conventions

Exemplars to imitate: `notes/agent_a2_24_composite.md` (the two `Z/24`
records), `notes/how_we_found_2220_examples.md` (the `[2,2,20]` record),
`notes/contact6_m36.md` (the `[6,6]` record).

- **The note is the project memory.** Write/update `notes/<route>.md` with:
  the chart and exact parameters; the integral `f`; the exact torsion
  invariants; the simplicity witness `(pp, chi)` (and how many primes);
  HOW it was found — script path, log path, run parameters, date; and the
  conclusion in one bold line. **Negative results are recorded deliberately**
  (counts, heights, filters) so nobody re-runs them — a killed scan with
  recorded counters is data (`running-torsion-searches` §9).
- **Keep the discovery log.** The greppable `HIT`/`TARGET` lines and the
  `SEARCH_DONE` counters of the discovering run stay in `data/` or
  `results/`; the note points at them.
- **Commit** code + note + log together. Message style: subject = the
  mathematical result; body = the key numbers (torsion, certificate prime,
  counts). Co-author line, verbatim:

  ```text
  Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
  ```

  **The user pushes** — agent `git push` fails on auth here; never attempt
  it, never treat an unpushed commit as published.
- **Retractions are explicit, never silent.** If a committed claim turns out
  wrong, add an **ERRATUM** section to the note stating exactly which claim
  is retracted and why, and keep the original text in place (struck or
  bracketed). The live example is `notes/agent_m18_416_R8_dA_quotients.md`,
  whose ERRATUM retracts invalid `d_A` "character quotient" conditions (the
  underlying error: a non-abelian `D4` tower — no `(Z/2)^3` character group);
  the same note's "Validation and a cautionary note" section separately
  retracts a spurious rank-0 quotient caused by the `Factorization` unit
  trap (`magma-lab-conventions` §1). Imitate its pattern: the ERRATUM names
  the false statement, the reason it is false, and what (if anything)
  survives.

## Pitfalls

- **Validating the logged `f` string instead of rebuilding from parameters.**
  The string can carry a scaling/transcription bug; the parameters + chart
  formulas are the ground truth. Rebuild, then compare.
- **Claiming `[2,24]` from `[24]` + 2-rank 2.** The extra factor must appear
  in `Invariants(TorsionSubgroup(J))`. 2-rank is necessary, not sufficient
  (`two-rank-and-factor-types`).
- **One-prime certificates.** A single witness prime is the *minimum*; the
  recorded standard is multiple ("D4 at 8 primes"). A reducible `chi` at one
  prime is not disproof either — see `simplicity-certificates` for both
  directions.
- **Family claims from one seed** — the `kappa=0` tale (step 7). This nearly
  produced a false "infinitely many simple Z/24" announcement. Two
  independent members, fully validated, or it is not a family.
- **Skipping the even-sextic / discriminant sanity block** because "the
  search already checked" — the search checked the *search's* model; you are
  validating the *rebuilt* one.
- **Silent retraction.** Deleting a wrong claim from a note destroys the
  project memory that prevents its rediscovery. ERRATUM, always.
- **Attempting `git push`** or treating local commits as published — the
  user pushes.

## See also

- `simplicity-certificates` — the certificate this checklist records (step 5).
- `two-rank-and-factor-types` — step 4, and the quintic accounting.
- `named-charts-reference` — the chart formulas for the independent rebuild
  (step 1).
- `magma-lab-conventions` — Mumford scaling (step 3), the Factorization unit
  trap behind the ERRATUM example, commit/push house rules.
- `running-torsion-searches` — the log-marker conventions the hit came from;
  recording killed scans.
- `component-boundary-analysis` — the proper extraction of a genuine cover
  (how the `kappa=0` branch was unmasked).
- `g2-torsion-lab` — hub; update its state-of-play when a target falls.
