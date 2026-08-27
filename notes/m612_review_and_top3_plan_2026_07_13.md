# `[6,12]` review and top-three execution plan

Date: 2026-07-13

Target: a geometrically simple genus-2 Jacobian over `Q` with exact rational
torsion `[6,12]`.

## Bottom line

No geometrically simple `[6,12]` example has yet been found.  The strongest
new structural advance is the contact-6 `P8` lane:

1. the missing independent rational `3`-class is represented by an explicit
   connected degree-`24` signed cover over the `e`-line;
2. modular reconstruction produces a compelling genus-4 quotient candidate
   `E8` covering a genus-2 curve `E4`; the reconstructed models pass exact
   internal identities, but an exact characteristic-zero map from the
   original signed cover has not yet been committed;
3. `E8` has two exact rational boundary places over `e=0`;
4. the Prym of `E8 -> E4` is represented by the genus-2 Jacobian

   ```text
   D: y^2 = -3*x^6 + 24*x^3 - 75,
   ```

   with reported exact rank `1`.

Once the quotient map and bigonal correspondence are certified exactly, the
global gatekeeper is in Prym-Chabauty range: `rank Prym = 1 < 2`.  This is a
highly promising reduction of the problem, but the current Coleman sweep is
not yet a proof and should not be continued unchanged.

## Review findings that affect the next run

### Critical precision defect in the current Coleman server

The pairings printed in `results/contact6_m612_E8_coleman3.log` are only known
to `O(37^5)`.  After their common factor of `37` is removed, the coefficients
of the annihilating differential are known only to `O(37^4)`.

In `code/contact6_m612_E8_coleman_server.m`, those coefficients are then
coerced into a parent field of precision `20`.  The Coleman library's
`zeros_on_disk` routine uses the precision of the *parent*, not the absolute
precision of the individual coefficients, and converts the resulting
coefficients to rational representatives modulo `37^20`.  Unknown digits are
therefore effectively filled in.  The printed zeros with `O(37^19)` are not
certified zeros of the true annihilating differential, and the partial sweep
cannot support a completeness claim.

There are three related cautions:

- the committed `p=37` log stops after disk `15/42` and has no total-zero or
  reconstruction section;
- the normalization uses `37^mv` even when the script is run with another
  prime such as `79`; this must be `p^mv`;
- rational reconstruction is a candidate detector, not a proof that an
  unreconstructed `p`-adic zero is irrational or that no large rational point
  exists.

The post-processing also discards the `inf` flag of each Tuitman point.  In an
infinity disk the stored coordinate is `1/e`, not `e`; such a zero is
currently tested on the wrong exact fiber.  The corrected driver must retain
the flag and invert the coordinate when appropriate.

The Prym rank-`1` output is recorded in a log, but the short Magma source that
produces the minimized model and `RankBounds = 1..1` was not committed with
that result.  It should be reconstructed and committed before the rank is
used as a publication-grade dependency.

Two theorem-level links also need promotion from evidence to exact
certificates:

- the current `E4/E8` reconstruction proves identities about the reconstructed
  polynomials, but not an exact quotient map from the original degree-`24`
  signed function field;
- the bigonal scripts prove useful norm and conic identities, while the Prym
  identification is checked only through several first Frobenius traces and
  a later `p`-adic transport.  The branch hypotheses and exact correspondence
  over `Q` should be recorded explicitly.

The Coleman checkout should also be pinned to the exact tested revision,
rather than cloned from a moving default branch.

The irreducible local quartic argument does rule out an elliptic quotient
defined over `Q(zeta_3)`.  It does not, as currently written, rule out every
possible splitting after an arbitrary extension.  That stronger sentence
should be removed unless absolute simplicity is proved separately.

### Results that look solid and useful

- The exact `P8` halving family and its `[2,12]` near miss.  The torsion and
  `R3` half are verified, but the cited polynomial
  `529*T^4+22*T^2+1` is even and therefore does **not** certify geometric
  simplicity: squaring a Frobenius root drops its degree to at most `2`.
  A different good prime/root-power certificate is required before the near
  miss may be called geometrically simple.
- The degree-`12` orthogonal support cover and connected degree-`24` signed
  lift, including exact recovery maps and bounded height sieves.
- The internally exact quotient-candidate models `E4` and `E8`, their genera,
  and the two exact boundary places of the reconstructed `E8` model.
- The bigonal formula for the Prym candidate, supported by exact identities
  and Frobenius-trace checks.  For a final proof, the correspondence and its
  divisor action should be written as an exact certificate rather than
  described only through trace agreement.
- The optimized HLP direction

  ```text
  G_A = 2+x-x^2+x^3+x^4+x^5+x^6
  ```

  and its exact finite-Jacobian masks.  Its affine finite-residue pass-density
  is ten times that of `1+x`; after the projective infinity class retained by
  a height sieve is included, the corresponding product ratio is `9/2`.
- The subsequent HLP geometry: the marked order-`3` branches lie on one
  irreducible degree-`40`, genus-`51` cover.  This makes the one-parameter
  `G_A` line a poor parametrization target.
- The corrected contact-6/TB/Richelot searches.  Their negative conclusions
  are bounded and correctly retain projective, saturation, and bad-reduction
  boundaries; they are not global nonexistence theorems.

### Additional contact-6 items to repair before reusing negative claims

- Two finite Richelot mask routines select a dual twist by matching only
  `#J(F_p)` and take the first sign that matches.  Equal cardinalities do not
  identify the correct twist: at `p=5`, `(a,b)=(3,3)` gives order `36` for the
  source and for both dual signs.  Since order `36` cannot contain `[6,12]`,
  the aggregate `p=5,7` obstruction survives; the class-by-class twist
  attribution and larger-prime counts should nevertheless be rerun using the
  full local Euler factor, testing both signs when necessary.
- The height-`30` TB solver receives saturation/`Variety` success flags but
  does not assert or report them.  A reproducible rerun must count every
  failure before calling all `81` fibers exact solves.
- The corrected height-`10` dual split-core driver omits the exceptional
  recovery chart `DR=NR=0` on `v=1`.  An independent elimination found no
  square-`M` rational candidate in the omitted height box, so the numerical
  conclusion appears unchanged; that check still needs to be committed.

### HLP reproducibility items before Attack 3

- The saved exhaustive-direction artifact begins with only the final `2000`
  records and omits the command-line bounds.  Save a transcript proving the
  actual `35136 -> 2000` staged run.
- Add the seed automorphism-group calculation showing that the three visible
  nonhyperelliptic involutions exhaust the Humbert-4 branches.
- Save the reproduced modulo-`101` degree-`40` irreducibility/genus-`51`
  transcript and fix the harmless unassigned `output_file` reference in the
  aggregate Magma runner.

## Attack 1: repair and complete Prym-Chabauty on `E8`

This has the highest leverage: after the quotient chain is certified, it can
either produce a nonboundary candidate fiber or globally close the most
promising contact-6 lane.

### Work

1. Certify the exact characteristic-zero quotient map from the degree-`24`
   signed cover to `E8`, not only the internally exact CRT-reconstructed
   minimal polynomials.
2. Add a self-contained verifier for
   `D: y^2=-3*x^6+24*x^3-75`, its rank bounds, the explicit infinite-order
   Mumford class, and the model transformations used by the bigonal map.
3. Make the bigonal divisor transport reproducible with exact checks of the
   two finite and two boundary contributions.
4. Pin the Coleman library revision and produce a precision-safe driver:
   - normalize by `p^mv`, not `37^mv`;
   - set the working precision from the minimum **absolute** precision of the
     annihilator coefficients and the Coleman integrals;
   - never coerce an `O(p^r)` coefficient into a computation at precision
     greater than `r`;
   - increase the initial Coleman precision until the normalized
     annihilator is known to a useful target precision (initial target:
     at least 12 reliable `p`-adic digits, adjusted after Newton-polygon
     diagnostics).
   - retain each zero's finite/infinity flag and use `e=x` or `e=1/x`
     accordingly.
5. Rerun a small set of disks at two increasing precisions.  Require stable
   root counts and compatible roots before launching all disks.
6. Complete all residue disks at `p=37`; repeat at `p=79` if it supplies
   independent information.

### Success criterion

- **Discovery:** an exactly reconstructed nonboundary rational point of `E8`,
  followed by exact lift tests through the degree-`24` cover and the `P8`
  family.
- **Closure input:** a certified list of `p`-adic zeros with precision and
  root-stability bounds.  Rational reconstruction alone will not be called a
  proof.

## Attack 2: a finite-prime Prym Mordell-Weil sieve

The first partial sweep suggests weak Chabauty may leave many mock zeros.  A
Mordell-Weil sieve should be built in parallel rather than relying on
reconstruction of individual `p`-adic approximations.

### Work

1. Make the Abel-Prym map of a rational point `P` explicit, either as
   `[P-iota(P)]` in the Prym or through the bigonal correspondence to `J(D)`.
2. For several good primes, compute:
   - the reduction of the rank-`1` Prym generator;
   - the Abel-Prym class of every point of `E8(F_p)`;
   - the congruence classes of `n` for which that class equals `n*G_p`.
3. Combine the congruences across primes.  Attach the surviving residue
   classes to the rigorously computed Coleman zeros from Attack 1.
4. If a nonboundary class survives, reconstruct and test it exactly on `E8`,
   then lift it to the full signed cover and finally to an exact genus-2
   Jacobian.

### Success criterion

- no compatible global generator multiple remains outside the two boundary
  points, proving `E8(Q)` is boundary-only; or
- a short, exact list of nonboundary fibers to test for `[6,12]`.

## Attack 3: replace the HLP line by a two-parameter section surface

The original HLP seed proves exact split `[6,12]`, and the tangent calculation
shows full-dimensional nonsplit directions.  The obstruction on `G_A` is the
choice of a generic line, not a lack of deformation space.

### Work

1. Work in the marked HLP incidence variety with two order-`3` contact blocks
   and one halving block.
2. Choose a two-parameter plane on which one contact block is forced to be a
   rational section.  Use the known tangent basis to ensure the plane is
   transverse to the visible Humbert-4 branches.
3. As a concrete first experiment, verify a second point in the diagonal HLP
   split family and test a plane/chord through two exact split `[6,12]`
   anchors.  The second seed and its marked classes must be verified in Magma
   before this construction is trusted.
4. Eliminate only the remaining independent contact and halving blocks.
   Reject planes whose modular covers are high-genus/irreducible in the same
   way as the `G_A` line; retain planes with a low-degree factor or rational
   section.
5. For every rational nonsplit specialization, require:
   - exact `TorsionSubgroup = [6,12]`, not mere containment;
   - explicit independence of the two order-`3` classes;
   - exact order `4` and no unwanted order `8` lift;
   - a geometric-simplicity certificate, preferably a `D4` Frobenius/root-
     power certificate at good reduction.

### Success criterion

A nonsplit rational specialization with exact torsion `[6,12]` and a rigorous
geometric-simplicity certificate.

## Deliberately deprioritized

- Finishing the present `E8` server run without fixing precision propagation.
- Treating rational reconstruction of Coleman zeros as a global proof.
- Extending the one-parameter HLP `G_A` resolvents: the contact cover already
  has genus `51`, while the primitive order-`4` cover is predicted to have
  genus `181` if connected.
- A wider blind affine TB/Richelot height search.  The existing masks show that
  unresolved points are concentrated in projective/saturation/bad-reduction
  charts; any return to this lane should derive those charts explicitly.

## Initial execution order

1. Commit the missing Prym rank/generator verifier and a precision audit.
2. Produce and validate the precision-safe `p=37` Coleman driver on a few
   disks.
3. In parallel, derive the finite-field Abel-Prym reduction map for the
   Mordell-Weil sieve.
4. Begin the two-anchor/two-parameter HLP experiment while the high-precision
   Coleman run is executing.

## Response to review (Fable, 2026-07-13)

Two review items closed immediately:

1. **Near-miss simplicity restored.**  The even-chi defect is real; the
   repaired certificate is at p=7 on the SOURCE curve of the near miss
   ((a,b)=(-409/200,0), torsion [2,2,6]):
   chi_7 = x^4 + x^3 + 5x^2 + 7x + 49, irreducible, a1 = -1 (odd part
   nonzero), no root-power degree drop through n=12.  Since the [2,12]
   dual is ISOGENOUS to the source (the R3-halving correspondence),
   geometric simplicity transfers to the dual.  The near miss stands.
   (code/contact6_m612_nearmiss_recertify.m, results/..._recertify.log)

2. **Prym rank verifier committed.**  Self-contained script asserting
   RankBounds(J(y^2=-3x^6+24x^3-75)) = 1..1, trivial torsion, the
   explicit generator (x^2+2x+4, 5x+5), and the exact model
   isomorphisms used by the bigonal transport.
   (code/contact6_m612_prym_rank_verifier.m, results/..._verifier.log)

Accepted without reservation: the precision-propagation defect in the
Coleman server driver (O(37^4) annihilator coerced into precision-20
parents => uncertified zeros), the 37^mv/p^mv bug, and the dropped inf
flag.  The server bundle is WITHDRAWN until the precision-safe driver
(review Attack 1.4-1.5) is committed; do not run the current
contact6_m612_E8_coleman_server.m.
