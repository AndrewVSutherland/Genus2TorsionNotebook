# Ranked unrealized torsion targets and realization plan

Date: 2026-07-17.

Scope: this note reviews the current project tree and ranks ten torsion groups
that the repository has not yet certified as occurring on geometrically
irreducible genus-2 Jacobians over `Q`, but which look most likely to occur
based on the accumulated scripts, data files, and notes.  I interpret
"unrealized" in the project-internal sense: no checked-in final certificate in
this repository currently gives a smooth genus-2 curve over `Q`, a verified
rational torsion subgroup of the listed isomorphism type, and a geometric
irreducibility / absolute-simplicity certificate for its Jacobian.

The ranking is deliberately constructive rather than purely speculative.  Each
entry is paired with the local evidence that makes it plausible and the shortest
plan I would run next.

## Ranking at a glance

| Rank | Target group | Confidence | Why it ranks here |
|---:|---|---|---|
| 1 | `[5,5]` | high | Open saturated finite-field points and smooth p-adic lifts already exist. |
| 2 | `[35]` | medium-high | The 7-contact/5-contact hybrid has narrowed to two coherent 3-adic central branches. |
| 3 | `[2,24]` | medium | Many tools exist; current low-height fibers are closed, but the route is mechanically clear. |
| 4 | `[48]` | medium | A16 infrastructure is mature and the missing 3-part now has an exact contact filter. |
| 5 | `[60]` | medium | Multiple independent order-30-to-60 and `M(12)`-plus-5 lanes exist, with known obstruction diagnostics. |
| 6 | `[8,8]` | medium-low | Shows up naturally in A(2,24) scans; likely needs a better 2-primary chart. |
| 7 | `[4,16]` | medium-low | Close to the successful `[4,8]` tangent-cover family and the Elkies `[32]` reconstruction. |
| 8 | `[2,2,20]` | medium-low | Simple specializations are strongly suggested, but a clean positive-dimensional/simple certificate is still missing. |
| 9 | `[72]` | low-medium | Contact-9 plus halving/root methods are available; no convincing search hit yet. |
| 10 | `[80]` | low-medium | Contact-5 order-40 infrastructure exists; doubling to 80 remains hard but well-posed. |

## Rank 1: `[5,5]`

### Why this is the best bet

The strongest unrealized target is the full independent 5-torsion problem.  The
`b2=0` full-norm branch is not merely a heuristic search region: after
saturation by the open product, the `h1=1, h2=0` slice has open finite-field
points.  The recorded `F_7` and `F_11` points are smooth in the five slice
variables, with unique p-adic lifting through the tested precisions.  Magma
also verifies that the contact class and the constructed `[U,V]` class are
independent 5-torsion classes in the sampled finite-field charts.

The failure so far is bounded rational reconstruction, not a local obstruction.
That makes `[5,5]` the most attractive target: the local geometry looks alive,
and the remaining problem is global point finding on a specific small chart.

### Plan

1. Rebuild the saturated `h1=1, h2=0` slice over `Q`, not only over `F_7` and
   `F_11`, and eliminate down to the lowest-degree projection in two variables.
2. Use the known smooth p-adic points as residue constraints for a lattice /
   Coppersmith-style rational reconstruction with larger height but fewer free
   coordinates.
3. Run a Mordell-Weil or elliptic/hyperelliptic quotient search if the
   eliminated curve has genus 0 or 1; otherwise compute maps to lower-genus
   quotients from the five-equation Jacobian-rank structure.
4. For each rational point, certify:
   - smooth genus-2 model over `Q`;
   - two independent 5-torsion divisor classes;
   - exact torsion lower bound containing `[5,5]`;
   - irreducible Frobenius polynomial at one good prime, followed by the
     12th-power/geometric-simplicity audit used elsewhere in the project.

## Rank 2: `[35]`

### Why it is plausible

The `[35]` search has already compressed dramatically.  The original
18 liftable first directions on the `b=0, r=1` pole chart are now reduced to two
central branches that survive through the deepest recorded 3-adic lifting.  The
central branches exhibit persistent affine-coset behavior controlled by two
left-obstruction rows.  This is a much better situation than an undifferentiated
height search: it suggests that either a symbolic obstruction exists or a
structured p-adic family remains to be parameterized.

### Plan

1. Replace raw lift enumeration with the compressed-state automaton proposed in
   the notes: state variables should be the two left-obstruction residuals, the
   affine constants, and the surviving central branch label.
2. Derive the first nonzero symbolic obstruction term along both central
   branches.  If it vanishes identically, extract a formal 3-adic parameter.
3. Use the formal parameter to create CRT targets combining 3-adic data with a
   small auxiliary prime where the 5- and 7-contact conditions are nonsingular.
4. Search rational reconstructions on the compressed parameter, not on all
   original coefficients.
5. Certify an exact order-35 class by separately checking the order-5 and
   order-7 multiples and excluding smaller order; then run the geometric
   simplicity audit.

## Rank 3: `[2,24]`

### Why it remains high despite recent closures

The recent A(2,24) work closed the four best low-degree fibers and showed that
the height-5 shell adds no low-degree rational branch.  That is negative
progress, but it is high-quality negative progress: the halving equations,
translation classes, boundary removal, and exact divisibility checks are now
well understood.  The target is still close to the repository's successful
`[2,12]` and `[4,8]` technology.

### Plan

1. Stop revisiting the four closed A(2,12) fibers and build a broader
   partitioned scan in a different A(2,12) chart or with a different height
   measure.
2. Prioritize fibers where the saturated halving projection has degree at most
   4, but do not require rational roots immediately; record genus-0 or genus-1
   components separately.
3. Add a modular prefilter for divisibility by 2 of the relevant order-12 class
   before exact `IsDivisibleBy` calls.
4. When a rational branch is found, certify `[2,24]` by combining the visible
   extra rational 2-torsion with a halved order-12 divisor and then performing
   the same Frobenius/geometric-simplicity audit.

## Rank 4: `[48]`

### Why it is close enough to keep prominent

The bounded A16 square-root route through RTHeight 4 is cold, but it produced a
useful lesson: the old point-count gate could be fooled by curves with torsion
`[16]` and no 3-part.  The new cubic-contact / exact-3 layer is the right
filter.  The A16 machinery is mature, and the missing condition is now stated
as an explicit cubic-contact system.

### Plan

1. Promote the cubic-contact exact-3 test from a diagnostic into the first-class
   production filter for all future A16 candidates.
2. Search new A16 candidates using parametrizations not exhausted by RTHeight 4,
   especially cubic-contact-native slices rather than square-root-first slices.
3. For each A16 candidate, require exact 3-torsion over `Q` before running any
   expensive `[48]` certification.
4. If exact 3-torsion appears, test whether the 2-primary generator and the
   3-primary class combine cyclically to order 48 or only give a product-type
   subgroup.
5. Finish with exact torsion and geometric-simplicity checks.

## Rank 5: `[60]`

### Why it is plausible

There are several routes toward `[60]`: halving order-30 classes, combining the
simultaneous contact-5/contact-6 family with an extra 2-primary condition, and
`M(12)`-plus-5 searches.  The repository records obstruction work for direct
order-30-to-60 halving, but the number of independent approaches keeps this
above the lower-confidence targets.

### Plan

1. Compare the existing order-30 family with the `M(12)` charts and identify
   the parameter loci where the order-30 class has the best chance of being
   divisible by 2.
2. Express the halving obstruction as a squareclass condition on the order-30
   parameter whenever possible.
3. Search for rational points on the resulting conic/elliptic covers before
   using larger height boxes.
4. Certify that the halved class is cyclic order 60, not merely a product of
   visible 2- and 30-torsion, and run exact/geometric checks.

## Rank 6: `[8,8]`

### Why it is not hopeless

The A(2,24) branch closure repeatedly produced degree-pattern rows involving
`[8,8]`.  Those rows did not yield rational points in the tested fibers, but
they indicate that the geometry naturally wants two independent 8-primary
classes.  This target likely needs a direct 2-primary normal form rather than
being treated only as a byproduct of A(2,24).

### Plan

1. Build a direct `A(8,8)` or `M_1(8,2,2)`-based cover using the explicit
   `M_1(8,2,2)` formulas in the inventory note.
2. Impose the second order-8 class via squareclass equations and saturate away
   full-split/product-locus boundary components.
3. Use finite-field point counts to reject components whose group orders are
   incompatible with `[8,8]`.
4. Certify independent 8-torsion classes and a geometrically simple Jacobian.

## Rank 7: `[4,16]`

### Why it belongs in the top ten

The repository already has a strong `[4,8]` tangent-cover family and a
reconstructed Elkies `[32]` family.  A mixed `[4,16]` target sits between these
successes: it asks for a second 4-primary class on a component that already
supports high cyclic 2-power torsion, or for one extra halving on a known
`[4,8]` family.

### Plan

1. Start from the `[4,8]` tangent-cover family and impose halving of the order-8
   generator or an independent order-4 translation, whichever gives lower-degree
   conditions.
2. In parallel, start from the Elkies `[32]` reconstruction and test for an
   independent rational 4-torsion class on special subloci.
3. Use Richelot/quadratic-pair decompositions to separate genuine simple
   components from split product components.
4. Certify exact `[4,16]` containment and geometric simplicity.

## Rank 8: `[2,2,20]`

### Why it is still unfinished

The inventory records a Lombardo-certified geometrically simple specialization
containing `[2,2,20]`, but also warns that the repository does not yet package a
positive-dimensional family with exact torsion in the way the main successful
rows are packaged.  If the task is to realize a subgroup at least once, this may
already be almost done; if the task is to produce a clean reproducible final
certificate, it remains a documentation and verification target.

### Plan

1. Locate the simple `[2,2,20]` specialization in the contact-5/order-20 extra-2
   data and write a single exact verification script for it.
2. Normalize the curve model and divisor representatives.
3. Check exact rational torsion contains `[2,2,20]` and record whether the full
   torsion is exactly `[2,2,20]` or larger.
4. Attach the Lombardo/12th-power geometric-simplicity certificate in a
   reproducible data file.

## Rank 9: `[72]`

### Why it is plausible but lower-ranked

The contact-9 family and rational-root order-18 subfamily are healthy simple
families.  A route to 72 would require adding substantial 2-primary structure,
probably repeated halving on an order-18 or order-36 class.  The project has
contact-9 and target-72 scripts, but the recorded status is still finite/root
search rather than a narrowed local branch.

### Plan

1. Reuse the contact-9 exact family and target the rational-root subfamily first
   so that the visible even part is explicit.
2. Formulate the first and second halving conditions for the order-18 class as
   squareclass covers.
3. Use local solubility at small primes to split the cover into viable residue
   classes before rational reconstruction.
4. Certify a cyclic order-72 class; if only `[8,9]` or `[2,36]` appears, record
   it separately rather than forcing a cyclic interpretation.

## Rank 10: `[80]`

### Why it makes the list

The contact-5/order-20 and order-40 infrastructure is substantial, including
explicit order-40 family work and plus-extra-torsion boundary searches.  Pushing
from 40 to 80 is a difficult additional halving problem, but the base family is
one of the best-understood odd/contact constructions in the repository.

### Plan

1. Start from the order-40 family, not the broader order-20 family.
2. Write the order-40 generator's halving obstruction as a norm/squareclass
   condition over the smallest possible parameter field.
3. Search low-degree components and rational points after saturating away the
   known boundary and fake full-split loci.
4. Certify an order-80 class and run exact torsion and geometric-simplicity
   checks.

## Cross-cutting verification checklist

Every candidate that survives to an exact curve should be promoted only after a
uniform certificate bundle is checked in:

1. the curve equation over `Q` and discriminant/non-singularity proof;
2. Mumford representatives for the claimed torsion generators;
3. exact multiplication checks proving the claimed subgroup;
4. a torsion upper-bound or enough reductions to rule out accidental ambiguity
   when exactness is claimed;
5. at least one irreducible degree-4 good-prime Frobenius polynomial;
6. the stronger geometric-simplicity / 12th-power transform audit where the
   project convention requires it;
7. a short note explaining whether the result is an isolated realization or a
   positive-dimensional family.

## Practical next work order

1. `[5,5]`: eliminate the smooth `h1=1,h2=0` slice over `Q` and attack rational
   reconstruction with p-adic residues built in.
2. `[35]`: implement the central-branch compressed automaton and symbolic
   obstruction calculation.
3. `[2,24]` and `[48]`: keep as parallel production searches, but only with the
   improved filters learned from the latest negative results.
4. `[60]`: re-express the order-30 halving obstruction as a low-genus cover.
5. Lower 2-primary and high-contact targets (`[8,8]`, `[4,16]`, `[72]`, `[80]`):
   run only after the top four lanes either produce a certificate or hit a new
   structural obstruction.
