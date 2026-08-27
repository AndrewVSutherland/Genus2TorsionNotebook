# Computation audit of the paper source (2026-08-22)

Goal (AVS request): a complete, correct, **minimal** list of the code needed to
justify every claim in the manuscript that is not supported by an argument
in the paper or a citation to the literature.  The existing "Reproducibility and
data availability" section is IGNORED here (8 of its 15 listed files do not
exist in the repo; it predates the current manuscript).  Target: a new public
repository with at most ~a dozen files.

Numbering below follows the manuscript as it will compile (theorem counter is
per-section; the intro is unnumbered, so Preliminaries = Section 1).

## A. Claims that need included code

| # | Paper item | Computational claim(s) | Status |
|---|-----------|------------------------|--------|
| A1 | **Table 1** (`tab:census`, 72 groups, geom. simple) | exact torsion of every displayed curve | **covered**: `code/verify_simple_torsion_table.m` (new, 2026-08-22; ALL 72 CHECKS PASSED, 2.0 s; machine-checked to match the paper table row-for-row) |
| A2 | **Table 1** | geometric simplicity of every displayed Jacobian; End = Z for the generic rows discussed in the text; RM for [31] and [2,22] | partially covered: `data/claude_endz_certificates.txt` + checker `code/claude_endz_certificates.m` (the 8 non-database curves + Elkies [2,2,2,10] member), `data/claude_z31_rm_certificate.txt` + `code/claude_z31_verify.m` ([31]).  Production-LMFDB rows citable (BSSVY/CMSV endomorphism data); **alpha-only rows currently rest on the unpublished extended DB** — recommend one script certifying simplicity for all 72 witnesses uniformly (root-power criterion, 1–2 good primes per curve, minutes of CPU) |
| A3 | **Table 2** (`tab:splitcensus`, 75 groups, split) | exact torsion of every displayed curve | **covered**: `code/verify_split_torsion_table.m` (regenerated 2026-08-22 after fixing the [2,2,4,4] witness drift, see D1; ALL 75 CHECKS PASSED, 3.3 s; row-for-row match with the paper table) |
| A4 | **Table 2** | geometric splitness of the 9 non-database witnesses ([45], [63], [70], [7,7], [8,8], [6,12], [2,2,24], [2,2,4,4], [2,2,4,8]) | needs one small script: even-sextic bielliptic quotients for [45]/[63]/[7,7], degenerate `RichelotIsogenousSurfaces` codomain for the (2,2)-gluings (verified interactively today for [2,2,4,4]: factors have conductors 210, 2310 = 210.c5 x 2310.o4).  DB/literature rows citable |
| A5 | **Theorem 3.2** (three points on S°) | exact torsion [2,2,2,12] for all three curves; pairwise non-isomorphism; geometric simplicity (paper cites Magma + Lombardo) | needs one consolidated script in the *paper's* coordinates (the three S°-points).  Existing session-era pieces: `code/verify_record_22212_order96.m` (curve 1), `code/claude_sib_curve2_verify.m` (curve 2, chart coordinates), `data/claude_ari_curve3.txt` (curve 3).  Note: repo certificates use the root-power/D4 criterion, not Lombardo's algorithm — align the citation with what the script does |
| A6 | **Section 3, "Geometry of S"** (prose) | S has exactly 36 singular points, all rational ordinary double points with coordinates in {0,1,-1}; G-orbit counts 12 + 24 | **no script exists anywhere in the repo** — needs a short new script (SingularSubscheme of the (2,4) complete intersection, node check, orbit enumeration).  The adjunction/general-type sentence is argued and fine |
| A7 | **Section 4, [2,2,20] subsection** | (i) Eq. (eq:2220) has exact torsion [2,2,20] and End = Z (verified today: eq:2220 is isomorphic to the Table 1 row, so torsion is covered by A1; End = Z is in the certificates file); (ii) the one-parameter family's order-5 and order-20 classes; (iii) the t(z) parameterization; (iv) the degree-3 cover normalizes to Y: y^2+(x^2+1)y = x^5-x, conductor 2528; (v) rank Jac(Y) = 1; (vi) **Chabauty at p = 3** proves z = -1/7, -7/9 are the only good points, hence uniqueness in the family | needs one script `verify_2220_family.m` in the paper's coordinates.  `code/richelot_2220_all_double_linear_seeds.m` does the analogous Chabauty on a different model of the auxiliary curve — redo against the paper's presentation |
| A8 | **Corollary 4.3** (infinitely many [2,2,4,4]) | "E_7 has rank 1 over Q" | one-line rank computation (verified today: rank 1).  Fold into A9's script |
| A9 | **Proposition 4.4** (rational curve on A(2,2,4,4)) + **Lemma 4.5** (E_u ≅ V_t; the commented-out proof literally reads "Magma.") | symbolic identities over Q(s)/Q(t): P'_1, P'_2 lie on E_u, x(P'_1+P'_2) formula, x(P'_2) is a square, the model isomorphism, and the resulting a,b,c,d formulas | needs a small symbolic script (seconds) |
| A10 | **Theorem 4.7** (M(2,2,2,8) K3, dense points) | Disc(Q) = 0 ⇔ the quartic (eq:2228-k3) (from Lemma 2.3(b)); the elliptic model A_t with section (1,1); the displayed rational curve a,b,c,d(t) satisfies the quartic | K3 identification is cited (Bertin–Lecacheux); the identities need a small symbolic check.  Can share a file with A9 |
| A11 | **[2,22] subsection** | both displayed curves have [2,22] in J(Q)_tors (exact for the first — covered by A1); both have End(J_Qbar) = Z[phi] (RM 5); order-22 class of the form P - W_i | small script (or fold into A2's uniform certificate script); "we know of only two" is a DB-scope statement — cite the extended DB and scope the sentence |

## B. Claims that are fine as argued / cited (no code needed)

Theorem 1.1 (M(6)) — proof in paper.  Lemma 2.1 — cites Stoll §4.
Lemma 2.3 + Corollary 2.4 — derived from 2.1/2.2 + Zarhin citation.
Theorem 3.1 — self-contained algebra.  Proposition 4.2 — argued (but see C3).
Corollary 4.6, Remark 4.1 — argued.  Section 5 (Special loci) — prose.
Open questions — no claims.  Abstract/intro superlatives ("largest documented
was 80", "31 largest known prime") — literature + tables.

## C. Claims that need a wording fix or citation rather than code

- **C1. Lemma 2.2** (descent formulas for 2-torsion points): stated with no
  proof and no citation.  Add "See [Stoll2001, §4]" or a two-line proof.
- **C2. "One can show that A(2,2,4,4) is a Calabi–Yau threefold"** (after
  Cor. 4.6): no proof, citation, or computation.  Prove, cite, or soften.
- **C3. Proposition 4.2** smoothness conditions: "tedious but straightforward
  computation" — fine if kept, or add to the A9 script.
- **C4. Remark 6.1**: "the displayed [2,2,14] curve is the smallest of ten"
  (generic witnesses) — a search result; either ship the ten curves as a data
  file with checks, or weaken to "one of ten found here".
- **C5. Abstract vs Table 1**: the abstract lists [2,4,4] among the new
  groups, but Table 1 sources [2,4,4] to the extended DB (BookerSutherland);
  the intro says "six other new groups" while 7 rows are marked new.  Also the
  inline author query at line 150: the table count 72 is correct (72 = all known
  groups incl. the new ones; 66 would be the pre-existing count minus... the
  6-vs-7 bookkeeping should be settled at the same time).
- **C6. "Smallest-conductor known example"** (both tables): DB-query claims;
  record the query date/scope (as torsion_realizations.tex already does).
- **C7. The birational map Z → S** at the end of "How we found S°": a
  checkable identity inside discovery narrative; verify in the A6 script or
  present as narrative.

## D. Drift found and fixed today

- **D1. [2,2,4,4] split witness drift**: the paper's Table 2 displays the small
  gluing model y^2+(x^2+x)y = 60x^5+1000x^4-671x^3-5657x^2+867x+4913
  (Filip's paste, commit 68cc88d), while paper/split_torsion_table.tex, the
  witness JSON, and verify_split_torsion_table.m still carried the old
  source-51 quintic (a different, non-isomorphic curve).  Both curves have
  exact torsion [2,2,4,4].  Fixed: witness JSON + new_split_witnesses.txt +
  generator preamble updated to the small model (with today's Richelot
  splitness certificate: factors of conductor 210 and 2310), table tex + PDF +
  verifier regenerated, 75/75 pass, row-for-row match with the paper table
  re-verified.  The provenance note for the small model
  (notes/split_census_2026_08_13.md) lives only on the unmerged branch
  agent/torsion-cover-followups.

## E. Proposed public-repo file list (~10 files)

1. `verify_simple_torsion_table.m` — Table 1 exact torsion (72 checks; exists).
2. `verify_split_torsion_table.m` — Table 2 exact torsion (75 checks; exists).
3. `verify_simple_certificates.m` — geometric simplicity for all 72 Table 1
   witnesses + End = Z / RM certificates used in the text (A2, A11; assemble
   from claude_endz_certificates + z31 certificate).
4. `verify_split_certificates.m` — splitness of the 9 non-DB Table 2
   witnesses (A4).
5. `verify_22212_theorem.m` — Theorem 3.2, in surface coordinates (A5).
6. `verify_22212_surface.m` — 36 nodes + orbits of S; optionally the Z → S
   map (A6, C7).
7. `verify_2220_family.m` — [2,2,20] family + moduli curve Y + rank +
   Chabauty uniqueness (A7).
8. `verify_moduli_identities.m` — Prop 4.4, Lemma 4.5, E_7 rank, the
   [2,2,2,8] identities (A8–A10).
9. `README.md` — file-by-claim map, Magma version, runtimes, expected output.
10. (optional) `data_2214_ten_generic.txt` + check, if C4 keeps "ten".

Every discovery search (C++ sieves, funnels, gluing sweeps) stays out: no
claim in the paper depends on search output, only on the certificates above.

## F. Execution record (2026-08-22, second session pass)

All scripts of section E were written and now live in
`paper/scripts_and_data/` (the intended public-repo content), each
self-contained and asserting its claims:

| script | result | time |
|---|---|---|
| verify_simple_torsion_table.m | ALL 72 CHECKS PASSED | 2.0 s |
| verify_split_torsion_table.m | ALL 75 CHECKS PASSED | 3.2 s |
| verify_simple_certificates.m (A2, A11, C-parts) | ALL CERTIFICATES VERIFIED (74 curves) | 14–30 s |
| verify_split_certificates.m (A4) | ALL SPLITNESS CERTIFICATES VERIFIED | 0.6 s |
| verify_22212_theorem.m (A5 / Thm 3.2) | VERIFIED | 0.7 s |
| verify_22212_surface.m (A6, C7) | VERIFIED | <1 s |
| verify_2220_family.m (A7, incl. Chabauty) | VERIFIED | 1.4 s |
| verify_moduli_identities.m (A8–A10, Lemma 2.3(b)(d), 4.5) | VERIFIED | 0.4 s |

Combined log: `results/verify_paper_all_2026_08_22.log`.  Row-for-row
machine check of both paper tables vs the two table scripts: 72/72 and
75/75.

**Two manuscript errors found and fixed in the paper source:**
1. "Geometry of S": the 12-point node orbit representative was printed as
   [1:1:0:0:0], which does not lie on S; corrected to [1:0:0:1:0]
   (verified: the 36 nodes are exactly the G-orbits of [1:0:0:1:0] (12)
   and [1:1:0:1:1] (24)).
2. Lemma 4.5 displayed the y-part of the E_u ~ V_t map as
   (t-1)^4(t+1) y0/(x0+2(t^2-1)), which is wrong (fails numerically and
   symbolically); corrected to 2(t+1) y0/(x0+u)^2, which is proved in
   verify_moduli_identities.m.

**Directory reorganization:** paper/scripts_and_data/ now contains exactly
the eight verification scripts + README.md.  Moved out: halving.m,
halving_all.m -> code/ (their generic §1 check is subsumed by
verify_moduli_identities.m part (1)); tor2244.txt -> data/tor2244_bank.txt,
tor2228.txt -> data/tor2228_bank.txt (renamed: data/ already had DIFFERENT
files with the same names), ten2248models_abcd.txt -> data/ (identical copy
already there; duplicate removed).  All code/ and product/code references
to the old paths updated.

## G. End = Z upgrade (2026-08-23, AVS suggestions)

Two rounds, both prompted by the owner:

1. First extension (conductor route): End = Z for the database rows via an
   odd prime with conductor exponent 1 -- f_p = 2u + t + delta = 1 forces
   semistable toric rank 1, and End^0 of a geometrically simple surface
   acts unitally on X_*(T) tensor Q = Q, so End^0 = Q.  (NB the argument
   does NOT prove simplicity by itself: Q x Q maps unitally onto Q, and
   split Jacobians with exponent-1 primes exist -- the split [2,2,4,4]
   witness has ord_11(N) = 1.)  All 70 non-RM rows certified; runtime 268 s.
2. Final version (Zywina route, now PRIMARY): the two-Frobenius-polynomial
   criterion of Zywina's "Determining monodromy groups of abelian
   varieties" in its unconditional form (SS 1.5-1.6, pattern of SS 1.8):
   ordinary good primes q, p with Phi_q, Phi_p free of rank 3 (freeness
   certified exactly via valuation vectors of the roots + the product
   formula; maximal Frobenius torus => Mumford-Tate by Larsen-Pink Thm 4.3;
   ordinariness removes Zywina's Conjecture 4.5 by Noot), deg Q(W_p) = 8,
   and Q(W_p) linearly disjoint from Q(W_q) => |W(MT)| = 8 => MT = GSp_4
   => End(Jac_Qbar) = Z AND geometric simplicity, in one certificate.
   ALL 70 non-RM rows certified by this route alone (no fallbacks fired);
   script runtime dropped to 53 s.  The conductor and
   disjoint-fields+torsion routes remain in the script as fallbacks only.

## H. RM rows sharpened (2026-08-23, AVS question)

Can RM be certified from L-polynomials?  NO in the lower-bound direction
-- L-polynomials bound monodromy from below / End^0 from above, and no
finite Frobenius data can prove an endomorphism exists (Zywina's
non-maximal outputs are Monte Carlo predictions under his Conjectures
4.3/4.5 with ineffective exceptional sets, Rmk 1.15; even conditionally
two L-polynomials cannot see the Brauer class, Rmk 1.14(ii)).  What IS
certifiable was added to Part C: End^0 in {Q, Q(sqrt 2)} for [31] and
{Q, Q(sqrt 5)} for [2,22] -- two strict primes with non-isomorphic
quartic Frobenius fields of Galois-closure degree 8 bound the center in
the unique real quadratic subfield (CLV), excluding quartic CM and
imaginary quadratic; #tors = 31, 44 > 18 excludes QM (LSSV).  The single
remaining cited statement per RM row is the lower bound "End^0 <> Q"
(DB-certified endomorphisms for [2,22]; the CEHJMPV modular construction
for [31]), supported by the ~40-prime discriminant-core evidence.

## I. Trivial group added (2026-08-23, AVS decision)

The trivial group is now row 1 of the census (73 rows).  Witness hunt:
the smallest-conductor geometrically simple curve with exact trivial
torsion in the PRODUCTION LMFDB is 461.a.461.2 (generic, same isogeny
class as the [7] row), but the extended DB contains 277.a3 -- conductor
277, generic, the third member of the class 277.a whose other two members
are the [15] and [5] rows -- and nothing smaller (all alpha trivial-torsion
curves below 461: 277.a3, 349.a2, 353.a2).  NB 249.a.249.1, the smallest-
conductor End=Z curve overall, has torsion [14] (it is the [14] row), so
it cannot witness this row.  Row added to claude_gen_torsion_table.py
(ALPHA, cite BSSVY since production realizes the group), regenerated:
torsion_realizations.tex title/preamble updated (73 rows), PDF rebuilt,
verify_simple_torsion_table.m (73 checks PASS), verify_simple_certificates
(73 curves).  The paper-table edit (73 rows + caption) is with AVS.
Split-side analogue, if Table 2 gets the same treatment: smallest
geometrically split trivial-torsion curve in the DB is 1083.b.390963.1
(Q x Q); Table 2's caption has the same "known to arise" literalism.

## J. Trivial group added to the split table too (2026-08-23, AVS)

Smallest known geometrically split trivial-torsion curve: the alpha curve
961.a2 (conductor 961 = 31^2, Q x Q) -- it was ALREADY in
product/data/split_min_witnesses.csv (the original DB query found it) and
both split generators were explicitly filtering trivial torsion out; the
filters are now removed.  The earliest example is production/BSSVY
(1083.b.390963.1, the curve AVS suggested), so the source is \cite{BSSVY}.
split_torsion_table.tex now 76 rows (title/preamble updated, PDF rebuilt);
verify_split_torsion_table.m 76 checks PASS (empty-invariant emission
fixed: iseq([]) = "[]", not the polynomial fallback "[0]").  The paper Table 2
edit (75 -> 76 + caption) is with AVS.

## K. CORRECTION to J: split trivial witness is 1083.b.390963.1 (AVS catch)

AVS flagged the 961.a2 pick: the alpha curve 961.a2 is ISOMORPHIC to
production 961.a.961.1 (explicit isomorphism computed) -- the J0(31)
isogeny class -- and is geometrically SIMPLE with RM by sqrt 5, NOT split.
Verified in Magma: root-power simplicity certificate at p=7
(chi = x^4+4x^3+13x^2+28x+49, strict), RM disc-core 5; the alpha row's
is_simple_geom=false / geom_end_alg='Q x Q' is an ALPHA DATABASE ERROR
(isolated: the other four members of alpha class 961.a are correctly
simple/RM; the torsion value [] of 961.a2 is correct).  The split trivial
row now displays production 1083.b.390963.1 (BSSVY) -- AVS's original
suggestion -- whose splitness we verified rigorously: degenerate Richelot
codomain with elliptic factors of conductors 57 and 19.  Alpha's
1083.a2 is the same curve as production 1083.b.390963.1 (alpha model is
the x -> -x flip; alpha class letters differ from production).
verify_split_torsion_table.m: 76 checks PASS on the production model.

## L. Compact models for the split [2,2,24] and [2,2,4,8] rows (AVS request)

Both witnesses remodeled as quintics (one Weierstrass point at infinity,
one at 0), verified isomorphic to the originals with the same exact
torsion:
- [2,2,24]:  y^2 = x(52316x-156025)(2500x+3969)(32400x^2-34360x+255881)
  (2-rank 3: three rational linear factors + irreducible quadratic;
  incidentally 156025=395^2, 2500=50^2, 3969=63^2, 32400=180^2,
  255881=41*79^2); found by scanning all 12 (infty,0) Weierstrass
  labelings x scalings, minimal factored height 2.6e5 vs 1.6e10 before.
- [2,2,4,8]: y^2 = x(x+336100^2)(x+835200^2)(x+841500^2)(x+877221^2)
  (full 2-rank square-branch form of Lemma 2.3(b); the unique
  square-class-aligned labeling among all 30, smallest of the three).
Witness JSON now carries f/h (verifier data) plus a "display" field
(factored tex); split generators updated; verify_split_torsion_table.m
(76 checks) and verify_split_certificates.m (Richelot splitness on the
new models) both PASS; table tex/pdf regenerated.

## M. Splitness certificates extended to all 76 rows (2026-08-23, AVS request)

verify_split_certificates.m now covers every Table 2 row (was: 9 non-DB
rows).  Exploration of the (2,2)-isogeny graph (depth <= 3 over Q; depth
<= 2 over all quadratic fields Q(sqrt d) with d supported on the bad
primes) gave:
- 56 rows: degenerate Richelot codomain over Q at depth 1 (SetCart product
  E x E', or CrvEll = E x E^sigma over a quadratic field);
- 2 rows (the conductor-256 curves, [2,2] and [2,10]): degenerate codomain
  over Q(sqrt 2) at depth <= 2;
- 18 rows: NO (2,2)-route exists -- these are glued along odd torsion, so
  every product in their isogeny class is odd-isogenous and unreachable by
  Richelot chains.  For them the script verifies the anti-simplicity
  signature (no root-power-strict prime among all good p < 200 -- one
  strict prime would PROVE simplicity, so a split Jacobian fails at every
  prime; NB [19] = J_1(13) is Q-simple, chi_p often irreducible, and only
  the power-drop is the geometric fingerprint -- the naive
  "chi_p reducible" signature is WRONG for Q-simple split curves and was
  caught by exactly this row) and cites splitness (certified LMFDB
  endomorphism data; Howe 2015 for [70]).  As with the RM lower bounds,
  "an endomorphism exists" is not certifiable from finitely many
  L-polynomials.
The script is now generated by code/claude_gen_split_table_check.py from
the same witness data as the table (template
code/claude_split_certificates_template.m), 76/76 PASS in ~56 s.

## N. Explicit elliptic covers for ALL odd-glued split rows (2026-08-24, AVS idea)

AVS suggested certifying the Richelot-resistant split rows by exhibiting a
map from the genus-2 curve to a genus-1 curve (possibly over an extension),
found by any means but certified independently.  Implemented in full:
code/claude_elliptic_cover_finder.m computes the analytic Jacobian and its
numerical endomorphism ring (van Wamelen intrinsics), extracts rank-2
rational idempotents from split-minimal-polynomial elements, forms the
EXACT image lattice via HNF, solves the analytic idempotent directly from
eA*Pi = Pi*e (the naive fs<->Generators pairing is CROSSED -- the decisive
bug), builds the quotient's g2/g3 and Weierstrass-P by q-series, recognizes
everything in a number field of degree <= 6 (LLL/PowerRelation with
verification, spurious-convergent-safe thresholds), fits u = p(x)/q(x)
EXACTLY over K from rational sample points (real points of C where
available; the fit is exact so node clustering is harmless), and verifies
the certificate identity N(p,q)q = g h^2 exactly.  Results (all 18):
- 13 degree-3 covers over Q: [9]->19a1, [14]->26b1, [18]->35a1, [20]->11a1,
  [10]->33a2, [28]->182a1, [36]->90c3, [2,2,6]->20a2, [2,14]->294b2,
  [2,18]->114a1, [2,20]->462f1, [2,2,10]->66c1, and [70]->858k1 (Howe's
  curve -- the citation is no longer load-bearing);
- [25]: degree-5 cover over a quadratic field (the (5,5)-glue);
- [27]: degree-3 over a quadratic field;
- [3], [21]: degree-2 covers over CUBIC fields (bielliptic over a cubic);
- [19] = J_1(13): degree-2 cover over an explicit SEXTIC field (bielliptic
  over degree 6; needed precision 1000 and degree-<=-12 recognition).
verify_split_certificates.m now certifies ALL 76 rows by explicit geometric
objects (58 Richelot + 18 covers), zero citations, in under a second; the
certificates are embedded in product/data/split_cover_certificates.json and
verified by pure polynomial identities, independent of the finder.
