# Skill-library verbatim-quote review — findings log

Reviewer: verification subagent. Method: every code/formula block that claims
to be verbatim/copied from a named repo file was diffed token-by-token against
the actual source; numeric claims checked against the cited notes; cited file
paths checked for existence. Severity: CRITICAL = math/formula differs from
source; MINOR = cosmetic/wording/abbreviation; OK = verified clean.

Context: a prior (crashed) review session already edited these files
(mtimes Jul 8): two-rank-and-factor-types 20:51, finite-prefilters 21:32,
magma-lab-conventions 21:38, running-torsion-searches 21:39,
validate-and-record-a-hit 22:20, component-boundary-analysis 22:20,
pell-cf-order 23:38. Everything below reviews the CURRENT (post-fix) text.

---

## 1. named-charts-reference — VERDICT: CLEAN

Verified OK:
- A(8) block (e, d, lambda, u, v, a, b, c, Q, q, f) vs `A8f` in
  code/agent_a2_24_composite8x3.m lines 35–49: token-level match
  (skill renames rv/pv/tv -> r/p/t consistently).
- g8 = x^2+u*x+v, L = r*x+(p-r^2), ellBase = -(q+Q*L),
  v8 = (-L*ellBase) mod g8, `fInt, L := IntModel(f)`, and the 8*D8=O /
  4*D8!=O order check all match the script.
- d=0 slice: d=0 <=> p = r(r+t)/2 (re-derived, correct; matches
  notes/agent_a2_24_d0_derivation.md), Q=x^2, f=q*(x^4+q), f(0)=c^2,
  curve B (1/3,-1/9,-1) on the slice, genus-1 rank-0 cover — all match.
- A(12) block (s, t, mu, lambda, T1, R, ell, Q, F, f) vs
  code/agent_a12_224_descent_setup.m lines 50–61: exact; F = Q^2 + R*ell^2
  is asserted at line 59 as claimed. P4/P6/P12 match (header lines 14–15 and
  numeric section); funnel quotes (u4/u6/L-scaled P4,P6; conic locus
  r=(v^2-p^2)/p; N1==N2==p(p+r)=v^2) match code/agent_a12_224_funnel.m.
- M(2,12) block (m, T, Y^2) vs notes/contact6_m36.md lines 661–663: exact.
  allowed_312 = 0 mod 5 and split/nonsimple small [3,12] points: match.
- M_1(8,4): t, A, B, f=x*A*B vs FamilyPolynomial in
  code/agent_m18_416_search_crt.m lines 334–340: exact. Y_R & Qfac vs
  YRational (343–344): exact. PlusDisc (347–350), MinusDisc (352–355):
  exact, coefficient by coefficient. FirstCoverPossible = IsSquareQ(Plus)
  or IsSquareQ(Minus): exact. T_x=[x,0], P_R=(-R,Y_R): match header line 10.
  Height-150 quiet second halving + aux primes 11,13,37,41: match
  notes/m18_416_component_and_smooth_strata_2026_07_02.md (§G, line 944).
  code/agent_m18_416_live_stratum_search.m exists.

Defects: none. Fixes applied: none needed. Open: none.

---

## 2. contact-torsion-constructions — VERDICT: 2 MINOR (no fix required)

Verified OK:
- Quintic-contact-5: h = 1+t*x+((t^2-1)/2)*x^2, f = h^2-((t+1)^4/4)*x^5 vs
  notes/how_we_found_2220_examples.md lines 11–12: exact.
  H = [x^2+2x/(t+1), (t+2)x+1], 2H = [x-1,0]: exact (lines 19, 25).
- Elementary version h = 1+a*x+b*x^2, f = h^2-(1+a+b)^2*x^5, P=(0,1),
  div(y-h)=5(P)-5(inf), x=1 forced root: exact vs
  notes/m10_quintic_contact5.md.
- Contact-6: h6, f = h6^2-(x-1)^6, and the 5-coefficient f-expansion vs
  notes/contact6_m36.md lines 18–22: exact (also re-derived independently —
  the expansion is arithmetically correct).
- Side conditions b != -3, h6(1)=a+b+2 != 0, divisor 6P-6inf, P=(1,a+b+2):
  exact. Factorization f = x*((b+3)x^2+(a-3)x+2)*(2x^2+(b-3)x+(a+3)):
  exact vs lines 496–497 and re-expanded (correct).
- Cubic-contact identity, lead(f)=lead(h3)^2+kappa, geometric origin: exact
  vs notes/agent_a2_24_contact.md. kappa = -m^2 reconciliation between the
  two notes: mathematically correct.
- B = c5*L^2+3U, Delta = 4*c4*L^2+12(U^2+v^2)-B^2, N = B/(2L), R = Delta/(8L),
  S = v^3/L, and all three [3,6] cover equations vs contact6_m36.md
  lines 68–83: exact, coefficient by coefficient.
- Extraction algorithm (g4 quotient, h3 = v+u*(sx+w), the deg-4 identity,
  R<S,W,K> := PolynomialRing(Q,3), Coefficient(expr,i) i in [0..4],
  Variety, V[1], verify residual & lead relation) vs
  code/agent_a2_24_contact_extract.m lines 51–73 and 103–107: exact.
- Validated data curve A (q3 = x^2-435/73x+2529/292, disc 4608/5329 =
  (48/73)^2*2, h3 = 1205280x^3-10769040x^2+31705560x-30713580,
  kappa = -1451998172160, lead check 701706240) and curve B
  (q3 = x(x-2/3), h3 = 13122x^3+26244x^2-34992x+15552, kappa = 1377495072,
  P1=(0,64/243), P2=(2/3,32/243)) vs notes/agent_a2_24_contact.md: exact.
  `curves := [ <5,-5/2,-9/2,...>, <1/3,-1/9,-1,...> ]` matches script line 87.
- Thin-set warning (Rabinowitsch ideal = (1); genus-1 rank-0 cover) vs
  notes/agent_a2_24_d0_derivation.md: exact.

MINOR-1: "(here `v` is the constant of `q = x^2 + U*x + v^2`)" — strictly
v^2 is the constant term, v its square root. Meaning is clear from the
displayed q; left as-is.

MINOR-2 (observation): the skill says the extraction script "writes
f = h3^2 - kappa*u^3 internally, then reports f = h3^2 + kappa*q3^3 after a
sign flip". In fact the script's COMPUTATION is consistently in the +
convention (identity at line 56 solves f - h3^2 = kappa*u^3; check at line
103 is f - (h3^2 + kappa*q3^3)); only the header comments (lines 7, 50) carry
the minus sign — i.e. the script's own comments are inconsistent, there is no
computational flip. The skill's warning is directionally right (both sign
conventions appear in that one file); left as-is.

Fixes applied: none. Open: none.

---

## 3. halving-and-doubling — VERDICT: CLEAN

Verified OK (all vs notes/contact6_m36.md):
- §1: u = x^2+alpha*x+beta, ell = l3*x^3+l2*x^2+l1*x+l0; ell^2-f = l3^2*u^2*q;
  ell+h3 == 0 mod q; six coefficients deg 0..5 of ell^2-f-l3^2*u^2*q; two
  coefficients of 8L*ell+(8x^3+4Bx^2+Delta*x+8v^3) mod q; sign only chooses
  E/-E; [6,6] check = two invariant factors divisible by 6; generated by
  code/contact6_m36_symbolic.m -> data/contact6_m36_symbolic_equations.txt.
  All exact (lines 94–132 of the note).
- §2: centered t = x-1; h6 = t^3+At^2+Bt+C with A=b+3, C=a+b+2, B=A+C-2;
  f = h6^2-t^6; u = t^2+pt+q, ell = rt^2+st-C; ell^2-f = -2At*u^2;
  ell(0)=-C -> residual (0,-C), 2H=D; H1..H4 exact token-for-token
  (note lines 637–640); p = (A^2+2B-r^2)/(4A), s = A*q^2/C-B exact
  (lines 646–647); symbolic script/data paths exist.
- §3: M(2,12) chart exact (lines 661–663); [3,12] search script path and
  filter description match (lines 672–679).
- §4 wall evidence: allowed_312 0 (line 686); 1234321 / 3837681 pairs
  (lines 695, 706); discriminant factorization
  r^6*(r+1)^21*z^2*(z+1)^11*(z-1)^11*H(z,r) (line 745); Rinf+Z0 (line 772);
  the three [3,12] cases z=-5/3,r=-3/5; z=-5/3,r=-9/35; z=-5/4,r=-32/35
  (lines 815–817); automorphism group order 4, two degree-2 subcovers
  (line 820). Known split specialization z=-5/3, r=-3/5, m=-10/9 (line 713).
- Pitfall "[3,6]/[6,6] finite filters coincide" matches line 156.
- Contact-7 corroboration ("p 5 ... h4_divisible_by_2 0", "NO OPEN H4
  HALVES") vs notes/contact7_family.md: VERIFIED (see session-2 addendum
  below).

Defects: none. Fixes applied: none.

---

## 4. simplicity-certificates — VERDICT: 3 MINOR (2 fixed this session)

Verified OK:
- SimplicityCertificate cited at lines 75–97 of
  code/agent_a2_24_composite8x3.m: correct span. Prime list, admissibility
  conditions (lead coeff, Numerator(dsc), Degree>=5 + IsSquarefree): exact.
- a1/a2/chi lines and the Step-3 certificate block: token-for-token match
  (lines 86–94).
- CountCurve description vs lines 67–73: exact, including the
  IsSquare(LeadingCoefficient) infinity term.
- [2,2,20] p=71: chi = X^4+2X^3+14X^2+142X+5041, 12th-power transform
  irreducible, cited lines 91–97 of notes/how_we_found_2220_examples.md: exact
  (5041 = 71^2 checked; functional equation consistent with a1 = -2).
- [6,6] p=23: L_p = 529*T^4-26*T^2+1 irreducible, curve/parameters
  a=133/39, b=-7/13, f coefficients: exact vs contact6_m36.md lines 570–576.
- Sage calls J.geometric_endomorphism_algebra_is_field(B=100) /
  J.geometric_endomorphism_ring_is_ZZ(B=100) vs
  code/m212_extra3_geom_simple_check.py: exact. Cited line numbers all
  correct: how_we_found 84–89 (True/True), contact5_order40_family.md
  341–342 (True/True), m212_three_torsion.md 162–163 (False).
- Anti-signature [2,6,6] examples and L_p = (11T^2+1)^2, (13T^2-2T+1)^2 vs
  contact6_m36.md lines 459–468: exact. Cremona 90c3/510g1 +
  Degree2Subcovers/AutomorphismGroup vs code/m212_extra3_split_certificate.m:
  match. Bielliptic skip idiom: exact. R25_4_SB_v4_certificate.m exists.

MINOR-1 (FIXED): the gloss "here `L_p(T)` is the reciprocal
`p^2 T^4 chi(1/T)` normalization" was off by a factor p^2 — the correct
reciprocal is `T^4 chi(1/T)` (leading coeff p^2, constant 1; check:
529T^4-26T^2+1 = T^4 chi(1/T) for chi = T^4-26T^2+529). Fixed in the skill.

MINOR-2 (FIXED): "whose 12th-power transform is irreducible — i.e. no degree
drop through n=12" overstated: irreducibility of the 12th-power transform
witnesses no drop at n=12 (hence at divisors of 12), not at every n <= 12.
Reworded to "at n = 12".

MINOR-3 (left): the skill says the certificate function "returns
`false, 0, 0`"; the code returns `false, 0, RT!0` (zero polynomial).
Semantically identical; left as-is.

Open: none.

---

## 5. local-obstructions — VERDICT: CLEAN

Verified OK:
- Cross-difference/squareclass criterion and the K3 surface
  (ab+ac+ad+bc+bd+cd)^2 = 4abcd vs notes/a2244_local_obstructions.md: exact.
- p=11: 240 good-reduction roots, 0 surviving; p=23: 5280, 0 — match both the
  table (line 71ff) and the mod-p counts section (lines 271–285).
- Signed-class rationale (-1 nonsquare at 11 and 23): matches note line 408.
- code/a2244_padic_signature_sieve.py: squareclass_qp_signed returns
  (parity, Legendre(unit)) or "0" (lines 87–92) — exact;
  partition_status_mod_pk killed/deep/resolved_ok semantics (lines 129–152)
  — exact; real_partitions_for_sorted_positive -> ("12|34",) (lines 187–190)
  — exact; CLI defaults --primes 11 23, --max-depth 3 — exact.
- Depth ladder numbers (118 tuples; k=1: 28/9; k=3: 14/4; k=4 + real: 0)
  vs note lines 425–433: exact.
- code/a2244_component_adic_analysis.py: verdict priority string
  "smooth_p3_resolved > modp_resolved > p2_resolved > deep_only > killed"
  is literally in the code (line 270); all 30 (component,partition) pairs
  smooth_p3_resolved at both primes vs note: exact.
- code/a2244_small_prime_congruences.py: default cases "2:5,3:3,5:2,7:1"
  match the (p=2 depth 5, p=3 depth 3, p=5 depth 2, p=7 depth 1) claim;
  "cannot support four distinct nonzero squared branch residues" quote: exact.
- [3,6] generic counts (4, 6, 40, 90, 98, 310, 224, 350, 1534) vs
  contact6_m36.md lines 273–281: all nine match.
- allowed_312 quote "p=5: allowed_312 0, bad 19, good rank counts [<1,6>]."
  vs line 686: exact.
- m18_416: "TwoDescent residual num_covers = 0",
  "MordellWeilShaInformation = [3, 3]", genus-5 square-condition cover vs
  notes/agent_m18_416_R25_4_SB_descent_followup.md lines 93–100: exact.
  R=-8: three V4 quotients rank 1, MW invariants [2,2,0] (= torsion [2,2] +
  rank 1), sieve does not close through p=89 (200001 CRT cap), no local
  obstruction to p=101, "useful negative information..." quote vs
  notes/agent_m18_416_R8_mwsieve_attempt.md: exact. R=-29/8 and R=-25/4
  rank-zero kills vs notes/agent_m18_416_other_ELS_v4_scan.md: match.
- Quote "the obstruction is not geometric emptiness ... thin arithmetic
  problem" vs contact6_m36.md lines 284–287: accurate ellipsis.

Defects: none. Fixes applied: none.

---

## 6. component-boundary-analysis — VERDICT: 2 MINOR (1 fixed this session)

(Reviewed post prior-session edit of Jul 8 22:20.)

Verified OK:
- §1 vs contact6_m36.md: D = 8*L^2*(v^3-1); G2/G3 numerator constructions;
  only common factor L^2; nonboundary product
  L*v*(v^3-1)*(U^2-4v^2)*(b+3)*(numerator of a+b+2); G2core deg 18 / 116
  terms, G3core deg 10 / 39 terms; v=1 branch H (deg 6, 10 terms) and Jcore
  (deg 18, 99 terms), gcd 1; closing quote — all exact.
- §2 vs notes/agent_a2_24_d0_derivation.md + code/agent_a2_24_d0_saturate.m:
  covering-condition explanation, X = {(r,t,m) : Eq4=Eq3=0, kappa != 0},
  genus-1 degree-12 plane curve (both projections), B at (1/3,-1,2/9): exact.
  Rabinowitsch block formulas (kap, bN, Eq4, Eq3, ideal) token-match the
  script (lines 11–14).
- §3 boundary analysis: discriminant factorization, 36 = 6 + 30 classes,
  Rinf+Z0, 48/45/0 counts, split certification: exact.
- Pitfalls: kappa=0 spurious branch -> torsion [8], rank-1 elliptic
  y^2 = x^3+x^2-4x (note lines 21–24): exact. GCD/resultant/factor-selection
  block vs code/agent_a2_24_d0_cover_ptsearch.m lines 23–24, 30: token match
  (the "// then select..." line is clearly reviewer annotation). 4-point
  result quote: matches note lines 43–45. EllipticCurve-slow quote: exact
  (note line 39). m18_416 d4 = (R-1)*4*F5 repeated deg-5 factor: exact
  (m18_416_component_and_smooth_strata note lines 54–99).

MINOR-1 (FIXED): quoted printf line was abbreviated —
skill had `printf "Ideal is (1)? %o\n", I eq ideal<S|1>;` where the source
line 17 is `printf "Ideal is (1) [no genuine soln generically]? %o\n", ...`.
Fixed to match source exactly.

MINOR-2 (left): the quoted 4-point block drops the source's parenthetical
"(h3 not a genuine cubic)" after "m=0, degenerate". Content unchanged;
left as-is.

Open: none.

---

## 7. two-rank-and-factor-types — VERDICT: CLEAN (pending 2 numeric checks)

(Reviewed post prior-session edit of Jul 8 20:51.)

Verified OK:
- TwoRank block vs code/agent_a2_24_wsplit_3tors.m lines 70–78:
  token-for-token. Identical copies confirmed in agent_a2_24_composite8x3.m
  (lines 57–65) and agent_a2_24_cf_search.m (lines 63–71).
- SplitWPvals block vs wsplit lines 50–62: token-for-token, including the
  leading comment. E(p) prose formula matches the script header (line 9).
- Sextic closed form and table: re-derived; consistent with the coded
  function (even-subset count 2^(k-1) if odd factor present, else 2^k).
  All 12 table rows arithmetically correct. Quintic accounting (append a
  degree-1 orbit) correct; "coded TwoRank returns 1 on the [1,2,2] quintic"
  consistent with k-2 = 1.
- CF family [1,1,2,2] claim matches agent_a2_24_cf_search.m header lines 4–5.
- CountCurve infinity-term quote: exact.
- Bielliptic skip quote: exact (composite8x3 line 141).

Still open at time of first crash-save:
- density 3e-5 / "31 per million", "~90% genuinely 2-rank 2", "H=26 run,
  ~1M curves" vs notes/agent_a2_24_composite.md — VERIFIED in session-2
  addendum below.
- Magma spot-checks of table rows — run in session 2, see Magma section.

---

## 8. magma-lab-conventions — VERDICT: in progress at crash-save

(Reviewed post prior-session edit of Jul 8 21:38.)

Verified OK so far:
- try/catch block vs code/agent_a2_24_wsplit_3tors.m lines 167–171:
  token-for-token (including `catch ee`).
- Parameter-guard idiom vs code/agent_a2_24_ztors_sample.m: both lines are
  verbatim (H := 18 at line 24; MemGB guard + SetMemoryLimit at lines 18–19).
  NOTE (cosmetic): in the file the MemGB/SetMemoryLimit lines precede the H
  line; the skill lists H first. Verbatim per-line, order differs — left.
- ClearToR3 coercion idiom vs code/agent_a2_24_d0_cover_ptsearch.m
  lines 14–21: faithful (skill splits one expression into two lines with a
  named temp; flagged as "pattern", acceptable).
- `okc, crat := RationalReconstruction(Integers(pk)!cc);` — exact
  (agent_a2_24_construct_lift.m line 176).
- TC_* list vs code/torsion_cover_lab_utils.m: exactly the 8 named functions
  exist (TC_SumInts 10, TC_MakeMonic 17, TC_HeightRationals 24,
  TC_FactorDegreeMults 37, TC_ContainsPoint 44,
  TC_GoodReductionPolynomial 53, TC_PointCountGate 67,
  TC_NormalizeTorsionInvariants 100), no TC_* omitted.
- SplitWPvals/Roots(Numerator(E)) and Variety idiom cross-refs: verified
  earlier. Rabinowitsch cross-ref: verified. Mumford scaling: integral-model
  `v8 := (-Lden*ellBase) mod g8` matches ztors_sample/wsplit; finite-field
  unscaled `v8 := (-ellBase) mod g8` confirmed at
  code/agent_a2_24_locus_geometry.m line 92.

Open at crash-save (completed in session-2 addendum below):
- AUTHORING.md: "at most 3 concurrent Magma jobs / OOM at 6 / MemGB:=3".
- ERRATUM section in notes/agent_m18_416_R8_dA_quotients.md and the
  three-Factorization-bugs claim (lost 2 in c4, lost -2 in C2).

---

## 9. pell-cf-order — VERDICT: in progress at crash-save

(Reviewed post prior-session edit of Jul 8 23:38.)

Verified OK so far:
- SqrtPolyPart + CFOrder blocks vs code/agent_a2_24_cf_search.m lines 36–61:
  token-for-token, including the `// exact D_infty order via CF; ...`
  comment and every guard (`Qi eq 0`, exact-division check, `i ge 1 and
  Degree(Qi) le 0 and Qi ne 0`).
- Family f = (x^2-1)(x^2+ax+b)(x^2+cx+d), factor type [1,1,2,2]:
  matches cf_search.m header lines 4–5. Validate:=true mode exists.

Open at crash-save (completed in session-2 addendum below):
- The `CFOrder(f, 40)` + `ord mod 24 eq 0` deployment claim (search loop).
- The three test vectors vs code/agent_a2_24_cf.m (f14 -> 14, f18 -> 18,
  f28-candidate -> 7 with degree sequence [3,1,2,1]).
- construct_seeds / construct_lift claims (IsUnit guard, degree pattern
  [3,1,1,...,1], ClosureG manual division, dead-end outcome).

---

(Sections 10–14 and the session-2 addenda follow as the review proceeds.)

---

## Session-2 addenda (resumed pass)

### Skill 3 addendum — contact-7 corroboration VERIFIED
notes/contact7_family.md lines 309–310: "p 5 checked 250 root_good 48
surface 4 exact_h4 4 h4_divisible_by_2 0" / "NO OPEN H4 HALVES" — the skill's
quote is accurate. Skill 3 verdict stands: CLEAN.

### Skill 4 addendum — fixes applied this session
- FIXED: "reciprocal `p^2 T^4 chi(1/T)`" -> "reciprocal `T^4 chi(1/T)` ...
  with leading coefficient `p^2`".
- FIXED: "no degree drop through `n=12`" -> "no degree drop at `n=12`
  (hence at every divisor of 12)".

### Skill 6 addendum — fix applied this session
- FIXED: Rabinowitsch block printf restored to the source's exact string
  ("Ideal is (1) [no genuine soln generically]? %o\n").

### Skill 7 addendum — remaining numeric claims VERIFIED -> CLEAN
- "~3e-5" blind 2-rank-2 density: grounded in
  code/agent_a2_24_wsplit_3tors.m header line 6 ("2-rank 2 is ~3e-5 rare").
- "31 per million": grounded in data/agent_a2_24_composite_h20_part0.log
  ("PROGRESS tested=1000000 ... rank2=31").
- "~90% genuinely 2-rank 2" and "~1M curves in the H=26 production run":
  grounded in notes/agent_a2_24_composite.md ("991,275 curves 2-rank-2
  (~90%!)") and data/agent_a2_24_wsplit_h26_part*.log (RH=TH=BH=26; part0
  ends with histogram [2,8]:1708).
Final verdict skill 7: CLEAN (Magma spot-checks below).

### Skill 8 addendum — remaining items VERIFIED; 1 fix applied -> VERDICT: 1 MINOR (fixed)
- AUTHORING.md lines 63–64: "at most 3 concurrent Magma jobs; the machine
  has OOM'd at 6" — matches the skill's house limits. MemGB:=3 default
  matches the search scripts.
- ERRATUM section exists: notes/agent_m18_416_R8_dA_quotients.md line 73.
  The three-unit-bugs claim is grounded: same note line 61 ("third
  unit-constant bug in this project (after c4 and C2)"), and the specific
  constants in notes/agent_m18_416_p7_blowup_notes.md lines 340–342
  ("`c4` carries a factor 2 and `C2` a factor `-2`; both bit us once").
- MINOR (FIXED): the skill pointed ONLY at the ERRATUM section for the
  spurious rank-0 retraction; in the source the unit-bug retraction is in
  the "Validation and a cautionary note" section (the ERRATUM retracts the
  related d_A character quotients). Fixed the pointer and added the
  p7_blowup_notes cite for the c4/C2 constants.
- Guard-idiom line order (MemGB before H in the file): cosmetic, left.

### Skill 9 addendum — remaining items VERIFIED -> VERDICT: CLEAN (pending Magma run)
- Test vectors vs code/agent_a2_24_cf.m lines 72–76: polynomials exact
  (f14 = (x^2+1)(x^4+5x^2+4x+4) expect 14; f18 = (x^2-x+1)(x^4-x^3+9x^2+8x-8)
  expect 18; "f28-cyclic?" = x^6+2x^5-5x^4-14x^3-3x^2+24x+28 with expect 0
  in the harness — the skill's "-> CFOrder = 7, degseq [3,1,2,1]" is a
  runtime claim, checked by Magma below).
- Deployment: `ord := CFOrder(f, 40)` (line 136) and `ord mod 24 ne 0 ->
  continue` (line 142) in code/agent_a2_24_cf_search.m; CF is the first
  nontrivial filter (only Discriminant(f) eq 0 precedes). Exact.
- construct_seeds: IsUnit(LeadingCoefficient(Qi)) guard (line 46), generic
  cell "[3,1,1,...,1]" (line 16). construct_lift: ClosureG with manual
  unit-lc division (lines 35–62), freezes (a,b) and lifts (c,d) (lines
  81–89), RationalReconstruction arity-2 (line 176). All match.
- "where the [2,24] effort went instead": consistent with
  notes/agent_a2_24_composite.md line 3.

---

## 10. running-torsion-searches — VERDICT: 3 defects (all FIXED)

Verified OK:
- Funnel stage order vs code/agent_a2_24_wsplit_3tors.m main loop (lines
  150–181): enumerate (r,t,beta) -> SplitWPvals -> Degree/Discriminant/
  even-sextic skips -> IntModel -> TwoRank ge 2 -> 3|#J prefilter (PreP has
  exactly 14 primes) -> D8 try/catch order check -> TorsionSubgroup ->
  SimplicityCertificate on hits only. Order exact.
- Sharding pattern (ridx/NParts/Part): token match (lines 136–139).
- Marker formats: PROGRESS field order matches printf line 147;
  TARGET_2_24 matches line 188; SEARCH_DONE matches line 196;
  TORSION_HISTOGRAM "[ 2, 8 ] : 1708" matches data/agent_a2_24_wsplit_h26_
  part0.log tail. The HIT24 sample line is verbatim-real (see fix 3).
- HeightRationals block: token match vs wsplit lines 122–128.
- "index ~7104 of 12175 at H=100": EXACT (recomputed in Python — 1/3 is
  1-based index 7104 of 12175).
- §8 budget numbers (~1M rank-2, 1709 prefilter survivors, [2,8]:1708):
  match notes/agent_a2_24_composite.md and the h26 logs.
- §9 "1 genuine point in 6.3M pairs": matches d0fast_h45_part*.log
  (2097524+2100033+2100033 tested, genuine_hits totals 1: GENUINE r=1/3
  t=-1 m0=2/9) and the d0 note.

DEFECT-1 (numeric, FIXED): "#HeightRationals(45) = 1263" — actual count is
2511 (recomputed; also 2511^2 = 6,305,121 ≈ the 6.3M pairs of the H=45 d0
scan, consistency check). Fixed to 2511 and tied to the 6.3M figure.

DEFECT-2 (numeric, FIXED): §7 said the d0 fast scan was "at H=15"; the
ground (notes/agent_a2_24_d0_derivation.md line 48 and
data/d0fast_h45_part*.log) says H=45. Fixed to H=45. (No H=15 run exists
in data/.)

DEFECT-3 (attribution, FIXED): the HIT24 sample line was attributed to
notes/agent_a2_24_composite.md; the literal line lives in
data/agent_a2_24_composite_h12_part1.log line 21 (the note records the same
curve in prose). Attribution corrected.

Open/unverifiable: §6's "~25k pairs/sec dropped to ~5k/sec at height ~100"
— no repo artifact records these rates (likely live-session observation).
Left in place, flagged here.

---

## 11. validate-and-record-a-hit — VERDICT: 2 MINOR (both FIXED)

(Reviewed post prior-session edit of Jul 8 22:20.)

Verified OK:
- The quoted HIT line "HIT24 r=5 p=-5/2 t=-9/2 torsion=[ 24 ] 2rank=1
  simple=true (q=17 chi=x^4 - 5*x^3 + 16*x^2 - 85*x + 289)" is verbatim-real:
  data/agent_a2_24_composite_h12_part1.log line 21.
- The certificate record "torsion [24]; SIMPLE (chi_17 =
  x^4-5x^3+16x^2-85x+289, irred; D4 at 8 primes)" is exact vs
  notes/agent_a2_24_composite.md line 32.
- ERRATUM reference: notes/agent_m18_416_R8_dA_quotients.md HAS an
  "## ERRATUM (same day)" section (line 73); it retracts the d_A character
  quotients (non-abelian D_4 quartic extension, no (Z/2)^3 character group),
  keeps the original text in place, names the false statement + reason +
  what survives (the C_A quartic route). Skill description accurate after
  fix 2 below.
- kappa=0 cautionary tale (rank-1 elliptic through B; every member except B
  torsion [8]; genuine cover genus-1 rank-0): exact vs
  notes/agent_a2_24_d0_derivation.md.
- Step-3 exact-order pattern: faithful positive-form of the wsplit check
  (labeled "verbatim shape", acceptable).
- Co-author line "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>":
  matches 20/20 recent commits in git log. "The user pushes (agent push
  fails on auth)": matches AUTHORING.md line 65.

MINOR-1 (FIXED): hit-line attribution pointed only at the note; now points
at the log for the verbatim line.

MINOR-2 (FIXED): the ERRATUM description bundled the Factorization unit trap
into the ERRATUM's underlying error; in the source the unit-trap retraction
is the separate "Validation and a cautionary note" section of the same note.
Reworded to separate the two retractions.

---

## 12. finite-prefilters — VERDICT: 1 MINOR (FIXED)

(Reviewed post prior-session edit of Jul 8 21:32.)

Verified OK:
- ThreeTorsionPrefilter + 14-prime PreP block vs
  code/agent_a2_24_wsplit_3tors.m lines 33, 112–120: token match after fix
  below. The 10-prime earlier list [7,11,13,17,19,23,29,31,37,41] matches
  code/agent_a2_24_composite8x3.m line 32.
- "bad reduction -> continue, never return false": matches the code (the
  Degree/IsSquarefree guard is a skip).
- Production numbers (~1M rank-2 -> 1709 prefilter survivors (~1 in 580:
  991275/1709 = 580.0) -> 0 genuine, histogram [2,8]:1708): match
  notes/agent_a2_24_composite.md "3-torsion is the wall" section.
- [3,12] subgroup-compatibility filter + allowed_312 pattern-B edge case:
  match contact6_m36.md (search script path exists).
- m2224_plus3 claims: "first genuine good-reduction obstruction at p=13"
  (note, Finite-field diagnostic section); "residue-driven p=13-boundary
  enumerator + point-count primes through 73 -> ~6 tuples per box, killed by
  cubic-contact filters" — grounded in the note's later residue-enumerator
  section ("point_survivors 6 ... The six tuples were ... those six are
  eliminated by the cubic-contact filters", kills at primes 23/29/37).
- m3222_plus3 claims: open finite condition EMPTY at p=7,11,13 (allowed3=0
  lines 56–58); "must reduce to boundary at all three" (line 61);
  "through p=73 eliminates every rational parameter of height at most 20"
  (line 126). Exact.
- Cost-ordering and "don't filter the chart's guarantee" sections:
  consistent with the wsplit funnel and A(8) chart facts.

MINOR-1 (FIXED): the "(verbatim)" ThreeTorsionPrefilter block had an added
inline comment "// bad reduction: skip" not present in the source; removed
(the explanation lives in the bullet below the block).

---

## 13. g2-torsion-lab (hub) — VERDICT: CLEAN

- Skill index: every indexed name resolves to an existing directory under
  skills/ — contact-torsion-constructions, halving-and-doubling,
  named-charts-reference, simplicity-certificates, two-rank-and-factor-types,
  validate-and-record-a-hit, magma-lab-conventions, running-torsion-searches,
  finite-prefilters, local-obstructions, component-boundary-analysis,
  pell-cf-order, target-playbook (+ the hub itself = 14 directories, exactly
  what exists; no referenced-but-missing skills, no unindexed skill dirs).
- Repo map: notes/*.md = 78 ("~80" OK); code = 258 .m + 32 .py ("~300 Magma
  + some Python" is a slight overcount of the .m files but within "~");
  paper/main.tex, results/, data/ all exist. TC_* helpers exist.
- Naming table: m2220, m36, m3222, m2224, m2226/m2228/m2248, m18_416,
  a2_24, Z35, Z48, z5x5, contactN all resolve to notes/ or code/ files.
  (m3222 = the M_1(8,2,2) [2,2,8] base whose worker line targets [2,2,16];
  the hub's "[2,2,16] worker line" phrasing matches repo usage, e.g.
  notes/agent_m18_416_p7_blowup_notes.md line 467.)
- State of play: [2,2,20], [6,6], cyclic Z/24 all documented in the cited
  notes (verified during skills 2/4/11); the six main_four_target_* notes
  exist; all listed frontier targets have notes.
- Operating rules: 3-job cap + OOM'd-at-6 (AUTHORING.md lines 63–64);
  -b buffering; user-pushes (AUTHORING.md line 65); co-author line matches
  20/20 recent commits.

Defects: none.

---

## 14. target-playbook — VERDICT: CLEAN

Spot-checked grounded claims, all verified:
- [2,24] "3-torsion is the wall": ~90% W-split density, absent 3-torsion
  through ~1M samples, "no local obstruction (F_13..F_29), sparse global
  point" — matches notes/agent_a2_24_composite.md; (r,p,t)=(5,-5/2,-9/2)
  with chi_17 matches the note's record.
- [6,6] example: a=133/39, b=-7/13, L=29/16, U=-9/4, v=5/2, torsion [6,6],
  certificate p=23, L_p=529T^4-26T^2+1 irreducible — exact vs
  contact6_m36.md lines 570–576; [1,2,2] core route and nonsimple
  extra-root caveat match.
- Z/5xZ/5 fixed-quartic claim: Phi(X) = X^5-(X-1)^5 =
  5X^4-10X^3+10X^2-5X+1, IsIrreducible true, "empty off the same-contact
  boundary" — exact (coefficient-by-coefficient) vs
  notes/agent_z5x5_contact5_contact5.md lines 55–95.
- (b) M(12) line a=(1-r)/4, good mod-7 residues r=3,4, combined residue
  sieve to height 300 leaving only singular r=-1: exact vs
  notes/m12_z12x2_halving.md (lines 6, 62, 105–111).
- [4,16]: height 150, six smooth Q_7 strata, P_R never halves, tangent
  sieve {17,23,29,47}: exact vs
  notes/m18_416_component_and_smooth_strata_2026_07_02.md (lines 700–944).
- Z/48: A16 candidates all exact-test to [16], has3=false: exact vs
  notes/agent_Z48_cubic_contact_route.md (lines 127–131).
- Fifth-pass priority list (Z/35 highest / Z/5x5 second / Z/48 bounded
  background / A(2,24) low, with those exact worker descriptions): exact vs
  notes/main_four_target_fifth_pass_2026_07_02.md ("Current expected
  priority" block).
- [2,2,20] route facts (z=-1/7, t=-8233/7225, type 1+1+2, p=71 Frobenius):
  match notes/how_we_found_2220_examples.md.
- m2224/m3222 "+3" claims: same grounds as skill 12 (verified above).
- notes/order60_attempts.md, notes/agent_Z35_next_route.md ("simultaneous
  contact7/contact5 point equations"), notes/m244_to_248_route.md,
  notes/contact5_order40_family.md: all exist and match their one-line
  descriptions.

Defects: none.

---

## Magma validation runs (single process, SetMemoryLimit(2*10^9), <2 min)

Script: scratchpad/review_checks.m; binary magma.
Results — 12/12 PASS:

```text
TwoRank [1,1,2,2]   = 2 (expect 2)  PASS
TwoRank [2,4]       = 1 (expect 1)  PASS
TwoRank [1,1,1,1,2] = 3 (expect 3)  PASS
TwoRank [2,2,2]     = 2 (expect 2)  PASS
TwoRank [6]         = 0 (expect 0)  PASS
CFOrder(f14) = 14 (expect 14)       PASS
CFOrder(f18) = 18 (expect 18)       PASS
CFOrder(f28) = 7  (expect 7)        PASS  degseq [3,1,2,1] as recorded
curve B: fInt = h3^2 + kappa*q3^3   PASS  (lead f = 1549681956)
curve B: lead(h3)^2+kappa = lead(f) PASS
curve A: fInt = h3^2 + kappa*q3^3   PASS  (lead f = 701706240)
curve A: lead(h3)^2+kappa = lead(f) PASS
```

The curve A/B checks rebuild f from the A(8) chart formulas (A8f + IntModel
exactly as in the skills) and verify the recorded (q3, h3, kappa) of
notes/agent_a2_24_contact.md on the integral models — this simultaneously
validates the named-charts A(8) block, the IntModel scaling convention, and
the contact-torsion validated data, end to end.

---

## SUMMARY

Scope: all 14 skills under skills/, every verbatim-quoted code/math
block diffed token-by-token against its cited source, numeric claims checked
against notes and data logs, cross-references and file paths resolved,
12 Magma spot-checks run (12/12 PASS).

Totals for THIS session (post prior-session fixes of Jul 8):
- CRITICAL (formula differs from source): 0 found anywhere.
- Defects found: 10 (2 wrong numbers, 1 wrong math gloss, 1 overstated
  gloss, 3 imprecise attributions/citations, 3 verbatim-block deviations).
- Fixes applied directly (10 edits, all logged per-skill above):
  1. simplicity-certificates: p^2*T^4*chi(1/T) -> T^4*chi(1/T) gloss.
  2. simplicity-certificates: "no degree drop through n=12" -> "at n=12".
  3. component-boundary-analysis: Rabinowitsch printf string restored.
  4. magma-lab-conventions: ERRATUM/cautionary-note attribution + c4/C2
     source cite.
  5. running-torsion-searches: #HeightRationals(45) 1263 -> 2511 (tied to
     the 6.3M pairs).
  6. running-torsion-searches: d0 fast scan H=15 -> H=45.
  7. running-torsion-searches: HIT24 line attribution -> the data log.
  8. validate-and-record-a-hit: HIT24 line attribution -> the data log.
  9. validate-and-record-a-hit: ERRATUM description untangled (D4 tower vs
     unit-trap retraction).
  10. finite-prefilters: removed non-source inline comment from the
      verbatim ThreeTorsionPrefilter block.
- Minor observations left as-is (5): "v is the constant of q" wording and
  the sign-flip characterization (contact-torsion-constructions);
  "false, 0, 0" vs RT!0 (simplicity-certificates); dropped source
  parenthetical in the 4-point quote (component-boundary-analysis);
  guard-idiom line order (magma-lab-conventions).
- Open item for the orchestrator (1): running-torsion-searches §6 claims
  "~25k pairs/sec dropping to ~5k/sec at height ~100" — no repo artifact
  records these rates (likely a live-session observation). Either keep as
  anecdote or ground it in a future measured run.

Per-skill verdicts: named-charts-reference CLEAN; contact-torsion-
constructions CLEAN (2 minor obs); halving-and-doubling CLEAN;
simplicity-certificates 3 minor (2 fixed); local-obstructions CLEAN;
component-boundary-analysis 2 minor (1 fixed); two-rank-and-factor-types
CLEAN (Magma-checked); magma-lab-conventions 1 minor (fixed);
pell-cf-order CLEAN (Magma-checked); running-torsion-searches 3 defects
(all fixed); validate-and-record-a-hit 2 minor (both fixed);
finite-prefilters 1 minor (fixed); g2-torsion-lab CLEAN (index complete,
14/14 directories); target-playbook CLEAN.
