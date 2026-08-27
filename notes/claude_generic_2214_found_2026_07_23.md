# Generic `[2,2,14]` realized: five End=Z curves (2026-07-23)

**Claim.** The torsion group `[2,2,14]` (order 56) is realized by geometrically
simple genus-2 Jacobians over Q with `End(J_Qbar) = Z`.  Five verified curves,
from the three-root route on the contact-7 chart
(`notes/claude_generic_222_2214_plan_2026_07_23.md`); this closes one of the
two `^RM` flags in `paper/torsion_realizations.tex`.  No torsion group
containing `(Z/2)^2 x Z/7` appears anywhere in the published literature for
ANY endomorphism type (literature dossier in the plan note), so these are new
beyond-RM and beyond-literature.

## Construction

Contact-7 chart: `h = 1-(7/2)x+ax^2+bx^3`, `f = (h^2+(x-1)^7)/x^2` (monic
quintic), marked order-7 class `[(1,h(1)) - inf]`.  Roots at `r_i = 1-v_i^2`
for `v_i in {s,t,u}` with the linear conditions `a + b(1-v_i^2) = A(v_i)`,
`A(v) = (2v^5+4v^4+6v^3+8v^2+10v+5)/(2(v+1)^2)`; `(s,t)` determine `(a,b)`,
and `u` is a root of the residual cubic `R(s,t,u)` (the three-root surface).
Factor type `[1,1,1,2]` (three finite rational Weierstrass points + infinity)
gives 2-rank 3; with the marked 7-class, torsion contains
`(Z/2)^3 x Z/7 = [2,2,14]`.

## The five generic curves (all: exact TorsionSubgroup = [2,2,14], marked order 7)

Verification: `code/claude_2214_threeroot_verify.m` (rebuilds each curve from
`(s,t,u)` from scratch, asserts the three roots, exact `TorsionSubgroup`,
marked-class order, End=Z certificate), log
`results/claude_2214_threeroot_verify.log`.  End=Z certificate = L-polynomial
at `p0` irreducible with all Frobenius-root powers of degree 4 (root-power =>
geometrically simple), second prime `q0` ALSO root-power strict (absolute
simplicity of the second reduction, so End^0 of its geometric reduction equals
Q(pi_{q0}) — upgraded after Codex's PR-8 review caught that an even chi like
x^4-2x^2+1681 at q0 leaves an M_2-type reduction where RM could hide) with
multiplicative splitting-field degrees (=> center Q, kills RM), torsion
`56 > 18` (=> not QM, by Laga–Schembri–Shnidman–Voight).  Cross-checks in PARI (independent):
`56 | #J(F_p)` at all good `p <= 59`, real-subfield-disc scatter (End=Z
signature) for all five.

```text
(s,t,u) = (-10, -1/2, -10/7)   p0=17, q0=31
  f = x^5 + 1136785381/11573604 x^4 - 75084797/964467 x^3 - 53745313/642978 x^2 + 151804/1701 x - 41701/2268
(s,t,u) = (-5, -15/8, -15/22)  p0=17, q0=83
  f = x^5 + 463110041/18593344 x^4 + 691213617/37186688 x^3 - 22161775159/297493504 x^2 + 824317/17248 x - 79143/8624
(s,t,u) = (-1/2, -15/8, -15/19) p0=17, q0=29
  f = x^5 + 1168230281/55472704 x^4 + 1039464297/110945408 x^3 - 62905712119/887563264 x^2 + 1587293/29792 x - 163047/14896
(s,t,u) = (-4/9, 4/17, -4/25)  p0=31, q0=43
  f = x^5 - 650775201299/140512522500 x^4 + 301045120409/35128130625 x^3 - 556345299313/70256261250 x^2 + 684788/187425 x - 504803/749700
(s,t,u) = (4/17, -5/18, -10/49) p0=43, q0=59
  f = x^5 - 4993712661259139/1085163373145124 x^4 + 176638001774783/20868526406637 x^3 - 24972622480333/3210542524098 x^2 + 58817786/16470909 x - 3320579/5067972
```

(Integral models `y^2 = den^2 f` and the certificate chi-pairs are in the
verification log.)

## Minimal models and conductors (same day; code/claude_2214_minimal_models.m + gp genus2red)

> **ERRATUM (2026-07-25).  Every `N` printed in the code block below is only
> the ODD PART of the conductor.**  They were computed with PARI `genus2red`,
> which is documented for `p > 2`; it returned the sentinel exponent `-1` at 2
> for all six curves and then omitted the 2-part from `N`
> (`results/claude_ov_erratum_sentinel.log`).  The corrected table is
> immediately below the block; the derivation, the validation of the
> replacement method and the full repository-wide sweep are in
> `notes/claude_ov_erratum_2026_07_25.md`.  The claim "all conductors are odd
> non-squares (no GL2-type smell)" was an ARTEFACT of the missing 2-part and
> has been struck.

All five re-verified (exact torsion + End=Z) on their reduced minimal
Weierstrass models:

```text
(-10,-1/2,-10/7):   y²+(x²+x)y = 9x⁶-1005x⁵+24243x⁴+304602x³-4244115x²-42264735x-37741275
                    odd(N) = 3²·5²·7²·13²·19·23·37·317 = 9550095752925   <- SMALLEST: displayed in torsion_realizations.tex
(-5,-15/8,-15/22):  y²+(x²+x)y = 900x⁶+6630x⁵-21900x⁴-167639x³+289767x²+1159725x-2056775
                    odd(N) = 3²·5²·7²·11²·19·59·67·179 = 17934836205825
(-1/2,-15/8,-15/19): y²+(x²+x)y = 4x⁶+1770x⁵+209208x⁴+1087252x³-357136500x²+5573377890x-16704297900
                    odd(N) = 3²·5²·7²·11²·13²·19²·31·59·97 = 14439206077220925
(-4/9,4/17,-4/25):  y²+(x²+1)y = 16x⁶-732x⁵-29172x⁴+989439x³+18286616x²-114265515x+10752131
                    odd(N) = 3²·5²·7²·17²·13·23·97·229 = 21161909161575
(4/17,-5/18,-10/49): y²+(x²+x)y = 16900x⁶+496470x⁵-2676360x⁴-102391553x³+295524693x²+3295748145x-6440348871
                    odd(N) = 3²·5²·7²·17²·13·31·37·61·71·157·1153 = 37247776656021702225
```

**Corrected conductors** (Magma `Conductor` on the reduced minimal Weierstrass
model; "Ogg" = whether `v_2(disc_min) < 12`, the regime in which Magma issues
no warning):

```text
(s,t,u)                v_2(disc_min)  Ogg   N
(-10,-1/2,-10/7)            10        YES   2^2 * 9550095752925        = 38200383011700
(-5,-15/8,-15/22)           14         no   2^2 * 17934836205825       = 71739344823300
(-1/2,-15/8,-15/19)         14         no   2^2 * 14439206077220925    = 57756824308883700
(-4/9,4/17,-4/25)           30         no   2^5 * 21161909161575       = 677181093170400
(4/17,-5/18,-10/49)         30         no   2-part NOT DETERMINED (Magma Conductor OOMs at 48 GB)
(-3,-3/4,-3/5) [RM]          7        YES   2^2 * 38025 = 152100 = 390^2
```

All six are EVEN, and the RM one is a perfect square (390^2) — the GL2-type
signature the original note said was absent.  The ordering is unchanged:
`(-10,-1/2,-10/7)` is still the smallest.  The [2,2,14]
certificate now heads `data/claude_endz_certificates.txt` (regenerated, and
carrying the corrected `N`);
`paper/torsion_realizations.tex` row + preamble updated (RM flag dropped,
displayed witness = the conductor-38200383011700 curve) and the PDF rebuilt;
generator `code/claude_gen_torsion_table.py` updated to match.

## The sixth curve: a certified RM(sqrt 2) witness at tiny height — of the KNOWN conductor

`(s,t,u) = (-3, -3/4, -3/5)`:
`f = x^5 + 4769/400 x^4 + 9009/400 x^3 - 104671/1600 x^2 + 1699/40 x - 42/5`,
torsion exactly `[2,2,14]`, marked order 7, but constant real-subfield disc
`{2}` at ~30 good primes (RM-by-Q(sqrt2) signature, like J_0(29)) and no
disjoint-splitting-fields pair.

> **CORRECTION (2026-07-25, Lane 7).**  This note called it "an RM witness far
> smaller than 152100.eb.2".  That is true of the HEIGHT and false of the
> CONDUCTOR.  The RM is now certified: `J` is isogenous over `Q` to the modular
> abelian surface `A_f` of the level-390 weight-2 newform with Hecke field
> `Q(sqrt 2)` (L-polynomial identity at 135 of 135 good primes `p < 800`), so
> `cond(J) = 390^2 = 152100` — **exactly the conductor of 152100.eb.2**, and
> the level-390 `Q(sqrt2)` orbit is a single Galois orbit, hence a single
> isogeny class.  This curve is therefore a small-height model of the KNOWN
> class, not a new RM witness.  The `38025` that the original note would have
> obtained from `genus2red` is the odd part only.  Full argument:
> `notes/claude_ov_lane7_2026_07_25.md` §4.

## Caveats / remaining hygiene

- The five curves are CURVES, not a family claim.  The three-root surface has
  6 known orbits at height <= 32 with common-numerator structure; fitting the
  rational curve(s) through them (M(24)-component playbook) is the natural
  next step toward an infinite generic family.
- DONE same day: minimal models + conductors (section above), certificates
  data file regenerated (with the strict-q0 upgrade), and the
  `torsion_realizations.tex` row/preamble + PDF update.  Still open:
  `paper/torsion_sources.tex` should eventually gain this note as the [2,2,14]
  generic source, and the RM(sqrt2) sixth curve needs formal RM certification.

## Files

- `code/claude_2214_threeroot_sweep.gp`, `code/claude_2214_threeroot_surface.gp` (construction + surface)
- `code/claude_2214_threeroot_verify.m` + `results/claude_2214_threeroot_verify.log` (verification)
- `results/claude_2214_threeroot_h16.log`, `_h32.log` (sweeps)
- plan + literature dossier: `notes/claude_generic_222_2214_plan_2026_07_23.md`
