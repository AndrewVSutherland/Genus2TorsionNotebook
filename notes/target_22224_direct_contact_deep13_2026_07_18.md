# Direct cubic-contact attack on `[2,2,2,24]` (2026-07-18)

## Scope and normalization

This is a fresh attack from rational 3-torsion, independent of the superseded
`Aaux` equations in the first full-family halving search.  Put

\[
 q=x^2+ux+t^2,\qquad
 h=x^3+\frac{1+3u}{2}x^2+wx+t^3,
 \qquad h^2-q^3=xF.
\]

Then `F` is monic quartic, and `[q,h mod q]` is the marked rational 3-torsion
class.  Writing

\[
 F=\prod_{z\in\{a,b,c,d\}}(x+z^2)
\]

gives four exact coefficient equations.  The four full-halving radicands used
throughout are

\[
 abcd,\quad a(a+b)(a+c)(a+d),\quad
 b(b+a)(b+c)(b+d),\quad c(c+a)(c+b)(c+d).
\]

For an open finite-field row the code also verifies, in the finite Jacobian,
that the marked classes have orders 3 and 4 and that the order-4 class belongs
to `2J(F_p)`.  Thus the finite masks are not merely coefficient tests.

## Two corrections that matter

1. The old `target_22224_full_family_halving.m` auxiliary formula was wrong.
   None of its masks or its old “Type I/II” deductions is used here.
2. The first version of the *fresh* p² lifting script accidentally used the
   forward difference `f(x+1)-f(x)` as a Jacobian.  That is not the formal
   derivative for these degree-2/4 equations.  It has been replaced by exact
   degree-four Lagrange differentiation.  The old 600-by-650 lift counts and
   the resulting 74,880,000-combination CRT run are therefore explicitly
   renamed `*_superseded_bad_gradient.*` and are not a certificate.

The finite and boundary enumerations evaluate the equations directly and were
not affected by the second bug.

## Fresh finite and boundary masks

The number of open target parameter triples `(u,t,w)` is:

| p | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| open triples | 0 | 0 | 1 | 2 | 4 | 3 | 8 | 23 | 54 | 76 |

The number of boundary parameter triples admitted by split roots and the four
cover radicands is:

| p | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| boundary triples | 22 | 37 | 46 | 61 | 83 | 112 | 132 | 199 | 237 | 283 |

At `p=17`, the corrected full-family computation maps to the two direct
triples `(3,14,2)` and `(13,2,6)`.  The first is in the open file; the second
has repeated `q` and is in the boundary file.  This is a useful independent
cross-check of the direct chart.

## Correct p² and explicit p⁵ lifting

After formal differentiation, the projected necessary p² masks are much
larger than first reported:

| p | signed boundary incidences | live incidences | base triples | projected p² triples |
|---:|---:|---:|---:|---:|
| 11 | 106 | 106 | 22 | 12,882 |
| 13 | 258 | 258 | 37 | 32,593 |

Consequently, the naive corrected cross with all open masks through `p=31`
would contain about 80.6 billion combinations and was not attempted.

The stronger script keeps explicit square-root variables
`(y0,y1,y2,y3)` and lifts all eight equations.  It exhibits genuine open,
off-rectangle branches modulo p⁵ at both forced primes:

- `p=11`, base `(a,b,c,d;u,t,w)=(0,5,2,8;4,9,10)`, rank 7, two zero
  radicands modulo 11, lifted modulo `11^5=161051`;
- `p=13`, base `(7,7,8,5;5,4,6)`, rank 8, one zero radicand modulo 13,
  lifted modulo `13^5=371293`.

For both certificates, all eight equations vanish modulo p⁵ and every tested
generic-fibre diagnostic is nonzero: `t`, `disc(q)`, branch products and pair
differences, `q(-a_i^2)`, all four square roots, and the three rectangle
differences.  Thus neither 11 nor 13 gives a local obstruction.  The rank-7
`p=11` component has a genuine persistent blow-up direction; superficially
open p² points can die at p³, so explicit roots and higher lifting are essential.

## Deep-branch reconstruction experiments

Eight independently sampled p⁵ centers were retained at 11 and twelve at 13;
the later centers differ already at earlier Hensel digits, rather than being
final-digit changes of the first four.  The p13 centers meet five different
rank-eight mod-13 base incidences.  Their final-digit tangent kernels contain
117,128 explicit open incidence states at 11 and 26,364 at 13,
projecting to respectively 968 and 144 distinct
parameter classes.  Incidence/root data is retained beside every projected
class.  This is exhaustive coverage of those final-digit tangent slices, **not**
coverage of every p⁵ branch.

Crossing these classes with every generic-open combination at
`p=17,19,23,29,31` gives

\[
968\cdot144\cdot(1\cdot2\cdot4\cdot3\cdot8)=26,763,264
\]

CRT combinations from the sampled tangent centers (included in the still
larger expanded run below).  The modulus is
`399365217381454753`, so classical coordinatewise reconstruction used
`|numerator|, denominator <= 446858600`.

The corrected-mask near miss `[50,528,-726,-891]` contributes a genuinely
different 11-adic component.  Its four radicands have exact square roots
`(130680,128180,148104,396396)`.  With the branch ratios fixed, the number of
contact parameter lifts at 11 is `1,11,121,121,1331` at exponents 1 through 5;
1,210 of the p⁵ classes are open.  Adding these gives 2,178 p11 classes in the
expanded bank.  Crossing them with the 144 p13 classes and all 192 open-good
combinations exhausts 60,217,344 CRT combinations: 4,369,680 reconstruct;
349,200 / 26,173 / 1,995 survive the p37 / p41 / p43 filters, and again
**none** of the 1,995 has an exact split quartic.

The same near miss is sharply obstructed at 13.  Introducing the necessary
projective scale gives the two mod-13 seeds `(lambda,u,t,w)=(6,2,1,10)` and
`(7,2,1,10)`, with rank 3.  The fixed-line counts are `2,26,338,338,0`
through exponents 1--5.  Moreover all 338 fixed p⁴ states fail the full
11-variable p⁵ linear lift, even when the final branch and square-root digits
move transversely.  Moving transversely one digit earlier, only 26 of the 338
p³ states reach p⁴, and 15,600 sampled five-dimensional directions produced
no p⁵ continuation.  The last sampling statement is heuristic, but the
fixed-line and fixed-p⁴ deaths are exhaustive.

A second open-prime near miss, `[144,697,-722,-1394]`, is easier to dispose
of.  Its radicands are exact squares, but modulo 13 its signed branches are
`(1,8,6,10)` and their squares are `(1,12,10,9)`, all nonzero and distinct.
It therefore has smooth reduction at 13.  The corrected finite enumeration
has no smooth target triple at 13, so this seed cannot generate a new
13-integral contact component.  Its mod-11 collision `(1,4,4,3)` is irrelevant
to that clean p13 exclusion.

Finally, coherent reconstruction was tested without discarding branch data.
The 2,178 p11 incidence rows and 144 p13 incidence rows give 313,632 deep
pairs modulo `11^5*13^5`; 2,403 reconstruct all seven coordinates
`(a,b,c,d,u,t,w)` within the classical bound 172,911, and none satisfies the
four exact contact coefficient identities.  Adding every coherent p17 cover
presentation raises the modulus to `1016550852031` and the bound to 712,934;
among 1,881,792 combinations, 11,594 reconstruct all seven coordinates and
again none satisfies the exact identities.  Adding all coherent p17 and p19
presentations tests 22,581,504 combinations at modulus `19314466188589` and
bound 3,107,608; 121,931 reconstruct all seven coordinates and **zero** obey
the exact contact identities.  This directly targets rational
split quartics and is complementary to reconstructing only `(u,t,w)`.

The later F4/F5/F7 boundary-looking seeds give no new p-adic component.  An
exhaustive search over projective scale and `(u,t,w)` finds no direct-contact
incidence even modulo 13 for any of
`[-126,28,49,-50]`, its partner with `d=-578`,
`[-112,14,49,-338]`, or `[-18,1,16,-1682]`.  They do have the expected two
scaled incidences modulo 11, but the mod-13 failure occurs before any Hensel
question.  The exact data is in
`results/target_22224_direct_contact_newseed_modp_diagnose.log`.

A simultaneous/common-denominator LLL search was also run:

- On the initial `484*52=25,168` deep pairs modulo `11^5*13^5`, 40 reduced-basis
  combinations per lattice gave 866,880 unique projective rational triples;
  none survived the masks through `p=31`.
- With the open `p=17,19` classes built into the lattice modulus, 50,336
  lattices gave 1,711,698 unique triples; none survived through `p=41`.

Earlier random CRT runs remain useful only as heuristic negative searches.
In particular 375,657 exact rational reconstructions had no split quartic,
but their p² mask source used the bad forward-difference Jacobian, so they must
not be presented as bounded completeness evidence.  The record-centered run,
which fixed the exact `(0,0,0)` classes at 11² and 13² and used masks through
43, produced 32,627 exact rational candidates and again no split quartic.

## Interpretation and next steps

The corrected computation rules out a substantial, precisely described set of
deep local tubes at a large reconstruction height, but it does not rule out
`[2,2,2,24]`.  The next useful expansion is to sample earlier Hensel digits on
the rank-7 11-adic blow-up and the rank-8 13-adic component, rather than merely
adding final digits to the same centers.  Denominator-nonunit charts at 11 or
13 are also outside the present affine reconstruction domain and require a
separate projective normalization.

Finally, the attractive rational cover surface
`(a,b,c,d)=(1,r,s,rs)` lies identically on the multiplicative-rectangle locus
`ad=bc`; its associated curves are bielliptic and therefore unsuitable for a
geometrically simple record.  Exact hits in any broader search should reject
all three pair-product equalities.

## Main files

- `code/target_22224_direct_contact_deep13_finite.m`
- `code/target_22224_direct_contact_deep13_boundary.m`
- `code/target_22224_direct_contact_deep13_lift.py`
- `code/target_22224_direct_contact_deep13_padic_lift.py`
- `code/target_22224_direct_contact_deep13_padic_tangent_masks.py`
- `code/target_22224_direct_contact_deep13_padic_tube_crt.py`
- `code/target_22224_direct_contact_deep13_padic_lll.m`
- `code/target_22224_direct_contact_deep13_coherent_crt.py`
- `code/target_22224_direct_contact_nearmiss_fixedbranch_lift.py`
- `code/target_22224_direct_contact_nearmiss_transverse_p13.py`
- `code/target_22224_direct_contact_deep13_verify.m`
