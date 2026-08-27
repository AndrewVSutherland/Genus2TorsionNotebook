# Repeated-fiber Mordell--Weil sieves for `[2,2,2,24]`

Date: 2026-07-18

## Setup

For a signed full-cover tuple with three fixed entries and a square-scaled
fourth entry,

```text
(z1,z2,z3,d0*T^2),
```

the three normalized square conditions are

```text
Ri(T) = (zi + d0*T^2)/(zi + d0),  i=1,2,3.
```

Each pair `Ri*Rj` is a genus-one quartic.  Its Mordell--Weil group supplies
all rational points satisfying two square conditions; the remaining square
condition is imposed locally and exactly.  At smooth good primes we also
require `3 | #J(F_p)`.  At `p=11,13`, a singular reduction is retained only
when its projective root multiset occurs in the corrected direct-contact
boundary bank.  Surviving coefficient classes are finally matched against
the explicit `p^5` tangent-tube banks.

The corrected projective root key quotients permutations, independent signs,
and common scaling of the four branches.  A mere repeated-root collision is
not treated as direct-contact incidence.

## Fiber and quotient ranks

The first seven new B=2000 repeated fibers have pairwise elliptic ranks:

| fiber | fixed triple; `d0`, second `T` | ranks `(12,13,23)` |
|---|---|---|
| F1 | `(-1470,-630,336); 25, 5` | `(2,2,3)` |
| F2 | `(-720,20,300); -363, 17/11` | `(2,2,2)` |
| F3 | `(-612,34,289); -338, 25/13` | `(3,3,3)` |
| F4 | `(-126,28,49); -50, 17/5` | `(2,2,2)` |
| F5 | `(-112,14,49); -50, 13/5` | `(2,2,2)` |
| F6 | `(-50,30,45); -48, 2` | `(2,2,2)` |
| F7 | `(-18,1,16); -50, 29/5` | `(3,2,2)` |

All 21 ranks are proved exactly.  The rank script and log are:

```text
code/target_22224_repeated_fibers_quotient_ranks.m
results/target_22224_repeated_fibers_quotient_ranks.log
```

The P8 shared-triple fiber is

```text
(528,-726,-891,50*T^2),
T = 1, 19/5.
```

All three elliptic quotients have rank 2 and torsion `[2]`.

## Completed coefficient-lattice exclusions

### P8

The direct rank-2 lattice walk through coefficient radius 20 tested 10,086
MW points and 5,396 distinct `T` values.  The only full-cover lifts were
`T=+-1,+-19/5`; neither survives the necessary 3-torsion gcd bound.

The corrected periodic sieve was subsequently pushed to `N=10000`, covering
800,080,002 raw coefficient/torsion triples.  It enumerates only 31,590
periodic points and leaves 29 p-level survivors.  Direct `Q_13` evaluation
gives 28 finite points and one infinity point; none matches a `13^5` tangent
tube.

```text
results/target_22224_p8_mw_lattice_B20.log
results/target_22224_P8_rank2_mwsieve_N10000.log
results/target_22224_P8_rank2_modular_N10000_p13deep.log
```

### F4

Fiber:

```text
(-126,28,49,-50*T^2), quotient <1,2>, rank 2.
```

At `p=13`, the 20 allowed MW coefficient residues all map to `T=infinity`;
there are no finite projective direct-contact classes.  The saturated box
`|m|,|n| <= 1001` contains 8,024,018 coefficient/torsion triples.  Nineteen
local profiles leave one class, whose exact lift fails the full cover.

```text
results/target_22224_F4_rank2_mwsieve_12_N1001.log
```

### N2

Fiber and quotient:

```text
(-25,-1,40,(4205/98)*T^2), T=1,49/29,
quotient <1,3>, rank 2, torsion [2].
```

The saturated `N=1001` box leaves two modular classes and no exact full-cover
lift.  Periodic CRT enumeration then covers `N=10000`, representing
800,080,002 raw coefficient/torsion triples.  It enumerates 2,930,697
periodic lattice points and leaves 182 local survivors.  Direct `Q_13`
evaluation gives 181 finite points (`T=1 mod 13`) and one infinity point;
none matches any of the 144 canonical `13^5` tangent profiles.

```text
results/target_22224_N2_rank2_mwsieve_N10000.log
results/target_22224_N2_rank2_modular_N10000_p13.log
```

### N9 (rank 3)

Fiber and quotient:

```text
(-46,47,49,(-1081/50)*T^2), T=1,5/7,
quotient <1,2>, rank 3, torsion [2].
```

The Magma `N=100` cube tested 16,241,202 coefficient/torsion triples; seven
survive 16 local profiles and all seven exact lifts fail the full cover.  A
dedicated C++ periodic engine reproduces exactly the same seven classes.

At `N=3000`, the raw cube has `2*6001^3` coefficient/torsion triples.  The
periodic engine enumerates 316,660,625 admissible `k` values; 219,487 survive
the p-level profiles.  Direct `Q_13` evaluation yields 219,486 finite points
and one infinity point.  None matches a `13^5` tangent tube.

```text
code/target_22224_N9_rank3_mwsieve.m
code/target_22224_N9_rank3_periodic.cpp
code/target_22224_N9_deep_padic_filter.m
results/target_22224_N9_rank3_periodic_N3000_p13deep.log
```

### M5

Fiber and quotient:

```text
(99,244,4026,6*T^2), T=1,29,
quotient <1,2>, rank 2, torsion [2].
```

The exact saturated `N=1001` box leaves two modular classes and no exact
full-cover lift.  Periodic `N=10000` covers 800,080,002 raw triples,
enumerates 3,572,323, and leaves 188 local survivors.  Of these, 187 are
finite at 13 and one is at infinity; none matches a `13^5` tube.

```text
results/target_22224_M5_rank2_mwsieve_12_N1001.log
results/target_22224_M5_rank2_mwsieve_N10000.log
results/target_22224_M5_rank2_modular_N10000_p13deep.log
```

### O5

Fiber and quotient:

```text
(-15,-10,16,(50/3)*T^2), T=1,5/3,
quotient <1,2>, rank 2, torsion [2].
```

This is a finite-boundary control rather than an infinity-only lane.  The
exact `N=1001` box has no full-cover lift among its final classes.  Periodic
`N=10000` enumerates 22,915,834 classes from 800,080,002 raw triples and
leaves 68 p-level survivors; all 68 are finite at 13 and none matches a
`13^5` tube.

```text
results/target_22224_O5_rank2_mwsieve_N10000.log
results/target_22224_O5_rank2_modular_N10000_p13deep.log
```

### O6 (torsion-coset control)

Fiber and quotient:

```text
(-15,-1,25,(605/27)*T^2), T=1,21/11,
quotient <1,2>, rank 2, torsion [4].
```

Periodic `N=10000` covers 1,600,160,004 raw coefficient/torsion triples
because there are four rational torsion cosets.  It enumerates 44,634,294
periodic points and leaves 308 p-level survivors.  Over `Q_13`, 161 are
finite and 147 lie in the infinity chart; none matches a `13^5` tube.

```text
results/target_22224_O6_rank2_mwsieve_N10000.log
results/target_22224_O6_rank2_modular_N10000_p13deep.log
```

### R1 (new rank-one quotient)

The live B=10000 stream at `a` about 779 produced

```text
(-1,2,25,(-529/338)*T^2), T=1,13/23.
```

Its pair `<1,2>` quotient is the first rank-one quotient in this stream:

```text
E: y^2 = x^3 - x^2 + 31*x + 33,
E(Q) = <(-1,0)> + <(1,8)>.
```

The full range `|m| <= 1,000,000` across both torsion cosets consists of
4,000,002 coefficient classes.  Twenty local profiles leave 13 classes.
Over `Q_13`, 12 are finite and one lies at infinity; none matches a `13^5`
tangent tube.

```text
results/target_22224_R1_rank2_mwsieve_12_N1000000.log
results/target_22224_R1_rank1_N1000000_p13deep.log
```

### R6

The same live stream produced

```text
(-3,5,75,(-405/121)*T^2), T=1,11/21.
```

The pair `<1,3>` quotient is

```text
E: y^2 = x^3 - 207*x + 8856,
rank 2, torsion [2].
```

The exact saturated `N=1001` box leaves one local class and no full-cover
lift.  A periodic `N=10000` run covers 800,080,002 raw coefficient/torsion
triples, enumerates 2,852,338 periodic lattice points, and leaves 34 modular
survivors.  The sole exactly reconstructed full-cover specialization in the
central range is `T=1`, whose prime-order gcd is 448 and hence has no
3-torsion.  Over `Q_13`, 33 survivors are finite and one lies at infinity;
none matches a corrected `13^5` direct-contact tube.

```text
results/target_22224_R6_rank2_mwsieve_13_N1001.log
results/target_22224_R6_rank2_mwsieve_N10000.log
results/target_22224_R6_rank2_modular_N10000_p13deep.log
```

### S1 (completed-stream rank-one quotient)

Completing the B=10000 stream through `a=1000` revealed three more repeated
fibers.  The most useful is

```text
(-4,9,30,(-2166/245)*T^2), T=1,7/19.
```

Its pair `<1,2>` quotient has proven rank one:

```text
E: y^2 = x^3 + x^2 + 3392*x + 62288,
E(Q) = <(-17,0)> + <(152,-2028)>.
```

The direct `|m| <= 1,000,000` sieve across both torsion cosets leaves 26
classes; none lifts to the corrected `13^5` contact bank.  A generalized-CRT
run then extends the coefficient range to `|m| <= 100,000,000`: from
400,000,002 raw coefficient/torsion classes it enumerates 2,051,283 periodic
points and leaves 3,412 modular survivors.  Over `Q_13`, 3,411 are finite and
one lies at infinity; again none matches a corrected `13^5` tube.  This is
the strongest bounded exclusion for any repeated fiber in this lane.

```text
results/target_22224_B10000_a1000_newfiber_ranks.log
results/target_22224_S1_rank1_mwsieve_N100000000.log
results/target_22224_S1_rank1_N100000000_p13deep.log
```

### S3 (completed-stream rank-three quotient)

The final genuinely new completed-stream fiber is

```text
(-9,35,65,(-3025/91)*T^2), T=1,5/11.
```

All three pair quotients have proven rank 3 and torsion `[2]`.  The pair
`<1,3>` model

```text
y^2 = x^3 - x^2 + 67007092*x + 1181772362562
```

gave the best local profile.  Nineteen usable prime profiles were exported
and fed to the compiled periodic rank-3 enumerator.  At `N=3000` this covers
more than 432 billion raw coefficient/torsion triples.  The CRT engine tests
134,139,967 periodic `k` representatives and leaves 979 modular survivors.
Over `Q_11`, 594 are finite and 385 lie in a pole/infinity chart; none matches
the corrected `11^5` direct-contact bank.

```text
results/target_22224_S3_rank3_fullprofile_13.log
results/target_22224_S3_rank3_periodic_N3000.tsv
results/target_22224_S3_rank3_periodic_N3000_p11deep.log
```

### U1/U2 (first fibers beyond `a=1000`)

The continuing B=10000 stream found two reciprocal presentations of the
same new repeated-fiber geometry:

```text
U1 = (-119,124,126,(-1054/9)*T^2), T=1,3/7,
U2 = (-1071,-1054,1116,1134*T^2), T=1,7/3.
```

Their quotient ranks are the same up to permutation: `(2,4,2)`.  For U2,
pair `<1,3>` is rank 2 with torsion `[2]` and minimal model

```text
y^2 = x^3 + x^2 + 1475392*x + 25086288.
```

At 13, only map-pole classes pass the projective contact mask.  The exact
`N=1001` box leaves two local classes and no full-cover lift.  Periodic
`N=10000` covers 800,080,002 raw triples, enumerates 7,145,002 lattice
points, and leaves 281 modular survivors.  A deterministic saturated basis
was fixed in the p-adic filter (Magma otherwise returned a different
unimodular basis for this model).  In that matching basis, 280 survivors are
finite over `Q_13` and one is a pole; none matches a corrected `13^5` tube.

**Infinity-chart resolution.**  A complete scan of both U1/U2 projective
`T`-lines first showed that the current 144-key tangent bank omits the
infinity component.  A dedicated lift then found the missing contact branch:
two smooth unit-`L` points modulo 13 lift through `13^5` to all 28,561
infinity `z=1/T` classes.  However, on the first blow-up the four full-cover
radicands are `z^6` times units reducing to `[9,7,11,4]`; the middle two are
nonsquares modulo 13.  Hence the actual A(2,2,2,8) cover has no U2 infinity
branch, independently of the sampled bank.  Exact recovery also showed that
all 281 MW survivors are finite (one apparent pole was p-adic precision loss,
with exact `T=664293/527543`).  This turns the earlier coverage warning into a
genuine local exclusion of the only projected contact lane.

```text
results/target_22224_B10000_a1075_newfiber_ranks.log
results/target_22224_U2_rank2_mwsieve_13_N1001.log
results/target_22224_U2_rank2_mwsieve_N10000_fast.log
results/target_22224_U2_rank2_modular_N10000_basisA_p13deep.log
results/target_22224_U2_global13_P1_scan.log
results/target_22224_U2_infinity_contact_lift.log
results/target_22224_U2_infinity_MW_check.log
```

### V1/V2 (continuing `a>1000` stream)

Two further repeated fibers appeared by `a=1524`:

```text
V1 = (-18,20,75,(-1470/121)*T^2), T=1,11/49,
V2 = (-9,14,49,(-2023/162)*T^2), T=1,3/17.
```

All six pair quotients have proven rank 2 and torsion `[2]`.  For V1 pair
`<1,2>`, the exact `N=1001` coefficient box (8,024,018 raw triples) is cut to
one local class by the prime profiles, and that class does not lift to the
full cover.  A periodic `N=10000` extension covers 800,080,002 raw classes
and leaves 33 finite-prime survivors.  A complete intrinsic depth-5 V1 mask
has 59,164 `T` keys: full cylinders over `T=+/-5` and 1,021 ramified keys in
each `T=+/-1` disk.  Exactly three MW classes meet this full cover/contact
mask.  Auxiliary prime 103 kills `(4900,-1960,ti=1)`, prime 101 kills the
identity, and `(-2,-1,1)` has exact
`T=-4574693/21848987` with a nonsquare third cover ratio.  Thus no
`N<=10000` V1 class remains, without any sampled-bank assumption.

For V2 pair `<1,3>`, six classes survive the primes through 97;
prime 101 kills five immediately, and exact reconstruction of the last
`(m,n,ti)=(3,-1,1)` gives
`T=-2063497452257/1698222432003`, for which the middle square condition
fails.  Thus both `N=1001` boxes are excluded without using the sampled
tangent bank.

```text
results/target_22224_B10000_a1524_newfiber_ranks.log
results/target_22224_V1_rank2_mwsieve_12_N1001.log
results/target_22224_V1_rank2_mwsieve_N10000.log
results/target_22224_V1_component_auxprime_filter.log
results/target_22224_V2_survivor_prime_filter.log
```

### W1/W2 (continuing `a>2000` stream)

By `a=2250` the next tranche produced

```text
W1 = (-2254,-2162,2303,4900*T^2), T=1,7/5,
W2 = (-8,9,25,(-2209/338)*T^2), T=1,13/47.
```

W1 has quotient ranks `(3,3,3)` with torsion `[2]`.  On its pair `<1,3>`
quotient, 33 usable prime profiles through 167 reduce the full rank-3
`N=1000` box (more than 16 billion raw triples) to the single coefficient
class `(2,-1,1,ti=2)`.  Exact reconstruction gives
`T=-6015613/4711271` and square flags `[true,false,true]`, so it does not lift
to the full cover.  W2 has ranks `(2,3,2)` and quotient torsion `[4]`, `[2]`,
`[2,2]`, respectively.  On the W2 pair
`<1,2>` quotient, the four torsion cosets make the `N=1001` box contain
16,048,036 raw triples.  Three classes survive the primes through 97; prime
101 kills two and prime 109 kills the last.  Thus this box is excluded using
only finite-prime necessary conditions, with no tangent-bank assumption.

```text
results/target_22224_B10000_a2250_newfiber_ranks.log
results/target_22224_W1_rank3_periodic_N1000_p167.tsv
results/target_22224_W1_exact_survivor.log
results/target_22224_W2_rank2_mwsieve_12_N1001.log
results/target_22224_W2_survivor_prime_filter.log
```

### X1

The same tranche next produced

```text
X1 = (-49,-21,81,(6069/121)*T^2), T=1,143/17.
```

All three pair quotients have proven rank 2 and torsion `[2]`.  Pair `<1,3>`
has minimal model `y^2=x^3+x^2+3884*x+43934`.  Its exact `N=1001` box is
reduced from 8,024,018 raw triples to one class by the finite-prime profiles,
and that class does not lift to the full cover.

```text
results/target_22224_B10000_a2650_newfiber_ranks.log
results/target_22224_X1_rank2_mwsieve_13_N1001.log
```

### Y1

The final new repeated fiber through `a=3000` is

```text
Y1 = (-15,28,50,(-5290/189)*T^2), T=1,3/23.
```

All quotient ranks are 2 with torsion `[2]`.  Pair `<1,3>` has minimal model
`y^2=x^3+4571700*x+5421458000`; its exact `N=1001` box again leaves one
finite-prime local class, which fails to lift to the full cover.

```text
results/target_22224_B10000_a3000_newfiber_ranks.log
results/target_22224_Y1_rank2_mwsieve_13_N1001.log
```

### Z1

The `3001<=a<=4000` tranche adds

```text
Z1 = (-169,-120,225,(5415/32)*T^2), T=1,26/19.
```

Its quotient ranks are `(2,2,3)`, with torsion `[4]`, `[2]`, `[2]`.
Pair `<1,3>` has minimal model
`y^2=x^3+x^2+1850420*x+244635950`; the exact `N=1001` box leaves one local
class and no full-cover lift.

```text
results/target_22224_B10000_a4000_newfiber_ranks.log
results/target_22224_Z1_rank2_mwsieve_13_N1001.log
```

### AA1

The first half-range endpoint adds

```text
AA1 = (-64,-40,65,(9800/117)*T^2), T=1,9/7.
```

All three quotient ranks are 2 with torsion `[2]`.  On pair `<1,3>`, five
classes survive the primes through 97 in the `N=1001` box.  Prime 101 kills
four; exact reconstruction of `(1,2,ti=2)` gives
`T=636151141061/155196346323` with square flags `[true,false,true]`.  Hence
the box contains no full-cover lift.

```text
results/target_22224_B10000_a5000_newfiber_ranks.log
results/target_22224_AA1_survivor_prime_filter.log
```

## Direct-box extension through `a=2000`

The B=10000 transverse search was extended from `a<=1000` to
`1001<=a<=2000`.  It tested 36,149,668,500 ordered triples, 538,190,783
squarefree-kernel `k` candidates, and 4,305,526,264 signed presentations.
There were 1,531 full-cover presentations and six non-rectangle local-contact
survivors, but no match in either sampled deep bank.  Exact verification of
all six survivors gave geometrically simple Jacobians with torsion exactly
`[2,2,2,8]` in every case.  Repeated-fiber mining increased from 42 to 48
square-ratio pairs; after removing two scaled known fibers, the genuinely new
families are precisely U1/U2 and V1/V2 above.

```text
results/target_22224_a2228_deep_p5_box_B10000_a1001_2000_allsigns.tsv
results/target_22224_a2228_deep_p5_box_B10000_a1001_2000_verify.log
results/target_22224_a2228_B10000_a2000_repeated_fibers.txt
```

The next tranche `2001<=a<=3000` tested a further 28,151,668,500 ordered
triples and 1,985,265,800 signed presentations.  Its four local-contact
survivors are all geometrically simple with exact torsion `[2,2,2,8]`.
Repeated-fiber mining adds W1, W2, X1, and Y1 (plus one scaled known family),
bringing the combined B=10000 stream to 53 square-ratio pairs.

```text
results/target_22224_a2228_deep_p5_box_B10000_a2001_3000_allsigns.tsv
results/target_22224_a2228_deep_p5_box_B10000_a2001_3000_verify.log
results/target_22224_a2228_B10000_a3000_repeated_fibers.txt
```

The `3001<=a<=4000` tranche tests 21,153,668,500 more ordered triples and
968,787,256 signed presentations.  Both local-contact survivors are simple
with exact torsion `[2,2,2,8]`.  It adds the single new repeated family Z1,
bringing the combined stream to 54 square-ratio pairs.

```text
results/target_22224_a2228_deep_p5_box_B10000_a3001_4000_allsigns.tsv
results/target_22224_a2228_deep_p5_box_B10000_a3001_4000_verify.log
results/target_22224_a2228_B10000_a4000_repeated_fibers.txt
```

The `4001<=a<=5000` tranche tests 15,155,668,500 ordered triples and
464,606,240 signed presentations.  Its sole local-contact survivor is simple
with exact torsion `[2,2,2,8]`.  AA1 is the only new repeated family, so the
half-range combined stream contains 55 square-ratio pairs.

```text
results/target_22224_a2228_deep_p5_box_B10000_a4001_5000_allsigns.tsv
results/target_22224_a2228_deep_p5_box_B10000_a4001_5000_verify.log
results/target_22224_a2228_B10000_a5000_repeated_fibers.txt
```

The `5001<=a<=6000` tranche tests 10,157,668,500 ordered triples,
26,100,279 squarefree-kernel candidates, and 208,802,232 signed
presentations.  Of 250 full prime-filter hits, the only non-rectangle local
survivor is `[5301,5808,-5952,-7600]`; exact Magma verification gives a
geometrically simple Jacobian with torsion `[2,2,2,8]`.  No new repeated
square-ratio fiber occurs, so the cumulative stream still contains 55 pairs.

```text
results/target_22224_a2228_deep_p5_box_B10000_a5001_6000_allsigns.tsv
results/target_22224_a2228_deep_p5_box_B10000_a5001_6000_verify.log
results/target_22224_a2228_B10000_a6000_repeated_fibers.txt
```

The `6001<=a<=7000` tranche tests 6,159,668,500 ordered triples,
10,337,938 squarefree-kernel candidates, and 82,703,504 signed
presentations.  There are 109 full prime-filter hits and three
off-rectangle hits, but no local-contact survivor, hence no exact torsion
candidate to verify.  Repeated-fiber mining remains unchanged at 55 pairs.

```text
results/target_22224_a2228_deep_p5_box_B10000_a6001_7000_allsigns.tsv
results/target_22224_a2228_deep_p5_box_B10000_a6001_7000_verify.log
results/target_22224_a2228_B10000_a7000_repeated_fibers.txt
```

The `7001<=a<=8000` tranche tests 3,161,668,500 ordered triples,
3,129,113 squarefree-kernel candidates, and 25,032,904 signed
presentations.  Its 59 full prime-filter hits include three off-rectangle
hits but no local-contact survivor.  Exact verification is therefore empty,
and the repeated-fiber count remains 55.

```text
results/target_22224_a2228_deep_p5_box_B10000_a7001_8000_allsigns.tsv
results/target_22224_a2228_deep_p5_box_B10000_a7001_8000_verify.log
results/target_22224_a2228_B10000_a8000_repeated_fibers.txt
```

The `8001<=a<=9000` tranche tests 1,163,668,500 ordered triples, 609,766
squarefree-kernel candidates, and 4,878,128 signed presentations.  Its 22
full prime-filter hits include one off-rectangle hit but no local-contact
survivor.  Exact verification is empty, and repeated-fiber mining is again
unchanged at 55 pairs.

```text
results/target_22224_a2228_deep_p5_box_B10000_a8001_9000_allsigns.tsv
results/target_22224_a2228_deep_p5_box_B10000_a8001_9000_verify.log
results/target_22224_a2228_B10000_a9000_repeated_fibers.txt
```

The final `9001<=a<=10000` tranche tests 165,668,499 ordered triples,
31,988 squarefree-kernel candidates, and 255,904 signed presentations.
Its four full prime-filter hits contain no off-rectangle or local-contact
hit.  The completed `B=10000` box contains 1,841 distinct primitive signed
tuples and exactly 55 repeated square-ratio pairs (52 distinct families
after removing aliases); no target occurs in the transverse stream.

```text
results/target_22224_a2228_deep_p5_box_B10000_a9001_10000_allsigns.tsv
results/target_22224_a2228_deep_p5_box_B10000_a9001_10000_verify.log
results/target_22224_a2228_B10000_complete_repeated_fibers.txt
```

## Whole-family local audit

A component-complete first-pass audit enumerates every finite full-cover
residue `T mod p` for each of the 52 distinct repeated families.  Whenever
all such reductions are smooth with `3` not dividing `#J(F_p)`, and the
`v_p(T)<0` leading cover units are not all squares, reduction injectivity
rules out rational 3-torsion for every rational specialization in that
family.  Through `p<=43` this gives a rigorous one-prime obstruction for 37
families.  Only

```text
F1 F4 F6 N3 N4 N10 N11 O4 O9 R3 R7 S1 S3 V1 W2
```

remain undecided by this inexpensive test.

```text
code/target_22224_repeated_family_local_audit.m
results/target_22224_repeated_family_local_audit.tsv
results/target_22224_repeated_family_local_audit.log
```

Extending the same complete audit through `p=199` kills four of those
fifteen: `F4` at 47, `N10` and `V1` at 73, and `S1` at 71.  The final
eleven families surviving every tested prime through 199 are

```text
F1 F6 N3 N4 N11 O4 O9 R3 R7 S3 W2
```

This is a local-viability shortlist, not evidence that a rational target
exists.  The `S1` and `W1` obstructions were also checked independently with
the corrected cubic-contact equations; their empty coefficient masks hold
at every Mordell--Weil height.

```text
results/target_22224_repeated_family_local_audit_live15_p47_199.tsv
results/target_22224_repeated_family_local_audit_live15_p47_199.log
notes/target_22224_S1_component_complete_2026_07_19.md
notes/target_22224_W1_component_complete_2026_07_19.md
```

## Reusable code

```text
code/target_22224_repeated_fiber_rank2_mwsieve.m
code/target_22224_p8_rank2_mwsieve.m
code/target_22224_rank2_deep_padic_filter.m
```

The first script performs direct bounded rank-1/rank-2 coefficient sieves.
The second performs generalized-CRT periodic enumeration.  The third applies
the corrected `p^5` tube bank directly over `Q_p`, avoiding enormous exact
rational coordinates for large MW indices.

## Current interpretation

Repeated full-cover fibers are abundant and their elliptic quotients nearly
always have rank at least 2.  The U2 analysis also shows that the current
deep tangent bank is sampled rather than component-complete: a zero bank
intersection is useful triage, but is not by itself a local obstruction.
The productive workflow is therefore:

1. identify repeated square-ratio fibers from the expanding transverse box;
2. prove quotient ranks and choose a saturated rank-1/rank-2 quotient when
   available;
3. impose corrected projective contact masks at 11 and 13 immediately;
4. build periodic coefficient lattices rather than raw boxes;
5. use the sampled `p^5` bank to prioritize high-index survivors, then derive
   a component-specific blow-up/lift for every projected contact lane not
   represented in that bank;
6. claim a local exclusion only from a complete component lift or a direct
   cover obstruction, and use exact rational reconstruction for the final
   small set of genuinely viable MW classes.
