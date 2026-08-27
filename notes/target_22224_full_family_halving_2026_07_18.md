# Full-family `[2,2,2,24]` halving attack

Date: 2026-07-18

> **Correction (2026-07-18, later audit).**  The finite contact formula used
> in the first version of this report has a factor-of-two error: if
> `M=L^2`, its auxiliary numerator must be
> `P=4*M*e1+12*(U^2+v^2)-(M+3*U)^2`.  The old code effectively replaced
> `P` by `4*M*e1+12*(U^2+v^2)-2*(M+3*U)^2`.  Consequently the finite masks,
> incidence ranks, and the claimed `13`-adic Type I/II lift conditions below
> are invalid and must not be used.  Independent direct point counts give
> **no smooth full-cover curve with 3-primary Jacobian order at either
> `p=11` or `p=13`**.  Correct target-presentation counts are `48,96,144,192,
> 336` at `p=17,19,23,29,31`.  The group-theoretic reduction in Section 1,
> the exact four-radicand cover, and the exact bank torsion scan remain valid.
> Corrected audit code is `code/target_22224_full_contact_audit.py` and
> `code/target_22224_full_cover_pointcount_audit.m`.

## Outcome

No rational realization was found, but this run leaves the old two-dimensional
`q=(x+t)^2` slice and replaces the degree-16 affine halving equations by the
correct **full three-dimensional** `A(2,2,2,8)+3` incidence cover.

The main concrete results are:

1. Halving the actually marked order-12 class is equivalent to halving its
   marked order-4 summand.  Thus the target is exactly the intersection of the
   full `A(2,2,2,8)` cover with rational 3-torsion.
2. The full cover has many smooth target residues at `p=11,17,19,23,29,31`;
   every recorded residue is a nonsingular Hensel seed of the expected rank.
3. At `p=13` the smooth target mask is empty.  A boundary computation and its
   first-order lift force every rational target into a genuinely deep
   `13`-adic chart.
4. All `3,675` rows of the supplied `tor2228.txt` bank satisfy the full four
   order-8 cover conditions.  A good-prime 3-primary scan leaves no survivor.
   The bank contains `619` distinct primitive tuples.
5. Only two primitive bank points reach the live mod-`169` boundary charts.
   Both have contact coordinates lifting through `13^5`, but exact rational
   torsion is only `[2,2,2,8]`; they are killed by good reduction at `23` and
   `31`, respectively.

The next search should generate new rational points on the full order-8 cover
directly inside the forced mod-`169` charts and combine them with the positive
target masks at `11,17,19,23,29,31`.  Enlarging the old q-square height box or
the point-on-the-curve K3 slice would miss most of this threefold.

## 1. Exact group-theoretic reduction

The direct construction in `target_22224_order12_halving.m` uses

```text
D3  = [q,h mod q],                         ord(D3)=3,
D4  = [g4,-x*L4 mod g4],                   ord(D4)=4,
D12 = D3+D4,                               ord(D12)=12.
```

This is not merely an abstract decomposition: `MarkedD12Data` constructs the
displayed `D4` and `D12`, and the exact record audit asserts

```text
D12 eq D3 + D4.
```

Since `D3` has odd order,

```text
D3 = 2*(2*D3).
```

Consequently,

```text
D4=2K       => D12=2*(K+2D3),
D12=2H      => D4 =2*(H-2D3).
```

Therefore

```text
D12 in 2J(Q)  <=>  D4 in 2J(Q).
```

Any half of `D12` has exact order `24`.  Since the square-branch model has
full rational `J[2]`, a rational solution contains `[2,2,2,24]` (with exact
torsion and geometric simplicity still to be certified at the end).

This also explains why the earlier square-quartic halving system was much
larger than necessary: the 3-primary summand contributes no 2-divisibility
obstruction.

## 2. The full three-dimensional chart

Use signed square-branch coordinates

```text
C: y^2=x(x+a^2)(x+b^2)(x+c^2)(x+d^2).
```

The marked order-4 class `D4` is divisible by 2 precisely on the finite cover
of projective `(a:b:c:d)` defined by the four square conditions

```text
r0^2 = a*b*c*d,
r1^2 = a*(a+b)*(a+c)*(a+d),
r2^2 = b*(b+a)*(b+c)*(b+d),
r3^2 = c*(c+a)*(c+b)*(c+d).
```

The analogous expression with leading factor `d` is dependent: the product
of all four branch-evaluation radicands is `abcd` times a square.

This is the full `A(2,2,2,8)` threefold.  It must not be confused with the
codimension-one K3 divisor

```text
s2(a,b,c,d)^2=4abcd,
```

where an order-8 half happens to be represented by a point of the curve.  The
earlier `M(2,2,2,8)+3` searches treated only that K3 divisor.

For the 3-primary part, put

```text
q=x^2+U*x+v^2,        L=1/m,        M=L^2,
```

and, for the elementary branch-square coefficients `e1,...,e4`, set

```text
Aaux = 2*M*e1 + 6*(U^2+v^2) - (M+3U)^2.
```

The eliminated cubic-contact equations are

```text
(M+3U)Aaux + 8v^3 - 4e2M - 4U^3 - 24Uv^2 = 0,
Aaux^2 + 16(M+3U)v^3 - 16e3M - 48(U^2v^2+v^4) = 0,
Aaux*v^3 - 2e4M - 6Uv^4 = 0.
```

They reconstruct a cubic `h` satisfying

```text
h^2-f=m^2*q^3.
```

There are eleven affine variables

```text
a,b,c,d,r0,r1,r2,r3,L,U,v
```

and seven equations.  Quotienting the common branch scale leaves dimension
three, as expected for a finite level-8-and-level-3 cover of genus-2 moduli.

## 3. Exhaustive finite masks and Hensel seeds

The finite program exhausts projective signed branch tuples, all four
order-8 square conditions, and the nondegenerate cubic-contact equations.
The counts below are labelled/projective counts; symmetry-equivalent rows are
deliberately retained in the TSV masks.

| `p` | full-cover bases | smooth cover bases | smooth curve keys | target bases | target curve keys | incidence ranks |
|---:|---:|---:|---:|---:|---:|---|
| 11 | 68 | 24 | 4 | 24 | 4 | all 7 |
| 13 | 104 | 48 | 6 | 0 | 0 | -- |
| 17 | 267 | 168 | 18 | 48 | 4 | all 7 |
| 19 | 376 | 192 | 28 | 48 | 4 | all 7 |
| 23 | 677 | 432 | 68 | 48 | 8 | all 7 |
| 29 | 1,368 | 912 | 108 | 144 | 24 | all 7 |
| 31 | 1,705 | 1,200 | 164 | 480 | 64 | all 7 |

Rank `7` is the full row rank of the four cover equations plus three contact
equations in eleven affine variables.  Hence every positive residue listed in
the mask files is a smooth Hensel seed with four affine tangent directions,
or three after projectivizing.

Prime `13` is qualitatively different: there are `48` smooth points of the
full order-8 cover, but none has a nondegenerate rational 3-contact.  Thus a
rational target must have bad/boundary reduction in these signed coordinates
at `13`.

## 4. The full `13`-adic boundary

Allowing zero radicands, singular branch reduction, repeated `q`, and
`gcd(q,f)>1`, the raw normalized mod-13 computation gives

```text
projective bases                         2380
full-cover bases including boundary       508
smooth full-cover bases                     48
bases with raw contact                     268
bases with contact-open                     44
contact-open presentations                  88
```

Every contact-open base is on a zero/collision intersection.  There is no
pure open, pure zero, or pure collision survivor.  Up to permutation, the
generic live signature is

```text
a=0 mod 13,       b+c=0 mod 13.
```

The first-order lift of the 88 contact-open presentations gives

```text
inconsistent mod 13^2                     32
liftable mod 13^2                          56
```

The 48 generic `Zi+Ejk-` presentations are all liftable, but every lift
forces both boundary equations to persist one level deeper:

```text
a=0 mod 169,      b+c=0 mod 169.
```

There is no first-order escape to a smooth base.  The remaining eight
liftable presentations lie on triple-zero/one-unit strata.  If all three
zero coordinates had valuation exactly one, `v_13(abcd)=3`, contradicting
the requirement that `abcd` be a rational square.  At least one of those
three coordinates must therefore be divisible by `169`.

An independent lift calculation which includes explicit square roots of all
four order-8 radicands finds `162` sheet-labelled contact-open incidence
points, `98` lifts through `13^2`, and `64` inconsistent points.  This is the
same boundary geometry with the cover sheets retained.

Thus the two practical target charts are, after permuting coordinates,

```text
Type I:   169 | a,       169 | (b+c);
Type II:  13 | a,b,c,    13 does not divide d,
          and at least one of a,b,c is divisible by 169.
```

## 5. Rational bounded searches

### Entire supplied full-cover bank

The footer of `paper/scripts_and_data/tor2228.txt` labels the generation bound
as `B=16384`.  The file contains `3,675` rows, or `619` distinct primitive
tuples after removing common scale.  Every row was checked against the full
four-radicand order-8 cover.

A good-prime 3-primary gcd scan eliminated the rows as follows:

```text
p=11: 974       p=13: 1654      p=17: 380
p=19: 266       p=23:   44      p=29: 205
p=31:  43       p=37:   76      p=53:  25
p=61:   7       p=67:    1
```

There were no reduction survivors and hence no rational target hit.

### Deep mod-169 survivors

Only four bank rows, representing two primitive tuples, lie on a
contact-open mod-13 component and satisfy the forced mod-169 conditions:

```text
[121,1919,3211,4949],
[1369,1711,2349,9971].
```

After projective normalization their residues mod `169` are

```text
(1,41,0,-1),       (1,41,-1,0).
```

For each fixed rational base, the contact coordinates have two Hensel lifts
at every level tested through `13^5`; the two lifts differ by the expected
sign of `L`.  They are genuine local positive controls for the deep chart.

Exact rational computation nevertheless gives

```text
tuple                         exact torsion     first 3-killing prime
[121,1919,3211,4949]          [2,2,2,8]         23
[1369,1711,2349,9971]         [2,2,2,8]         31
```

Both Jacobians are strictly geometrically simple by root-power certificates:
the certificate primes are `53` and `71`, respectively.  Thus failure is
specifically the absence of rational 3-primary torsion, not split geometry.

## 6. Recommended continuation

The highest-value next computation is a resumable generator for rational
points of the **full four-radicand cover** subject from the outset to either
Type I or Type II mod-169 conditions.  Each generated point should then be
tested against the exact target masks at `11,17,19,23,29,31` before any
Jacobian computation.

The two deep bank points show that these 13-adic strata really contain simple
rational `A(2,2,2,8)` points and p-adic 3-contact lifts.  What remains is to
move along the three-dimensional cover until the local 3-classes at the other
good primes become globally compatible.  This is a substantially sharper
search target than either:

* the old two-dimensional q-square slice; or
* the point-on-the-curve `M(2,2,2,8)` K3 divisor.

## 7. Artifacts

```text
code/target_22224_full_family_halving.m
code/target_22224_full_family_halving_padic13.py
code/target_22224_full_boundary13.py
code/target_22224_full_boundary13_lift.py

results/target_22224_full_family_halving_bank.log
results/target_22224_full_family_halving_deep_seeds.log
results/target_22224_full_family_halving_finite.log
results/target_22224_full_family_halving_finite_p11.tsv
results/target_22224_full_family_halving_finite_p13.tsv
results/target_22224_full_family_halving_finite_p17.tsv
results/target_22224_full_family_halving_finite_p19.tsv
results/target_22224_full_family_halving_finite_p23.tsv
results/target_22224_full_family_halving_finite_p29.tsv
results/target_22224_full_family_halving_finite_p31.tsv
results/target_22224_full_family_halving_padic13.txt
```

Reproduction from `torsion_jac`:

```bash
magma -b mode:=finite PrimeList:=11,13,17,19,23,29,31 \
  log_file:=results/target_22224_full_family_halving_finite.log \
  code/target_22224_full_family_halving.m

magma -b mode:=bank \
  log_file:=results/target_22224_full_family_halving_bank.log \
  code/target_22224_full_family_halving.m

magma -b mode:=deep \
  log_file:=results/target_22224_full_family_halving_deep_seeds.log \
  code/target_22224_full_family_halving.m

python3 code/target_22224_full_family_halving_padic13.py \
  --output results/target_22224_full_family_halving_padic13.txt \
  --bank paper/scripts_and_data/tor2228.txt
```
