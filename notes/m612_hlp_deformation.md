# Unrestricted marked deformation of the split HLP `[6,12]` seed

Date: 2026-07-10.

## Outcome

There is a genuine infinitesimal deformation of the explicit HLP `[6,12]`
fiber that preserves both independent rational order-`3` classes and the
marked rational order-`4` class, while leaving every `(2,2)`-split/Humbert-4
branch visible at the symmetric seed.

The unrestricted incidence system has

```text
28 variables, 21 coefficient equations,
Jacobian rank 21 at the HLP seed,
tangent dimension 7.
```

Each of its three `7 x 7` auxiliary blocks is invertible.  Hence projection
to the seven sextic coefficients has full tangent rank `7`.  After quotienting
by `PGL_2` and the equation scaling, this is the expected three-dimensional
fine-moduli tangent space.

An explicit rational transverse direction is

```text
delta F = 1+x.
```

It has a fully rational lift to all marked-contact variables.  Intersecting
the incidence system with

```text
F_t = F_HLP+t*(1+x)
```

gives a one-dimensional algebraic slice over `Q`, smooth at the rational
point `t=0`, whose tangent leaves the split locus.

What is not yet proved is that this slice has a second rational point.  The
implicit branch is a rational **formal** arc; it is not yet a rationally
parametrized algebraic curve.  Thus no nonsplit exact specialization was sent
to Magma in this bounded first pass.

Arithmetic continuation on the exact slice now gives a strong bounded
negative result: projective finite-support masks through `p=89`, followed by
a conservative CRT enumeration of every primitive `t=n/d` with
`|n|,d <= 1000000`, leave only

```text
t=0,
```

the original split HLP fiber.  No nonzero rational specialization survives
to the exact-solving stage.


## The seed

Put

```text
F = 183*(x^2+1)*(32*x^2+61*x+32)*(32*x^2-61*x+32)
  = 187392*x^6-118767*x^4-118767*x^2+187392.
```

The exact rational torsion is `[6,12]`.  The extraction file

```text
code/m612_hlp_marked_torsion_extract.m
```

prints two independent order-`3` Mumford classes and an order-`4` class:

```text
qA = x^2-61/8*x+1,
vA = 197213/8*x-3721;

qB = x^2-13/48,
vB = 18605/48;

uG = x^2-32/29,
vG = 11163/29*x,
2G = Q = [x^2+1,0].
```


## Exact global identities

The printed `vA,vB` are reduced Mumford remainders, not yet the global
contact polynomials.  If `q` is a Mumford support and `v` its remainder, write

```text
H = v+q*(s*x+t).
```

The condition that `H^2-F` be divisible by `q^3` first gives two linear
equations for `(s,t)`: reduce

```text
(F-v^2)/q-2*v*(s*x+t)
```

modulo `q`.  The resulting lifts are as follows.

### First order-`3` direction

Here

```text
s = 3904/9,       t = 29585/9,
```

and

```text
HA = -3904/9 + 61/3*x - 61/3*x^2 + 3904/9*x^3.
```

The exact identity is

```text
HA^2-F = 62464/81*(x^2-61/8*x+1)^3.                  (A)
```

Moreover `HA mod qA = vA`, so this is the extracted Mumford class.

### Second order-`3` direction

Here `(s,t)=(0,-183)`, and

```text
HB = 2623/6-183*x^2.
```

The exact identity is

```text
HB^2-F = -187392*(x^2-13/48)^3.                     (B)
```

Again `HB mod qB = vB`.

### The order-`4` half

Let

```text
q0  = x^2+1,
uG  = x^2-32/29,
ell = 183*x*(x^2+1).
```

Then

```text
ell mod uG = 11163/29*x = vG
```

and

```text
ell^2-F = -153903*(x^2-32/29)^2*(x^2+1),            (G)
153903 = 183*29^2.
```

Because `ell` is divisible by `q0`, identity (G) makes `[q0,0]` a rational
`2`-class and `[uG,ell mod uG]` one of its halves.  Thus it preserves the
`[2,4]` primary structure without imposing a contact-6/rational-Weierstrass
normalization.

All identities are checked exactly in

```text
code/m612_hlp_deformation_tangent.py
```

using only rational arithmetic.


## The unrestricted simultaneous system

Let

```text
F(x) = f0+f1*x+...+f6*x^6.
```

For the two independent `3`-directions introduce

```text
qA = x^2+a1*x+a0,       HA = ha3*x^3+ha2*x^2+ha1*x+ha0,
qB = x^2+b1*x+b0,       HB = hb3*x^3+hb2*x^2+hb1*x+hb0,
```

and constants `kA,kB`.  For the `2`-halving introduce

```text
q0  = x^2+c1*x+c0,
u   = x^2+u1*x+u0,
L   = l1*x+l0,
ell = q0*L,
```

and `kG`.

The complete system is the `21` coefficient equations

```text
coeff_i(HA^2-F-kA*qA^3)       = 0,   i=0,...,6,       (SA)
coeff_i(HB^2-F-kB*qB^3)       = 0,   i=0,...,6,       (SB)
coeff_i(ell^2-F-kG*u^2*q0)    = 0,   i=0,...,6.       (SG)
```

There are `28` variables:

```text
7 curve coefficients
+ 7 first-contact variables
+ 7 second-contact variables
+ 7 halving variables.
```

The seed coordinates are

```text
(f0,...,f6) = (187392,0,-118767,0,-118767,0,187392),

(a0,a1,ha0,ha1,ha2,ha3,kA)
 = (1,-61/8,-3904/9,61/3,-61/3,3904/9,62464/81),

(b0,b1,hb0,hb1,hb2,hb3,kB)
 = (-13/48,0,2623/6,0,-183,0,-187392),

(c0,c1,u0,u1,l0,l1,kG)
 = (1,0,-32/29,0,0,183,-153903).
```

The same equations are prepared in Magma in

```text
code/m612_hlp_deformation_system.m
```

but that file was not launched while the shared Magma slots were occupied.


## Jacobian calculation

For fixed `F`, each contact block has seven auxiliary variables and seven
equations.  At the seed their exact ranks are

```text
rank J_A = 7,
rank J_B = 7,
rank J_G = 7.
```

The full `21 x 28` Jacobian therefore has

```text
rank = 21,       kernel dimension = 7.
```

Since all three auxiliary blocks are invertible, every first-order sextic
deformation has a unique lift to the marked incidence space.  In particular,
the marked system is smooth at the HLP point and the projection to sextic
coefficient space is etale there.

This is the key difference from the codimension-one contact-6 construction:
the two order-`3` contacts plus the half of a `2`-class do not force a local
equation on the genus-2 moduli.  They give a fine level cover of the full
three-dimensional moduli space.


## Leaving the split locus

The seed sextic is even and reciprocal.  Its three visible reduced elliptic
involutions are

```text
x -> -x,       x -> 1/x,       x -> -1/x.
```

For each marked involution, the tangent to its Humbert-4 branch is the sum of
the invariant binary-sextic subspace and the infinitesimal `PGL_2` orbit of
`F`.  Exact linear algebra gives the following primitive normals, written on
`(delta f0,...,delta f6)`:

```text
N_{-x}    = (0,1298,0,2423,0,1298,0),

N_{1/x}   = (-649,0,-3072,0,3072,0,649),

N_{-1/x}  = (-649,0,-3072,0,3072,0,649).
```

The last two branches are tangent to one another at this extra-symmetric
fiber.  Now take

```text
delta F = 1+x.
```

Then

```text
N_{-x}(delta F)   = 1298,
N_{1/x}(delta F)  = -649,
N_{-1/x}(delta F) = -649.
```

Hence this direction is outside all three split-branch tangent hyperplanes.
It is not a coordinate change of the HLP curve and is a visible first-order
direction away from the `(2,2)`-split locus.

Escaping Humbert `H_4` does not by itself certify absolute simplicity of a
future specialization: another Humbert surface could still be met.  An
irreducible good-prime Frobenius polynomial remains the intended final
simplicity certificate.


## A smooth algebraic slice over `Q`

Impose the six linear equations

```text
f0-f1 = 187392,
f2 = -118767,
f3 = 0,
f4 = -118767,
f5 = 0,
f6 = 187392.
```

Equivalently,

```text
F_t = F+t*(1+x),       t=f1.
```

Together with `(SA),(SB),(SG)`, this gives `27` equations in `28` variables.
At the HLP point its exact Jacobian has

```text
rank 27,       local dimension 1.
```

Thus this is an exact algebraic curve over `Q`, smooth at the rational seed,
with `t` a local parameter and transverse tangent `1+x`.  The formal implicit
function theorem gives a unique branch over `Q[[t]]`; the script prints the
complete rational first-order lift.

For example the three auxiliary tangent blocks begin

```text
A: da0=4279/697285632,
   da1=-1703383/17432140800, ...

B: db0=db1=9/22698100, ...

G: dc0=dc1=1/680943,
   du1=-32/19747347, ...
```

The complete values are intentionally left in the reproducibility script
rather than copied here.


## Exact finite-support masks on the transverse slice

Rather than computing large symbolic resultants, the three finite covers can
be tested exactly over each finite field.  This is an equivalent finite
resolvent test and retains the block structure of the incidence system.

The implementation is

```text
code/m612_hlp_slice_finite_masks.py
```

and works with the homogeneous/projective contact identities, not merely
finite Jacobian orders.

### The two cubic-contact blocks

For a monic affine quadratic

```text
q=x^2+U*x+V,
h=h3*x^3+h2*x^2+h1*x+h0,
```

the `x^6,x^5,x^4,x^3` coefficients of

```text
h^2-F_t=k*q^3
```

determine `k,h2,h1,h0` successively after choosing `(U,V,h3)`.  The three
low coefficients are then checked exactly.

To avoid losing a rational class whose support meets infinity after
reduction, the script also enumerates the remaining projective charts

```text
q=x+V,       q=1.
```

Thus it enumerates all projective binary-quadratic contact supports.  A
rational `[3,3]` subgroup supplies at least two distinct supports at every
good prime; distinct supports represent distinct order-`3` lines.

### The halving block

For every monic quadratic divisor `q0` of `F_t`, divide

```text
F_t=q0*R.
```

The halving identity becomes

```text
q0*L^2-R=k*u^2.                                       (H)
```

After choosing the linear polynomial `L`, the left side is a quartic.  The
script tests directly whether it is a nonzero scalar times the square of a
projective quadratic.  It includes the affine-degree `2`, `1`, and `0`
charts for `u`, so halves specializing toward infinity are retained.

For a smooth degree-six fiber, a quadratic branch factor `q0` cannot itself
contain infinity; hence the monic `q0` chart is complete.  Singular curve
residues are never rejected.

### Masks

The exact block masks are in

```text
data/m612_hlp_slice_masks_p43.txt
data/m612_hlp_slice_masks_p47_89.txt
```

For each prime the stored allowed affine set is

```text
bad/singular t residues
union
smooth t residues with >=2 contact supports and >=1 halving support.
```

The projective parameter `t=infinity` is added conservatively by the height
sieve.  Prime `61`, where the displayed sextic loses degree because
`61 | 187392`, is omitted entirely.

Some especially sharp masks are

```text
p=11: allowed affine {0}, no bad residue;
p=13: allowed affine {0}, no bad residue;
p=31: allowed affine {0}, no bad residue;
p=43: allowed affine {0}, no bad residue;
p=71: allowed affine {0}, no bad residue.
```

At `p=5`, `t=0` is singular and every smooth residue fails the simultaneous
target.  These block masks agree with the independent aggregate finite-group
test in `code/m612_hlp_transverse_finite.m` on the decisive residues; the
aggregate test is occasionally stronger because it checks the full invariant
factor pattern.


## Projective height sieve through one million

The height enumerator is

```text
code/m612_hlp_slice_height_sieve.py
```

It combines the singleton masks at `11,13,31,43` by CRT before visiting a
numerator.  For every prime dividing the denominator `d`, it retains the
projective infinity class and permits every numerator residue compatible
with `gcd(n,d)=1`.  Thus no candidate is discarded merely because a marked
support or the parameter has a pole.

The completed run was

```text
python3 code/m612_hlp_slice_height_sieve.py \
  --height 1000000 \
  --progress 100000 \
  --masks data/m612_hlp_slice_masks_p43.txt \
          data/m612_hlp_slice_masks_p47_89.txt \
  --out data/m612_hlp_slice_candidates_h1000000.txt
```

Its exact summary is

```text
SURVIVOR n=0 d=1 t=0/1

DONE height 1000000
     leaves 12645732
     bounded 132995974
     primitive 82753525
     survivors 1

kill_counts:
p=71  80822104
p=83   1877992
p=89     51622
p=59      1745
p=53        59
p=79         2
```

All later masks had nothing left to kill.  The candidate file contains only
the original split point `t=0`.  Consequently there was no nonzero survivor
on which to solve the three rational blocks or run exact Magma torsion and
simplicity tests.

The compact run transcript is

```text
data/m612_hlp_slice_height_h1000000_summary.txt
```


## Current arithmetic status and next move

The positive conclusion is local and geometric: an unrestricted simple
deformation direction is present, and the marked incidence is nonsingular.
The remaining obstruction is arithmetic rationality on the transverse slice.

The slice is a fiber product over the `t`-line of three finite covers:

1. the selected `qA,HA,kA` contact cover;
2. the selected `qB,HB,kB` contact cover;
3. the selected `q0,u,L,kG` halving cover.

A rational value of `t` must split the selected local branch in all three
covers simultaneously.  The height-one-million computation shows that this
does not happen at accessible naive height on the particular linear slice
`F+t*(1+x)`.

The next deformation computation should therefore change the algebraic slice
rather than enlarge this same height box immediately.  Natural options are:

1. choose another transverse curve direction, preferably with a low-degree
   rational parametrization in one contact block;
2. allow a two-parameter transverse plane and search rational sections of
   one block before intersecting the other two;
3. compute symbolic block resolvents only after such a lower-degree slice is
   identified;
4. for a hit, verify exact torsion `[6,12]`, independence of the two
   order-`3` classes, and an irreducible local Frobenius polynomial.


## Small transverse-direction optimization

A deterministic follow-up searched all primitive sign-normalized directions

```text
G=g0+g1*x+...+g6*x^6,       gi in {-2,-1,0,1,2},
```

subject to both independent Humbert-4 normal pairings being nonzero.  This
leaves `35,136` directions.  A fast projective-block prefilter at a smooth
finite fiber imposed

```text
rational 2-rank >= 2,
at least two projective cubic-contact supports,
at least one rational half of a nonzero rational 2-class.
```

These conditions are necessary for the finite Jacobian to contain `[6,12]`,
but not sufficient on every auxiliary boundary: degenerate support data can
give false positives.  Consequently every reported final mask below was
recomputed independently from the full finite Jacobian invariant factors in
Magma.  All projective support charts are included in the prefilter.  A
smooth degree-five fiber is first replaced by the birational degree-six model

```text
z^6*f(a+1/z),       f(a) != 0,
```

so a `2`-class involving the branch point at infinity is not lost.

The block-prefilter staging was exhaustive for the requested core primes:

```text
all 35,136 directions at p=7,11,13;
best 2,000 core-prime directions at p=17,19,23.
```

Prime `5` was scored separately.  In fact no smooth finite fiber in the
box-two search even has `72 | #J(F_5)`; `t=0` is singular modulo `5`.
Finite bad residues and the projective parameter `t=infinity` were kept
separate throughout.

The scripts and transcripts are

```text
code/m612_hlp_direction_order_prefilter.py
code/m612_hlp_direction_exact_block_search.py
code/m612_hlp_direction_exact_staged.py
code/m612_hlp_direction_search.m
data/m612_hlp_direction_exact_staged_b2.txt
data/m612_hlp_direction_exact_exhaustive_core_b2.txt
data/m612_hlp_direction_aggregate_finalists_magma.txt
```

### Best broad direction

Up to the symmetries `x -> -x` and coefficient reversal induced by
`x -> 1/x`, the best direction in the exhaustive box-two block search, and
also the best after aggregate invariant-factor verification of the leading
finalists, is

```text
G_A = 2+x-x^2+x^3+x^4+x^5+x^6.
```

Its two independent normal pairings are

```text
(N_minus,N_plus)=(5019,5495),
```

so it is genuinely transverse.  The exact masks are:

```text
p    smooth allowed             bad/singular
5    {}                         {0,1}
7    {0}                        {}
11   {0,6,8}                    {2,4}
13   {0,8}                      {11}
17   {0}                        {1}
19   {0}                        {}
23   {0,8}                      {}
```

There are no unresolved residues.  Magma gives the following exact invariant
factors at the four nonzero smooth target fibers:

```text
p=11,t=6 : [12,12]
p=11,t=8 : [12,12]
p=13,t=8 : [12,12]
p=23,t=8 : [24,24].
```

If bad residues are passed conservatively for a rational height sieve, the
product of the seven finite pass densities is exactly `10` times the product
for `G=1+x`.  Thus `G_A` is the preferred slice for auxiliary-block
resolvents: it improves genuine smooth splitting at three primes, not merely
by accumulating boundary fibers.

### Clean second direction

The next aggregate-verified direction with genuine smooth target fibers at
both middle core primes is

```text
G_C = 2+x-2*x^2+x^3+x^4+2*x^5,
(N_minus,N_plus)=(6317,7918).
```

Its exact masks are

```text
p    smooth allowed             bad/singular
5    {}                         {0}
7    {0}                        {}
11   {0,7}                      {}
13   {0,6}                      {8,10}
17   {0}                        {}
19   {0}                        {}
23   {0}                        {}
```

The two nonzero fibers have invariant factors `[12,12]` at `p=11,t=7` and
`[6,36]` at `p=13,t=6`.  This is a cleaner but narrower second resolvent
experiment; `G_A` has the broader smooth-target support.

One important diagnostic is that the block prefilter alone falsely retained
nonzero fibers for

```text
2+x+x^2-2*x^3-x^4-2*x^5
```

at `p=7,17`; the aggregate Magma test kills both.  This is why the final
ranking is based on invariant factors rather than the auxiliary block counts.

The exhaustive core-prime block search and aggregate checks of the leading
finalists did not find a line with a nonzero smooth target residue at every
tested prime.  The gain is therefore substantial but not a
local-unobstructedness theorem for one linear slice.

No exact nonsplit specialization was constructed.  The combined result is:
the simple formal deformation is visible; the original `1+x` slice has no
nonzero rational candidate through height one million; and `G_A` is the
strictly better next linear slice on which to compute the three auxiliary
block resolvents.

### Follow-up on the `G_A` auxiliary covers

The resolvent computation has now been carried out far enough to decide the
low-genus question.  The two marked order-3 contacts are two points above
`t=0` on the same degree-40 support cover.  A good reduction modulo `101`
has an irreducible degree-40 resolvent over `F_101(t)`, certifying that this
cover is irreducible in characteristic zero.  The pencil has ten simple
nodal fibers and smooth fiber at infinity; Picard--Lefschetz plus
Riemann--Hurwitz gives genus `51`.  Thus neither marked order-3 branch lies
on a rational or elliptic component of this linear slice.

The primitive order-4 support cover has degree `120`; if connected its genus
is `181`.  Connectivity of the marked halving component has deliberately
not been asserted without its own modular factorization.  Full details and
reproducibility files are in
`notes/m612_hlp_GA_probe_2026_07_10.md`.

Accordingly `G_A` remains useful for finite-support sieving, but the next
geometric deformation should use a two-parameter surface on which one
contact block may acquire a section, not attempt to parametrize a marked
branch on this line.
