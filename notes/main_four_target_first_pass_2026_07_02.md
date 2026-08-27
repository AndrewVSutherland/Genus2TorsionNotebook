# Four target first pass: `Z/5 x Z/5`, `Z/2 x Z/24`, `Z/48`, `Z/35`

Date: 2026-07-02.

This records the first coordinated pass on the four routes selected from
`ideas.txt`.

## `Z/5 x Z/5`

Detailed files:

```text
notes/agent_z5x5_contact5_contact5.md
code/agent_z5x5_contact5_contact5_symbolic.m
code/agent_z5x5_contact5_contact5_search.m
```

The literal two-rational-contact construction is empty over `Q`, not merely
cold in a search box.  After normalizing the two contact points to `0` and
`1`,

```text
f = h0^2 - K*x^5 = h1^2 - K*(x-1)^5
```

forces

```text
(h0-h1)*(h0+h1) = K*(x^5-(x-1)^5).
```

But

```text
x^5-(x-1)^5 = 5*x^4 - 10*x^3 + 10*x^2 - 5*x + 1
```

is irreducible over `Q`, while both `h0-h1` and `h0+h1` have degree at most
`2`.  The same-contact branch is exactly the boundary `r=0`; `K=0` is square
and not a smooth genus-2 curve.

Finite-field smoke tests confirm the model: over primes where the quartic
splits into quadratic factors, independent `F_p` 5-contact pairs appear.

Next move: keep one rational 5-contact, but make the second 5-torsion class a
degree-2 Mumford/contact condition, or use a quadratic conjugate 5-contact pair
and descend the Galois-invariant class.

## `Z/2 x Z/24`

Detailed files:

```text
notes/agent_A2_24_halving_cover.md
code/agent_A2_24_halving_cover_probe.m
```

The translated order-12 halving problem is the standard square-quartic system.
For `D_T=[u_T,v_T]`,

```text
ell_T = v_T + u_T*(M*x+N),
S_T = (ell_T^2-f)/u_T
    = s4*x^4+s3*x^3+s2*x^2+s1*x+s0.
```

On `s4 != 0`, halving is cut out by

```text
E1 = 8*s4^2*s1 - s3*(4*s4*s2 - s3^2),
E0 = 64*s4^3*s0 - (4*s4*s2 - s3^2)^2.
```

For the known height-10 local survivor, the raw degree-32 resultants are
misleading: the rational linear factors are exactly the boundary

```text
s4 = M^2 - 81/15625 = 0.
```

At those rational boundary points `deg(S)=2`, so they are not halves.  After
removing `s4=0`, each of the four translated order-12 classes has affine
degree `16`, split as two irreducible degree-8 components, with no rational
half.  The curve has exact torsion `[2,2,12]`.

Next move: do not build Chabauty/descent from the unsaturated degree-32
resultants.  If this lane continues, form the saturated global cover over a
chosen `A(2,12)` resolvent branch with `s4=0` removed first.

## `Z/48`

Detailed files:

```text
notes/agent_Z48_A16_plus3.md
code/agent_Z48_A16_plus3_probe.m
code/main_Z48_known_A16_3filter.m
code/main_Z48_A8_plain_prefilter.m
```

For a curve with a rational point of order `16`, producing `Z/48` is exactly
the same as adding nonzero rational `3`-torsion.  Hence every good prime
`p != 3` must satisfy

```text
3 | #J(F_p).
```

All known small certified A(16) curves fail this immediately.  The documented
simple `[16]` curves are killed by small good primes:

```text
(r,t,p)=(3,1/3,2):        p=5,  #J(F_5)=16
(r,t,p)=(3,1/3,34/9):     p=13, #J(F_13)=128
```

The split sanity-check `[4,16]` and `[2,2,16]` curves also fail the 3-part
test.

I also added a small A(8) scout that first enforces `48 | #J(F_p)` at two good
primes and only then attempts the exact order-16 halving.  Results:

```text
RH=PH=TH=4: filterPass=324,  halvable=0, z48Hits=0
RH=PH=TH=5: filterPass=2151, halvable=0, z48Hits=0
```

Next move: future A16 searches for `Z/48` should include the `3 | #J(F_p)`
gate before any expensive halving/cubic-contact work.  The current known A16
points do not lead to `Z/48`; a serious continuation needs a simultaneous
filtered A16-plus-3 construction rather than mining the existing hits.

## `Z/35`

Detailed files:

```text
notes/agent_Z35_contact7_p3_boundary.md
code/agent_Z35_contact7_boundary_probe.m
```

The existing `simple_35_attempt.md` had already shown that the open contact-7
plus 5 route is obstructed at `p=3`, and that the broad `p=3` boundary search
through height `30` has no survivors.

The remaining plausible boundary branches were the degenerate `h(1)=0`
classes `(a,b)=(0,1)` and `(2,2)`.  On the blow-up chart

```text
a = a0 + 3*u,
b = 5/2 - a + 3*s,
h(1)=3*s,
```

the first layer `v3(s)=0` has node thickness `7`.  Thus the contact-7 class is
not killed by the blow-up; direct Magma samples confirm exact order `7`.

However, the necessary `5 | #J(F_p)` filters away from `p=3` kill all tested
local lifts:

```text
(0,1), first layer, height 18:       59994 smooth, 0 survivors
(2,2), first layer, height 16:       36498 smooth, 0 survivors
(0,1), layers v3(s)=0,1, height 14:  33930 smooth, 0 survivors
```

Next move: direct contact-7 plus 5 now looks weak.  If revisited, it should be
via a reason that forces high 3-adic depth or a different simultaneous
contact-7/contact-5 construction, not another broad height increase on the
same boundary chart.

