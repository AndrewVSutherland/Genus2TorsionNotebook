# Four target second pass

Date: 2026-07-02.

This continues `main_four_target_first_pass_2026_07_02.md`.

## `Z/5 x Z/5`

Detailed files:

```text
notes/agent_z5x5_degree2_contact.md
code/agent_z5x5_degree2_contact_symbolic.m
code/agent_z5x5_degree2_contact_search.m
```

After the rational contact-5/contact-5 obstruction, the next smaller test was:
one rational 5-contact at `x=0`, plus a degree-2 divisor forced by a normalized
`B=1` contact identity

```text
f = H^2 - U^5 = h^2 - K*x^5,      U=x^2+s*x+t.
```

This is also too small over `Q`.  The exact elimination leaves only:

```text
t = s^2/4        double-root boundary for U,
t = 5*s^2/4      q4=0 boundary; tail is not a rational square.
```

This is a new obstruction, not the previous irreducible quartic
`x^5-(x-1)^5`: the finite-field smoke tests do produce genuine independent
pairs at `p=11,31,41`.

Next object: the full degree-2 Mumford norm system on the one-contact family,

```text
f = h^2 - K*x^5,
A(x)^2 - B(x)^2*f = Lambda*U(x)^5,
U=x^2+s*x+t,     deg(A)<=5, deg(B)<=2.
```

Equivalently, compute `[5][U,V]=0` by Cantor on the contact-5 family.  The
normalized `B=1` contact slice cannot contain the answer.

## `Z/2 x Z/24`

Detailed files:

```text
notes/agent_A2_24_saturated_global_cover.md
code/agent_A2_24_saturated_global_cover.m
```

The first-pass local survivor showed why `s4=0` must be removed before reading
the degree-32 resultants.  This pass tested an actual split A(2,12) branch:

```text
p=-5/3, z=1, r=2/3,
D = P12 + T_extra.
```

For the saturated square-quartic halving equations:

```text
s4 = M^2 + 225/1331,
resultant_N degree = 32,
boundary contribution = <2,8>,
saturated affine degree = 16,
affine factor degrees = [<16,1>].
```

So this branch has the expected degree-16 halving fiber, but it is one
irreducible degree-16 orbit on the rational split point.  Finite-field counts
on the oriented split slice are sparse after removing `s4=0`.

Decision: this particular branch is not an attractive Chabauty/descent target
unless a later quotient or special subbranch lowers the degree.

## `Z/48`

Detailed files:

```text
notes/agent_Z48_simultaneous_A16_plus3.md
code/agent_Z48_simultaneous_A16_plus3.m
code/main_Z48_A8_plain_prefilter.m
code/main_Z48_A8_wsplit_prefilter.m
```

The simultaneous A16 script enumerates rational square-root A16 candidates on
fixed `(r,t)` slices, applies the good-prime gate

```text
48 | #J(F_l)
```

before exact halving or torsion, and only then would run exact certification.
The medium run

```text
RTHeight=3, SearchBound=10, PrimeBound=43, MinGood=3
```

found 37 rational A16-equation roots.  33 were singular, and the four
nonsingular roots were all killed by the point-count gate at `p=5`.  Thus no
candidate reached exact certification.

I also added two broader A8 scouts:

```text
main_Z48_A8_plain_prefilter.m
main_Z48_A8_wsplit_prefilter.m
```

These first require `48 | #J(F_p)` at good primes, then attempt exact halving
of the visible order-8 class.  Results so far:

```text
plain A8, RH=PH=TH=6, MinGood=2:
  filterPass=3860, halveTried=3860, halvable=0, z48Hits=0

W-split A8, RH=TH=BH=6, MinGood=2:
  filterPass=410, halveTried=410, halvable=0, z48Hits=0
```

A stricter plain A8 run also completed:

```text
plain A8, RH=PH=TH=7, MinGood=3:
  tested=347900, smooth=340345,
  filterPass=4570, halveTried=4570, halvable=0, z48Hits=0
```

Next object: either a larger partitioned A16 square-root scan with the
`48 | #J(F_l)` gate, or the actual cubic-contact equations

```text
h^2 - f = Lambda*q^3
```

added after the point-count gate.  The current small A16/A8 boxes are cold.

## `Z/35`

Detailed files:

```text
notes/agent_Z35_next_route.md
code/agent_Z35_next_route_probe.m
```

The direct contact-7 plus 5 and shallow `p=3` boundary branches were already
cold.  This pass returned to the simultaneous point-contact equations:

```text
q^2 - f = const*(x-r)^5,      q=c0+c1*x+c2*x^2.
```

Writing

```text
d=c2-b,       e=c2+b
```

and eliminating `c1` from the two residual equations gives

```text
Res_c1(N0,N1) = d^3 e^3 (d-e)^8 (d+e)^4 Phi38(d,e).
```

Modulo `3`, the nondegenerate residual chart has no points.  All residual
solutions are forced onto five centers:

```text
d=0, e=0, d=-e, and two d=e,r=1 centers.
```

The best next targets are the two `d=e, r=1` centers, i.e. the `b=0` pole
charts lost by the old division formulas.  Then check the `d=-e` / `c2=0`
pole chart.  The `d=0` and `e=0` centers land back on the already-cold
`(a,b)=(2,2), h(1)=0` branch.

Next object: blow up the original coefficient equations near `b=0, r=1`
without dividing by `b`, divide by the correct powers of `3`, and run a
finite-field viability table on the transformed chart.
