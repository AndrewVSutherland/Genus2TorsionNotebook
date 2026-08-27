# Agent Z/48: simultaneous A(16) plus 3-part gate

Date: 2026-07-02.

This continues the `Z/48` lane after the known simple A16 examples failed the
3-part point-count gate.

## Context

I read the requested A8/A16 notes and scanners, including:

- `notes/agent_Z48_A16_plus3.md`
- `code/main_Z48_A8_plain_prefilter.m`
- `a8_order16_cover_status.md`
- `search_A8_plain_halve16.m`
- `scan_A8_order16_target_rt.m`
- the existing plus-3/contact scripts
  `m2228_three_torsion_search.m`, `m3222_plus3_search.m`,
  and `m212_extra3_search.m`.

The main-thread scouts `main_Z48_A8_plain_prefilter.m` and
`main_Z48_A8_wsplit_prefilter.m` already cover raw A8 parameter boxes.  The
new script below instead works on the sign-aware A16 square-root slice
coordinates, so it is not a repeat of those boxes.

## A16-Slice Gate

On the A8 chart, write the visible order-8 class as before.  To halve it, set

```text
ell = ell8base + g8*(mu*x + N)
S = (ell^2 - f)/g8
```

and impose the square-root ansatz

```text
S = C*(x^2 + y*x + z)^2,
C = mu^2 - 2*r*mu + r/t.
```

For fixed `(r,t)`, the coefficient equations are sequential:

```text
x^3 equation: solve N,
x^2 equation: solve z,
x^1,x^0 equations: F1(p)=F0(p)=0.
```

The script enumerates rational `(mu,y)`, takes `gcd(F1,F0)` in `p`, and for
each rational root applies the Z/48 reduction gate before any exact Jacobian
halving:

```text
48 | #J(F_l)
```

at `MinGood` good primes `l != 3`.  The script also supports `GateMode:=3`,
but `GateMode:=48` is the natural default because a certified A16 point would
already force the 16-part at every good prime.

New file:

```text
code/agent_Z48_simultaneous_A16_plus3.m
```

## Cubic-Contact Formulation

The exact rational 3-torsion condition can be appended after this gate.  For a
degree-5 or degree-6 A8/A16 model `C: y^2=f(x)`, a nonzero rational 3-torsion
class represented by a degree-2 Mumford divisor is equivalent to a cubic
contact identity

```text
h(x)^2 - f(x) = Lambda*q(x)^3,
q = x^2 + U*x + V,
h = H3*x^3 + H2*x^2 + H1*x + H0.
```

For odd degree with monic normalization this specializes to the familiar
`Lambda=H3^2` used in `m2228_three_torsion_search.m`.  For a sextic A8 model,
the leading coefficient equation gives `Lambda = H3^2 - lc(f)`.  Thus the
fully simultaneous system is the A16 square-root system above plus the six
coefficient equations from this cubic-contact identity, with the usual
nonboundary conditions `q` squarefree enough and `gcd(q,f)=1`.

In practice the point-count gate is much cheaper, so the script only sends
gate survivors to exact halving/torsion checks.

## Runs

Smoke:

```text
magma -b RTHeight:=2 SearchBound:=6 PrimeBound:=31 \
  MinGood:=2 MaxExact:=50 \
  code/agent_Z48_simultaneous_A16_plus3.m
```

Result:

```text
slices=30, tested=66270, commonRootPairs=8, rationalRoots=8,
singular=8, pointGateReject=0, pointGatePass=0,
exactTried=0, certified=0, z48Hits=0.
```

Medium:

```text
magma -b RTHeight:=3 SearchBound:=10 PrimeBound:=43 \
  MinGood:=3 MaxExact:=100 \
  code/agent_Z48_simultaneous_A16_plus3.m
```

Result:

```text
slices=182, tested=2935478, commonRootPairs=37, rationalRoots=37,
singular=33, pointGateReject=4, pointGatePass=0,
exactTried=0, certified=0, z48Hits=0.

FIRST_POINT_GATE_KILLS:
  p=5 : 4
```

## Outcome

No candidates survived.  In the medium A16 slice box, every nonsingular
rational A16-equation root was killed by the `48 | #J(F_5)` gate before exact
halving, exact torsion, or cubic-contact work.  This improves on checking the
known A16 examples because the gate is now built into the A16 candidate
generation stage itself.

