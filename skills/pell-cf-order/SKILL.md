---
name: pell-cf-order
description: Compute the EXACT order of the class at infinity D_inf = inf+ - inf- on J(y^2=f), deg f = 6, via the polynomial continued fraction of sqrt(f) (Platonov-Petrunin / polynomial Pell) - a cheap exact order filter with NO Jacobian arithmetic, plus the mod-p seed/Newton-lift machinery and its recorded dead end. WHEN hunting [2,2n] targets on 2-rank-2 families, when you need an exact torsion order cheaply, or when tempted to construct curves backward from a prescribed CF period. Trigger words - continued fraction, Pell, D_infinity, quasi-period, CFOrder, partial quotient, order at infinity.
---

# Exact order of D_infinity via polynomial continued fractions

## When to use this

Load this when you want the exact order of the specific torsion class
`D_inf = inf+ - inf-` on `J(y^2 = f)` for a **monic squarefree sextic** `f` —
for example as a cheap exact filter over a 2-rank-2 family when hunting
`[2,2n]` targets. The computation is pure polynomial arithmetic (no Jacobian,
no point counting), so it costs microseconds per curve and is **exact**, not
a necessary condition.

## The mathematical statement

For monic squarefree `f` of degree 6, expand `sqrt(f)` as a polynomial
continued fraction. The class `D_inf` is torsion **iff** the expansion is
(quasi-)periodic, and its exact order is

```text
N  =  SUM of Degree(a_i) over one quasi-period,  INCLUDING  deg(a_0) = 3.
```

(the Platonov–Petrunin polynomial-Pell correspondence: a period gives a
function with divisor `N*(inf+ - inf-)`, the polynomial Pell solution). A
return of `0` from the tool means "no period within the step budget" — that
is **either** non-torsion **or** budget exceeded; `0` is NOT a proof of
non-torsion.

## The algorithm (verbatim, `code/agent_a2_24_cf_search.m`)

```magma
function SqrtPolyPart(f)
    s := x^3;
    for k in [1..3] do
        d := f - s^2;
        if Degree(d) le 2 then break; end if;
        s := s + (Coefficient(d, 6-k)/(2*Coefficient(s,3)))*x^(3-k);
    end for;
    return s;
end function;

// exact D_infty order via CF; returns 0 if not periodic within maxsteps
function CFOrder(f, maxsteps)
    s := SqrtPolyPart(f);
    Pi := P!0; Qi := P!1; total := 0;
    for i in [0..maxsteps] do
        if Qi eq 0 then return 0; end if;
        ai := (Pi + s) div Qi;
        total +:= Degree(ai);
        Pn := ai*Qi - Pi;
        if (f - Pn^2) mod Qi ne 0 then return 0; end if;
        Qn := (f - Pn^2) div Qi;
        Pi := Pn; Qi := Qn;
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then return total; end if;
    end for;
    return 0;
end function;
```

The guards, each load-bearing:

- `s = SqrtPolyPart(f)` is the polynomial part of `sqrt(f)` — degree 3,
  computed by matching the top coefficients. Its degree 3 is what puts
  `deg(a_0) = 3` into the total.
- `(f - Pn^2) mod Qi ne 0 -> return 0`: the division must be exact at every
  step; a nonzero remainder means the recursion has left the CF orbit
  (non-generic input, e.g. non-monic or wrong degree) — abort, don't trust.
- `i ge 1 and Degree(Qi) le 0 and Qi ne 0 -> return total`: a **constant
  nonzero** `Qi` marks the end of the quasi-period (the `i ge 1` excludes the
  start). This is quasi-periodicity — the constant need not be 1.

## Validated test vectors

`code/agent_a2_24_cf.m` is the standalone validated tool with a self-test
(run it: `magma -b code/agent_a2_24_cf.m`). Its recorded vectors:

```text
f14 = (x^2+1)*(x^4+5*x^2+4*x+4)          -> CFOrder = 14   (expected, verified)
f18 = (x^2-x+1)*(x^4-x^3+9*x^2+8*x-8)    -> CFOrder = 18   (expected, verified)
"f28-cyclic?" = x^6+2*x^5-5*x^4-14*x^3-3*x^2+24*x+28   -> CFOrder = 7
                 (verified by running the self-test: the naive "order 28"
                 guess is wrong — the actual D_inf order is the proper
                 divisor 7, degree sequence [3,1,2,1]. THE cautionary vector.)
```

Any modification to the CF code MUST re-pass `f14 -> 14` and `f18 -> 18`
before being used in a search.

## Search deployment

`code/agent_a2_24_cf_search.m` deploys `CFOrder` over the 2-rank-2 family

```text
f = (x^2 - 1)*(x^2 + a*x + b)*(x^2 + c*x + d)      // factor type [1,1,2,2] -> 2-rank 2
```

as the FIRST nontrivial filter (before any Jacobian arithmetic): compute
`ord := CFOrder(f, 40)`; only `ord mod 24 eq 0` candidates proceed to exact
`TorsionSubgroup` + 2-rank + simplicity. It also has a `Validate:=true` mode
(run that after any edit). The honest recorded outcome of the `[2,24]` hunt on
this family: scans found large CF orders but never a verified simple order-24
`D_inf` at accessible height — the deployment is validated machinery on a cold
family; see `notes/agent_a2_24_composite.md` for where the `[2,24]` effort
went instead (the composite `8x3` route).

## Mod-p seeds and the backward construction (recorded DEAD END)

The tempting inverse problem — *construct* `(a,b,c,d)` with CF order exactly
24 by working backward from the period structure — was tried and is recorded
here so nobody re-walks it:

- `code/agent_a2_24_construct_seeds.m`: scans `(a,b,c,d)` in `F_p^4`,
  computes CF orders over `F_p` (guard: `IsUnit(LeadingCoefficient(Qi))`,
  abort otherwise — over `F_p` a non-unit leading coefficient breaks the
  division), and records order-24 seeds with their **degree pattern**; the
  generic order-24 cell is `[3,1,1,...,1]` (twenty-one degree-1 partial
  quotients after `a_0`).
- `code/agent_a2_24_construct_lift.m`: for each generic seed, freezes `(a,b)`,
  Newton-lifts `(c,d)` p-adically along the two closure conditions (the
  `ClosureG` function computes the closure coefficients over `Z/p^k` with
  MANUAL division by unit-leading-coefficient divisors — generic `div` does
  not exist there), then attempts `RationalReconstruction` (arity 2 — see
  `magma-lab-conventions`) and exact verification over `Q`.
- **Outcome:** lifts converge p-adically but rational reconstructions do not
  verify over `Q` at any reasonable height — the order-24 CF locus has
  algebraic degree too high for this route. The machinery (mod-p^k CF,
  finite-difference Newton) is preserved and reusable; the route as a way to
  reach `[2,24]` is dead. Do not restart it without a new structural idea.

## Pitfalls

- **Forgetting `deg(a_0) = 3` in the total.** The order includes the degree of
  the integer part. Omitting it reports `N - 3` and every downstream
  conclusion is wrong. The validated vectors catch this instantly — run them.
- **Reading `return 0` as "not torsion".** `0` = no period within
  `maxsteps`. Raise the budget (the searches use 40–60) before concluding
  anything; and even then, `0` only means "no period found", full stop.
- **Non-monic or degree-5 input.** `SqrtPolyPart` assumes a monic sextic
  (`s := x^3` and division by `2*Coefficient(s,3)`). For a degree-5 model
  there is only ONE point at infinity and `D_inf` is not even the right
  object. Normalize to a monic sextic model first, or don't use this tool.
- **Skipping the exactness check** `(f - Pn^2) mod Qi ne 0`. Without it the
  recursion silently produces garbage on non-generic input instead of
  aborting.
- **Quasi-period vs period.** The stop condition is `Qi` constant **nonzero**
  (quasi-period). Testing `Qi eq 1` (strict period) misses valid torsion and
  loops to the budget.
- **Over `Z/p^k`, dividing by a non-unit.** Both the mod-p seed scan and the
  lift must guard `IsUnit(LeadingCoefficient(...))`; a non-unit leading
  coefficient means that seed/step is outside the generic cell — skip it, do
  not force it.
- **`D_inf`'s order is not the exponent of the torsion group.** It is the
  order of ONE class. A curve can have large torsion with small `D_inf` order
  and vice versa; combining `D_inf` order with independence of other classes
  still requires exact `TorsionSubgroup` (`validate-and-record-a-hit`).

## See also

- `two-rank-and-factor-types` — why the `(x^2-1)(...)(...)` family has 2-rank
  2 built in.
- `running-torsion-searches` — where a CF filter slots into the funnel (very
  early: it is exact AND cheap).
- `magma-lab-conventions` — `RationalReconstruction` arity; manual division
  idioms; batch-mode discipline for the seed scans.
- `component-boundary-analysis` — the covering-condition lens that explains
  *why* backward constructions like the CF lift hit high algebraic degree.
- `local-obstructions` — deciding whether a cold family is obstructed or thin.
- `g2-torsion-lab` — hub.
