# Geometry triage for the reduced `[8,8]` and corrected `[64]` covers

This is a bounded symbolic pass only; it performs no further height search. The reproducible driver is `code/task5_cover_geometry.sage`, with output in `data/task5_cover_geometry.txt`.

## Corrected order-64 cover

Sage confirms that the projective closure of `C32(z,r)=0` is irreducible of degree 8, arithmetic genus 21, and geometric genus 0. Factoring Sage's bulky normalization and changing the source coordinate gives

```text
r = 1/(t^4-1),
z = 2*(t^4+t^3+t^2+t+1)/(t^2*(t^2+t+1)).
```

The printed point `(z,r)=(11/6,1/15)` is `t=-2`. Using the corrected scale `p=4*Pnum/Pden` gives

```text
p = t^5*(t^2+1)*(t^2+t+1)^2/(t^4+t^3+t^2+t+1),
a = 2*t^3*(t^2+1)*(t^2+t+1).
```

The explicit `b,c` are in the data file. Eliminating `n` from the generic point-halving equations gives the expected boundary factor `f5^4` and one degree-16 factor in `m`. (The other sign is isomorphic via `(m,n)->(-m,-n)`.) After specialization to this `t`-model, that degree-16 factor is irreducible in `QQ(t)[m]`. Thus the good-open corrected order-64 cover is a connected degree-16 cover of `P1_t`.

Clearing denominators gives bidegree `(176,16)` in `(t,m)`, 1712 terms, and Newton polygon with vertices `(0,0),(176,0),(64,16),(0,16)`. Its 1785 interior lattice points give a Baker/Newton upper bound `g <= 1785` (equality would require nondegeneracy). Sage's exact `FunctionField.genus()` integral-basis computation timed out at 180 seconds; a direct discriminant computation timed out at 90 seconds, so no exact genus is claimed.

The next exact arithmetic object is the irreducible degree-16 extension `QQ(t)[m]/(H64)`, constructed as `h` in `analyze_64_extension()`. The next useful calculation is a local integral-basis/different computation at the finite boundary factors visible in `a,b,c`, not another parameter-height scan.

## Reduced `[8,8]` cover

The global target is a surface, not a curve: seven reduced coordinates `(R,s,U,M,N,z,tau)` satisfy three first-half and two second-half equations and map finitely to `(R,s)`. Hence it has dimension 2 and no single genus.

On the good square subcover, the first sign-eliminant has two nonboundary quadratic factors in `U`; adjoining compatible `(M,N)` signs gives generically eight first-half sheets over `(R,s)`. For one first-half sheet and fixed `eta`, `G(z)=Res_tau(E1,E2)` is an irreducible 64-term polynomial of degree 8 in `z`; `E2` is linear in `tau`, so `tau` is unique. The two `eta` branches give degree 16 per first-half sheet, hence generic degree 128 over `(R,s)` in either isomorphic collapsed `(delta,eps)` presentation. Further decomposition after imposing a selected quadratic `U` branch remains unproved.

The next exact object is one normalized component formed by adjoining a root of one quadratic `U` branch, compatible `(M,N)`, then a root of the compact degree-8 `G(z)`. To obtain a curve and genus, take a controlled one-coordinate fiber (for example fixed rational `R`) and normalize that degree-8 cover; the global target itself remains two-dimensional.

## Reproduction

```text
timeout 60s sage code/task5_cover_geometry.sage 64 88
```

No Magma computation was used, and all timed processes terminated.
