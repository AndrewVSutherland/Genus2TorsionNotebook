# M(2,2,4,8): HPL normalization and direct cover search

This records the next large-torsion pass after the negative `5`-torsion and section-family searches.

## Why this route

The current best large-simple target remains a `2`-power lift toward `M(2,2,4,8)`.  The direct K3 height searches are strongly obstructed at `11` and `23`, and the low-degree AT/BL section-family searches did not produce witnesses.  The concrete remaining lead in the notes is the Howe--Poonen--Leprevost split `(2,2,4,8)` example, known to be non-simple but useful as a base point for the equations.

## HPL normalization check

The script

```text
code/m2248_hpl_normalization_check.m
```

tests all ordered choices of two Weierstrass points to move to `0` and `infinity`.  This matters because a raw tuple `[a,b,c,d]` for

```text
y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2)
```

need not already be in the normalized `M2248` chart.

For both HPL tuples appearing in `paper/NotesAndTodo.tex`, the `use_squares=true` interpretation gives:

```text
square_normalizations 6
intermediate_sources 4
intermediate_witnesses 64
full_sources 4
full_witnesses 32
```

The `use_squares=false` interpretation gives no normalized square model.  Thus the HPL examples are recognized by the current full `M(2,2,4,8)` equations once normalized correctly.  Feeding the raw tuple directly to `m2248_sieve.m` gives zero only because it has not chosen the correct `0` and `infinity` branch points.

A compact output summary is in

```text
data/m2248_hpl_normalization_summary.txt
```

## Section-family searches

The best symbolic AT section families were `n=1` with torsion translates `T2` and `T3`.  The best BL family was `n=1`, torsion translate `T`.  Exact searches gave no intermediate witnesses:

```text
magma -b search_output_file:=data/m2248_at_n1_intermediate_h200.txt   max_multiple:=1 height:=200 full_check:=false   code/m2248_section_multiples.m

Total sources: 0
Total witnesses: 0
```

```text
magma -b bl_search_output_file:=data/m2248_bl_n1_intermediate_h200.txt   max_multiple:=1 height:=200 full_check:=false   code/m2248_bl_multiples.m

BL normalized tuples: 97848
BL total sources: 0
BL total witnesses: 0
```

## Direct rho-sigma-tau search

The HPL witnesses suggest searching directly in the cover coordinates `(rho,sigma,tau)`.  The first script is

```text
code/m2248_rst_direct_search.m
```

It sets `B=1`, `A=rho^2`, and reconstructs

```text
C = (sigma^2 - A)/(sigma^2 - 1),
D = (tau^2   - A)/(tau^2   - 1).
```

It then requires `C,D` and the four full-cover expressions `F1..F4` to be rational squares.  Through height `40` it found cover-square triples but no cases with both `C,D` square:

```text
height 40:
checked 13953516
cover_square 68
cd_square 0
hits 0
```

The improved script

```text
code/m2248_rst_param_search.m
```

parametrizes the `C=c^2` and `D=d^2` conditions automatically.  Through height `12` it tested `5860440` valid reconstructed tuples and found no full-cover square hit:

```text
height 12:
checked 5910632
valid_param 5860440
cover_square 0
hits 0
```

The summary is in

```text
data/m2248_rst_direct_search_summary.txt
```

## Conclusion

The useful positive result is that the known HPL `(2,2,4,8)` points are now correctly located in the normalized `M2248` equations.  The negative result is that small-height direct searches in `(rho,sigma,tau)` do not find new points, simple or otherwise.

The next serious step should be algebraic: use one HPL witness as a base point and parametrize or elliptically fiber the four full-cover square conditions `F1..F4`.  Continuing naive boxes in either K3 tuples or `(rho,sigma,tau)` is unlikely to be efficient.


## HPL Slice Attempts

I then tried to turn the HPL point into an actual search surface by slicing through it.

The first slice keeps the full-cover condition

```text
F1 = (1+rho)(1+sigma)(1+tau)
```

square by construction.  It writes

```text
sigma = sigma0*(1 + ds*u),
tau   = tau0  *(1 + dt*u),
q1    = q10   *(1 + dq*u),
rho   = q1^2/((1+sigma)(1+tau)) - 1.
```

The script is

```text
code/m2248_hpl_f1_slice_search.m
```

A run with `dir_bound=5`, `height=30` gave

```text
checked 1476300
valid_slice 1472980
branch_square 0
cover_square 0
hits 0
```

The stronger slice also parametrizes `C=c^2` automatically using the conic formula for `(rho,c)` in terms of `(sigma,m)`, while still forcing `F1` square.  The script is

```text
code/m2248_hpl_c_f1_slice_search.m
```

The same bounds gave

```text
checked 1476300
valid_slice 1472960
d_square 0
cover_square 0
hits 0
```

Finally, I diagnosed the genus of the first remaining square condition `D=d^2` on these HPL line slices:

```text
code/m2248_hpl_slice_genus_diagnostic.m
```

For all directions with `dir_bound=5`, the lowest squarefree degree was `6`, i.e. genus `2`:

```text
rank 1 degree 6 genus 2 dir [ 0, 0, -5 ]
...
rank 10 degree 6 genus 2 dir [ 0, 0, 5 ]
rank 11 degree 12 genus 5 dir [ 0, -5, 0 ]
```

So there is no genus-one line slice of this type through HPL in the tested small directions.  This explains the negative numeric slice searches: even after forcing `C=c^2` and `F1` square, the first remaining square condition is already genus `2` in the best line slices.

The practical conclusion is that a simple line/tangent slice through HPL is not the promised elliptic fibration.  The next algebraic move would need a richer fibration that makes at least two of `D=d^2`, `F2`, `F3`, `F4`, and `rho*sigma*tau` automatic before searching rational points.


## Fixed-rho HPL elliptic fibration

I then corrected the attempted `D=d^2` parametrization.  For fixed `rho`, the condition

```text
D = d^2 = (tau^2-rho^2)/(tau^2-1)
```

is not a conic in `tau`.  After clearing denominators it is the genus-one quartic

```text
Y^2 = (tau^2-rho^2)(tau^2-1),   Y = d*(tau^2-1).
```

At the normalized HPL value

```text
rho0  = 58466134224/53109477625
sigma0= 719363573659505664/749082246897952705
tau0  = 307598400/352612321
c     = 58466134224/30294861575
d     = 72946054224/53109477625
```

Magma identifies this quartic as an elliptic curve with

```text
rank_bounds        2 2
torsion_invariants [ 2, 4 ]
```

This is the first genuinely promising structure found in the HPL neighborhood: it is a rank-2 elliptic fibration through the known `(2,2,4,8)` point, not a blind height box.

The script

```text
code/m2248_hpl_fixedrho_elliptic_search.m
```

uses the HPL point as the identity on the quartic, searches bounded combinations of the obvious rational points, and keeps

```text
sigma*tau/rho0 = sigma0*tau0/rho0
```

fixed so the intermediate condition `F0=rho*sigma*tau` is automatic.  With `coeff_bound=20` it produced:

```text
checked            17000
distinct_E         417
affine_pullbacks   413
valid_cover_coords 413
c_square           3
cover_square       1
exact_sources      1
simple_hits        0
hits               1
```

The one exact full-cover source is precisely the known HPL point again, with nonsimple Jacobian.

I also used the simpler 2-isogenous model

```text
V^2 = U(U-1)(U-rho0^2),   U=tau^2, V=tau*Y.
```

The script

```text
code/m2248_hpl_fixedrho_legendre_lift_search.m
```

searches multiples of the HPL image on this Legendre model and lifts those with `U` square.  The model again has

```text
rank_bounds        2 2
torsion_invariants [ 2, 2 ]
P0_torsion_order_le_32 0
```

so the HPL image is a real non-torsion direction.  Through `max_multiple=40` the run gave:

```text
checked       323
distinct_E0   323
u_square      322
valid_lifts   640
c_square      4
cover_square  1
exact_sources 1
simple_hits   0
hits          1
```

Again, the only full-cover hit is HPL.  A higher exploratory run to `max_multiple=200` became dominated by very large rational arithmetic and printed no additional `C_SQUARE` event before being stopped, so the saved data file records the clean `max_multiple=40` run.

The conclusion is mixed but useful.  The fixed-`rho` HPL surface really does contain a rank-2 elliptic object, so the route is mathematically sensible.  However, the residual condition `C=c^2` is already very restrictive on the tested non-torsion lifts, and the remaining `F1..F4` full-cover squares cut it back to the original nonsimple point.  The next improvement should be to find an independent low-height rank direction on the Legendre model, or to vary `rho` in a second elliptic direction, rather than pushing single-generator multiples to enormous height.


## Independent rank direction from 2-descent

I next tried to get a second rank direction on the Legendre model instead of pushing multiples of the HPL point.  A direct call to `Generators(E0)` did not finish within a 120 second timeout:

```text
code/m2248_hpl_legendre_mw_probe.m
data/m2248_hpl_legendre_mw_probe.txt
```

The useful method was Magma `TwoDescent` with torsion and the known HPL point removed:

```text
TwoDescent(E0 : RemoveTorsion := true, RemoveGens := {P0})
```

This completed quickly and produced seven 2-coverings.  A search to bound `1000` found rational points on cover 2, mapping to

```text
P1 = (185261445034489541/75208331264762500,
      2019033491444667568388406032021/950637494391639082835996875000)
```

on

```text
E0: V^2 = U(U-1)(U-rho0^2).
```

Since this point comes from the descent quotient after removing `P0` and torsion, it gives the desired independent rank direction.  The probe output is in

```text
code/m2248_hpl_legendre_twodescent_probe.m
data/m2248_hpl_legendre_twodescent_probe.txt
```

I then searched the rank-2 lattice

```text
a*P0 + b*P1 + T,   T in E0[2],
```

keeping points with `U` square so they lift to the fixed-rho quartic, and then applying the same residual `C=c^2` and full-cover tests.  The script is

```text
code/m2248_hpl_fixedrho_legendre_rank2_search.m
```

The run with `coeff_bound=25` gave:

```text
checked       10400
distinct_E0   10400
u_square      5096
valid_lifts   10192
c_square      8
cover_square  1
exact_sources 1
simple_hits   0
hits          1
```

The only exact full-cover source was again the known HPL point.  This is an important negative refinement: even after using both rank directions on the fixed-rho elliptic fibration, the residual `C=c^2` condition remains extremely restrictive and no new simple `(2,2,4,8)` Jacobian appears in the tested rank-2 box.

The next logical move is therefore not larger coefficients on this fixed-rho lattice.  It is to let `rho` vary, ideally by building a two-parameter fibration where the fixed-rho rank directions are paired with a second square condition, or by running a modular obstruction analysis on the rank-2 lattice to explain why `C=c^2` collapses to HPL.


## Modular sieve on the fixed-rho rank-2 lattice

I then ran the proposed modular obstruction analysis on the rank-2 lattice

```text
Q(a,b,T) = a*P0 + b*P1 + T,   T in E0[2].
```

The script is

```text
code/m2248_hpl_rank2_mod_sieve.m
```

It reduces the Legendre model modulo good primes, enumerates residue classes for `(a,b,T)` modulo the orders of `P0` and `P1`, and applies the necessary finite-field conditions:

```text
U square,
C = (sigma^2-rho0^2)/(sigma^2-1) square,
F1,F2,F3,F4 square.
```

An important sanity check is needed: primes where the known HPL point does not survive the nonzero-square chart must not be used as obstructions.  The run identified

```text
p = 11, 41, 61
```

as bad for this chart in that sense: they kill even the known HPL class, so they are skipped in the CRT combination.

A clean run through `p <= 211` is in

```text
data/m2248_hpl_rank2_mod_sieve_p211.txt
```

The final combined result is:

```text
used_primes [19, 29, 43, 59, 67, 71, 73, 79, 83, 89, 97, 101,
             103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
             157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211]
final_mod_a 136936800
final_mod_b 12448800
final_full_survivors 7904
final_truncated false
```

Thus small good primes do not force failure of the fixed-rho rank-2 lattice.  They thin the search very strongly, but they leave many CRT-compatible full-cover residue classes beyond the HPL classes.  A longer run through `p <= 300` is in

```text
data/m2248_hpl_rank2_mod_sieve_p300.txt
```

That run remains nonempty but hits the `max_crt_classes` cap after `p=223`, because the coefficient moduli jump and many compatible residue classes lift.

Conclusion: the negative rational searches on fixed `rho` are not explained by a simple small-prime local obstruction on the rank-2 lattice.  The obstruction is either genuinely global, height-related, or tied to the residual rational square condition in a way the finite-field square tests do not see.  This makes a blind larger coefficient search unattractive, but it also means we should not claim the fixed-rho fibration is locally impossible.


## Rho-varying D+F0 rational surface

The fixed-rho searches were not locally impossible, but they also did not produce new rational points.  I therefore moved to a rho-varying surface that still makes two conditions automatic.

The key observation is that `D=d^2` can be parametrized rationally if `d` is treated as a parameter.  For fixed `d`,

```text
rho^2 = (1-d^2)*tau^2 + d^2
```

is a conic in `(rho,tau)` with boundary point `(1,1)`.  Taking the line

```text
rho = 1 + n*(tau-1)
```

gives a two-parameter parametrization in `(d,n)`.  Then setting

```text
sigma = q^2*rho/tau
```

makes

```text
F0 = rho*sigma*tau = (q*rho)^2
```

automatic.  The remaining tests are `C=c^2` and the four full-cover square conditions `F1..F4`.

The absolute-height script is

```text
code/m2248_df0_surface_search.m
```

It verifies the HPL point in these coordinates, where the parameter heights are already very large:

```text
HPL_param_heights d 72946054224 n 53109477625 q 352612321
```

The runs were:

```text
height 6:
  checked 87308
  valid_surface 1898
  c_square 3820
  cover_square 0
  hits 0

height 10:
  checked 1926036
  valid_surface 15286
  c_square 30824
  cover_square 0
  hits 0

height 12:
  checked 5875324
  valid_surface 32282
  c_square 64992
  cover_square 0
  hits 0
```

The data files are:

```text
data/m2248_df0_surface_h6.txt
data/m2248_df0_surface_h10.txt
data/m2248_df0_surface_h12.txt
```

This rho-varying surface finds many `C=c^2` branches but no full-cover point at low absolute parameter height.  That confirms the bottleneck has moved to `F1..F4`, not to `C`.

Since the HPL parameters are far outside small height, I also searched multiplicative line slices through HPL in the `(d,n,q)` coordinates:

```text
d = d0*(1 + ed*u)
n = n0*(1 + en*u)
q = q0*(1 + eq*u)
```

The script is

```text
code/m2248_hpl_df0_slice_search.m
```

The runs were:

```text
dir_bound 3, height 30:
  checked 379620
  valid_slice 378858
  c_square 18
  cover_square 18
  exact_sources 1
  simple_hits 0
  hits 1

dir_bound 5, height 30:
  checked 1476300
  valid_slice 1472990
  c_square 30
  cover_square 30
  exact_sources 1
  simple_hits 0
  hits 1
```

The data files are:

```text
data/m2248_hpl_df0_slice_d3_h30.txt
data/m2248_hpl_df0_slice_d5_h30.txt
```

Every reported `C_SQUARE` in these HPL-centered slices is just the original HPL point recovered by sign-return directions, for example changing `d` or `q` by a factor `-1` while leaving the actual `rho,sigma,tau` unchanged.  No genuinely varied point appeared.

Conclusion: allowing `rho` to vary while forcing `D` and `F0` is the right algebraic move, but the tested rational surface still does not produce new full-cover points.  The full-cover expressions `F1..F4` are now clearly the main obstruction.  A next serious step would need to force one of `F1..F4` by construction on this surface, rather than searching it as a residual square condition.


## Forcing F1 on the rho-varying D+F0 surface

The previous rho-varying surface made `D=d^2` and `F0` automatic, but still tested all of `F1..F4` as residual square conditions.  I next forced one of those conditions, namely `F1`, by using its conic structure.

For fixed `(d,n)`, hence fixed `(rho,tau)`, and with

```text
sigma = q^2*rho/tau,
```

the condition

```text
F1 = (1+rho)(1+sigma)(1+tau)
```

becomes the conic

```text
y^2 = A*(1+B*q^2),
A = (1+rho)(1+tau),  B = rho/tau.
```

The script first searches for a small seed point `(q0,y0)` on this conic, then parametrizes the conic by lines of slope `m` through that seed.  This makes `D`, `F0`, and `F1` automatic before testing `C,F2,F3,F4`.

The script is

```text
code/m2248_df0_f1_conic_search.m
```

Runs:

```text
dn_height 6, seed_height 5, slope_height 6:
  dn_checked 2116
  valid_surface 1898
  f1_seed_dns 76
  f1_seeds 152
  f1_param_points 7128
  c_square 36
  f234_square 0
  hits 0

dn_height 8, seed_height 6, slope_height 8:
  dn_checked 7396
  valid_surface 6994
  f1_seed_dns 202
  f1_seeds 404
  f1_param_points 35108
  c_square 54
  f234_square 0
  hits 0

dn_height 10, seed_height 8, slope_height 10:
  dn_checked 15876
  valid_surface 15286
  f1_seed_dns 382
  f1_seeds 764
  f1_param_points 96968
  c_square 128
  f234_square 0
  hits 0
```

The data files are:

```text
data/m2248_df0_f1_conic_h6_s5_m6.txt
data/m2248_df0_f1_conic_h8_s6_m8.txt
data/m2248_df0_f1_conic_h10_s8_m10.txt
```

Conclusion: forcing `F1` is computationally effective and produces many rational points on the `D+F0+F1` locus, but the remaining conditions `F2,F3,F4` still eliminate every tested candidate.  This is a stronger negative result than the earlier `D+F0` surface search: the obstruction is not just that `F1` was rare.  At least one of `F2,F3,F4` must be built in next, or the search should pivot away from this `M(2,2,4,8)` HPL neighborhood.


## Forcing F2 and F4 on the rho-varying D+F0 surface

After forcing `F1`, I repeated the same conic construction for `F2` and `F4`.  With

```text
sigma = q^2*rho/tau,
```

and fixed `(rho,tau)`, the conditions `F2` and `F4` also become conics in `q` after removing automatic square factors:

```text
F2/rho^2 = (1+rho)(rho+tau) * (1 + q^2/tau),
F4/tau^2 = (1+tau)(rho+tau) * (1 + rho*q^2/tau^2).
```

The target-selectable script is

```text
code/m2248_df0_force_fi_conic_search.m
```

It supports `target := "F1"`, `"F2"`, or `"F4"`.  I ran the new `F2` and `F4` targets with the same broad bounds used for the forced-`F1` comparison:

```text
target F2, dn_height 10, seed_height 8, slope_height 10:
  dn_checked 15876
  valid_surface 15286
  target_seed_dns 356
  target_seeds 712
  target_param_points 90266
  c_square 208
  residual_square 0
  hits 0

target F4, dn_height 10, seed_height 8, slope_height 10:
  dn_checked 15876
  valid_surface 15286
  target_seed_dns 1220
  target_seeds 2440
  target_param_points 309732
  c_square 248
  residual_square 0
  hits 0
```

The data files are:

```text
data/m2248_df0_force_f2_conic_h10_s8_m10.txt
data/m2248_df0_force_f4_conic_h10_s8_m10.txt
```

This is now a fairly strong negative result for the HPL-centered `M(2,2,4,8)` strategy.  On the rho-varying `D+F0` surface, each individually conic full-cover condition (`F1`, `F2`, or `F4`) can be forced effectively.  In all three cases, many rational points and many `C=c^2` branches are found, but the remaining cover equations eliminate every tested candidate.  The difficult obstruction is not isolated in one of `F1`, `F2`, or `F4`; it is in their simultaneous compatibility, very likely involving `F3` or the genus-one/fiber-product intersection of these conditions.

The practical conclusion is that another blind box in this HPL neighborhood is not the best use of time.  The next large-torsion move should either study the simultaneous `F1/F2` or `F2/F4` genus-one fibers symbolically, or pivot to a different high-torsion family rather than continuing ad hoc searches near the known nonsimple HPL point.

## Simultaneous pair-fiber study

I next studied the simultaneous `Fi/Fj` fibers on the rho-varying `D+F0` surface.  For fixed `(d,n)`, hence fixed `(rho,tau)`, write

```text
sigma = q^2*rho/tau.
```

Each of `F1,F2,F4` has the conic form

```text
y^2 = A_i*(1+B_i*q^2).
```

Given a seed point `(q0,y0)` on the first conic and `K=A_1*B_1`, parametrizing by slope `m` gives

```text
q(m) = (q0*m^2 - 2*y0*m + q0*K)/(m^2 - K).
```

Imposing a second conic gives the quartic

```text
W^2 = A_2*((m^2-K)^2 + B_2*(q0*m^2 - 2*y0*m + q0*K)^2).
```

This quartic is not arbitrary.  It is reciprocal:

```text
a0 = a4*K^2,  a1 = a3*K.
```

Equivalently, with `z=m+K/m` and `W=m*Z`, the pair fiber is the fiber product of two conics:

```text
Z^2 = A_2*((z^2 - 4*K) + B_2*(q0*z - 2*y0)^2),
u^2 = z^2 - 4*K.
```

Thus the simultaneous `Fi/Fj` locus is generically genus one, not a rational surface.  The script for this diagnostic is

```text
code/m2248_df0_pair_fiber_diagnostic.m
```

At the HPL point, all three pair fibers are genuine genus-one quartics:

```text
HPL F1/F2: degree 4, squarefree degree 4, genus 1
HPL F1/F4: degree 4, squarefree degree 4, genus 1
HPL F2/F4: degree 4, squarefree degree 4, genus 1
```

The HPL symbolic check is in

```text
code/m2248_pair_fiber_symbolic_hpl.m
data/m2248_pair_fiber_symbolic_hpl.txt
```

Magma converts all three HPL pair fibers to elliptic curves, each with torsion invariants `[2,2]`.  A trial `RankBounds` run on these large models was interrupted after about a minute, so I did not use rank data here.

The broader diagnostic run

```text
data/m2248_df0_pair_fiber_diag_h10_s8_m10.txt
```

gave:

```text
F1/F2:
  dn_checked 15876
  valid_surface 15286
  seeds 382
  degree 3 fibers 2
  degree 4 fibers 380
  genus 1 fibers 382
  pair_slope_points 100
  c_square_unique 22
  all_cover_square 0

F1/F4:
  seeds 382
  degree 4 fibers 382
  genus 1 fibers 382
  pair_slope_points 72
  c_square_unique 0
  all_cover_square 0

F2/F4:
  seeds 356
  degree 4 fibers 356
  genus 1 fibers 356
  pair_slope_points 80
  c_square_unique 0
  all_cover_square 0
```

So the simultaneous fibers really are genus-one fibers in the sampled range.  The only pair that produced any `C=c^2` points was `F1/F2`, and every reported `C`-square point had `q=1`.

This led to a useful symbolic explanation.  On the diagonal `q^2=1`, one has

```text
sigma = rho/tau,
F0 = rho^2,
C = rho^2/d^2,
F2 = rho^2*F1.
```

Thus the `F1/F2/C` successes in the diagnostic are dependent: `C` is automatic on the `D` surface, and `F2` is the same square class as `F1`.  The remaining independent full-cover conditions reduce to `F1` and `F4`.

In the `(d,n)` parametrization of the `D` surface,

```text
rho = -(2*d^2*n - d^2 + n^2 - 2*n + 1)/(d^2+n^2-1),
tau = -((d-n+1)*(d+n-1))/(d^2+n^2-1).
```

On `q^2=1`, the two remaining cover expressions factor as

```text
F1 = -8*d^2*n^2*(d-1)*(d+1)*(n-1)^2
     /((d-n+1)*(d+n-1)*(d^2+n^2-1)^2),

F4 = 8*d^2*n^2*(n-1)^2
     *(d^4 + d^2*n^2 + 2*d^2*n - 2*d^2 + n^2 - 2*n + 1)
     /(d^2+n^2-1)^4.
```

The targeted diagonal script is

```text
code/m2248_q1_diagonal_search.m
data/m2248_q1_diagonal_h30.txt
```

At `dn_height 30` it found many individual square events but no simultaneous full-cover point:

```text
dn_checked 1232100
valid_surface 1227002
f1_square 916
f4_square 1410
both_square 0
unique_full 0
simple_certified 0
```

I also saved a small affine good-reduction residue check:

```text
code/m2248_q1_diagonal_local.py
data/m2248_q1_diagonal_local.txt
```

For the simplified diagonal square classes, the affine good chart has no simultaneous residues at `p=3,5,13`; other tested primes do have residues.  This is not a complete p-adic proof because rational points could reduce to bad chart/boundary residues at those primes, but it explains why the diagonal is arithmetically tight.

Conclusion: the symbolic picture supports the negative computational evidence.  Pairing two of `F1,F2,F4` generically puts us on a genus-one fiber.  The only low-height route that also passes `C` collapses to the special `q^2=1` diagonal, where `F1` and `F2` are dependent and the remaining `F4` square class appears incompatible with `F1` at good affine residues for small primes.  This makes the HPL-centered `M(2,2,4,8)` path look increasingly obstructed rather than merely under-searched.

