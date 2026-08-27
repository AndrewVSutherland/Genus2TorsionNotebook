# V1/V2 component-aware local analysis (2026-07-18)

> **Follow-up:** the fixed-(T) strict-transform Jacobian was subsequently
> proved invertible and the complete V1 (13^5) mask was intersected with
> the (N=10000) MW search.  See
> `notes/target_22224_V1_depth5_MW10000_2026_07_18.md`.  That follow-up
> supersedes the cautious finite-depth discussion of the exceptional
> (T\equiv\pm1) branches below and excludes every (N=10000) survivor.

## Outcome

The two new repeated fibers behave very differently.

- **V2 is ruled out completely at (p=11).**  This is an intrinsic local
  obstruction, not a conclusion from the sampled 144-key contact bank.
  Every finite (T)-disk allowed by the full (A(2,2,2,8)) cover has good
  reduction and
  
  \[
  \#J(\mathbf F_{11})=128\not\equiv0\pmod 3.
  \]
  
  Hence none can carry rational 3-torsion.  The entire nonintegral/infinity
  disk is independently killed by nonsquare normalized cover units.

- **V1 remains genuinely live at (p=13).**  Its two finite disk types have
  open simultaneous contact-plus-cover lifts through (13^5).  Therefore the
  old projected mask cannot be upgraded to a local obstruction.  The next
  search should concentrate on V1, not V2.

The scaled projective representatives used throughout are

\[
\begin{aligned}
V1&=[-2178,2420,9075,-1470T^2],\\
V2&=[-1458,2268,7938,-2023T^2].
\end{aligned}
\]

The common scalings are fourth powers in the halving radicands, so they do
not change the square conditions at the primes under consideration.

## Separating the sampled bank from the intrinsic equations

The earlier projected computation asked whether the projective root
multiset appeared among a finite sample of direct-contact tangent keys.  It
was useful for triage, but absence from that bank is not a theorem.

The present computation instead uses the corrected cleared contact
equations themselves.  Put (M=L^2), let (e_i) be the elementary
symmetric functions of the four squared branch parameters, and set

\[
 P=4Me_1+12(U^2+v^2)-(M+3U)^2.
\]

The three equations are

\[
\begin{aligned}
F_1&=(M+3U)P+16v^3-8U^3-48Uv^2-8Me_2,\\
F_2&=P^2+64(M+3U)v^3-192(U^2v^2+v^4)-64Me_3,\\
F_3&=Pv^3-12Uv^4-4Me_4.
\end{aligned}
\]

The full cover is imposed simultaneously with square-root variables for

\[
\begin{aligned}
R_0&=abcd,\\
R_1&=a(a+b)(a+c)(a+d),\\
R_2&=b(b+a)(b+c)(b+d),\\
R_3&=c(c+a)(c+b)(c+d).
\end{aligned}
\]

Thus the reported lifts are solutions of the intrinsic contact and cover
systems, not matches against sampled keys.

## Resolving the universal (L=0) component

At (L=0), the cleared equations contain the universal component

\[
M=0,\qquad 2v-U=0,
\]

independently of the curve.  Counting these points as contact would be a
false positive.  Write

\[
 2v-U=Ls
\]

and divide (F_1,F_2,F_3) by (L^2).  On the exceptional divisor (L=0),
the strict transform is

\[
\begin{aligned}
0&=3s^2U-12U^2+12Ue_1-8e_2,\\
0&=3s^2U^2-8U^3+6U^2e_1-8e_3,\\
0&=3s^2U^3-6U^4+4U^3e_1-32e_4.
\end{aligned}
\]

These equations were derived symbolically and the full polynomial strict
transform is what was used for lifting.  The identity was also checked
directly after substituting (M=L^2, v=(U+Ls)/2).

## Complete infinity-disk tests

For (z=1/T), the infinity chart is

\[
[Az^2,Bz^2,Cz^2,D].
\]

After dividing each cover radicand by (z^6), the unit residues are:

| fiber, prime | normalized units | square pattern |
|---|---|---|
| V1, (p=13) | ((1,2,4,5)) | yes, **no**, yes, **no** |
| V2, (p=11) | ((4,2,5,8)) | yes, **no**, yes, **no** |
| V2, (p=13) | ((1,7,3,6)) | yes, **no**, yes, **no** |

Consequently all three infinity disks are impossible on the full cover.
This is complete for every (z\in p\mathbf Z_p\setminus\{0\}), and is
independent of the contact variables.

## V2: definitive (p=11) obstruction

The full cover permits only

\[
T\equiv\pm1,\ \pm5\pmod {11}.
\]

The corresponding radicands are

| (T) type | ((R_0,R_1,R_2,R_3)\bmod11) |
|---|---|
| (pm1) | ((4,1,4,9)) |
| (pm5) | ((1,5,3,3)) |

All are nonzero squares.  Direct enumeration gives no unit-(L) contact
point in any of these four disks.  More importantly, the first blow-up of
the universal (L=0) component has no exceptional point over any of them.

There is also a chart-independent proof.  The reduced branch tuples are

\[
[5,2,7,1]\quad(T=\pm1),\qquad
[5,2,7,3]\quad(T=\pm5),
\]

their squares are pairwise distinct and nonzero, and Magma gives

```text
T = +/-1:  #J(F_11) = 128
T = +/-5:  #J(F_11) = 128
```

Since (3\nmid128), prime-to-11 torsion injectivity rules out rational
3-torsion throughout all four residue disks.  Together with the infinity
unit obstruction this exhausts (\mathbf P^1(\mathbf Q_{11})).  V2 cannot
produce ([2,2,2,24]).

The apparent V2 (p=13) contact branches are therefore diagnostically
interesting but globally irrelevant.  They do lift to (13^2); the
(p=11) argument kills the fiber before any deeper (p=13) work is useful.

## V1: viable (p=13) components

The full cover permits

\[
T\equiv\pm1,\ \pm5\pmod {13}.
\]

The two disk types are

| (T) type | radicands mod 13 | unit-(L) contact | exceptional points |
|---|---|---|---|
| (pm1) | ((1,3,9,0)) | ((M,U,v)=(9,2,1)) | ((s,U)=(6,2),(7,2)) |
| (pm5) | ((12,12,1,3)) | ((M,U,v)=(9,2,1)) | ((s,U)=(6,2),(7,2)) |

The exact involutions (T\mapsto-T) and independent signs of the four
cover roots allow one representative of each sign orbit to be lifted.

### Exhaustive (13^2) layer

| mode | (T_0) | starts mod 13 | nodes mod (13^2) | (T\)-classes mod (13^2) |
|---|---:|---:|---:|---:|
| unit (L) | 1 | 8 | 1,352 | 1 |
| unit (L) | 5 | 16 | 2,704 | 13 |
| exceptional, each (s) | 1 | 8 | 104 | 1 |
| exceptional, each (s) | 5 | 16 | 208 | 13 |

For (T_0=1), the zero cover radicand forces one (T)-class at the first
lift.  For (T_0=5), every next (T)-digit occurs.

### Open (13^5) certificates

The unit-(L) components over both (T_0=1) and (T_0=5) have explicit
simultaneous contact-plus-cover solutions modulo (13^5) which have left
the curve-collision and gcd boundaries.  The exceptional components over
(T_0=5) also move off (L=0) and have open certificates through (13^5).

For the exceptional (T_0=1) points, exhaustive deterministic lifting did
not find a point that had left all boundaries by (13^5).  This does not
matter for local viability because the unit-(L) (T_0=1) component is
already open.  It should not be promoted to a formal impossibility without
a further blow-up.

These are finite-depth certificates on rank-deficient special fibers, not
by themselves a proof of an infinite (\mathbf Q_{13})-branch.  They are,
however, enough to show that the (p=13) projected mask is not an
obstruction to V1.

## Mordell--Weil comparison

The completed V1 rank-two sieve at (N=1001) left the single modular row

```text
(m,n,torsion_coset)=(-2,-1,1).
```

Exact map evaluation gives

\[
T=\frac{317596969778}{493651296635}\equiv-1\pmod {13},
\]

so it lands in a locally viable V1 disk.  Nevertheless all three exact
remaining cover ratios are nonsquares; it is not a global full-cover hit.

The six recorded V2 modular rows are all congruent to the same class
((m,n,ti)=(3,3,1)) modulo the (p=11) generator orders ((8,4)).  The
small representative ((3,-1,1)) maps to

\[
T=-\frac{2063497452257}{1698222432003}\equiv5\pmod {11}.
\]

Thus the generic sieve's conservative “map infinity” label at (p=11) is
a removable-map presentation issue: the exact point is finite.  It lies in
one of the disks with (#J(\mathbf F_{11})=128), and one exact cover ratio
already fails to be a square.

## Recommended next move

Stop all V2 work.  For V1, use the intrinsic (p=13) components as a
necessary local filter inside a larger Mordell--Weil search, but do not
expect (p=13) alone to cut strongly: the (T\equiv\pm5) components admit
all second digits.  The useful next computation is therefore an expanded
V1 MW lattice search (or an auxiliary-prime refinement) with exact cover
tests, prioritizing classes mapping to the four viable (p=13) disks.

## Files

```text
code/target_22224_V_components_analysis.py
code/target_22224_V_components_symbolic.m
code/target_22224_V_components_finitefield.m
code/target_22224_V_components_MW_compare.m
results/target_22224_V_components_analysis.log
results/target_22224_V_components_p2_components.tsv
results/target_22224_V_components_padic_certificates.tsv
results/target_22224_V_components_finitefield.log
results/target_22224_V_components_MW_compare.log
```
