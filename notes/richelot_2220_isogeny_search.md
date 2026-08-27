# Richelot search around the simple `[2,2,20]` example

Seed curve: the contact-5 specialization

```text
t = -8233/7225,
h = 1 + t*x + ((t^2 - 1)/2)*x^2,
f = h^2 - ((t+1)^4/4)*x^5.
```

This is the known geometrically simple example with

```text
J(Q)_tors = [2,2,20].
```

The search script is

```text
code/richelot_2220_isogeny_search.m
```

and the saved output is

```text
data/richelot_2220_isogeny_search.txt
```

For exact torsion, the script square-scales each rational-coefficient
hyperelliptic model `y^2=f(x)` to an integral model `Y^2=d^2 f(x)`, which is
`Q`-isomorphic via `Y=d*y`.

## Immediate Richelot neighbors

Magma finds three rational Richelot codomains, all Jacobians:

```text
richelot_1: J(Q)_tors = [2,20],   branch factor degrees [2,2,2]
richelot_2: J(Q)_tors = [2,20],   branch factor degrees [2,2,2]
richelot_3: J(Q)_tors = [2,2,10], branch factor degrees [1,1,1,1,2]
```

So no direct Richelot neighbor improves on the seed torsion order `80`.

## Full 2-power traversal

I then used Magma's built-in

```text
TwoPowerIsogenies(J)
```

which traverses Richelot and double-Richelot steps over `Q`.

It returned:

```text
twopower_jacobians          5
twopower_products           0
twopower_weil_restrictions  0
```

The exact rational torsion on those five Jacobians is:

```text
twopower_1: [2,20],   order 40
twopower_2: [2,20],   order 40
twopower_3: [2,2,10], order 40
twopower_4: [2,10],   order 20
twopower_5: [2,10],   order 20
```

Thus the largest torsion order in the computed rational 2-power
Richelot/double-Richelot neighborhood is `40`, below the seed's order `80`.

Because every Jacobian found is `Q`-isogenous to the geometrically simple seed,
geometric simplicity is inherited for these Jacobian codomains.  However, none
of them gives larger rational torsion.

## Conclusion

The Richelot/isogeny-class strategy around this particular simple `[2,2,20]`
example does not find `[2,2,2,20]`, `[2,4,20]`, `[2,40]`, or anything larger.
In fact, all computed 2-power Richelot neighbors have smaller rational torsion.

This is useful negative information: the known simple `[2,2,20]` point appears
to be locally maximal inside its rational 2-power Richelot component.


## All double-linear seeds

I then checked whether there are other `[2,2,20]` seeds in the same
contact-5/order-20 double-linear locus before repeating the Richelot search.
The script is

```text
code/richelot_2220_all_double_linear_seeds.m
```

with output

```text
data/richelot_2220_all_double_linear_seeds.txt
```

The double-linear locus is controlled by

```text
Y^2 = (r+1)(r^2+2r+2)(r^3-r^2-4r+2).
```

For this auxiliary genus-2 curve, Magma gives

```text
RankBounds(J) = [1,1]
J(Q)_tors     = [4].
```

Using the infinite-order divisor

```text
(x^2 - x/3, -14*x/9 + 2, 2)
```

Chabauty returns exactly nine rational points:

```text
two points at infinity,
r=-2 with Y=+-2,
r=-1 with Y=0,
r=0  with Y=+-2,
r=1/3 with Y=+-40/27.
```

Only `r=1/3` is nondegenerate.  It gives

```text
z = -1/7,  w = -7/9,
t = -8233/7225,
J(Q)_tors = [2,2,20].
```

The other points are boundary, collision, or degenerate:

```text
r=-2: boundary
r=-1: z=w collision
r=0:  t=-1 degenerate
```

Thus there are no other `[2,2,20]` seeds in this double-linear locus to feed
into the Richelot search.  Re-running the Richelot/`TwoPowerIsogenies` test on
the complete seed list again leaves the original seed as the best object:

```text
best_order_including_seed 80
best_labels [seed_1_seed]
```


## Richelot sweep over many `[2,20]` sources

I next broadened the Richelot idea from the unique `[2,2,20]` seed to many
ordinary `[2,20]` sources in the same contact-5/order-20 family.  The script is

```text
code/richelot_contact5_extra2_sweep.m
```

It enumerates the two parametrized extra-2 loci:

```text
linear 1+3:
t = -(z^4 + 4*z + 4)/(z^4 + 4*z^3 + 8*z^2 + 8*z + 4)

quadratic-quadratic 2+2:
t = -(r^6 - 2*r^5 + 2*r^4 - 4*r^3 + 4*r^2 - 8*r + 8)
     /((r^2 - 2)^2*(r^2 - 2*r + 2)).
```

For each smooth source that exact-checks to torsion order at least `40`, it
computes:

```text
RichelotIsogenousSurfaces(J)
TwoPowerIsogenies(J)
```

and exact torsion on every Jacobian codomain.

### Complete height 7

Run:

```text
magma -b height:=7 max_sources_per_label:=1000 \
    code/richelot_contact5_extra2_sweep.m \
    > data/richelot_contact5_extra2_sweep_h7_complete.txt
```

Summary:

```text
checked_linear 71
checked_qq     71
tested_linear  69
tested_qq      46
unique_t       115
best_order     80
interesting_count 1
```

The only order-`80` object was again the known source:

```text
linear_63, z = -1/7, t = -8233/7225,
J(Q)_tors = [2,2,20].
```

No Richelot or double-Richelot codomain had order `80`, exponent greater than
`20`, or any new larger structure.

Aggregate torsion pattern:

```text
linear sources:
  [2,20]     67
  [2,2,20]    1

linear Richelot/TwoPower codomains:
  mostly [10], plus the known seed's old [2,20], [2,2,10], [2,10] neighbors

qq sources:
  [2,20]     46

qq Richelot/TwoPower codomains:
  Richelot:  [2,2,10] once per source
  TwoPower:  [2,2,10] once and [2,10] twice per source
```

### Capped height 20

Run:

```text
magma -b height:=20 max_sources_per_label:=100 \
    code/richelot_contact5_extra2_sweep.m \
    > data/richelot_contact5_extra2_sweep_h20_m100.txt
```

Summary:

```text
checked_linear 103
checked_qq     147
tested_linear  100
tested_qq      100
unique_t       200
best_order     80
interesting_count 1
```

Again, the only order-`80` object was the known source

```text
linear_63, z = -1/7, t = -8233/7225,
J(Q)_tors = [2,2,20].
```

All other tested sources had exact torsion `[2,20]`.  None of their Richelot or
`TwoPowerIsogenies` codomains improved torsion; the codomains had torsion among

```text
[10], [2,10], [2,2,10], [2,20].
```

I also attempted a complete height-`10` run, but it was interrupted at
`linear_107` because `TwoPowerIsogenies` on that larger-coefficient source was
slow.  The partial file is

```text
data/richelot_contact5_extra2_sweep_h10_complete.txt
```

and before interruption it had found no new interesting object beyond the same
known `[2,2,20]` source.

Conclusion: the broader Richelot sweep over many `[2,20]` sources did not find
a source whose rational Richelot or 2-power-isogenous neighbor gains torsion.
The observed behavior is instead torsion loss: `[2,20]` sources usually map to
codomains with exponent `10`, and the unique order-`80` source remains locally
maximal in the tested 2-power isogeny neighborhoods.
