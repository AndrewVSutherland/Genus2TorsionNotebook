# A genuine contact-7 route toward cyclic 49

For the contact family
```
h = 1 - (7/2)x + alpha*x^2 + beta*x^3,
f = (h^2 + (x-1)^7)/x^2,
D7 = [(x-1), h(1)],
```
the marked class (D_7) has exact order 7 on the smooth generic fiber.

A generic degree-two Mumford class (Q=[U,V]) with (7Q=\pm D_7)
is encoded by the norm equation
```
A^2 - B^2*f = -(x-1)*U^7,
```
where (U=x^2+s x+t), (deg A\leq 7), and (B) is monic of
degree 5.  Recover (V=-A/B\pmod U).  This gives 15 coefficient
equations in 17 variables; keeping those equations uneliminated is
small (about 20 MB), whereas triangular symbolic substitution causes
severe expression swell.  The repeated-root specialization
(U=(x-p)^2) is the sublocus (Q=2P-K) mentioned in the motivating
example and is therefore retained by the norm equation.

An exact finite-group screen at (p=3) finds no smooth residue pair
for which (D_7\in 7J(\mathbf F_3)).  Thus a rational solution in
this normalized family must meet a bad or nonintegral (3)-adic
boundary.

For the two ordinary one-node branches
((\alpha,\beta)\equiv(0,1),(2,2)\pmod3), the component-group
thickness is
```
n = 7*v_3(h(1)).
```
The order-7 component of (D_7) is divisible by 7 only if (49\mid n).
The first possible layer is consequently (v_3(h(1))=7), parameterized
by
```
alpha = alpha0 + 3*u,
beta  = 5/2 - alpha + 3^7*s,   v_3(s)=0,
alpha0 in {0,2}.
```
The exact away-from-3 sieve through rational parameter height 12 tested
11,610 smooth pairs on each branch and left no survivor (23,220 total);
maximum resident memory was below 37 MB.

This is a rigorous exclusion only for that bounded, first ordinary-node
layer.  Remaining charts are deeper layers (v_3(h(1))=14,21,\ldots),
the (Q_5=0) branch at ((1,1)), the more singular intersection
((1,0)), and nonintegral/infinity parameter charts.

Reproduce with:
```
magma -b mode:="symbolic" code/contact7_to49_division.m
magma -b mode:="finite" primes:="3,5,11,13" code/contact7_to49_division.m
magma -b branch_a0:=0 height:=12 first_layer_only:=1 code/contact7_to49_boundary3_search.m
magma -b branch_a0:=2 height:=12 first_layer_only:=1 code/contact7_to49_boundary3_search.m
```
