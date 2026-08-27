# New height-80 target-54 search in the contact-9 root family

This is a genuinely new specialization sieve for rational `[3,18]` torsion.
The earlier height-80 runs used targets 36 and 72; neither is the necessary
order-54 condition.

## Files and commands

```text
data/audit_contact9_root_target54_h80_all_exact.txt
code/audit_contact9_target54_chart_counts.py
data/audit_contact9_target54_chart_counts_h80.txt
```

The recorded computations were:

```bash
timeout 600s magma -b mode:=root height:=80 target:=54 \
  progress_interval:=2000 max_exact:=20000 \
  code/contact9_family_search.m \
  > data/audit_contact9_root_target54_h80_all_exact.txt

python3 code/audit_contact9_target54_chart_counts.py --height 80 \
  > data/audit_contact9_target54_chart_counts_h80.txt
```

Both completed with exit status zero.  The Magma run took less than two
minutes on the current machine.

## Good-prime sieve

`PassesTargetReduction` only uses a prime after `f` reduces to a degree-five
squarefree polynomial.  Undefined or singular displayed reductions are
skipped.  At a good prime `p`, it requires the prime-to-`p` part of 54 to
divide `#J(F_p)`.  Hence every rejection is a valid good-reduction
obstruction.  Divisibility by 54 is only a necessary condition for
`Z/3 x Z/18`, so this sieve cannot create false hits.

The exact search box contains 7863 reduced values `s` of height at most 80
and both signs `eps`.  Its totals were

```text
checked             15726
smooth              15718
marked order 9      15718
target survivors        0
exact tests              0
hits                     0
```

The exact-test count is zero because no candidate survives the necessary
good-prime sieve; thus every survivor was exact-tested vacuously.  The first
obstructing-prime distribution begins

```text
p=7: 1958,  p=11: 5678,  p=13: 4706,  p=17: 2184,
p=19: 474,  p=23: 312,
```

and the remaining candidates are killed by primes through 71.  No rational
parameter reaches `TorsionSubgroup`.

## Mod-7 chart coverage

The effective q-cover parameter is `t=eps*s`.  Among the 15718 smooth marked
order-9 candidates, the exact chart counts are

```text
denominator divisible by 7: 2000
t mod 7 = 0: 2000
t mod 7 = 1: 1944
t mod 7 = 2: 1958
t mod 7 = 3: 1956
t mod 7 = 4: 1958
t mod 7 = 5: 1958
t mod 7 = 6: 1944
```

The 1958 candidates at `t=2 mod 7` are exactly the candidates first killed
at the only smooth displayed mod-7 fiber.  The sieve correctly skips the
singular or undefined charts at the other residues.

In particular, the search includes and eventually kills at later good
primes all

```text
1958 candidates in the live t=4 ordinary-node chart,
1944 candidates in the t=-1 mod-7 pole chart,
2000 candidates whose t denominator is divisible by 7.
```

Thus neither of the two locally viable mod-7 boundary charts produces a
rational `[3,18]` specialization within height 80.  This does not rule out
larger height or rational points on the generic degree-12/27 q-cover
components, but it closes the previously untested target-54 height box.
