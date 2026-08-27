# A(2,12) -> A(2,24) alternate scout

Date: 2026-07-02.

This is a bounded height-6 shell scout beyond the now-closed four height-4
fibers and the cold height-5 split scan.  It reuses the saturated
square-quartic halving projection, but defaults to

```text
Height:=6, ShellOnly:=true
```

so the active search is only the height-6 shell, not another pass through the
height-5 box.

New code:

```text
code/agent_A2_24_alt_scout.m
```

Logs:

```text
results/A2_24_alt_syntax.log
results/A2_24_alt_h6_shell_full.log
```

## Commands

Syntax/counting smoke:

```text
magma -b Height:=6 ShellOnly:=true MaxChecked:=200 Progress:=0 \
  code/agent_A2_24_alt_scout.m \
  > results/A2_24_alt_syntax.log 2>&1
```

Full height-6 shell:

```text
magma -b Height:=6 ShellOnly:=true Progress:=5000 \
  code/agent_A2_24_alt_scout.m \
  > results/A2_24_alt_h6_shell_full.log 2>&1
```

## Shell size

The height-6 box has

```text
full_enumerated_triples=99452
```

and the exact height-6 shell has

```text
full_active_triples=43136
```

The other `56316` triples are exactly the height-5 box and were skipped by
`ShellOnly:=true`.  Since the full shell was only 43,136 active triples, no
smaller partition was needed; the script still supports `PStart/PStop` and
prints partition estimates.

## Result

The full shell found only four split/order-12 fibers:

```text
p=-1/5, z=-6/5, r= 5/4
p=-1/5, z= 6/5, r= 5/4
p= 1/5, z=-6/5, r=-5/4
p= 1/5, z= 6/5, r=-5/4
```

All four have residual quartic factor type

```text
[ <1,1>, <1,1>, <2,1> ]
```

and eight rational 2-torsion translations.  The 32 translated order-12 rows
had saturated affine `M`-degree distribution:

| saturated affine `M` degrees | translated rows |
|---:|---:|
| `[16]` | 24 |
| `[8, 8]` | 8 |

Equivalently, the minimum affine factor degree distribution was:

| minimum degree | translated rows |
|---:|---:|
| 16 | 24 |
| 8 | 8 |

Final counters:

```text
active_checked=43136
split_fibers=4
order12_split_fibers=4
translated_order12_rows=32
low_rows_le_4=0
degree4_or_less_rows=0
rational_M_root_rows=0
rational_M_root_factors=0
rational_point_rows=0
exact_divisible_rows=0
torsion_cert_rows=0
positive_dim_rows=0
errors=0
```

## Verdict

No new height-6 shell branch has saturated degree `<=4`, and no saturated
affine branch of any degree has a rational `M`-root.  Therefore the height-6
split A(2,12) shell gives no new A(2,24) candidate and triggers no `N`-lift or
exact halving/torsion certification.

Recommendation: do not spend more time on naive split height-box expansion
unless it is heavily partitioned for background coverage.  The better next
move is a different chart or a descent/resolvent diagnostic that avoids
waiting for rare split fibers to expose low-degree saturated factors.
