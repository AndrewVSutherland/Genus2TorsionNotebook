# R = -8 MW-sieve attempt for the remaining ELS fiber

This records the MW-sieve style attack on the remaining `R=-8` ELS fiber after
the corrected `S_B` V4 scan killed `R=-25/4` and `R=-29/8`.

## Objects

The corrected `S_B` parameterization for `R=-8` gives

```text
X(lambda)  = (-2016*lambda - 1824)/(lambda^2 - 1)
S(lambda)  = (-1008*lambda^2 - 1824*lambda - 1008)/(lambda^2 - 1)
Y_B^2      = (1143072*lambda^3 + 3229632*lambda^2
              + 3002400*lambda + 919296)/(lambda^3-lambda^2-lambda+1).
```

All three V4 quotients have rank 1.  In particular the `m` quotient is

```text
E_m : y^2 = x^3 - 1263*x + 17138,
MW invariants [2, 2, 0].
```

The original C2 elliptic fiber is also the same minimal elliptic curve:

```text
g^2 = 15876*m^4 - 28957824*m^2 + 16131032064,
Emin : y^2 = x^3 - 1263*x + 17138,
MW invariants [2, 2, 0].
```

The map from the minimal model to the C2 `m`-coordinate is

```text
m = -28*(6*x + y - 24)/(2*x + y - 128).
```

## Scripts added

- `code/agent_m18_416_R8_fullstage_mwsieve.m`
  - finite MW sieve on the `E_m` lambda quotient;
  - tests the full A/B component conditions over finite fields.

- `code/agent_m18_416_R8_C2_fullstage_mwsieve.m`
  - finite MW sieve on the original C2 elliptic fiber;
  - uses the actual finite-field `m` and `sqrt(G)` data.

- `code/agent_m18_416_R8_local_stage_scan.m`
  - fixed-fiber local scan for `R=-8`, testing `m` classes mod `p^2`.

- `code/agent_m18_416_R8_C2_padic_residue_sieve.m`
  - exploratory p-adic residue sieve on the C2 MW group.

- `code/agent_m18_416_R8_qp_map_probe.m`
  - small probe used to extract/check the p-adic map formula above.

## Results

The larger local scan found no local obstruction up to `p=101`.
Typical witnesses:

```text
p=29 soluble, witness m=1
p=43 soluble, witness m=7
p=89 soluble, witness m=13
p=101 soluble, witness m=3
```

The conservative finite MW sieve on `E_m` does not close.  With primes through
`89`, the CRT list exceeds the `200001` cap:

```text
R=-8 full-stage MW sieve on E_m quotient
BranchMode=conservative MaxPrime=101
...
p=83 ... classes 27216 -> 136080
p=89 ... classes 136080 -> 200001 TRUNCATED
```

The conservative finite MW sieve on the original C2 fiber also does not close.
With primes through `83`, it likewise exceeds the cap:

```text
R=-8 C2 full-stage MW sieve
BranchMode=conservative MaxPrime=101
...
p=79 ... classes 45360 -> 45360
p=83 ... classes 45360 -> 200001 TRUNCATED
```

The p-adic residue probes show the same high-density behavior rather than a
hidden quick contradiction.  For example:

```text
p=11, Depth=2: tested=2904, survivors=2904
p=17, Depth=2: tested=11560, survivors=11560
p=29, Depth=2: tested=60552, survivors=40032
p=31, Depth=2: tested=53816, survivors=46002
p=43, Depth=2: tested=147920, survivors=103197
```

The p-adic script is conservative around non-unit chart denominators, so these
numbers should be read as a sieve diagnostic, not as a Chabauty certificate.

## Conclusion

The straightforward MW-sieve style argument has been carried out in the two
natural rank-one ways:

1. sieve on the `S_B`/`m` quotient `E_m`;
2. sieve on the original C2 elliptic fiber, retaining the actual `sqrt(G)`.

Neither closes the `R=-8` fiber.  This is useful negative information: the
remaining case is not a small-prime or plain finite-MW-sieve artifact.  A proof
for `R=-8` will need an additional global input, most likely a genuine
Chabauty/cover-descent computation on the full second-stage cover, or a more
structured use of the A/B fiber product than independent residue sieving.

