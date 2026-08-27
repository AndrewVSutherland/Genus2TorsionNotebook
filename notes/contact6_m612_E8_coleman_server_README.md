# Server run: Prym-Chabauty sweep on E8 (the [6,12] gatekeeper)

One self-contained Magma script; needs Tuitman's Coleman library.

## Setup

```bash
git clone https://github.com/jtuitman/Coleman.git
cp torsion_jac/code/contact6_m612_E8_coleman_server.m Coleman/
cd Coleman
```

## Run

```bash
# main run (verified prime):
nohup magma -b p:=37 contact6_m612_E8_coleman_server.m > e8_p37.log 2>&1 &

# optional, in parallel (the only other verified prime; ~2-4x slower):
nohup magma -b p:=79 contact6_m612_E8_coleman_server.m > e8_p79.log 2>&1 &
```

Single-core each, < 8 GB.  p=37: roughly 12-24 h (the long poles are one
boundary-disk integral at expansion length e:=300, then 42 residue
disks at ~10-20 min each, with per-disk progress lines).

## What the script does

1. coleman_data on E8's integral plane model (genus 4).
2. Finds the 2 rational (boundary) points; transports the rank-1 Prym
   generator through the bigonal correspondence (all steps re-verified
   in-run: model-error valuations ~19-20, s^2 = 2A+2u consistency).
3. Pairs the two anti-invariant (Prym) differentials against the
   generator; forms the annihilating differential omega_ann.
4. Sweeps every residue disk for zeros of Int_b^P omega_ann (every
   rational point of E8 must be such a zero).
5. Post-processing: rationally reconstructs each zero and EXACTLY tests
   its fiber for rational points over Q.

## Reading the output

- `PAIRING basis[3] = -205372*37 + O(37^5)` and
  `PAIRING basis[4] = 473735*37 + O(37^5)` — at p=37 these are
  deterministic; matching them validates the run.
- `disk k/42 swept ...: n zero(s)` — progress.
- `TOTAL ZEROS ...` then the reconstruction:
  - `!!! RATIONAL POINT(S) ON E8 ... [6,12] CANDIDATE !!!` — the
    jackpot line: a nonboundary rational point of E8, i.e. a live
    [6,12] candidate fiber.  Send the log immediately.
  - otherwise every zero is classified `(boundary e=0)` or
    `(... mock zero)` — then E8(Q) = {2 boundary points} is supported
    at this prime, with any residual mock zeros to be killed by the
    second prime / Mordell-Weil sieve.

Known caveat: at p=79 the boundary integral may demand a longer local
expansion; if the run errors with "e is too small: try greater than X",
raise the one constant `ee := 300;` in the script above X and rerun.

Context and derivation: notes/contact6_m612_relative3_s3_quotients_2026_07_11.md.
