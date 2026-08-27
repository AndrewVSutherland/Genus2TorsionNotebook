# Infinitude census: which torsion groups have provably infinitely many realizations (2026-07-31)

Ari's question: for which of the 72 groups in paper/torsion_realizations.tex can we PROVE
infinitely many realizations? Levels: **L1** = infinitely many geometrically simple
genus-2 Jacobians /Q with torsion CONTAINING G; **L2** = with torsion EXACTLY G;
**L3** = L2 with End(Jac_Qbar) = Z. Statuses: **PROVEN** (written theorem, verified),
**PMS** = provable-modulo-standard (a committed family + committed exact/simple member
exist; what remains is assembling the standard specialization argument, with the exact
missing lemma named in the workflow records), **OPEN** (sporadic witnesses only, or a
genuine mathematical gap). Produced by an 8-agent workflow (framework + six group batches
+ adversarial verifier); every PROVEN claim was independently re-verified (chart
identities recomputed, commits and logs checked). Full per-group proof sketches and
gaps: the workflow journal (session artifact); framework reference sheet R1-R11
summarised at the end.

## Headline

- **PROVEN at full strength (L2+L3): exactly one group — [2,2,14]**, by the lab's own
  Theorem C/C+ (commit d36bb18: 7898/7920 congruence classes, each an infinite set of
  members with exact torsion [2,2,14] AND End = Z, strict disjoint-prime certificates
  (43,307) classwide). This is currently the ONLY torsion group with an unconditional
  infinitude theorem at the exact+generic level.
- **PROVEN at L1: precisely the subgroups of [2,2,14]** — [2], [2,2], [2,2,2], [7],
  [14], [2,14], [2,2,14] — all inherited from Theorem C by containment.
- **PMS at L2/L3** (the realistic to-write list; each needs only the R1-sandwich +
  Cadoret-Tamagawa R3a assembly over a committed chart with a committed exact member):
  small cyclic [3]-[14], [18], [20], [22], [23], [30], [32], and
  [2,4], [2,8], [2,12], [2,14], [2,20], [3,3], [4,4], [2,2,4], [2,2,8],
  [2,2,2,2], [2,2,2,4], [2,2,2,10] (Elkies 2024).
- **OPEN everywhere** (only finitely many witnesses known): [17], [19], [21] (Leprevost's
  classical families fail simplicity scrutiny), [24]-[29], **[31]** (one RM witness),
  [33]-[40], [2,16], [2,18], **[2,22]** (two RM moduli points; finiteness is the live
  hypothesis), [2,26], [2,28], [3,6], [3,9], [4,8], [2,2,6]+, [2,2,10]+, [2,2,12],
  [2,2,20], [2,4,8], [6,6], [2,4,4], [2,2,2,6], [2,2,2,8], [2,2,4,4], [2,2,2,12]
  (provably stuck at three known curves as of the orbit census).
  (+ = L1 is PMS there; exactness is the gap.)

## Full table (verifier-stamped)

group | L1 | L2 | L3 | confidence
---|---|---|---|---
[2] | PROVEN | PMS | PMS | high
[3] | PMS | PMS | PMS | high
[4] | PMS | PMS | PMS | high
[5] | PMS | PMS | PMS | high
[6] | PMS | PMS | PMS | high
[7] | PROVEN | PMS | PMS | high
[8] | PMS | PMS | PMS | medium
[9] | PMS | PMS | PMS | high
[10] | PMS | PMS | PMS | high
[11] | PMS | PMS | PMS | high
[12] | PMS | PMS | PMS | high
[13] | PMS | PMS | PMS | high
[14] | PROVEN | PMS | PMS | high
[15] | PMS | OPEN | OPEN | medium
[16] | PMS | OPEN | OPEN | medium
[17] | OPEN | OPEN | OPEN | high
[18] | PMS | PMS | PMS | high
[19] | OPEN | OPEN | OPEN | high
[20] | PMS | PMS | PMS | high
[21] | OPEN | OPEN | OPEN | high
[22] | PMS | PMS | PMS | high
[23] | PMS | PMS | PMS | high
[24] | OPEN | OPEN | OPEN | high
[25] | OPEN | OPEN | OPEN | high
[26] | OPEN | OPEN | OPEN | high
[27] | OPEN | OPEN | OPEN | high
[28] | OPEN | OPEN | OPEN | medium
[29] | OPEN | OPEN | OPEN | high
[30] | PMS | PMS | PMS | high
[31] | OPEN | OPEN | OPEN | high
[32] | PMS | PMS | PMS | high
[33] | OPEN | OPEN | OPEN | high
[34] | OPEN | OPEN | OPEN | high
[36] | OPEN | OPEN | OPEN | high
[39] | OPEN | OPEN | OPEN | high
[40] | OPEN | OPEN | OPEN | high
[2,2] | PROVEN | PMS | PMS | high
[2,4] | PMS | PMS | PMS | medium
[2,6] | PMS | OPEN | OPEN | medium
[2,8] | PMS | PMS | PMS | high
[2,10] | PMS | OPEN | OPEN | high
[2,12] | PMS | PMS | PMS | high
[2,14] | PROVEN | OPEN | OPEN | high
[2,16] | OPEN | OPEN | OPEN | high
[2,18] | OPEN | OPEN | OPEN | high
[2,20] | PMS | PMS | PMS | high
[2,22] | OPEN | OPEN | OPEN | high
[2,26] | OPEN | OPEN | OPEN | high
[2,28] | OPEN | OPEN | OPEN | high
[3,3] | PMS | PMS | PMS | high
[3,6] | OPEN | OPEN | OPEN | high
[3,9] | OPEN | OPEN | OPEN | high
[4,4] | PMS | PMS | PMS | high
[4,8] | OPEN | OPEN | OPEN | high
[6,6] | OPEN | OPEN | OPEN | high
[2,2,2] | PROVEN | PMS | PMS | high
[2,2,4] | PMS | PMS | PMS | medium
[2,2,6] | PMS | OPEN | OPEN | high
[2,2,8] | OPEN | OPEN | OPEN | high
[2,2,10] | PMS | OPEN | OPEN | high
[2,2,12] | OPEN | OPEN | OPEN | high
[2,2,14] | PROVEN | PROVEN | PROVEN | high
[2,2,20] | OPEN | OPEN | OPEN | high
[2,4,4] | PMS | PMS | PMS | medium
[2,4,8] | OPEN | OPEN | OPEN | high
[2,2,2,2] | PMS | PMS | PMS | high
[2,2,2,4] | PMS | PMS | PMS | high
[2,2,2,6] | PMS | PMS | PMS | medium
[2,2,2,8] | PMS | PMS | PMS | medium
[2,2,2,10] | PMS | PMS | PMS | high
[2,2,2,12] | OPEN | OPEN | OPEN | high
[2,2,4,4] | PMS | PMS | PMS | high

## Verifier outcomes

All 20 PROVEN-level claims CONFIRMED on independent re-verification (safety-review pass
was degraded for the verifier agent; its verdicts were reviewed by the orchestrator).
Two conservative gradings overturned upward:
- **[2,14] L2/L3 -> PMS**: the committed contact7 two-root-surface note already contains
  the parametrization, generic-exactness statement and simplicity certificate the batch
  agent thought missing (second strict disjoint prime verified at q=11).
- **[7] L2**: the generic-exactness input is derivable from committed material alone
  (R1 injection into the committed [28]-member), no uncommitted session member needed.

## Framework in one paragraph

R1 (torsion specialization injectivity, char 0, with the SANDWICH corollary: containment
identity + one exact member => generic exactness); R2 (endomorphism specialization: one
certified End=Z member => generic End=Z); R3/R3a (Cadoret-Tamagawa I+II: in a 1-dim
non-isotrivial family, members of bounded degree with non-generic endomorphisms are
FINITE — the correct replacement for the useless countable-union arguments; Serre's
thin-set lemma as the cheap version); plus verified readings of Flynn 1991, Leprevost
1995, Elkies 2024, and the lab's Theorem C/C+. Full reference sheet with precise
statements, references and scope warnings: workflow wf_22fc9102-89e journal (framework
agent result).
