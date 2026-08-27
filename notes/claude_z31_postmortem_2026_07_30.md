# Post-mortem: why the lab never targeted Z/31 (2026-07-30)

Question (Ari): the notes are heavily 2-/3-torsion focused — was that anchoring us
away from [31]? Answer, from a full audit of all four top-10 editions, the dossier,
the playbook, and the paper: **yes, structurally**. Details with file:line evidence
in the session's research-lane report; summary here.

## The five mechanisms

1. **No prime target was ever ranked.** Across the four top-10 editions (07-17,
   07-19, 07-25, 07-26): ~40 ranked slots, zero prime-order targets; the largest
   prime factor of any ranked group is 11 ([2,22]). The 07-19 edition is entirely
   {2,3}-smooth.

2. **The ranking criteria structurally excluded primes.** Criterion A (top10_ranking
   13-18) scored targets by *split-Jacobian existence evidence* — but a split
   Jacobian /Q cannot carry rational prime torsion >= 11 (Mazur caps elliptic
   torsion; odd torsion transfers across the gluing isogeny), so every prime >= 11
   sat permanently in the "nowhere" tier. The HLP 2000 split table ("the strongest
   existence ledger in the project") contains only elliptic-gluing composites,
   max prime factor 7. Criterion C ("machinery in this repo") requires a route the
   toolkit can express — and the construction alphabet (contact-5/6/7/9, coprime
   composition, halving, +3, 2-rank engineering) cannot express any prime beyond 7
   (11 only via the CF/D_inf ansatz). The target-playbook's decomposition step
   returns the empty set for a bare prime: no factorization, no contact chart above
   9, nothing to halve.

3. **Z/31 was named once and lost.** claude_review_gpt56_plan.md:96-99 (2026-07-17)
   corrected "(3,12) is the smallest open group" to note (5,5), Z/31, Z/35 are
   smaller; Z/31 landed on a near-list as "no known construction anywhere"
   (top10_ranking:80) and was silently absent from every later gap table — the
   07-25 "conspicuous gaps" table jumps 25 -> 35 despite its own frontier list
   making 31 the first cyclic gap.

4. **The dossier knew and filed it under traps.** Both dossier records touching 31
   are framed as things to avoid: trap 4 ("Elkies' 31 is a subgroup without rational
   generator — not a realization") and trap 7 (the gl2tors conjectural list
   "not Jacobian realizations"). The paper mentions 31 once — to justify that it
   has no row.

5. **Modularity/RM was contamination, never a source.** The RM screen existed purely
   to *discard* RM curves from claims; no note ever proposed scanning modular
   abelian surfaces as a torsion source, despite the dossier tabulating that prime
   orders 5/7/11 were historically first realized by cuspidal classes on modular
   Jacobians. Compounding this: every lab construction forces *visible* (rational or
   Weierstrass-supported) classes; the 1830 generator is invisible
   (conjugate-quadratic support), unreachable by design. The one tool with unbounded
   order reach (CF/D_inf) is also visible-by-construction, and the literature's CF
   school stalled at <= 33 for the same reason.

## Corrective actions taken this session

- The Eisenstein sweep (`notes/claude_z31_censuses_and_eisenstein_2026_07_30.md`) is
  the missing "modularity as a source" scan — run over all 80,387 dim-2 forms; it
  retro-finds 1830.2.a.q deterministically and outputs the next candidates
  (2190.2.a.v = Z/37).
- Generic-pool censuses with no construction bias (isogeny-completion 618k;
  Bruin–Stoll 196k) — both negative, meaning the miss was NOT hiding in the
  generic pools we could reach; it was in the GL2-type world the toolkit excluded.
- Target policy for future lists: track open groups by ORDER gaps (31, 35, 37, ...)
  regardless of whether the toolkit has a route; a group with no route is a reason
  to look for a new mechanism (Eisenstein, invisible classes), not to unrank it.
