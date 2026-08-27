# Split-gap campaign, 2026-08-26: derivation-chain scripts

Session provenance for the campaign that produced the first split-Jacobian
realization of **[11]** (see `notes/claude_split_gaps_campaign_2026_08_26.md`
for the narrative and `product/code/` for the productionized lanes).  These
are verbatim working scripts, kept in discovery order — including the dead
ends, which carry the structural lessons.  Some attach machine-local paths
(`/home/claude/torsion_jac/genus2.m`, Drew's Magma spec) or Sam Frengley's
N-congruences package (`github.com/SamFrengley/N-congruences`); adjust to
your checkout when replaying.  The clean end-to-end reproduction of the [11]
hit is `product/code/lane_qglue.m` (N:=11, ~1 s); the clean verification is
`product/code/verify_split11.m`.

Chain, in order:

| script | what it established |
|---|---|
| `probe_x19_disc.m` | square-discriminant (cyclic-cubic self-glue) loci of X1(7)/X1(9): genus 2/3 curves with only degenerate small points — [9,9] self-glue route parked |
| `probe_frengley.m` | Frengley's ZNrModuli works: Z(3,2) points -> genuine anti-3-congruent pairs; blind sweeps are torsion-thin |
| `probe_c_fiber.m` | the Z(3,2) fiber over fixed j1 is an irreducible degree-12 plane curve of geometric genus 0 (the X_E(3;2)-twist) |
| `probe_contrav.m` | **the contravariant-pencil discovery**: Jacobian(P + lam*Q) for (P,Q) = Contravariants(GenusOneModel(3,E)) is the anti-3-congruent family with correct twists (glues for every member); the Hesse pencil is the symplectic twin (never glues) |
| `probe_x111.m` | X1(11) raw model from psi_11(0): r^2 - r(s^3-3s^2+4s) + s = 0 (deg 2 in r => rational-s fibers are quadratic points) |
| `probe_x111b.m` | **the scan that found the hit**: among 1109 quadratic points of X1(11), exactly one (s0 = 4/5, K = Q(sqrt 11)) has E[2] ~ E^sigma[2]; no sigma-3-matches, no Q-curves |
| `probe_s45.m` | deep dive at s0 = 4/5: E(K)tors = [11], no isogeny E -> E^sigma of degree <= 13, 2-division field S3 over K (=> the equivariant iso is unique and sigma-descent is automatic) |
| `probe_s45glue.m` | the HLP/BHLS (2,2)-glue reimplemented over K (genus2.m's version is Q-specific); the glued sextic lands in K[x]; Igusa-Clebsch invariants are rational; Mestre reconstruction over Q "succeeds" |
| `s45diag.m` | **the Mestre trap**: the reconstructed curve is NOT any quadratic twist of the true descent (bielliptic => Aut V4 => Q-forms exceed quadratic twists); per-prime P(±1) mod 11 diagnostic |
| `s45check.m` | the K-glue itself is right: #Jac(CK)(k_p) = #E(k_p1) #E(k_p2) at six split primes |
| `s45descend.m` | honest Weil descent: sigma acts on the even glue sextic by x -> 1/x; matrix Hilbert-90 (N = A + M A^sigma) + scalar rescale gives an explicit rational model |
| `s45final.m` | twist sieve (11-divisibility of #J(F_p)): d = 2 gives **exact torsion [11]**; d = -5 is the twin descent with torsion [] |
| `s45validate.m` | independent validation from integer coefficients alone: minimal model, exact [11], Weil-restriction count identities at split and inert primes, Q-simplicity |
| `probe_x113.m` | X1(13) raw model (genus 2, hyperelliptic) |
| `dbg13b.m` | X1(13) fibral quadratic points have E^sigma = quadratic twist of E (hyperelliptic involution = diamond); x0 = 2 recovers the classical Q(sqrt 17) point; the unique glue datum is the degree-1 twist iso => Kani-degenerate => 2-glue dead on the family |
| `glue_lmfdb11.m` | none of the LMFDB's five quadratic-11-torsion conjugate pairs (conductor norms 46..4338) is sigma-2-congruent |
| `lmfdb11_c3.m` | ...nor sigma-3-congruent (22-32 trace mismatches each) — the Q(sqrt 11) configuration is unique among all known curves |
