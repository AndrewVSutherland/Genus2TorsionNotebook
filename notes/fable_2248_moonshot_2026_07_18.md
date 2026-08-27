# Fable moonshot: [2,2,4,8] (order 128) Richelot-walk campaign — 2026-07-18

Target: geometrically simple genus-2 Jacobian /Q with exact J(Q)_tors =
[2,2,4,8], order 128.  This note records a complete, validated, NEGATIVE
sampling campaign over the isogeny-neighborhood flank, complementing (not
duplicating) the same-day prod-07 closures.

## Where the target stood before this campaign

1. prod-07 sign-reduction theorem + tier-1 audit: on the all-squares model
   y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2), NO tuple with 0<a<b<c<d<=65535 has
   torsion containing (2,2,4,8), through ANY component or twist; twisted
   enumeration independently 0 through d<=100000 (claude_prod_07_2248.md).
2. The genus-3 split-locus curve of S' is completely closed (third elliptic
   quotient E3: y^2=(x-1)(x+1)(x+3) has rank 0 unconditionally).
3. The simple [2,4,8] seed's own rational Richelot component: closed, 20
   nodes, no target (target_2248_attack_2026_07_18.md).

The unexplored flank: rational Richelot codomains of the (2,2,4,4)/(2,2,2,8)
banks and of the high-height certified curves.  Torsion is NOT pulled back
through an isogeny, so codomains can have strictly larger torsion — the
record's own component contains order-24 [2,12] leaves one step from the
order-96 source.  Simplicity is isogeny-invariant, so hits inherit it when
the source is certified (and get their own certificate regardless).

## Campaign design

Funnel per codomain: normalize to an integral y^2=F model -> gcd of
#J(F_p) over up to 4 good primes from {13,...,41} must be divisible by 96
-> exact TorsionSubgroup.  [2,2,4,8] forces 128 | #J(F_p) at every good p,
so the 96 gate is safe for the target and also catches any new order>=96
group.

Validation (Stage 0): walking from the record's first [2,12] Richelot leaf,
the funnel re-discovers an order-96 [2,2,2,12] vertex (VALIDATE_OK in
results/fable_2248_bank_walk.log).  The pipeline detects a jump when one
exists.

Control (Stage 1): the depth-1 fan of the known SPLIT [2,2,4,8] tuple
(ten2248models_abcd.txt row 1) has NO codomain passing even the 96 gate —
neighbors of a [2,2,4,8] curve do not retain order 96, so jumps are
directional (you must land ON the special fiber, not near it).

## Runs and outcomes (all code/fable_2248_*, results/fable_2248_*)

1. fable_2248_richelot_bfs.m — depth-2 BFS from the two certified seeds
   (2244 row 26629 = (36,57,64,132); 2228 row 3 = (1,55,99,125)):
   48 distinct vertices; every vertex torsion <= 64; the only rationally
   halvable 4-classes are those inside existing 8-chains; the [2,2,4,4]
   seed's 48 order-4 classes have NO rational halves (consistent with the
   sign-reduction theorem).  Codomain inventory: [2,8], [2,2,4], [4,4],
   [2,2,2,2], down to [2,2] — Richelot steps TRADE 2-structure, they do
   not accumulate it.
2. fable_2248_bank_walk.m — depth-1 fans over 170 stride-sampled tor2244
   rows (spanning rows 1..26654, mostly >1300 = never explored) plus 170
   stride-sampled tor2228 rows (never explored): ~5,000 codomains, 277
   gate-passers exact-tested, torsion inventory of gate-passers:
   101 x [4,4], 95 x [2,2,4], 60 x [2,8], 19 x [2,2,2,2], 2 x [2,2,2,4].
   ZERO jackpot, ZERO order >= 96.
3. fable_2248_highseeds.m — depth-1 fans of the six certified-simple
   off-rectangle [2,2,2,8] curves at heights 1e13-1e16 (already outside
   the audit box) and the twisted (29,121,125,145) curve; model convention
   verified per source (64 | #J(F_p)).  105 codomains, ZERO gate passes.

Also staged: fable_2248_validate.m — fresh-session validator (exact
torsion + two-prime root-power certificate) for any future hit.

## Conclusion

[2,2,4,8] is now boxed from three independent sides: the low-height model
box is closed by theorem; the split locus is closed unconditionally; and
the rational isogeny neighborhoods of every accessible simple order-64
curve (both banks, both certified seeds, all six high-height curves, the
[2,4,8] seed) contain no order >= 96 vertex at depth <= 2 / wide depth-1.
Order-96+ torsion does not arise by drifting near it: the record itself
fell to a structural pencil point, not to sampling.  [2,2,4,8] will fall
the same way or not at all.

## Next-step ladder (ranked)

1. STRUCTURAL: the HLP two-parameter section-surface program (two contact
   blocks + one halving block, plane through two split [2,2,4,8] anchors,
   transverse to Humbert branches) — the exact analogue of the T5-pencil
   that produced the record.  The ten split anchors are in
   paper/scripts_and_data/ten2248models_abcd.txt.
2. ENUMERATIVE: extend identity2248_prod past d=65535 on the tor2244-side
   4-square conditions (prod-07's stated route; heavy compute, cluster
   job, not this machine).
3. INCREMENTAL: depth-2 frontier expansion on the 277 gate-passers
   (cheap; expected yield low given (1) directionality).
4. STANDING RULE: depth-1 fan every NEW simple curve of torsion order
   >= 48 the moment it is found (cost ~seconds; the record's component
   shows a single lucky fiber can hide a 4x jump).
