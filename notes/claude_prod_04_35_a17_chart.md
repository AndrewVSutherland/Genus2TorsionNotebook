# Z/35 Route B: the A_1(7) threefold — reference identified, equations NOT yet accessible (2026-07-31)

## CORRECTION + RESOLUTION (2026-08-01)

Two updates from the follow-up retrieval session:

1. The published Elkies paper (ConM 796, 2024,
   doi:10.1090/conm/796/16001) had already been consulted by 2026-07-01
   (this note's "not web-accessible" claim was wrong).
2. **The [NDE 2003] mystery is resolved and it is NOT a retrievable
   document**: the published paper has no [NDE 2003] bibliography entry;
   footnote 13 (extract line 1299) says "Theorems 2.3 and 2.5 date back to
   2003" — i.e. the slides' bracket cites Elkies' own unpublished 2003
   notes, whose 5-torsion content is now published as Theorems 2.3/2.5.
   The published text contains NO 7-torsion formulas (the slides' remark
   about "similar formulas for 6-, 7-, or 8-torsion" refers to the same
   unpublished notes).  So A_1(7) has exactly two sources: email Elkies for
   the 2003 notes, or re-derive with the paper's now-fully-published method
   (Theorem 2.3's proof is the template: express 2D+i*(D'), D+2D' as
   divisors cut by cubics, coefficient-match; for 7-torsion the analogous
   divisor identities involve K+aT multiples with a length-7 relation).

Lane 2A of the 2026-07-31 campaign.  Route B for Z/35
(`notes/claude_prod_04_35.md`, `notes/claude_top10_10_35.md`): sieve the full
rational moduli threefold `A_1(7)` of pairs (C, T7) for `5 | #J`, mirroring
the executed `A_1(5)` sweep (which is clean to height 88).

## What "NDE 2003" is (resolved)

Elkies' LuCaNT 2023 slides ("Families of genus-2 curves with 5-torsion",
ICERM 2023-07-14, published Contemp. Math. 796 (2024); slides at
app.icerm.brown.edu/assets/381/5442/5442_4022_Elkies_071420231100_Slides.pdf,
local text extract in the session tmp) state, after deriving the universal
5-torsion model `y^2 + (L'Q' - LQ)y = Q^2 Q'`:

> "[This is already [NDE 2003], with similar formulas for genus-2 Jacobians
>  with a 6-, 7-, or 8-torsion point.]"

So the A_1(7) equations EXIST in an Elkies 2003 work.  Attempts to retrieve
it this session:

- `people.math.harvard.edu/~elkies/g2_tors.html` — fetched; contains the
  high-order examples (incl. the N=31 curve) but NOT the 7-torsion family.
- AMS full text of Contemp. Math. 796 art. 16001 (whose bibliography would
  pin the exact citation) — HTTP 403, not open-access from this machine.
- NMBRTHRY archive search for a 2003 Elkies posting — not found via web
  search (archives poorly indexed).

## Decision (per the session plan's rule)

A0 could not be completed in one working session => Route B is DEMOTED
behind the [8,8] primary lane.  No compute was spent.

## Unblocking actions (cheap, for Filip)

1. Pull the published Elkies paper (Contemp. Math. 796, article 16001,
   "Families of genus-2 curves with 5-torsion") through the university
   library / MathSciNet and read off the [NDE 2003] bibliography entry —
   with the equations either in that reference or reconstructible from the
   paper's method in an afternoon.
2. Or email Elkies directly for the 2003 notes (the 6/7/8-torsion formulas).
3. Fallback (research session, not a lookup): re-derive A_1(7) with the
   LuCaNT method — model a genus-2 curve with marked divisor T, impose
   7T ~ 0 via two function identities generalizing the 5-torsion
   `y^2 + (L'Q' - LQ)y = Q^2 Q'` trick.  The 5-torsion derivation is fully
   spelled out in the slides (pages 9-14) and is the template.

Once the chart is in hand, the prepared next steps are A1 (order-7
validation gate, 10/10 members), A2 (F_p^3 5-divisibility scans at
p=3,11,13,17), A3 (C sweep cloned from `claude_prod35_sweep.c`) — see the
session plan.

## Bonus extract from the same slides (relevant elsewhere)

- The "atypical" (C,T) locus recovers the Boxall–Grant–Leprévost family
  `y^2 + (x^3+a1x^2+a2x)y = a5x`, and full 2-level + 5-torsion has moduli
  open in the Clebsch–Klein cubic surface (sum r_i = sum r_i^3 = 0);
  Elkies' order-80 record [2,2,2,10] curve
  `y^2 = x(x+1)(x-1)(3x-7)(8x-13)(24x+25)` comes from the C-K point
  (1:5:-7:-8:9).  Potentially useful chart for future [2,2,2,2n]+5 work.
- RM5 + atypical intersection: the 1-parameter hypergeometric family
  `y^2 = (x^3+5x^2+5x)^2 + 4x`.
