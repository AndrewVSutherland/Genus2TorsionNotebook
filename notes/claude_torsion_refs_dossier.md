# Literature dossier: earliest realizations of torsion subgroups on geometrically simple genus-2 Jacobians /Q

Compiled 2026-07-20 (literature-mining pass over nicholls.txt [line-verified],
gl2tors.txt, qmtors.txt, Elkies' g2_tors page, and the primary sources), as the
reference basis for paper/torsion_realizations.tex. Simplicity grades:
AS = proven absolutely/geometrically simple; QS = Q-simple only; U = unaddressed;
DB = certified via LMFDB endomorphism data.

## Source key
- Fly90 = Flynn, J. Number Theory 36 (1990) 257-265 (orders 6, 9, 10, 11, 13; U).
- Fly91 = Flynn, Invent. math. 106 (1991) 433-442.
- Lep91a = Leprevost, C.R. 313 (1991) 451-454 (order 13; U).
- Lep91b = Leprevost, C.R. 313 (1991) 771-774 (orders 15, 17, 19, 21; U).
- Lep95 = Leprevost, JTNB 7 (1995) 283-306 (AS via his D4 criterion: 22, 23, 24,
  26, 29; his C21 NOT Q-simple pp. 300-301; C25 QS+geometrically split (per
  PZP13); C27 split pp. 301-302; [3,9] family QS only pp. 303-304).
- Lep97 = Leprevost, Manuscripta 92 (1997) 47-63 (order 20 = 2g(2g+1); U).
- Oga94 = Ogawa, Proc. Japan Acad. 70 (1994) 295-298 (order 23 family; U).
- BGL01 = Boxall-Grant-Leprevost, JLMS 64 (2001) (5-torsion families; U).
- Elk02 = Elkies g2_tors page (2001-02): AS realizations 32 (1-param family),
  34 (two curves), 39, 40; order-31 SUBGROUP without rational generator.
- PP12/PZP13 = Platonov-Petrunin Dokl. 85/86 (2012), Platonov-Zhgun-Petrunin
  Dokl. 87 (2013): simple 14, 18, 28, 33 (AS); their 36, 48 are split.
- How15 = Howe, BLMS 47 (2015) 127-135: AS 27, 39; order-70 split.
- LMFDB16/BSSVY = Booker-Sijsling-Sutherland-Voight-Yasaki (2016); geometric
  simplicity certification Costa-Mascot-Sijsling-Voight, Math. Comp. 88 (2019).
- Nic18 = Nicholls thesis (2018): first proven-simple 19, 21, 25 (first ever),
  30 (first published; Leprevost had an unpublished 4/1996 order-30 preprint).
- DS18 = Daowsud-Schmidt, JNT 189 (2018) (CF order-11 families; U).
- LSSV24 = Laga-Schembri-Shnidman-Voight, Forum Math. Sigma 12 (2024): PQM
  torsion; exact (Z/3)^2 AS curve (first published simple [3,3]); (Z/2)^3 and
  Z/2x(Z/3)^2 impossible for PQM Jacobians.
- Elk24 = Elkies, LuCaNT, Contemp. Math. 796 (2024) 165-185: the 2-dim
  (Z/2)^4+Z/5 family, generic torsion [2,2,2,10] order 80, generically End = Z;
  simplest member y^2 = x(x+1)(x-1)(3x-7)(8x-13)(24x+25).
- SS24 = Sadek-Suluyer arXiv:2410.14455 (genus-2: 16, 18, 23-family; QS only).

## Attribution decisions used in the table
Cyclic: [5] BGL01+dag; [6],[9],[10],[11] Fly90+dag; [13] Fly90+Lep91a+dag;
[14],[18],[28] PP12+PZP13; [15],[17],[19],[21] Lep91b+dag ([21]: first
proven-simple is Nic18); [20] Lep97+dag; [22],[24],[26],[29] Lep95;
[23] Oga94+dag (AS first by Lep95); [25],[30] Nic18; [27] How15;
[32],[34],[39],[40] Elk02; [33] PZP13; [2]-[4],[8],[12],[16] BSSVY.
Non-cyclic: [3,3] LSSV24; [3,9] Lep95+dag (QS family; geom-simple witness only
in extended DB); [2,2,2,10] Elk24; production rows BSSVY; alpha-only rows
Booker-Sutherland (in prep); project rows ([2,2,20],[2,4,8],[6,6],[2,4,4],
[2,2,4,4],[2,2,2,12]) "new".

## Traps recorded (do not re-trip)
1. Ogawa 1994 is order TWENTY-THREE, not 13.
2. J1(13) (Z/19 torsion) and J1(16) ([2,10]) are geometrically SPLIT — never
   cite them as simple realizations.
3. Leprevost's C21/C25/C27 are non-simple/split; the 1991 families may contain
   simple members but nobody proved it before Nic18/How15.
4. Elkies' 31 is a subgroup without rational generator — not a realization.
5. [36] simple: NO published source (PZP13's 36 is split); extended-DB only.
6. Unverified leads: Stoll C.R. 321 (1995) (simple Jacobians, rank >= 19 —
   torsion group unknown, could predate Elk24 for [2,2,2,2] if full 2-torsion);
   Platonov-Petrunin Dokl. 91 (2015) "new curves" (paywalled, orders unknown);
   Platonov-Fedorov Dokl. 98 (2018) infinite order-28 family (simplicity
   unverified); Leprevost unpublished 4/1996 order-30 family.
7. gl2tors (Alessandri-Coppola) contains conjectural ORDER lists for GL2-type
   surfaces, not Jacobian realizations; its g=2 conjectural set is
   {1..24, 28, 31, 37, 44, 56}.

## End(Jac_Qbar) = Z refinement (2026-07-20, second pass)

- Production census: all 35 minimal-conductor simple representatives already have
  geom_end_alg = 'Q' (generic); non-generic simple curves in production are rare
  (116 RM, 6 CM, 3 QM) and never displace a minimal generic witness.
- Zywina-style certificates (two irreducible L-polys, linearly disjoint splitting
  fields; Zywina RNT 2022, Costa-Lombardo-Voight RNT 2021) computed for all 14
  non-production curves in the table: the six project curves, all three
  (2,2,2,12) curves, the [40]/[28] jumps, Elkies' 34/39/40, and the Elkies-2024
  [2,2,2,10] member. All PASS: results/claude_end_z_certificates.log.
- [3,3] caveat: LSSV24's realization is PQM (non-generic); flagged with a
  double-dagger in the table; generic witness = extended-DB 5100.a.1 (pending
  Booker's endomorphism filter).
- DISCOVERY: the repo's contact-5 t=-1/3 order-40 curve is Q-ISOMORPHIC to
  Elkies' 2001 order-40 curve y^2 = (3x+4)(x^4+5x^3+8x^2+(19/4)x+1) — the
  contact-5 halving family rediscovers Elkies' record curve exactly.

## Corrections from Sol's review of torsion_realizations.tex (2026-07-21, applied)

- [6],[9],[10]: Flynn citations REMOVED (Flynn's genus-2 constructions are orders 11
  and 13 only; the nicholls.txt table lists [Fly90a] for 6/9/10 in the UNVERIFIED
  column — misread by the first literature pass). Rows now BSSVY.
- [5]: now Boxall-Grant, Trans. AMS 352 (2000) (y^2-y=x^5, proven absolutely simple;
  CM, hence double-dagger non-generic flag like [3,3]). Coleman-era work underlies
  the torsion bound but no clean construction citation was pinned; BGL 2001 dropped.
- [21]: Lep91b(dag) + Nicholls2018 (Leprevost's C21 proven non-simple).
- [23]: Ogawa1994(dag) + Leprevost1995 (Ogawa constructs, Leprevost proves AS).
- [27]: Howe2015(dag) + Nicholls2018 (Howe does not discuss simplicity).
- [3,9]: Lep95(dag) + BookerSutherland.
- Column relabeled "Source(s)"; preamble reworded (dagger = torsion-construction
  only; certification via linked DB or co-cited source).
- Bib fixes: PP2012 correct 2012 title ("New orders of torsion points in Jacobians
  of curves of genus 2 over the rational number field"); PZP2013 pages 318-321;
  BSSVY vol 19(A); Elkies page "2001-2002, updated 2010".
