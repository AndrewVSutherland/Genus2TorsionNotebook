# Compact rewrites of six non-database rows in `torsion_realizations.tex`

Date: 2026-07-23.  Follow-up to Ari's `[2,2,2,12]` model swap
(`notes/ari_22212_split_model_2026_07_23.md`): a systematic audit of the
remaining non-database rows for nicer presentations.  Unlike the
`[2,2,2,12]` change, **every rewrite here is the same curve** — either the
literal same polynomial re-typeset, or the isomorphic model
`Y^2 = h^2 + 4 f_0` (`Y = 2y + h`) of a row displayed as
`y^2 + h y = f_0`.  All six identities verified exactly in Magma
(scratchpad `rewrite_check.m`, all `true`); the certificate generator
`code/claude_endz_certificates.m` was rewritten in the same factored forms
and rerun — `data/claude_endz_certificates.txt` came out **byte-identical**
(log `results/claude_endz_certificates_regen_20260723b.log`), confirming
both the identities and that no certificate is affected.

## Pure re-typesetting (identical polynomial)

- `[2,2,2,8]`: `x(x+1)(x+3025)(x+9801)(x+15625)`
  `= x(x+1)(x+55^2)(x+99^2)(x+125^2)`  (and `125^2 = 5^6`).
- `[2,2,4,4]`: `x(x+1296)(x+3249)(x+4096)(x+17424)`
  `= x(x+36^2)(x+57^2)(x+64^2)(x+132^2)`.
  Bonus structure: many root differences are `7·square`:
  `64^2-57^2 = 7·11^2`, `64^2-36^2 = 7·20^2`, `132^2-57^2 = 7·45^2`,
  `132^2-36^2 = 7·48^2`.

## Completed square + factorization (isomorphic model, same curve)

- `[2,2,14]` (was `y^2+(x^2+x)y = 9x^6-1005x^5+...-37741275`):
  `y^2 = (x+1)(x+9)(2x-115)(6x+55)(3x^2-220x+2652)`.
  Factor type `[1,1,1,1,2]` displays the 2-rank 3 (2-part `[2,2,2]`).
- `[2,2,20]` (was `y^2+(x^2+x)y = -391671x^6+...+254016`):
  `y^2 = -(x-1)(6x+1)(2x+7)(6217x+1008)(21x^2-161x+144)`.
  Same `[1,1,1,1,2]` shape, 2-rank 3 (2-part `[2,2,4]`).
- `[2,4,8]` (was `y^2+(x^2+x)y = 60x^5+...-239760x`):
  `y^2 = x(3x-5)(5x+333)(16x^2+23x+576)`.
  Quintic type `[1,1,1,2]` (+ infinity), 2-rank 3.

## Same-model factorization

- `[6,6]`: `1872x^5-3000x^4+6969x^3-1691x^2+4875x`
  `= x(39x^2-69x+125)(48x^2+8x+39)`.
  Quintic type `[1,2,2]` displays the 2-rank 2 (2-part `[2,2]`); note the
  constant term of the second quadratic equals the leading coefficient of
  the first (39).

## Left alone

- `[2,4,4]` (user decision: keep `x(x+16)(x+(644/799)^2)(x^2+x+16)`).
- `[2,2,2,10]` (already Elkies's compact form).
- All LMFDB-linked rows: their equations are the stable identifiers per
  the paper's own conventions paragraph, so they are not re-typeset.
