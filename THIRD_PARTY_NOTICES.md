# Third-party material and external dependencies

The MIT license in `LICENSE` covers the original content of this
repository (the authors' and their AI assistants' notes, code, data, and
reports). The following items have their own provenance:

## Third-party code included in this repository

* **`product/code/gluing.m`** — the `BHLS3` block (the (3,3)-gluing
  algorithm) is copied, as marked inline, from Everett Howe's
  `Genus2.magma` (<https://ewhowe.com/Magma/Genus2.magma>), implementing
  Bröker–Howe–Lauter–Stevenhagen (MR3349314); the (2,2)-gluing follows
  Howe–Leprévost–Poonen (MR2514865). The surrounding file is extracted
  from A. Sutherland's `genus2.m` package.

## External software and data this repository depends on or was built with

None of the following is included here; the code refers to them and the
notes/logs record their use.

* **Magma** (computational algebra system; most `.m` files), **PARI/GP**,
  **SageMath**, **Singular** — the computer algebra systems used
  throughout.
* **CHIMP** (E. Costa et al., <https://github.com/edgarcosta/CHIMP>) and
  **ModularAbelianSurfaces** (Costa–Elkies–Hashimoto–Jha–Martin–Poonen–
  Voight, <https://github.com/edgarcosta/ModularAbelianSurfaces>) — used
  by the genus-2 curve reconstruction lane (`code/claude_g2rec_run.m`);
  the `[31]` witness curve derives from the latter's public dataset.
* **N-congruences** (S. Frengley,
  <https://github.com/SamFrengley/N-congruences>, Zenodo DOI
  10.5281/zenodo.21704448) — an external dependency of
  `product/code/lane_z3.m`, which contains clone instructions.
* **The LMFDB** (<https://www.lmfdb.org>) and the extended
  Booker–Sutherland genus-2 database served on the LMFDB alpha site —
  the source of the database-derived curve data and labels appearing in
  notes, data files, and the paper's tables (dated snapshot identifiers;
  see the paper for citations).
* **Elkies–Kumar ancillary data** (arXiv:1209.3527) — two files from the
  arXiv source tarball were used and are cited, not included; see
  `data/claude_z31_ek8_source/SOURCE.md`.
* Published papers and web tables (Elkies, Howe–Leprévost–Poonen,
  Leprévost, Nicholls, Platonov–Petrunin, Alessandrì–Coppola, and
  others) are cited where used; no copyrighted article text is
  reproduced in this repository.

## A note on running the code

This is frozen research code, not maintained software. Scripts assume
their inputs are the trusted data files of this repository (several use
`eval`-style parsing on those inputs by design); do not point them at
untrusted data, and prefer running them in an isolated environment. See
also the warning in `AGENTS.md`.
