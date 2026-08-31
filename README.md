# Genus2TorsionNotebook

This is the shared **lab notebook** for the project behind the paper [Rational Torsion on Simple Genus Two Jacobians](http://arxiv.org/abs/2608.28543) by Jennifer S. Balakrishnan, Filip Najman, Ari Shnidman, and Andrew V. Sutherland.


It is not a polished software package. It is the working record — largely
unplanned, only loosely coordinated, and deliberately preserved warts and
all — of several research mathematicians and a diverse collection of AI
coding agents (Anthropic's Claude, OpenAI's GPT models and codex, and
others) attacking one problem in a shared mono-repo: **which finite abelian
groups arise as the rational torsion subgroup of a genus-2 Jacobian over
Q**, on both the geometrically simple and the geometrically split side.

Half-baked ideas, abandoned routes, duplicated effort, failed searches, and
the occasional wild goose chase are all here on purpose: some of the
project's breakthroughs came when one agent picked up another's dead end
and viewed it from a different angle, and we would like both future humans
and future agents to be able to do the same.

## Who wrote this

The repository contains human-generated content contributed by **Jennifer
Balakrishnan, Filip Najman, Ari Shnidman, and Andrew Sutherland** (the
paper's authors), but it consists very largely of **AI-generated content**:
the design notes, the search and verification code, the run logs, and most
of the session write-ups were produced by AI models operating as
computational assistants under the authors' direction. To the best of our
recollection the models involved were OpenAI's **GPT‑5.2, GPT‑5.4, GPT‑5.5,
and GPT‑5.6** (the last both as the interactive assistant nicknamed "Sol"
and as the "codex" pull-request reviewer), Anthropic's **Claude Opus 4.8,
Claude Opus 5, and Claude Fable 5** (mostly via Claude Code), and Google's
**Gemini 3.1 Pro**. Nothing in the paper rests on an assertion by any of these systems:
every result claimed in the paper carries an independently checkable
certificate (see `paper/scripts_and_data/`).

## What's in here

A snapshot of roughly **300 MB** across **~6,850 files**. Approximate
counts by kind (files / lines):

| Kind | Files | Lines |
|---|---:|---:|
| Magma (`.m`) | 1,082 | 180,000 |
| Python (`.py`) | 148 | 29,000 |
| C / C++ (`.c`, `.cpp`) | 57 | 17,000 |
| SageMath (`.sage`) | 29 | 5,600 |
| PARI/GP (`.gp`) | 74 | 5,000 |
| Shell (`.sh`) | 17 | 540 |
| Markdown notes (`.md`) | 347 | 84,000 |
| TeX (`.tex`) | 12 | 2,400 |
| HTML reports (`.html`) | 9 | 5,500 |
| Data (`.txt`, `.tsv`, `.csv`) | 1,620 | 5,500,000 |
| Run logs (`.log`) | 3,390 | 820,000 |

The bulk of the byte count is captured computation — data tables and run
logs — not source. The human-readable narrative lives in the `.md` notes
and the `.html` reports.

## An example conversation

The AI agents and the authors held many conversations over the course of
this project. Most are not preserved here, but we have chosen to include
**one** as a representative example: the full transcript of the overnight
Claude Code session in which the order‑96 group (Z/2)³×Z/12 — a record for
a geometrically simple genus‑2 Jacobian over Q — was found. It is a fair
illustration of how the human/AI collaboration actually proceeded. Read it
as
[`reports/order96-discovery-session/transcript.md`](reports/order96-discovery-session/transcript.md)
(renders directly on GitHub, with the agent's tool calls in expandable
blocks); a styled standalone version is in `index.html` in the same folder.

Start with:

* **[HUMANS.md](HUMANS.md)** — orientation for human mathematicians: what
  the problem is, what was found, how the repo is organized, where the
  interesting unfinished threads are.
* **[AGENTS.md](AGENTS.md)** — orientation for AI agents: how to navigate,
  the naming and workflow conventions the notebook uses, the distilled
  methodology in `skills/`, and the open problems.

The paper's formal verification scripts (self-contained Magma scripts that
certify every computational claim in the paper) live in their own minimal
repository, [AndrewVSutherland/Genus2Torsion](https://github.com/AndrewVSutherland/Genus2Torsion);
a copy is included here under `paper/scripts_and_data/`.

This repository is released under the MIT License (see `LICENSE`);
third-party material and external dependencies are itemized in
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

Everything in this notebook that matters mathematically was re-verified
before being claimed in the paper. Everything else is a snapshot of work in
progress and should be treated accordingly: notes record what we believed
at the date they carry, including beliefs that later turned out to be
wrong — corrections and errata are themselves part of the record, and
[`STATUS_AND_ERRATA.md`](STATUS_AND_ERRATA.md) maps the significant
corrected claims to their corrections.

A final practical note: this is frozen research code, not maintained
software. The scripts assume their inputs are the trusted data files of
this repository (several parse input with `eval`-style constructions by
design); do not point them at untrusted data, and prefer running them in
an isolated environment.
