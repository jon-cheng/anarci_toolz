# CLAUDE.md — anarci_toolz

## Ground rules

- Check in before anything consequential: new dependencies, deleting/moving/renaming files, rewriting an existing public function's signature or behavior, force-pushing, `git clean`/`reset --hard`, or committing/pushing. Propose the change and wait rather than just doing it, unless explicitly told to proceed autonomously for a given task.
- For every non-trivial decision, state the reasoning, not just the result: why this library/dependency, why this test design or mocking strategy, why this Docker base image, why this API/MCP shape. A one- or two-line "because X, which trades off Y" beats silent action. If there were real alternatives, name them and say why they lost.
- If a spec here is ambiguous, don't silently guess — state the assumption being made and why, so it's easy to correct. Prefer asking a short clarifying question over quietly picking a default when the choice actually matters (e.g., FastAPI response shape, MCP tool argument design).
- `origin/main` on GitHub (https://github.com/jon-cheng/anarci_toolz) is the source of truth, not necessarily the local working tree. Check `git status` / diff against `origin/main` before describing "current state" — the local tree has had stale scratch files and uncommitted edits sitting around before.
- Small, unrelated fixes (a missing dependency, a typo, a dead-code duplication) should be their own small commit, not bundled into a big refactor branch/PR.

## What this project is

`anarci_toolz` wraps ANARCI and AbNumber to do batch antibody sequence numbering/annotation from tabular (CSV/DataFrame) input, using multiprocessing for throughput. Built by Jonathan during his tenure at Bristol Myers Squibb; already deployed in production there. Currently: an installable CLI + importable Python library. No tests, no Docker, no API layer, no MCP server yet — that's what this refactor adds.

**Hard environment constraint, shapes everything below:** ANARCI and AbNumber are not pip-installable. They require conda (bioconda) plus HMMER, and HMMER does not run natively on Apple Silicon (arm64) — it needs Rosetta / an osx-64 conda subdir (documented in README.md's addendum). Any code path that actually calls into ANARCI or AbNumber can only run somewhere those are genuinely installed — this is why testing and Docker are sequenced together below.

## Refactor plan, in order

**0. Testing, especially integration tests**
- Two tiers. Tier 1: pure-logic unit tests with no bio deps, runs anywhere/in plain CI — covers `utils.py`, `validation.py`, `imgt_index_ordering.py`, the formatting helpers in `numbering.py`, the DataFrame-shuffling helpers in `abnumber_tool.py`, the parsing-only pieces of `anarci_tool.py`, and `pipeline.file_or_dir_input`/`validate_aa_translation`. Tier 2: bio-dependent tests that call real ANARCI/AbNumber — mark with a pytest marker (e.g. `@pytest.mark.bio`), skip cleanly when the deps aren't present, and run for real inside the Docker image from step 1.
- Top-level integration test: call `run_anarci_toolz()` end-to-end against `test_files/therasabdab_sample.csv` (the README's own 100-row worked example, column `sequence_aa`). Spot-check known values already documented in the README's output table (e.g. Timigutuzumab -> v_gene IGHV3-66*01) rather than full-dataframe golden-file equality — `e_value`/`bitscore`/`bias` can drift slightly across HMMER versions, so exact-match there would be flaky.

**1. Dockerize**
- Needs HMMER + ANARCI + AbNumber baked into the image. README.md's addendum already has example Dockerfile snippets for installing both from source — reuse/adapt those rather than re-deriving the install steps.
- This image becomes the actual environment Tier 2 tests run in, both locally and in CI.

**2. Build a FastAPI layer**
- Wraps `run_anarci_toolz` as a service. Open questions to raise explicitly rather than silently deciding: sync vs. async given the multiprocessing-heavy ANARCI/AbNumber calls, file-upload vs. JSON-records request shape, how much of the existing CLI argument surface (scheme, allowed_species, retain_indices, etc.) maps directly to request fields, and — when the time comes — whether to support batch upload (multiple CSVs in a single request, mirroring the CLI's directory-of-files mode) or single-file-per-request only.
- This is also a good place to put throttling/rate-limiting logic if this ends up serving a public demo — caps request rate or sequence-count-per-request to keep compute costs bounded. Not needed for internal/local use; only wire it in if the API is actually exposed publicly.

**3. Build an MCP server**
- Expose the tool for agentic use on top of the FastAPI layer.

## Known issues already identified — verify current state before assuming these are unfixed

- `pipeline.py` defines its own `validate_aa_translation`, nearly identical to the one in `validation.py`; only the `pipeline.py` copy is actually called anywhere. The `validation.py` copy is dead code.
- `anarci_tool.run_parallel_anarci` hardcodes `Pool(processes=cpu_count())`, ignoring the `num_cpu` parameter entirely; `pipeline.py` doesn't even forward `num_cpu` to it in the first place (unlike `abnumber_tool.run_parallel_abnumber`, which does respect `num_cpu` via `get_num_cpu`).
- `requirements.txt` was missing `biopython` despite `pipeline.py` doing `from Bio.Seq import Seq` — a fresh `pip install` would break DNA-mode with an `ImportError`. Check whether this has been committed yet.

## Test fixture

`test_files/therasabdab_sample.csv` — 100-row TheraSabDab sample, `sequence_aa` column. This is the README's own worked example; reuse it rather than inventing new fixtures for the integration test.
