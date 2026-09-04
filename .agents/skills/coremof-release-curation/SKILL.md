---
name: coremof-release-curation
description: Run and audit the transferred CoRE-MOF v26.0.2 strict CR/NCR ML benchmark on a separate GPU machine. Use for verifying the restricted benchmark handoff, building target-independent whole-block assignments, constructing or attaching the combined as-of-cutoff adsorption targets after assignment freeze, or reporting benchmark provenance and coverage from a CoRE-MOF-Tools clone.
---

# CoRE-MOF v26.0.2 ML benchmark handoff

This repository-scoped skill is deliberately limited to the portable ML
benchmark. It contains no original-host paths, scheduler state, licensed-node
configuration, or release-production instructions.

## Start from the two exact bindings

1. Treat the Git clone as the source root. Read `ML_BENCHMARK_HANDOFF.md`
   completely before running the benchmark.
2. Obtain the separately transferred restricted handoff from the project
   coordinator. Do not expect release tables or targets in Git.
3. Verify the archive SHA-256 before extraction, then verify the extracted
   top-level `SHA256SUMS` and every nested ledger named by
   `TRANSFER_MANIFEST.json`.
4. Compare the manifest's exact repository URL, branch, commit, package
   version, and source-file hashes with the checkout. Stop on any mismatch;
   do not select a similarly named commit or data directory by date.
5. Use Python 3.9--3.11 and install `.[benchmark]`. Run `coremof doctor` and
   require the exact five-package backend recorded by the guide and manifest.
6. Keep the extracted handoff and all derived assignments outside Git.

The transferred snapshot is intentionally current-finished rather than a
completed target campaign. Preserve its coverage/error/null accounting and all
false completion, promotion, publication, and official flags. Never fill a
missing target from another structure, campaign, or unit convention.

## Freeze target-blind assignments

Strict computation-ready (CR) means MOFClassifier, original MOFChecker,
Chen--Manz, MOSAEC, and SETC-GAT are all available and PASS. Strict
non-computation-ready (NCR) means those same five results are all available and
FAIL. Any missing, timed-out, failed, unsupported, or otherwise
`NOT_AVAILABLE` result is a non-vote that makes the row UNCHECKED; it never
becomes FAIL or NCR. A complete mixture of PASS and FAIL is AMBIGUOUS.

For `group_criteria="priority_main"`, `priority_main` is the conflict-aware
explanatory hierarchy: exact available 264-value depth-5 RAC5 groups seed
components, then exact complete MOFid-v2 and MOFid-v1 groups attach unresolved
rows. A weaker group touching multiple stronger components records
`PARENT_METHOD_CONFLICT` and does not merge them; missing evidence remains a
singleton. It excludes targets, Zeo++, CrystalNets, source identifiers, CIF
hashes, and StructureMatcher.

The separate `main_union` leakage guard is the transitive connected-component
union constructed over the complete release before any label, source, variant,
metal, ID, or target filtering. Its direct edges come from exact full CIF SHA-256,
database-namespaced source siblings, and available release-authorized exact
RAC5, MOFid-v2, and MOFid-v1 edges. It is a partition guard, not an explanatory
parent or identity claim. The benchmark adds the requested criterion edges and
takes connected-component closure; every resulting effective leakage block is
indivisible across train, validation, and test.

The checksum-bound v26.0.2 integration has 6,294 raw strict CR and 2,299 raw
strict NCR rows. Some share an effective block with another checker label, so
the default must fail closed. Explicitly request
`complete_release_label_pure_effective_blocks`: it excludes 1,601 CR and 572
NCR rows and leaves eligible pools C=4,693 and M=1,727. This is a sensitivity
cohort, not relabeling and not the complete strict population. At NCR-pool
fraction q, select `round_half_up(q*M)` NCR and
`C-round_half_up(q*M)` CR rows; q=1 is 1,727 NCR plus 2,966 CR, not a 100%-NCR
cohort. Recompute and verify all counts from the bound release.

Run the exact command from `ML_BENCHMARK_HANDOFF.md`. The representative
diversity profile may read only complete finite RAC5, otherwise the declared
complete Zeo++ vector, otherwise the explicit no-numeric tier. It must never
read adsorption values or target availability. It requires exact NumPy
1.26.4, scikit-learn 1.5.0, SciPy 1.13.1, joblib 1.5.3, and threadpoolctl
3.6.0; numerical libraries run with a one-thread limit and the receipt records
their non-path runtime identity. Cross-architecture bit identity is not
guaranteed, so compare receipts and retain the frozen assignment digest.

Require all 30 requested runs, one identical whole-block pure-CR test across
ratios and seeds, zero crossed effective blocks, q=1 containing every eligible
NCR row, zero partially selected blocks, and `official_split=false`. The
supplementary `full_cr_diagnostic` covers the complete raw strict-CR pool,
including policy-excluded rows; it is not the independent paper test.

## Attach targets only after freeze

Verify the benchmark `SHA256SUMS`, suite receipt, release binding, and frozen
assignment digest before opening the target table. Prefer the independently
audited combined as-of-cutoff snapshot, not the completion-only current-results
view. Build it with `examples/build_combined_target_dataset.py`, require exact
release/source hashes and expected counts, rebuild independently, and audit it
with `examples/audit_combined_target_dataset.py --comparison-dataset` before
promotion. Read `COMBINED_TARGET_DATASET.md` for the current count and null
contract.

At cutoff `2026-09-04T05:43:23Z`, the combined exact-ID left join spans all
42,574 published structures and has 28,979 finite CH4, 28,974 finite H2, and
28,944 finite raw CO2/N2 Widom-ratio labels. The earlier 2,335 / 3,744 / 14,167
counts are new current-finished evidence only, not total target availability.
`HISTORICAL_SCIENTIFIC_NULL` means the frozen source status is `EXISTING` with
`explicit_null=true` and a nonempty diagnostic; it remains eligible but
unavailable/native-null and cannot be filled because only `MISSING` keys accept
current results. The current combined snapshot has one historical Widom
`ZERO_DENOMINATOR`, so the finite count is one below the raw-existing-plus-new
arithmetic.

Use the target config from the restricted snapshot with `missing="keep"` unless
the analysis contract explicitly requires another policy. `keep` preserves
every assigned ID and a native null for unavailable targets. `drop` makes only
a derived filtered view; it must not refill, rebalance, or resplit.

The CLI attaches one `runs/<run_key>.csv` at a time using the suite
`receipt.json`; `membership_manifest.csv` repeats IDs across runs and is not
an attachment manifest. For all runs together, use
`suite.attach_targets(...)` in Python before serialization or loop over the
per-run CSVs. Target hashes never enter or alter the original assignment
receipt, and every derived output keeps `official_split=false`.

## Publication and reporting boundary

- Keep release tables, checker evidence, structure-resolved targets, and
  derived manifests in the approved restricted transfer area. Do not commit
  them, add them to Git LFS, or upload them to an unapproved public cloud.
- The clone contains legacy tracked SI archives and is not a sanitized
  data-free repository. Do not add any new restricted asset.
- Confirm recipient rights, institutional CSD/CCDC entitlement, and
  asset-specific redistribution terms before moving CSD, SI, MOSAEC, or other
  licensed-derived content.
- Report raw, excluded, eligible, missing/null/error, and endpoint coverage
  separately. Do not describe the current target snapshot, sensitivity cohort,
  exploratory assignment, or full-CR diagnostic as complete, official, or a
  publication benchmark.
