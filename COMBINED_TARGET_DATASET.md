# Constructing combined target datasets

Simulation targets are separate from CoRE-MOF-COD core metadata. Loading the
core release does not load targets. Use the optional target supplement only
when an analysis requests it, with an independently verified manifest hash:
see [explicit target attachment](docs/source/target_supplements.rst).

This guide documents the reproducible, fill-only target builder. It does not
declare a particular target snapshot to be the current release. Read coverage,
endpoint definitions and source hashes from the selected dataset's own receipt
and `coverage_summary.json`, not from example counts in software documentation.

## Deterministic construction

`examples/build_combined_target_dataset.py` accepts the release CIF manifest,
frozen base/additions completion manifests and summaries, and an independently
audited saved-result evidence table. It combines accepted values by exact CoRE
ID without changing groups or assignments. It fails on an unexpected input
hash, duplicate assignment, endpoint or release mismatch, non-finite success,
overwrite attempt, inconsistent exclusion or requested-count mismatch. It
writes transactionally and refuses to overwrite an existing snapshot.

Only source keys marked `MISSING` may receive a new value. An accepted scientific
null is not an empty slot: `HISTORICAL_SCIENTIFIC_NULL` means the source status
is `EXISTING`, with `explicit_null=true` and a nonempty diagnostic. Preserve it
as a native null, including a `ZERO_DENOMINATOR` result when applicable.

`examples/audit_combined_target_dataset.py` independently reloads inputs and
outputs, recomputes assignments and coverage, validates `targets.json` and the
receipt, and compares builds byte for byte with `--comparison-dataset`.

Run `--help` on both scripts for the full argument contract. A production build
should supply every `--expected-input-sha`, all three `--expected-finite`
counts, and `--expected-finite-total`. Independently build the same declared
snapshot twice and audit one against the other before adopting it.

The output directory contains:

- `targets_for_attachment.csv`: one exact-ID row per release structure.
- `target_assignments.csv.gz`: source, status and value by structure and endpoint.
- `targets.json`: typed target-attachment configuration.
- `coverage_summary.json`: finite, null, excluded and intersecting populations.
- `BUILD_RECEIPT.json` and `SHA256SUMS`: input, policy and output bindings.

The builder retains its exact source endpoint identifiers for compatibility.
Do not rename a stored response merely to change a plot label. The current
optional supplement describes CH4 and H2 uptake in mol/kg-framework and
dimensionless CO2/N2 Henry selectivity. Interpretation of a saved
Rosenbluth-weight ratio as Henry selectivity requires its documented simulation
and thermodynamic conditions. Do not infer individual Henry constants or
volumetric uptake from it.

## Select the intended analysis workflow

For a **new target-complete benchmark**, attach the declared targets explicitly,
select structures with all required finite endpoints, then construct cohorts
and partitions using full-release grouping and checker evidence. Zero is a
finite target. Missing-target structures must still participate in full-release
group construction so that filtering cannot hide connecting relationships.
Target magnitudes never determine groups, diversity or assignments. Use the
[target-first guide](docs/source/target_first_benchmark.rst).

For **existing frozen assignments**, verify their receipt before attachment.
`missing="keep"` preserves every assigned structure and null target.
`missing="error"` checks completeness. `missing="drop"` creates a derived
filtered view only, without refill, rebalancing or resplitting. Attaching targets
does not change the original assignment digest. A newer target snapshot does
not replace the target data of a completed experiment without defining a
separate analysis. See [benchmark handoff](ML_BENCHMARK_HANDOFF.md).

## Dataset identity and distribution

Both builder and auditor require `--identity-contract PATH` and
`--identity-contract-sha256 VERIFIED_SHA256`, received independently with the
trusted input hashes. The JSON specifies the exact final `release_version`
and `phase`/`version` for `base` and `additions`. For example:

```json
{
  "release_version": "CoRE-MOF-COD",
  "base": {"phase": "CoRE-MOF-COD-base", "version": "CoRE-MOF-COD-base"},
  "additions": {"phase": "CoRE-MOF-COD-additions", "version": "CoRE-MOF-COD"}
}
```

Use the exact manifest identities. The final version must match the additions
edition, and the two cohort phases and source versions must be distinct. Source
and saved-result evidence rows must match the declared identity. Receipts bind
the contract hash. Existing immutable receipts can be audited with a separately
pinned identity contract without changing their files or hashes.

Structure-resolved targets and their restricted transfer manifests are not
included in this code repository. Apply the source-specific access and
redistribution conditions described in [database access](README_DATABASE_ACCESS.md).
A successful audit is not permission to publish data.
