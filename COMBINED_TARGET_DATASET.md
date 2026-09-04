# Combined CoRE-MOF v26.0.2 target dataset

The canonical ML attachment dataset is a target-only, exact-ID left join over
all 42,574 published v26.0.2 structures. It combines accepted historical
values with collector-validated current calculations without changing parent
groups, leakage blocks, cohorts, or train/validation/test assignments.

“Final as of the cutoff” means the immutable union of the accepted, hash-bound
inputs available at `2026-09-04T05:43:23Z`. It does not mean that the remaining
calculation campaign is complete or that the structure-resolved dataset is
authorized for public redistribution.

## Validated coverage

| Endpoint | Finite unique IDs | Published-release coverage |
|---|---:|---:|
| CH4 absolute loading, 298 K / 65 bar | 28,979 | 68.0674% |
| H2 absolute loading, 77 K / 100 bar | 28,974 | 68.0556% |
| Raw CO2/N2 Widom weight ratio, 298 K / 1 bar | 28,944 | 67.9852% |

There are 86,897 finite structure-endpoint assignments, 28,983 structures
with at least one finite target, and 28,937 with all three. The independent
audit verifies 42,574 wide rows, 127,722 long rows, exact source bindings,
fill-only behavior, native nulls, complete coverage accounting, checksum
integrity, and byte-identical double builds. The machine-readable public-safe
counts are in `V2602_COMBINED_TARGET_COVERAGE_20260904.json`.

One accepted historical Widom source record is an explicit scientific null
with diagnostic `ZERO_DENOMINATOR`. `HISTORICAL_SCIENTIFIC_NULL` means that the
frozen source status is `EXISTING` with `explicit_null=true` and a nonempty
diagnostic. The structure remains eligible, but the endpoint is unavailable
and native-null. It cannot be replaced because current results may fill only
source keys marked `MISSING`.

## Deterministic construction

`examples/build_combined_target_dataset.py` accepts the published release CIF
manifest, frozen base/additions completion manifests and summaries, and an
independently audited current-result evidence table. It fails on an unexpected
input hash, duplicate assignment, endpoint or release mismatch, non-finite
success, overwrite attempt, inconsistent exclusion, or requested-count
mismatch. It writes transactionally and refuses to overwrite an existing
snapshot.

`examples/audit_combined_target_dataset.py` independently reloads the inputs
and outputs, recomputes every assignment and coverage field, validates the
complete `targets.json` contract and receipt, and can compare two builds byte
for byte with `--comparison-dataset`.

Run `--help` on both scripts for the full argument contract. A production run
should supply every `--expected-input-sha`, all three `--expected-finite`
counts, and `--expected-finite-total`. Generate two snapshots with the same
snapshot ID and cutoff in separate directories, then audit one against the
other before promotion.

The output directory contains:

- `targets_for_attachment.csv`: one exact-ID row per published structure;
- `target_assignments.csv.gz`: one status/provenance row per structure and
  endpoint;
- `targets.json`: typed deferred-target attachment configuration;
- `coverage_summary.json`: finite, null, exclusion, source, and intersection
  counts;
- `BUILD_RECEIPT.json` and `SHA256SUMS`: input, policy, output, and checksum
  bindings.

Use `targets.json` only after target-independent assignments are frozen.
`missing="keep"` preserves the cohort and its null targets. `missing="drop"`
creates a derived filtered view only; it must never refill, rebalance, or
resplit. Target hashes do not enter the original split receipt.

## Redistribution boundary

The structure-resolved target files and their restricted checksum manifests
remain outside Git. They may contain CSD-derived values and require the
appropriate institutional rights and an approved transfer channel. This
repository contains the reproducible code, tests, workflow documentation, and
aggregate coverage only.
