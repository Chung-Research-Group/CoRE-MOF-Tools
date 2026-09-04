# CoRE-MOF v26.0.2 GPU benchmark handoff

This guide is the minimal transfer-and-run contract for the target-independent
strict five-checker benchmark. Use the exact source commit and restricted data
digests supplied by the project coordinator; do not select similarly named
files by timestamp.

## 1. Freeze the benchmark inputs

Strict computation-ready (CR) means that MOFClassifier, original MOFChecker,
Chen–Manz, MOSAEC, and SETC-GAT are all available and PASS. Strict
non-computation-ready (NCR) means that the same five results are all available
and FAIL. A missing, failed, timed-out, or otherwise `NOT_AVAILABLE` checker
result makes the row UNCHECKED and never becomes NCR.

For `group_criteria="priority_main"`, `priority_main` is the conflict-aware
explanatory hierarchy in which exact available 264-value depth-5 revised
autocorrelation (RAC5) groups seed components, then exact MOFid-v2 and MOFid-v1
groups attach unresolved rows. A weaker group spanning multiple stronger
components records `PARENT_METHOD_CONFLICT` and does not merge them; missing
rows remain singletons. It excludes Zeo++, CrystalNets, source ID, common name,
CIF hash, StructureMatcher, and optional reference relations.

The separate `main_union` leakage guard is built over the complete release
before filtering. It takes transitive components over exact full CIF SHA-256,
database-namespaced source siblings, and available release-authorized RAC5,
MOFid-v2, and MOFid-v1 edges; it is not a parent or identity claim. The
benchmark adds the selected criterion edges and takes connected-component
closure. Each resulting effective leakage block is indivisible across train,
validation, and test.

The checksum-bound v26.0.2 integration run has:

| Pool | Raw strict rows | Label-pure eligible rows | Excluded with a mixed-label block |
|---|---:|---:|---:|
| CR | 6,294 | 4,693 | 1,601 |
| NCR | 2,299 | 1,727 | 572 |

The raw pools cannot be selected in full without crossing effective leakage
blocks. Therefore the v26.0.2 sensitivity ladder must explicitly request:

```text
complete_release_label_pure_effective_blocks
```

That policy excludes every strict CR or NCR row whose complete-release
effective block contains another checker label. It is a declared sensitivity
cohort, not a relabeling of excluded rows. Omitting the option fails closed.
Every real run recomputes these counts from its checksum-bound release.

After this policy, `full_cr` uses eligible `C=4,693` and `M=1,727`.
For NCR-pool fraction `q`, it selects
`round_half_up(q*M)` NCR and `C-round_half_up(q*M)` CR structures. Thus
`q=1` contains all 1,727 eligible NCR structures and 2,966 eligible CR
structures; it is not a 100%-NCR cohort.

## 2. Install the exact source and numerical environment

Clone or update this repository, check out the exact commit recorded in the
handoff manifest, and install the benchmark extra:

```bash
git fetch origin
git checkout COMMIT_SHA_FROM_HANDOFF_MANIFEST
python -m pip install -e ".[benchmark]"
coremof doctor
```

Use Python 3.9, 3.10, or 3.11. The representative diversity backend requires
exactly NumPy 1.26.4, scikit-learn 1.5.0, SciPy 1.13.1, joblib 1.5.3, and
threadpoolctl 3.6.0. Require the doctor report to contain
`[OK] Representative benchmark`; missing packages or version drift fail
explicitly.
Numerical libraries run with a one-thread limit and their non-path runtime
identity is recorded. The resulting assignment digest is frozen, but
cross-architecture bit identity is not guaranteed; compare receipts rather
than assuming it.

The repository also carries the project workflow at
`.agents/skills/coremof-release-curation/`. Codex scans that repository-scoped
location automatically when launched anywhere inside the clone; invoke it as
`$coremof-release-curation`. If a newly pulled update is not visible, restart
Codex. See the [official Codex skill documentation](https://developers.openai.com/codex/build-skills).

The compact handoff's release root is loader-complete for
`verify_cif_files=False`: it contains `dataset_info.json`, the checksum-bound
metadata, parents, feature tables, and CIF manifest, but intentionally omits
the CIF bytes. Do not pass `--verify-cifs` with that compact root. If modelling
requires structure files, obtain the separate restricted full archive and
verify its own receipt and SHA-256 before use. A flattened share archive is not
automatically a loadable release root; never mix v26.0.1 and v26.0.2 files.

## 3. Build assignments without targets

```bash
coremof benchmark-cr-ncr /secure/path/to/coremof_v26.0.2 \
  --group-criteria priority_main \
  --cohort-eligibility complete_release_label_pure_effective_blocks \
  --ncr-pool-fractions 0.0 0.2 0.4 0.6 0.8 1.0 \
  --seeds 42 43 44 45 46 \
  --fractions 0.8 0.1 0.1 \
  --output-directory /secure/coremof-ml-work/benchmark_outputs
```

Equivalent Python:

```python
from CoREMOF.dataset import CoREMOFDataset

dataset = CoREMOFDataset.from_release(
    "/secure/path/to/coremof_v26.0.2"
)
classified = dataset.classify("5checker")
suite = classified.build_cr_ncr_benchmark(
    ncr_pool_fractions=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    seeds=(42, 43, 44, 45, 46),
    total_size="full_cr",
    train=0.8,
    val=0.1,
    test=0.1,
    group_criteria="priority_main",
    cohort_eligibility="complete_release_label_pure_effective_blocks",
    diversity="representative",
    test_policy="fixed_pure_cr",
    include_full_cr_diagnostic=True,
)
suite.write("/secure/coremof-ml-work/benchmark_outputs")
```

The common `fixed_pure_cr` test uses whole label-pure effective blocks and is
shared across every ratio and seed. The writer records raw counts, eligible
counts, exclusions, requested and achieved partition counts, assignment
digests, and `official_split=false`. No audited official assignment manifest
currently exists. The supplementary `full_cr_diagnostic` instead covers the
complete raw strict-CR pool, including rows excluded by the label-pure
sensitivity policy; it is not the independent test.

## 4. Attach targets only after assignment

Do not expose adsorption values, availability, or target-derived features to
cohort construction, diversity balancing, or partition assignment. First
verify the frozen suite receipt and assignment digest; then attach the
checksum-bound target table:

```bash
coremof attach-targets /secure/path/to/coremof_v26.0.2 \
  --manifest /secure/coremof-ml-work/benchmark_outputs/coremof_cr_ncr_benchmark/runs/seed42_q0p0.csv \
  --receipt /secure/coremof-ml-work/benchmark_outputs/coremof_cr_ncr_benchmark/receipt.json \
  --config /secure/path/to/targets.json \
  --missing keep \
  --output-directory /secure/coremof-ml-work/attached_targets
```

`keep` preserves every assigned ID and represents unavailable targets as
null. `drop` creates only a derived filtered view; it never refills,
rebalances, or resplits. Target hashes do not change the frozen assignment
receipt, and attached outputs retain `official_split=false`. The command
attaches one run CSV at a time; loop over `runs/*.csv` for persisted suites, or
use `suite.attach_targets(...)` in Python before serialization to attach all 30
runs together. `membership_manifest.csv` is an accounting table with repeated
IDs and is not an attachment manifest.

## 5. Restricted-data transfer boundary

This Git update adds only source code, public documentation, dependency
metadata, and the repository-scoped workflow skill. The hash and transfer
manifests belong to the separately transferred restricted archive, not the Git
repository. The repository still contains legacy tracked SI archives and
therefore is not a sanitized data-free clone; do not add release tables, CIF archives,
structure-resolved target data, checker findings, or derived benchmark
manifests to GitHub, Git LFS, or an unapproved public-cloud service.

Transfer restricted release and target archives separately through an approved
institutional channel, SFTP, or SSH-based `rsync`, and verify their supplied
SHA-256 ledgers before extraction and use. Confirm the recipient's project
access and institutional CSD/CCDC entitlement first. CSD CIFs and
structure-resolved CSD-derived outputs remain licence-gated unless written
permission covers them; SI assets remain rights-pending asset by asset; COD
material retains its applicable licence and attribution requirements; and
MOSAEC/CCDC-derived findings require their own redistribution clearance.

Treat newly generated benchmark and target-attached outputs as restricted until
their row-level redistribution rights are reviewed. A modelling source filter
does not itself sanitize the complete-release leakage graph or grant
redistribution permission.
