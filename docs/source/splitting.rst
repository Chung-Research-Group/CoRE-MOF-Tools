Classification and dataset splitting
====================================

The splitting API loads an existing CoRE MOF release, independently
recomputes its checker-consensus labels, and assigns whole parent groups to
train, validation, or test. It does not run a checker or calculate a missing
RAC, MOFid, or Zeo++ result.

.. warning::

   This is an initial development API. Splits made from the current v26 parent
   tables are reproducible exploratory splits, not an official CoRE MOF
   benchmark. MOFid-dependent parent relations and the recommended main policy
   remain provisional until the pinned MOFid results and rebuilt parent tables
   are promoted. An official split will require a release-provided, audited
   split manifest.

Project-defined split identifiers
---------------------------------

``priority_main``, ``main_union``, ``auto``, and ``parent_only`` are
CoREMOF-tools API identifiers, not community-standard crystallographic terms.
``parent_only`` means that the selected parent method's resolved groups alone
form the split blocks; it adds no relationship from any other method.

At first use here, **release-authorized and available** mean that a criterion's
``status/group/size`` triad is declared in ``parent_group_methods.json``,
stored in ``parent_groups.csv``, and validated by the loader. ``MATCHED``
means available with an observed group size of at least two; ``UNMATCHED``
means available with size one; ``NOT_AVAILABLE`` supplies no edge and its
size-one group only preserves the table shape. The size must equal the number
of release rows carrying that group. Unavailable rows never match and become
unique singletons or explicit exclusions.

At first use here, **exact RAC5** means equality of all 264 ordered finite
values after IEEE-754 binary64 parsing, mapping ``-0.0`` to ``+0.0``, and
``float.hex()`` serialization, with ``rtol=atol=0`` and no scaling, deletion,
imputation, or rounding. **Exact MOFid** means equality of the complete
published string after converting to text; collapsing each Unicode-whitespace
run to one ASCII space; trimming; case-insensitively rejecting an empty field
or ``-``, ``nan``, ``none``, ``null``, ``n/a``, ``na``, ``unknown``,
``missing``, ``timeout``, ``timed out``, ``error``, ``failed``, ``fail``,
``fail process``, ``failed process``, or ``process failed``; applying Unicode
NFKC; and case-folding, in that order. It is not partial/fuzzy matching and
does not change CIF bytes, atoms, bonds, occupancies, coordinates, chemistry,
or unit cells.

``priority_main`` is the explanatory parent hierarchy. It is computed over the
complete release from ``rac_status/rac_group/rac_size``,
``mofid2_status/mofid2_group/mofid2_size``, and
``mofid1_status/mofid1_group/mofid1_size`` in
``parent_groups/parent_groups.csv``. It seeds components from every available
release-authorized exact RAC5 group, then processes exact MOFid-v2 groups, and
finally processes exact MOFid-v1 groups. At each lower step, a group touching no stronger
component creates a component; one touching exactly one stronger component
attaches only its unresolved rows; and one touching two or more stronger
components never merges them, records ``PARENT_METHOD_CONFLICT``, and leaves
lower-only rows unresolved. A structure with none of the three inputs becomes
a unique singleton unless explicit exclusion is requested. This is not a
row-wise first-nonmissing fallback. Zeo++, CrystalNets, source ID, common name,
CIF hash, provisional source-ID/MOFid transitive groups, and
StructureMatcher do not enter this hierarchy.

``main_union`` is a separate conservative leakage guard, not an explanatory
parent claim. It reads full hashes from ``manifests/cif_manifest.csv`` and the
source/RAC5/MOFid status, group, and size columns from the parent table. It forms
full-release transitive connected components from exact full CIF SHA-256,
source-sibling, RAC5, MOFid-v2, and MOFid-v1 edges before any experiment
filter. A missing CIF SHA-256 fails closed; missing optional evidence adds no
edge; nulls never match; and all available edges are unioned without
precedence. The source group compares the ordered ``(source_database,
source_id)`` pair after applying the text procedure above to each field, so
equal IDs in different databases never match. Normal release loading requires
that triad and never uses metadata fallback. A manually constructed object
with no source-criterion fields instead strips and uppercases the database
(``UNKNOWN`` only if blank), and strips, case-folds,
Unicode-whitespace-splits, and ASCII-space-joins the source ID; that fallback
does not apply NFKC or reject nonblank placeholders, and explicit
``NOT_AVAILABLE`` never falls back. ``leakage_guard="auto"`` is only a selector: it chooses
``main_union`` for ``priority_main`` and ``parent_only`` for an explicitly
selected direct/reference method. ``parent_only`` uses only that selected
explanatory grouping as blocks and adds no cross-method edge. An unresolved
``priority_main`` conflict is excluded under ``parent_only`` but may remain
safely assigned and diagnosed under ``main_union``. Consequently, two
different explanatory groups may share one ``main_union`` block without being
labelled the same parent.

Quick example
-------------

Load a release and build a COD/SI CR-versus-NCR split:

.. code-block:: python

   from CoREMOF.dataset import CoREMOFDataset

   ds = CoREMOFDataset.from_release("/path/to/coremof_v26.0.2")
   classified = ds.classify(checkers="5checker")

   # Direct CR/NCR membership for screening or export:
   cr_ids = classified.cr_ids
   ncr_ids = classified.ncr_ids

   split = classified.train_valid_test_split(
       parent_method="priority_main",
       fractions=(0.8, 0.1, 0.1),
       random_state=42,
       labels=("CR", "NCR"),
       sources=("COD", "SI"),
       leakage_guard="auto",
   )

The result keeps an indivisible leakage block in exactly one partition. The
achieved fractions can differ from the requested fractions when a block is
large; this is expected and is reported in the result summary.

The COD/SI filter above is applied *after* leakage blocks are constructed from
the complete COD+CSD+SI release universe. It therefore requires the authorized
full-release tables and manifests. A standalone open COD/SI bundle cannot
reproduce a bridge mediated by a hidden CSD row unless the release supplies an
audited projection of the full-universe block IDs.

.. code-block:: python

   print(split.counts)
   train_idx, valid_idx, test_idx = split
   split.write("model_splits", stem="cod_si_5checker_seed42")

Use the explicit ``train_ids``, ``validation_ids``, and ``test_ids`` fields
when saving or joining a split. ``train_indices``, ``validation_indices``, and
``test_indices`` are zero-based positions in the classified view passed to the
splitter. The corresponding ``*_release_indices`` fields refer to the original
``metadata.csv`` order. In the assignment CSV, ``view_index`` is blank for a
row outside the classified view, whereas ``release_index`` is always present.
Neither index is persistent; ``structure_id`` is the stable cross-table key.

A compact one-call form is also available:

.. code-block:: python

   from CoREMOF.splitters import split_release

   split = split_release(
       "/path/to/coremof_v26.0.2",
       checkers="3checker",
       parent_method="rac5",
       fractions=(0.8, 0.1, 0.1),
       random_state=42,
       labels=("CR", "NCR"),
       sources=("COD", "SI"),
   )

Five separate choices
---------------------

``checkers``, ``parent_method``, ``leakage_guard``, the row filters, and
``required_targets`` answer different questions. They are intentionally not
combined into one opaque "split mode".

``checkers``
   Decides how the classification label is recomputed. It does not describe
   structural relatedness.

``parent_method``
   Selects the scientific parent or duplicate relation to analyse. It does not
   change a CR/NCR label.

``leakage_guard``
   Decides which known relations must remain together during partitioning. It
   can be broader than the selected explanatory parent method.

Filters
   Select eligible labels, source databases, structure variants, metals, or
   explicit structure IDs. Filtering does not relabel a row and does not erase
   a relation that was found in the complete release universe.

``required_targets``
   Names the target endpoints that an eligible row must actually contain.
   Target availability is checked before assignment so target-null rows are
   excluded from train/validation/test, while the leakage graph still uses the
   complete release universe and therefore keeps hidden bridge relations.

Checker views
-------------

The three official views are:

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - View
     - Included checkers
   * - ``3checker``
     - MOFClassifier, original MOFChecker, and Chen-Manz
   * - ``4checker``
     - The 3-checker view plus MOSAEC
   * - ``5checker``
     - The 4-checker view plus SETC-GAT

Every view uses strict consensus:

* all included checkers ``PASS`` gives ``CR``;
* all included checkers ``FAIL`` gives ``NCR``;
* a complete PASS/FAIL mixture gives ``AMBIGUOUS``; and
* any unavailable or non-voting result gives ``UNCHECKED``.

There is no majority vote. A timeout, execution error, unsupported input, or
missing upstream output can never become ``NCR``.

An explicit ordered sequence of canonical checker names is accepted for a
user experiment:

.. code-block:: python

   custom = ds.classify(
       checkers=("MOFClassifier", "MOFChecker", "MOSAEC")
   )

The same strict rule applies, but the receipt identifies this as
``USER_DEFINED`` with a stable ``custom:...`` identifier. It is not presented
as an official release label.

When ``split_release()`` receives an already classified dataset, omitting
``checkers`` preserves that dataset's view. Supplying ``checkers`` explicitly
requires an exact match and raises on a conflict; an explicit request is never
silently ignored. A release path or an unclassified dataset still defaults to
the 5-checker view.

Parent methods
--------------

Use the exact API keys below.

Here **canonicalized identifier text** has a narrow release-defined meaning:
convert the published value to text, collapse every Unicode-whitespace run to
one ASCII space, trim it, reject a whole-field missing/execution
placeholder, apply Unicode NFKC, and apply Unicode case-folding, in that
order. Rejected case-insensitive whole fields are empty text, ``-``, ``nan``,
``none``, ``null``, ``n/a``, ``na``, ``unknown``, ``missing``,
``timeout``, ``timed out``, ``error``, ``failed``, ``fail``, ``fail process``,
``failed process``, and ``process failed``. Inherited v26.0.1 cleanup
additionally rejects a whole-field ``*`` and leading ``MOFid-v1.NA...`` and
``MOFid-v2.NA...`` placeholder families. A source key is the exact ordered pair of separately canonicalized
``source_database`` and ``source_id``. A MOFid key is the complete
canonicalized published MOFid string, never a prefix or fuzzy match.
This text processing does not change or modify CIF bytes, atoms, bonds,
occupancies, coordinates, chemistry, or unit cells.

The audited v26.0.1 identity components retained in v26.0.2 used a narrow v11
cleanup. For a v11 refcode, it converts to text; removes the first semicolon
and everything after it; changes backslashes to ``/``; removes the path and
one terminal ``.cif``/``.cif.gz``; trims; deletes every Unicode whitespace
character; case-folds without NFKC; rejects literal placeholders; removes one
terminal ``_ASR_pacman``, ``_FSR_pacman``, ``_ION_pacman``,
``_ION_ASR_pacman``, or ``_ION_FSR_pacman``; and rechecks the result. Those
preserved v11 refcode edges did not add a database namespace; new v26.0.2
source-ID edges do.

For a v11 MOFid, the seed removes the first semicolon record suffix and trims;
rewrites only literal ``MOFidv2.`` to ``MOFid-v2.``; collapses whitespace,
trims, and case-folds without NFKC; removes terminal ``.no_ref`` only when a
MOFid-v1 marker is present; then rejects literal execution placeholders and
the leading ``MOFid-v1.NA...``/``MOFid-v2.NA...`` families. No step modifies
an atom, bond, occupancy, coordinate, unit cell, general chemical punctuation,
or CIF byte. The splitter consumes the release-authorized group/status/size
columns and does not rerun this cleanup.

The **RAC5 fingerprint** is exact equality of all 264 finite descriptor
columns listed in order at ``criteria.rac5.ordered_descriptors`` in
``parent_group_methods.json``. It uses no scaling, feature deletion,
imputation, or tolerance: every value is parsed as finite binary64, ``-0.0`` is
mapped to ``+0.0``, and ``float.hex()`` tokens are compared in order with
``rtol=atol=0``. The **selected Zeo++ fingerprint** is exact equality
of ``density_g_cm3``, largest cavity, pore-limiting, and largest-free-path
diameters; accessible and nonaccessible N2 surface areas in m2/cm3 and m2/g;
accessible and nonaccessible N2 volumes in cm3/g and as fractions; and
``he_void_fraction``, plus exact N2 channel dimension and available
bonded-framework periodic dimension. Its 13 numbers use the same finite
binary64, negative-zero, and ``float.hex()`` procedure; its numeric tolerance
is zero and N2/He probe radii are 1.655/1.32 A. One invalid number or hard gate
makes the fingerprint unavailable. Unit-cell-extensive areas/volumes, OMS,
topology labels, and framework-component counts are excluded.

The **current CrystalNets scientific fingerprint**, where current means the
topology evidence authorized by the loaded release rather than a runtime
search for newer output, requires ``SUCCESS``,
``topology_available=true``, ``error=null``, and nonempty complete SingleNodes
and AllNodes subnets. It includes network dimension, interpenetrated-subnet
count, catenation degree, subnet count, top-level single-node/all-node nets and
agreement, and, for every canonically sorted subnet, single/all agreement plus
each node view's status, dimension, topology key, topology name, and
topological genome. Canonical sorted-key JSON is used; complete subnet
projections are sorted, duplicates are retained, and the result is SHA-256
hashed. Runtime, CIF paths/hashes, diagnostics, software boilerplate, and
original subnet order are excluded. Both count fields must be integers at
least one and equal the subnet count; subnet indices must be unique/contiguous
before sorting. Top-level dimension, nets, and agreement may be null for
heterogeneous subnets and remain part of the fingerprint. Each node view must
have ``SUCCESS``, dimension 0--3, and a nonblank topology key; topology name
and genome may be null. An incomplete/unsuccessful result adds no
edge; this is exact fingerprint equality, not a topology-similarity tolerance.

.. list-table::
   :header-rows: 1
   :widths: 22 20 58

   * - ``parent_method``
     - Role
     - Relation
   * - ``priority_main``
     - Recommended
     - Conflict-aware RAC5, then MOFid v2, then MOFid v1 evidence hierarchy
   * - ``rac5``
     - Main/direct
     - Exact equality of every validated depth-5 RAC descriptor
   * - ``mofid_v2``
     - Main/direct
     - Exact equality of complete canonicalized release-authorized MOFid v2 text; missing values never match; current v26 input remains provisional until pinned-result promotion
   * - ``mofid_v1``
     - Main/direct
     - Exact equality of complete canonicalized MOFid v1 text; missing values never match
   * - ``rac5_zeo``
     - Reference
     - Exact combined RAC5 and release-selected Zeo++ fingerprint
   * - ``rac5_topology``
     - Optional reference
     - ``RT-`` group: exact 264-value RAC5 plus a complete successful current CrystalNets scientific fingerprint
   * - ``mofid_v2_topology``
     - Optional reference
     - ``M2T-`` group: exact complete canonicalized MOFid v2 plus that CrystalNets fingerprint; provisional whenever its MOFid-v2 input is provisional
   * - ``structure_matcher_strict``
     - Optional reference
     - ``SM-`` connected component of exhaustive composition-compatible pairs whose forward/reverse pymatgen 2024.2.8 ElementComparator fits pass with ltol=stol=0.001, angle_tol=0.01, primitive_cell/attempt_supercell true, scale/allow_subset false, supercell_size=num_sites, and no ignored species
   * - ``zeo``
     - Reference
     - Exact release-selected Zeo++ fingerprint
   * - ``source_id``
     - Reference
     - Exact canonicalized ``(source_database, source_id)`` sibling key; no cross-database match
   * - ``common_name``
     - Reference
     - Exact NFKC/whitespace/case-folded common name; sparse and non-unique
   * - ``identity_union``
     - Reference
     - Provisional source-ID/MOFid transitive groups: preserve audited v26.0.1 components, then union exact canonicalized source-ID, MOFid-v2, or MOFid-v1 edges; no structural descriptor input
   * - ``none``
     - Control
     - One independent singleton per structure

The recommended direct sensitivity analyses are ``rac5``, then
``mofid_v2``, then ``mofid_v1``. The other methods are kept as references;
they should not be interpreted as interchangeable proofs of parentage.

Exact RAC or RAC/Zeo equality is fingerprint equivalence, not proof of a
shared synthetic parent. The project-defined ``identity_union`` key means
exactly the provisional source-ID/MOFid transitive groups stated above. Each
group, and therefore each value counted by ``identity_size``, is one
transitive connected component of structures joined by those
identifier-equality edges; it is not a count of edges or identifiers. All
three edge types are unioned with equal status and no priority or conflict
rule, so one
can bridge components made by another. A null or rejected placeholder adds no
edge and missing identifiers never match. The relation uses no RAC5, Zeo++,
CrystalNets topology, CIF hash, common name, or StructureMatcher input; it is
not proof of chemical or geometric identity and is excluded from both
``priority_main`` and ``main_union``. Current v26 identity groups remain
provisional until the pinned MOFid bundle is promoted and the parent table is
rebuilt.

``RT-``, ``M2T-``, and ``SM-`` are criterion prefixes, not scientific
results. They mean an exact RAC5-plus-current-CrystalNets group, an exact
MOFid-v2-plus-current-CrystalNets group, and a connected component of strict
direct StructureMatcher edges, respectively. The following compact digest is
only a group label, not a topology name, MOFid, RMSD, or similarity score.
The release builder hashes the criterion name and its complete comparison key
with length-delimited UTF-8 SHA-256, starts with eight uppercase hexadecimal
characters after the prefix, and extends every actually colliding prefix one
character at a time until unique in the v26.0.2 superset.
Missing or failed RAC5/CrystalNets/MOFid input creates no ``RT-`` or ``M2T-``
evidence. Directional fit disagreement or an unresolved pair creates no
strict edge; an incomplete strict component is projected to structure-specific
``NOT_AVAILABLE`` singletons.

``structure_matcher_strict`` is accepted only when the release declares its
audited ``sm_*`` columns. The frozen method uses Python 3.9, pymatgen 2024.2.8,
and NumPy 1.26.4. It blocks candidates by the parsed structure's
``ElementComparator`` fractional-composition hash and attempts every pair in
an equal block; routine formula, DOI, source, RAC5, topology, site count, and
volume are not exclusionary blocks. Each forward/reverse call uses
``fit(..., symmetric=True)`` with ``ltol=stol=0.001``,
``angle_tol=0.01``, ``primitive_cell=true``, ``scale=false``,
``attempt_supercell=true``, ``allow_subset=false``,
``supercell_size=num_sites``, and no ignored species. An undirected direct edge
requires both calls to succeed.

The pinned direct ``CifParser`` expands declared symmetry, merges generated
sites at ``site_tolerance=1e-4``, rounds coordinates near 1/3 or 2/3 at
``frac_tolerance=1e-4``, checks occupancy, sorts the periodic structure, and
preserves disorder. It performs no manual repair, occupancy choice, atom
deletion, or chemistry edit. Parser, timeout, OOM, and matcher failures are
``NOT_AVAILABLE`` rather than nonmatches; directional disagreement is
``ASYMMETRIC_RESULT`` and creates no edge.

The group is a convenience connected-component view over authoritative direct
edges. Since the tolerance relation can be non-transitive, a component does
not assert that every member pair directly matches; degree, edge count,
possible-edge count, completeness, clique status, and the direct ledger must
be consulted. Forward and reverse normalized RMS/maximum displacement divide
periodic site displacement by the corresponding direction's
``(V/Nsites)^(1/3)``. They are dimensionless, are not angstrom RMSD, and do not
come from the separate ``charnley/rmsd`` molecular Kabsch program. This
optional method changes neither ``priority_main`` nor ``main_union``; the old
``ltol=0.2``, ``stol=0.3``, ``angle_tol=5``, ``scale=true`` profile is
documentation-only and is neither executed nor exposed.

The loader accepts ``sm_*`` columns only when
``criteria.structure_matcher_strict`` declares the exact
``pymatgen_structure_matcher_strict_v2`` optional-reference contract. The
contract requires direct symmetric strict-match edges, identifies connected
components as a non-transitive convenience view, and keeps the criterion out
of both default hierarchies. Incomplete components must be represented by
unique ``NOT_AVAILABLE`` singletons, never ``PARTIAL`` public rows.

The method declaration binds
``parent_groups/structure_matcher_strict_evidence_receipt.json`` by a full
lowercase SHA-256. That passing release-adapter receipt is itself bound to the
exact parent table and strict pair-ledger hash and records pair, edge,
structure, and unavailable counts. Both the declaration and receipt must say
that historical relaxed evidence was neither executed nor exposed. The receipt hash is included
among the dataset and split input hashes.

The topology-combined methods are accepted only when the release declares
their validated parent columns. They remain reference/sensitivity criteria
and do not change ``priority_main``. Incomplete or unsuccessful CrystalNets
results never match. A release with provisional MOFid v2 evidence must also
mark ``mofid_v2_topology`` provisional in its method contract.

Missing parent evidence
~~~~~~~~~~~~~~~~~~~~~~~

The default is ``missing_parent="singleton"``. Each unavailable criterion gets
its own structure-specific group and therefore never matches another missing
value. Set ``missing_parent="exclude"`` only when the experiment deliberately
requires complete parent evidence; exclusions are reported and are never
silent.

A conflict in ``priority_main`` is different from missing evidence. When a
lower-priority MOFid relation would join two stronger RAC-anchored components,
the explanatory resolver does not silently merge those anchors. With the
default ``auto``/``main_union`` guard, the row remains assigned inside its full
leakage block and is marked ``PARENT_METHOD_CONFLICT`` in
``parent_diagnostics``; it is not excluded. If ``priority_main`` is explicitly
paired with ``parent_only``, the unresolved row is excluded.

Every lower-priority group that spans multiple stronger components is also
retained in the structured ``parent_conflicts`` ledger, including conflicts
that contain only already anchored rows and therefore require no exclusion.
Each entry records the lower method/group, stronger components, component
members, all affected members, and unresolved members. The row-oriented
``parent_diagnostics`` and group-oriented ledger answer different questions.

Leakage guards
--------------

The accepted guards are:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Guard
     - Behavior
   * - ``auto``
     - Resolve to ``main_union`` for ``priority_main`` and to ``parent_only``
       for an explicit direct or reference method
   * - ``main_union``
     - Keep every transitive component formed by exact CIF identity,
       source-sibling, available RAC5, available MOFid v2, and available MOFid
       v1 edges in one partition
   * - ``parent_only``
     - Keep only the selected method's published grouping together

``priority_main`` is not implemented as a naive row-by-row fallback. Such a
fallback can choose RAC5 for one row and MOFid for another even though a valid
lower-level or transitive relation connects them. With ``leakage_guard="auto"``,
the explanatory priority remains RAC5 > MOFid v2 > MOFid v1, but assignment is
constrained by the complete ``main_union`` component.

Leakage blocks are built on the complete release before applying ``labels``,
``sources``, ``variants``, ``metals``, or explicit-ID filters. A filtered-out
row can therefore remain an invisible bridge between two selected rows. This
prevents a convenient filter from creating train/test leakage.

Target and feature joins
------------------------

Join external target results before splitting so structures without the
required endpoints are excluded before partition assignment. CSV, JSON, and
JSONL target files can be combined in one modelling table. Units, conditions,
and CSV type conversions are explicit; the package does not infer or impute
them.

.. code-block:: python

   from CoREMOF.dataset import CoREMOFDataset
   from CoREMOF.targets import AliasRegistry, TargetSource

   ds = CoREMOFDataset.from_release("/path/to/coremof_v26.0.2")
   sources = (
       TargetSource(
           "/path/to/uptake.csv",
           id_column="structure_name",
           target_columns=("uptake",),
           target_names={"uptake": "xe_uptake"},
           value_types={"xe_uptake": "float"},
           units={"xe_uptake": "mol/kg"},
           conditions={"xe_uptake": {"temperature_K": 298, "pressure_bar": 1}},
       ),
       TargetSource(
           "/path/to/selectivity.jsonl",
           target_columns=("selectivity",),
           units={"selectivity": "dimensionless"},
           conditions={"selectivity": {"temperature_K": 298}},
       ),
   )
   modelling = ds.merge_targets(
       sources,
       feature_tables=("rac5", "zeo", "topology"),
   )
   classified = modelling.classify("5checker")
   split = classified.train_valid_test_split(
       parent_method="priority_main",
       required_targets=("xe_uptake", "selectivity"),
       required_target_mode="all",
       random_state=42,
   )

Exact current public IDs are accepted directly. Earlier/v11 identifiers are
accepted only through an explicit hash-bound registry, with every allowed
alias column named by the caller:

.. code-block:: python

   aliases = AliasRegistry(
       "/path/to/audited_name_registry.csv",
       current_id_column="structure_id",
       alias_columns=("previous_public_id", "v11_structure_id"),
   )
   modelling = ds.merge_targets(sources, alias_registry=aliases)

There is no fuzzy matching. Ambiguous aliases, unknown IDs, conflicting
duplicate values, or conflicting definitions for one target name fail the
complete merge. Null values remain null. Output rows contain current IDs only.
Use ``target_names`` to map a generic input column such as ``uptake`` to a
canonical endpoint such as ``xe_uptake`` without editing a source file.
Units, conditions, and value types are declared against canonical names;
provenance retains the original input column.

Scientific meaning is also preserved literally. A Rosenbluth weight should be
named, for example, ``rosenbluth_weight`` and declared ``dimensionless``; it
must not be relabelled as a Henry coefficient or uptake without a separate,
documented conversion and a new canonical target.

``TargetMergedDataset.write()`` emits a merged CSV, observation-level
provenance JSONL, and hash-bound receipt. The receipt is also embedded in a
subsequent split receipt. A required-target exclusion uses
``MISSING_REQUIRED_TARGET``. Crucially, the ``main_union`` leakage graph is
still constructed from the complete release, including structures with no
target data.

For repeatable command-line workflows, put the target sources, definitions,
optional alias registry, and requested feature tables in a JSON configuration:

.. code-block:: bash

   coremof merge-targets /path/to/coremof_v26.0.2 \
     --config targets.json \
     --output-directory model_inputs

   coremof split /path/to/coremof_v26.0.2 \
     --target-config targets.json \
     --require-target xe_uptake \
     --require-target selectivity \
     --output-directory model_splits

High-throughput screening
-------------------------

``examples/screen_candidates.py`` validates, filters, and ranks a release or
target-merged table in one unattended command:

.. code-block:: bash

   python examples/screen_candidates.py /path/to/coremof_v26.0.2 \
     --target-config targets.json \
     --rank-by xe_uptake \
     --require-target xe_uptake \
     --checkers 5checker \
     --label CR \
     --source COD \
     --order descending \
     --limit 1000 \
     --split \
     --parent-method priority_main \
     --leakage-guard auto \
     --output-directory screening/xe_uptake

Checker, identity, and required-target filters run before ranking. Null,
non-numeric, and non-finite values are excluded without imputation, and ties
use ascending ``structure_id``. Integer values and decimal text are ranked
without conversion through binary64, so values above ``2**53`` remain exact;
native floats retain their IEEE-754 values. Booleans, nulls, NaN, and
infinities are not rankable. Required targets are flat CSV columns. The receipt
hashes the script, exact package source closure, release and embedded
target/config/alias/feature inputs, ranked CSV, and optional split outputs.
An optional split uses the complete release leakage graph and remains
exploratory rather than official.

Both standalone split receipts and the nested split section of screening
receipts retain the selected parent method, its machine-readable definition,
the requested guard (for example ``auto``), the resolved guard (for example
``main_union``), and the resolved guard's machine-readable definition. Thus a
receipt does not require a reader to infer policy semantics from a coined
identifier.

The historical local ``henry.txt`` files contain dimensionless average
Rosenbluth weights, not Henry coefficients, uptake, or selectivity. Ranking
them in descending order demonstrates the workflow only.

Filters do not change labels
----------------------------

Classification retains all four labels unless a filter is requested:

.. code-block:: python

   all_labels = ds.classify(checkers="3checker")
   modelling = all_labels.filter(
       labels=("CR", "NCR"),
       sources=("COD",),
       variants=("ASR", "FSR"),
       metals=("Cu", "Zn"),
   )

The train/validation/test convenience methods default to
``labels=("CR", "NCR")``. AMBIGUOUS and UNCHECKED structures are not lost:
they remain explicit exclusions in the assignment table and receipt. Change
``labels`` deliberately to include them. A structure passes the metal filter
when its published metal list contains at least one requested element.

Stratification is independent from filtering. The default
``stratify_by=("label",)`` tries to preserve the selected label distribution
while assigning complete blocks. Additional public metadata fields can be
requested, for example:

.. code-block:: python

   split = modelling.train_valid_test_split(
       parent_method="priority_main",
       stratify_by=("label", "source_database"),
       random_state=42,
   )

More or smaller strata can make exact balance harder. The grouping constraint
always wins over a requested stratum fraction.

Compare parent assumptions
--------------------------

Run the direct methods separately when assessing how sensitive a model is to
the parent definition:

.. code-block:: python

   results = {}
   for method in ("rac5", "mofid_v2", "mofid_v1"):
       results[method] = classified.train_valid_test_split(
           parent_method=method,
           leakage_guard="parent_only",
           fractions=(0.8, 0.1, 0.1),
           random_state=42,
           labels=("CR", "NCR"),
       )

These are three different scientific assumptions, not three replicas of the
same split. The result receipt records the method, resolved guard, missing-data
policy, selected IDs, requested and achieved fractions, label counts, seed,
algorithm version, input fingerprints, and the parent-conflict ledger.
Implementation versions and source hashes are frozen when the result is
constructed, so a later code edit cannot silently change the receipt for an
existing in-memory split.

Reproducibility and balance
---------------------------

Assignment is deterministic for the same validated release, parameters, and
``random_state``. Metadata row order does not change the
``structure_id``-to-partition mapping or the assignment digest. Positional
indices and complete receipt bytes may change because they intentionally bind
the supplied release/view order and input files. Each guard block is
indivisible, so exact 80/10/10 counts may be impossible. The splitter optimizes
total-size and per-label balance at block level and reports the achieved
fractions, per-label counts, largest block, and deviations.

Changing ``random_state`` may change which complete blocks enter a partition;
it never changes labels, parent groups, or block membership.

Provisional and official splits
-------------------------------

A split generated by the current API is exploratory and emits
``official_split=false``. Passing ``official=True`` fails closed because the
current releases contain no audited official assignment manifest. Loading such
a manifest is a planned feature, not an inferred mode.

When a future release adds structures to a frozen base split, the package will
not move base rows silently. If an addition connects parent groups already
frozen in different partitions, the addition is reported as
``SPLIT_BRIDGE_CONFLICT`` and excluded from the modelling partitions while the
scientific relation is preserved. This frozen-extension logic is not yet
implemented. In the current independent exploratory runs, 124 of 8,593 shared
5-checker CR/NCR IDs moved between separately generated v26.0.1 and v26.0.2
splits; that is acceptable only because both are explicitly provisional.

The current public v26.0.2 parent table also authorizes fewer RAC parent edges
than the RAC feature table contains: 25,618 rows have an available RAC parent
status, while 29,891 rows have complete RAC features. This is intentional at
the package boundary. The splitter trusts published parent statuses and never
constructs an unapproved edge from feature values.

Data and licensing
------------------

Loading and splitting use release tables and manifests by default. They do not
need a CCDC
licence or installations of MOSAEC, SETC-GAT, Zeo++, molSimplify, MOFid, or
CrystalNets. The splitter does not copy CIFs.

Default validation checks manifest membership, paths, declared sizes, and full
SHA-256 syntax but does not rehash the CIF bytes. Use
``CoREMOFDataset.from_release(path, verify_cif_files=True)`` or the one-call
``split_release(..., verify_cif_files=True)`` form (CLI:
``--verify-cifs``) for byte-level size and hash verification. ``main_union``
fails closed if any structure lacks a full manifest SHA-256.

The ability to select CSD or SI rows does not grant redistribution rights.
COD is the clearly open source, but the API applies no source filter when
``sources=None``. CSD structure files remain licence-gated, and SI assets
require an asset-level rights check before redistribution.
