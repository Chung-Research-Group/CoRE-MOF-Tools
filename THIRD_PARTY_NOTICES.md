# Third-party components and data

The root LICENSE records the project's CC-BY-4.0 declaration. Third-party
components retain their own terms. This notice preserves attribution and does
not relicense an upstream work or establish permission for every database asset.

## External checker results, not bundled checker implementations

The maintained distribution excludes the MOSAEC implementation and reference
tables, and the MOFChecker/Chen–Manz/MOSAEC/SETC-GAT execution workers. Small
compatibility modules only report that execution is unavailable. They contain
no checker algorithm. Earlier private versions and their applicable licence
notices are retained separately for historical reproduction, not redistributed
as part of this package.

Precomputed checker findings retain their method attribution. For new
calculations use the original authors' software, including
[MOSAEC](https://github.com/uowoolab/MOSAEC) and
[MOFChecker](https://github.com/Au-4/mofchecker_2.0), under their own terms.
The CSD Python API is a separate licensed dependency for CSD retrieval and is
not supplied here. Removing checker software does not establish blanket
permission to redistribute checker outputs or source structures.

## Heat-capacity implementation

`CoREMOF/models/cp_app/` adapts
[tools-cp-porousmat](https://github.com/SeyedMohamadMoosavi/tools-cp-porousmat).
Copyright (c) 2022 Seyed Mohamad Moosavi. Retain the full MIT license in
`licenses/cp-app-MIT.txt`. The inspected upstream revision is
`50bf59b55002177fe2d01b79c7b321ae90564c42`. Local adaptations include input
validation, dependency compatibility, isolated temporary files and integration.
Method: https://doi.org/10.1038/s41563-022-01374-3.

The 300 heat-capacity model files in a full local checkout are not bundled in
the wheel or source distribution. Their integrity manifest identifies local
assets, not proof of their original publication or redistribution permissions.
The predictor does not download missing models automatically.

## Historical stability models

The separately retained local assets `final_model_T_few_epochs.h5` and
`final_model_flag_few_epochs.h5` match the corresponding
[MOFSimplify](https://github.com/hjkgrp/MOFSimplify) files byte for byte at
commit `5693968b3e9b9e26eab3bdb1db908ae2877d4bb7`.
Copyright (c) 2023 Kulik Group. Retain the MIT license in
`licenses/MOFSimplify-MIT.txt`. Method: https://doi.org/10.1021/jacs.1c07217.

The local scalers and water model are separate historical assets. Their
hashes are checked before use, but the exact original export records have not
been established by the package audit. The related WS24 archive
(https://doi.org/10.5281/zenodo.12110918) declares CC-BY-4.0. Its files do not
match these local water-model/scaler bytes, so that record alone is not an
asset-identity certificate. Method: https://doi.org/10.1021/jacs.4c05879.
Do not substitute the later CoREMOF-COD benchmark models for these assets.

## Database tables and node structures

The legacy `CoREMOF/data/CR.json` and `NCR.json` contain structure-resolved
information and are not the current CoREMOF-COD release. The bundled
`CoREMOF/data/mofid/nodes.zip` contains 1,182 node structures. These files are
excluded from the code-only contribution. The node archive
matches the project upstream copy, but its asset-specific source permissions
need separate documentation. A code license is not blanket clearance for
CSD-derived, SI-derived or other third-party structure information.

The current distributions are private review artifacts until the outstanding
data and combined-software redistribution questions are resolved. Scientific
citations supplement, rather than replace, applicable license notices.
