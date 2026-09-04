"""CoRE MOF database utilities and property-prediction tools.

The package root eagerly exports its lightweight input, splitting, benchmark,
and deferred-target APIs. Optional numerical and domain-specific scientific
backends remain lazily imported by the feature that needs them; for example,
``from CoREMOF.calculation import Zeopp`` accesses the Zeo++ integration.
"""

from .inputs import collect_cifs, resolve_cif_inputs

__version__ = "0.4.0.dev0"

from .attachments import (
    TargetAttachedBenchmarkSuite,
    TargetAttachedView,
    TargetAttachmentError,
    attach_targets,
)
from .benchmarks import (
    BenchmarkDependencyError,
    BenchmarkError,
    BenchmarkFeasibilityError,
    CRNCRBenchmarkRun,
    CRNCRBenchmarkSuite,
    CRNCRCohort,
    CRNCRCohortSuite,
    DataSplitResult,
    DiversityIndex,
    available_group_criteria,
    build_cr_ncr_benchmark,
    build_cr_ncr_cohorts,
    data_split,
)

__all__ = [
    "__version__",
    "BenchmarkDependencyError",
    "BenchmarkError",
    "BenchmarkFeasibilityError",
    "CRNCRBenchmarkRun",
    "CRNCRBenchmarkSuite",
    "CRNCRCohort",
    "CRNCRCohortSuite",
    "DataSplitResult",
    "DiversityIndex",
    "TargetAttachedBenchmarkSuite",
    "TargetAttachedView",
    "TargetAttachmentError",
    "attach_targets",
    "available_group_criteria",
    "build_cr_ncr_benchmark",
    "build_cr_ncr_cohorts",
    "collect_cifs",
    "data_split",
    "resolve_cif_inputs",
]
