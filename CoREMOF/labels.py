"""Strict CoRE-MOF checker consensus labels.

This module deliberately has no scientific-Python dependencies.  It implements
the public three-, four-, and five-checker definitions used by the v26 release
metadata.  Execution failures and other known non-votes are never interpreted
as scientific failures.
"""

from collections.abc import Sequence as SequenceABC
from types import MappingProxyType
from typing import Iterable, Mapping, Sequence, Tuple, Union


CHECKER_COLUMNS = MappingProxyType(
    {
        "MOFClassifier": "mofclassifier_status",
        "MOFChecker": "mofchecker_status",
        "Chen-Manz": "chen_manz_status",
        "MOSAEC": "mosaec_status",
        "SETC-GAT": "setc_gat_status",
    }
)

CHECKER_PRESETS = MappingProxyType(
    {
        "3checker": ("MOFClassifier", "MOFChecker", "Chen-Manz"),
        "4checker": (
            "MOFClassifier",
            "MOFChecker",
            "Chen-Manz",
            "MOSAEC",
        ),
        "5checker": (
            "MOFClassifier",
            "MOFChecker",
            "Chen-Manz",
            "MOSAEC",
            "SETC-GAT",
        ),
    }
)

LABELS = frozenset({"CR", "NCR", "AMBIGUOUS", "UNCHECKED"})
VOTE_STATUSES = frozenset({"PASS", "FAIL"})

# These are execution/data-availability states seen in CoRE-MOF curation
# evidence.  They are explicit non-votes.  Keeping this list closed is
# intentional: a new or misspelled token must not silently change a label.
NONVOTE_STATUSES = frozenset(
    {
        "NOT_AVAILABLE",
        "UNCHECKED",
        "ERROR",
        "FAILED",
        "INPUT_ERROR",
        "PARSE_ERROR",
        "PROCESS_ERROR",
        "TIMEOUT",
        "NOT_APPLICABLE",
        "NOT_APPLICABLE_NO_METAL",
        "MISSING",
        "NOT_RUN",
        "PLANNED",
        "PENDING",
        "QUEUED",
        "RUNNING",
        "UNKNOWN",
        "SKIPPED",
        "UNAVAILABLE",
        "CANCELLED",
        "BLOCKED",
        "BLOCKED_CCDC",
        "BLOCKED_CCDC_MOSAEC",
    }
)


class UnknownCheckerStatusError(ValueError):
    """Raised when a checker status is not part of the closed vocabulary."""


def resolve_checker_preset(
    preset: Union[int, str]
) -> Tuple[str, Tuple[str, ...]]:
    """Return the canonical preset name and ordered checker names.

    Only the scientifically defined 3-, 4-, and 5-checker views are accepted.
    """

    if isinstance(preset, bool):
        raise ValueError("checker preset must be 3, 4, 5, or '<n>checker'")
    if isinstance(preset, int):
        name = "{}checker".format(preset)
    elif isinstance(preset, str):
        name = preset
    else:
        raise TypeError("checker preset must be an integer or string")
    if name not in CHECKER_PRESETS:
        raise ValueError(
            "unknown checker preset {!r}; choose one of {}".format(
                preset, ", ".join(CHECKER_PRESETS)
            )
        )
    return name, CHECKER_PRESETS[name]


def resolve_checker_view(
    view: Union[int, str, Sequence[str]]
) -> Tuple[str, Tuple[str, ...], bool]:
    """Resolve an official preset or an explicit ordered custom checker view.

    Custom views must use canonical checker names, may not contain duplicates,
    and receive a ``custom:`` identifier.  The boolean return value indicates
    whether the view is an official published preset.
    """

    if isinstance(view, (int, str)) and not isinstance(view, bool):
        name, checkers = resolve_checker_preset(view)
        return name, checkers, True
    if isinstance(view, bool):
        resolve_checker_preset(view)
    if not isinstance(view, SequenceABC):
        raise TypeError(
            "checker view must be an official preset or an ordered checker sequence"
        )
    checkers = tuple(view)
    if not checkers:
        raise ValueError("a custom checker view may not be empty")
    if any(not isinstance(checker, str) for checker in checkers):
        raise TypeError("custom checker names must be strings")
    unknown = [checker for checker in checkers if checker not in CHECKER_COLUMNS]
    if unknown:
        raise ValueError(
            "unknown checker(s) in custom view: {}".format(
                ", ".join(repr(checker) for checker in unknown)
            )
        )
    if len(set(checkers)) != len(checkers):
        raise ValueError("a custom checker view may not contain duplicates")
    identifier = "custom:{}".format("+".join(checkers))
    return identifier, checkers, False


def consensus_label(statuses: Iterable[str]) -> str:
    """Classify an ordered collection of checker statuses.

    ``PASS`` and ``FAIL`` are the only scientific votes.  All PASS is ``CR``;
    all FAIL is ``NCR``; a complete PASS/FAIL mixture is ``AMBIGUOUS``.  Any
    known non-vote makes the result ``UNCHECKED``.  Unknown tokens raise instead
    of being guessed, which prevents misspellings from becoming labels.
    """

    values = tuple(statuses)
    if not values:
        raise ValueError("at least one checker status is required")

    unknown = []
    for value in values:
        if not isinstance(value, str):
            unknown.append(repr(value))
        elif value not in VOTE_STATUSES and value not in NONVOTE_STATUSES:
            unknown.append(repr(value))
    if unknown:
        raise UnknownCheckerStatusError(
            "unknown checker status token(s): {}".format(", ".join(unknown))
        )

    if any(value in NONVOTE_STATUSES for value in values):
        return "UNCHECKED"
    if all(value == "PASS" for value in values):
        return "CR"
    if all(value == "FAIL" for value in values):
        return "NCR"
    return "AMBIGUOUS"


def classify_checker_row(
    row: Mapping[str, str],
    preset: Union[int, str, Sequence[str]] = "5checker",
) -> str:
    """Recompute one release row's consensus label for ``preset``."""

    preset_name, checkers, _ = resolve_checker_view(preset)
    statuses = []
    for checker in checkers:
        column = CHECKER_COLUMNS[checker]
        if column not in row:
            raise KeyError(
                "missing checker column {!r} for {}".format(column, preset_name)
            )
        statuses.append(row[column])
    return consensus_label(statuses)


__all__ = [
    "CHECKER_COLUMNS",
    "CHECKER_PRESETS",
    "LABELS",
    "NONVOTE_STATUSES",
    "UnknownCheckerStatusError",
    "classify_checker_row",
    "consensus_label",
    "resolve_checker_preset",
    "resolve_checker_view",
]
