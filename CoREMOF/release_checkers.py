"""Compatibility notice for retired Chen-Manz and MOFChecker execution."""
from ._checker_execution import (
    CheckerExecutionUnavailableError as ReleaseCheckersError,
    results_only,
    retired_cli,
)

__all__ = ["ReleaseCheckersError", "calculate_release_checkers", "main"]


def calculate_release_checkers(*args, **kwargs):
    return results_only("Chen-Manz and MOFChecker")


def main(argv=None):
    return retired_cli("Chen-Manz and MOFChecker", argv)


if __name__ == "__main__":
    raise SystemExit(main())
