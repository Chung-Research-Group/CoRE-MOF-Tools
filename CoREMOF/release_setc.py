"""Compatibility notice for the retired external SETC-GAT replay interface."""
from ._checker_execution import (
    CheckerExecutionUnavailableError as ReleaseSETCError,
    results_only,
    retired_cli,
)

__all__ = ["ReleaseSETCError", "calculate_release_setc", "main"]


def calculate_release_setc(*args, **kwargs):
    return results_only("SETC-GAT")


def main(argv=None):
    return retired_cli("SETC-GAT", argv)


if __name__ == "__main__":
    raise SystemExit(main())
