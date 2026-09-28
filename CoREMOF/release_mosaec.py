"""Compatibility notice for the retired external MOSAEC replay interface."""
from ._checker_execution import (
    CheckerExecutionUnavailableError as ReleaseMOSAECError,
    results_only,
    retired_cli,
)

__all__ = ["ReleaseMOSAECError", "calculate_release_mosaec", "main"]


def calculate_release_mosaec(*args, **kwargs):
    return results_only("MOSAEC")


def main(argv=None):
    return retired_cli("MOSAEC", argv)


if __name__ == "__main__":
    raise SystemExit(main())
