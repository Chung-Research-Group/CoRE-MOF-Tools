"""Compatibility notice. The third-party MOSAEC engine is not bundled."""
from ._checker_execution import results_only


def run(*args, **kwargs):
    """Reject retired execution without importing software or writing files."""
    return results_only("MOSAEC")


def check(*args, **kwargs):
    return results_only("MOSAEC")
