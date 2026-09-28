"""Migration errors for retired third-party checker execution interfaces.

This module contains no checker algorithm, model, data table or dependency
loader. Precomputed evidence is read by the dataset and classification APIs.
"""


class CheckerExecutionUnavailableError(RuntimeError):
    """This distribution provides external checker results, not execution."""


def results_only(checker):
    raise CheckerExecutionUnavailableError(
        "{} execution is not distributed with CoRE-MOF-Tools. "
        "Load existing release results with CoREMOFDataset.from_release(...), "
        "then use dataset.classify(checkers=(...)) to select a CR/NCR policy. "
        "Detailed findings are in metadata/checker_findings.csv or .jsonl. "
        "For new calculations, use the original checker's software separately. "
        "No calculation was run and no result was written.".format(checker)
    )


def retired_cli(checker, argv=None):
    import argparse
    parser = argparse.ArgumentParser(
        description="{} execution has been retired from this distribution. "
        "Use the precomputed release results. See README.md.".format(checker)
    )
    parser.parse_args(argv)
    parser.error("Results-only distribution: no checker calculation is available")
