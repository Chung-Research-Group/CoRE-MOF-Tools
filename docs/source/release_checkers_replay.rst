Precomputed checker results
============================================

The public package distributes result-reading and CR/NCR-selection interfaces,
not third-party checker implementations. The former Chen–Manz/MOFChecker
calculation helper now raises an explanatory error before any file access,
output creation, dependency loading or subprocess launch.

Read a release and choose a checker combination:

.. code-block:: python

   from CoREMOF.dataset import CoREMOFDataset

   dataset = CoREMOFDataset.from_release("path/to/coremof_release")
   view = dataset.classify(checkers=("MOFChecker", "Chen-Manz"))
   print(dict(view.label_counts()))

All selected PASS gives CR, all selected FAIL gives NCR, mixed completed
votes give AMBIGUOUS, and any NOT_AVAILABLE gives UNCHECKED. This operation
combines saved votes. It does not execute a checker or change the saved results.

Detailed findings and raw scores are in ``metadata/checker_findings.csv`` and
``.jsonl``. See ``examples/read_checker_results.py`` and the repository README.
For new calculations, use the original authors' software separately. Existing
results and frozen benchmark assignments are unchanged. Historical execution
code is retained outside the distributable package for local reproduction.
