Recorded MOFClassifier ensemble
========================================

``CoREMOF.release_mofclassifier.calculate_release_mofclassifier`` runs the
recorded MOFClassifier 0.1.1 core ensemble on an isolated CIF copy. It uses
all 100 saved models, positive class 1, CPU inference with one thread and
the recorded random seed. Its CL score is ``math.fsum(bag_scores) / 100``.
A completed score at or above 0.6 gives PASS, otherwise FAIL. CL is a
learned crystal-likeness score, not a diagnosis of a particular chemical error.
Missing models, invalid scores, parser failures and timeouts never give FAIL.

.. code-block:: python

   from CoREMOF.release_mofclassifier import calculate_release_mofclassifier

   result = calculate_release_mofclassifier(
       "structure.cif",
       structure_id="2016[Co][sqc27]3[FSR]3",
       output_dir="new_classifier_replay",
       python="/path/to/recorded/env/bin/python",
       model_root="/path/to/recorded/env/lib/python3.9/site-packages/MOFClassifier",
       timeout_seconds=300,
       memory_limit_mb=8192,
   )
   print(result["execution_status"], result["mean_score"])

Use ``python -m CoREMOF.release_mofclassifier --help`` or the runnable
``examples/replay_release_mofclassifier.py`` example for the CLI.

Requirements and output
-----------------------

This explicit profile checks the recorded Python executable, package versions,
worker/configuration bytes, atom embeddings and every model checkpoint by hash.
It requires the already-installed Python 3.9.23 environment with Torch
2.7.0+cu118, NumPy 1.26.4, ASE 3.23.0 and pymatgen 2024.8.9. No software,
weights or data are downloaded and no model is trained. Install the upstream
assets separately before invoking it. The base package does not distribute
those model weights or import Torch on loading this API.

The upstream fallback can reformat its private CIF. This is recorded by
``private_copy_rewritten`` and the model-input hash. The original input is
checked before and after inference and is never passed to that parser.
The output directory must be new. It contains the summary, complete raw model
record, all 100 scores when available, stdout/stderr and runtime/asset receipts.
Resource failure retains unavailable results and diagnostics, not fabricated
scores or votes. Existing release metadata is never updated by this call.

Numerical reproducibility
-------------------------

The historical runtime record does not establish the byte identity of every
numerical library or the original CPU. Two COD replay controls reproduce their
PASS/FAIL decisions and source/model bindings but show small differences in
individual float32 probabilities and their means. Exact equality to those
historical scores is therefore **not established**. The original results remain
authoritative and are not replaced by replay values. The receipt explicitly
sets ``historical_full_runtime_byte_identity_proven=false``.

PyTorch documents that mathematically equivalent floating-point operations
need not be bitwise equal across platforms or implementations, even with
controlled randomness. See its `numerical accuracy notes
<https://docs.pytorch.org/docs/2.7/notes/numerical_accuracy.html>`_ and
`reproducibility notes <https://docs.pytorch.org/docs/2.7/notes/randomness.html>`_.
This explains a general limitation, not a proven cause of these specific
historical differences. No tolerance, rounding or replacement value is applied
to the saved scores by this API.

Ordinary batch prediction
-------------------------

``curate.run_mofclassifier`` remains the generic upstream batch wrapper and
retains its model choices, batch-size default, return mapping and upstream
mean. It now supplies only private CIF copies, validates complete finite
100-model outputs, imports MOFClassifier lazily and requires assets to exist
before import. An existing result requires ``overwrite=True``. The generic
wrapper is not the recorded CPU protocol and must not be presented as an exact
release-score replay.

The model and checker method should be cited as described in the upstream
`MOFClassifier project <https://github.com/Chung-Research-Group/MOFClassifier>`_
and its `associated paper <https://doi.org/10.1021/jacs.5c10126>`_.
