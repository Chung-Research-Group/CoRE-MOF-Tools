SETC-GAT results and retired execution
============================================

The SETC-GAT calculation protocol and execution worker are not distributed.
``release_setc.calculate_release_setc`` now raises an explanatory migration
error without loading a model, launching a process or creating output.

Saved SETC-GAT results remain part of the checker metadata and can be selected
with ``dataset.classify(checkers=("SETC-GAT",))``. See
:doc:`release_checkers_replay` for the result-reading workflow. An unavailable
result remains NOT_AVAILABLE, not FAIL. Historical calculations and their
recorded inputs are preserved separately. New calculations require the
original software and its specified input representation, not a substitute
checker or guessed feature vector.
