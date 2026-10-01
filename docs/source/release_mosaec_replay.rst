MOSAEC results and retired execution
============================================

The MOSAEC implementation, its reference tables and execution worker are not
included in this distribution. ``curate.run_MOSAEC``, ``mosaec.run`` and
``release_mosaec.calculate_release_mosaec`` now raise migration errors without
running calculations. The module CLI provides a help notice only.

Existing MOSAEC votes, findings and recorded diagnostics are still available
from the release metadata. Read them through the workflow in
:doc:`release_checkers_replay`. Reading results requires no CCDC installation.
For new calculations use `MOSAEC <https://github.com/uowoolab/MOSAEC>`_ separately
under the relevant software and data terms. Do not treat removing its code as
blanket redistribution permission for its inputs or outputs.
