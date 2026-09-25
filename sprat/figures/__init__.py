"""The seven figures of the paper, rebuilt from the records and the analysis outputs.

Modules: ``style`` (house style: type, palette, sizes, panel letters, saving), ``inputs`` (access
to the records and to ``analysis.json`` / ``spectra.json``), ``schematic`` (figure 1),
``campaign`` (figures 2, 3, 5, 6 and 7), ``mechanism`` (figure 4) and ``build``
(``build_all``, the ``sprat figures`` command).  Needs Matplotlib (``pip install sprat[figures]``).
"""

__all__ = ["build", "campaign", "inputs", "mechanism", "schematic", "style"]
