"""SPRAT: Side-coupled Photonic-crystal Resonator Analysis Toolkit.

Two-dimensional FDTD (Meep), MPB and plane-wave-expansion workflow for a point-defect
microcavity side-coupled to a W1 waveguide in a square lattice of dielectric rods immersed
in a liquid analyte.  The structure is read from a ``.phc`` text file, the run parameters
from a ``.par`` text file; runs are executed in parallel on a workstation; every run writes
one self-describing JSON record.  The analysis layer and the figure scripts reproduce the
tables and figures of the accompanying paper from the records alone.

Units: the lattice constant ``a`` and the vacuum speed of light are 1; frequencies are
``f = a / lambda``.
"""

__version__ = "1.2.1"
__all__ = ["__version__", "SCHEMA_VERSION"]

SCHEMA_VERSION = "sprat-record-1.0"
