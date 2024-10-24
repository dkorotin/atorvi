"""
atorvi - ATomic ORbitals VIsualization 
a library for visualizing individual atomic orbitals and their various linear combinations.
The result is a file in .xsf format, which can be opened and visualized using software like XCrysDen or VESTA.
"""

__author__ = "Dmitry Korotin"
__author_email__ = "dmitry@korotin.name"
__version__ = "0.1.1"
__license__ = "MIT"

from .atomic_orbitals import (
    radial_part,
    get_orbital,
    supported_orbitals,
    p_orbitals,
    d_orbitals,
    f_orbitals,
    get_atomic_number
)

from .atorvi import (
    OrbitalFile,
    main
)

__all__ = [
    'OrbitalFile',
    'radial_part',
    'get_orbital',
    'supported_orbitals',
    'p_orbitals',
    'd_orbitals',
    'f_orbitals',
    'get_atomic_number'
]
