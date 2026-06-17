__author__ = "Dmitry Korotin"
__author_email__ = "dmitry@korotin.name"

import re

import numpy as np

from .atomic_orbitals import d_orbitals, f_orbitals, p_orbitals


BOHR_TO_ANGSTROM = 0.529177210903

ORBITALS_BY_L = {
    0: ["s"],
    1: p_orbitals,
    2: d_orbitals,
    3: f_orbitals,
}

SPINS = {"up": 1, "down": 2}

QE_ORBITAL_FIX = {
    # Quantum ESPRESSO builds atomic wavefunctions with real spherical harmonics
    # from upflib/ylmr2.f90, in the order m=0,+1,-1,+2,-2,...
    # Modules/atomic_wfc_mod.f90 then multiplies the whole l shell by i**l in
    # reciprocal space. These signs map QE's ylmr2 angular convention onto the
    # real harmonics implemented in atomic_orbitals.py.
    "s": 1,
    "p_z": 1,
    "p_x": -1,
    "p_y": -1,
    "d_{3z^2-r^2}": 1,
    "d_{xz}": -1,
    "d_{yz}": -1,
    "d_{x^2-y^2}": 1,
    "d_{xy}": 1,
    "f_{z^3}": 1,
    "f_{xz^2}": -1,
    "f_{yz^2}": -1,
    "f_{z(x^2-y^2)}": 1,
    "f_{xyz}": 1,
    "f_{x(x^2-3y^2)}": -1,
    "f_{y(3x^2-y^2)}": -1,
}


def parse_qe_structure(lines):
    alat = _parse_last_alat(lines)
    cell = _parse_last_cell_parameters(lines, alat)
    fallback_cell = _parse_initial_crystal_axes(lines, alat)
    lattice_for_positions = cell if cell is not None else fallback_cell
    atoms = _parse_last_atomic_positions(lines, lattice_for_positions, alat)

    if cell is not None and atoms is not None:
        return cell, atoms

    if fallback_cell is not None and atoms is not None:
        return fallback_cell, atoms

    fallback_atoms = _parse_initial_positions(lines, fallback_cell, alat)

    if fallback_cell is None or fallback_atoms is None:
        raise ValueError("Could not parse final crystal structure from QE output")

    return fallback_cell, fallback_atoms


def parse_qe_hubbard_occupations(lines):
    occupations = {}
    section_start = 0
    for index, line in enumerate(lines):
        if "HUBBARD OCCUPATIONS" in line:
            section_start = index

    for index, line in enumerate(lines[section_start:], start=section_start):
        match = re.match(r"\s*-+\s*ATOM\s+(\d+)\s*-+", line)
        if not match:
            continue

        atom_number = int(match.group(1))
        cursor = index + 1
        spins = {}
        while cursor < len(lines):
            if re.match(r"\s*-+\s*ATOM\s+\d+\s*-+", lines[cursor]):
                break
            if "Number of occupied Hubbard levels" in lines[cursor]:
                break

            spin_match = re.match(r"\s*SPIN\s+(\d+)", lines[cursor])
            if not spin_match:
                cursor += 1
                continue

            spin_number = int(spin_match.group(1))
            eigenvalues, eigenvectors, cursor = _parse_qe_spin_occupation(
                lines, cursor + 1
            )
            dimension = len(eigenvalues)
            if eigenvectors.shape != (dimension, dimension):
                raise ValueError(
                    f"Malformed Hubbard eigenvectors for atom {atom_number}, "
                    f"spin {spin_number}"
                )
            spins[spin_number] = {
                "eigenvalues": np.array(eigenvalues, dtype=float),
                "eigenvectors": eigenvectors,
            }

        if spins:
            dimension = len(next(iter(spins.values()))["eigenvalues"])
            occupations[atom_number] = {
                "l": _l_from_dimension(dimension),
                "spins": spins,
            }

    return occupations


def select_qe_hubbard_atoms(occupations, atoms, l, atom_count):
    if atoms is None:
        selected = [
            atom_number
            for atom_number, atom_data in occupations.items()
            if atom_data["l"] == l
        ]
        if not selected:
            raise ValueError(f"No Hubbard atoms found for l={l}")
        return selected

    selected = list(atoms)
    for atom_number in selected:
        if not isinstance(atom_number, int):
            raise ValueError("atoms must contain QE 1-based integer atom numbers")
        if atom_number < 1 or atom_number > atom_count:
            raise ValueError(
                f"Atom {atom_number} is outside QE atom range 1..{atom_count}"
            )
        if atom_number not in occupations:
            raise ValueError(f"Atom {atom_number} has no Hubbard occupation data")
        if occupations[atom_number]["l"] != l:
            raise ValueError(
                f"Atom {atom_number} has Hubbard l={occupations[atom_number]['l']}, "
                f"not requested l={l}"
            )

    return selected


def _parse_last_alat(lines):
    alat = None
    for line in lines:
        match = re.search(r"lattice parameter \(alat\)\s*=\s*([-+0-9.Ee]+)", line)
        if match:
            alat = float(match.group(1)) * BOHR_TO_ANGSTROM
            continue

        match = re.search(r"celldm\(1\)\s*=\s*([-+0-9.Ee]+)", line)
        if match:
            alat = float(match.group(1)) * BOHR_TO_ANGSTROM

    return alat


def _parse_last_cell_parameters(lines, alat):
    cell = None
    for index, line in enumerate(lines):
        if not line.lstrip().startswith("CELL_PARAMETERS"):
            continue

        unit_text = _qe_unit_text(line)
        matrix = np.array(
            [_float_values(lines[index + offset])[:3] for offset in range(1, 4)],
            dtype=float,
        )
        cell = _convert_vectors_to_angstrom(matrix, unit_text, alat)

    return cell


def _parse_last_atomic_positions(lines, lattice_matrix, alat):
    atoms = None
    for index, line in enumerate(lines):
        if not line.lstrip().startswith("ATOMIC_POSITIONS"):
            continue

        unit_text = _qe_unit_text(line)
        parsed_atoms = []
        cursor = index + 1
        while cursor < len(lines):
            parts = lines[cursor].split()
            if len(parts) < 4:
                break
            try:
                coords = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
            except ValueError:
                break
            parsed_atoms.append((parts[0], coords))
            cursor += 1

        atoms = _convert_positions_to_angstrom(
            parsed_atoms, unit_text, lattice_matrix, alat
        )

    return atoms


def _parse_initial_crystal_axes(lines, alat):
    if alat is None:
        return None

    cell = None
    for index, line in enumerate(lines):
        if "crystal axes:" not in line:
            continue
        rows = []
        for offset in range(1, 4):
            values = _parenthesized_floats(lines[index + offset])
            if len(values) != 3:
                rows = []
                break
            rows.append(values)
        if rows:
            cell = np.array(rows, dtype=float) * alat

    return cell


def _parse_initial_positions(lines, lattice_matrix, alat):
    if alat is None:
        return None

    atoms = None
    for index, line in enumerate(lines):
        if "site n." not in line or "positions (alat units)" not in line:
            continue
        parsed_atoms = []
        cursor = index + 1
        while cursor < len(lines):
            match = re.match(
                r"\s*\d+\s+([A-Za-z][A-Za-z0-9_]*)\s+tau\(\s*\d+\)\s*=\s*"
                r"\(\s*([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*\)",
                lines[cursor],
            )
            if not match:
                break
            parsed_atoms.append(
                (
                    _clean_qe_element(match.group(1)),
                    np.array(
                        [
                            float(match.group(2)),
                            float(match.group(3)),
                            float(match.group(4)),
                        ]
                    )
                    * alat,
                )
            )
            cursor += 1
        if parsed_atoms:
            atoms = parsed_atoms

    return atoms


def _parse_qe_spin_occupation(lines, cursor):
    while cursor < len(lines) and "eigenvalues:" not in lines[cursor]:
        cursor += 1
    if cursor >= len(lines):
        raise ValueError("Malformed Hubbard occupations: missing eigenvalues")
    cursor += 1

    eigenvalues = []
    while cursor < len(lines) and "eigenvectors" not in lines[cursor]:
        eigenvalues.extend(_float_values(lines[cursor]))
        cursor += 1

    if cursor >= len(lines):
        raise ValueError("Malformed Hubbard occupations: missing eigenvectors")
    cursor += 1

    rows = []
    while cursor < len(lines) and "occupation matrix" not in lines[cursor]:
        values = _float_values(lines[cursor])
        if values:
            rows.append(values)
        cursor += 1

    return eigenvalues, np.array(rows, dtype=float), cursor


def _l_from_dimension(dimension):
    if dimension % 2 != 1:
        raise ValueError(f"Unsupported Hubbard matrix dimension {dimension}")
    l = (dimension - 1) // 2
    if l not in ORBITALS_BY_L:
        raise ValueError(f"Unsupported Hubbard matrix dimension {dimension}")
    return l


def _qe_unit_text(line):
    match = re.search(r"\(([^)]*)\)", line)
    return match.group(1).strip().lower() if match else "alat"


def _convert_vectors_to_angstrom(vectors, unit_text, alat):
    scale = _qe_length_scale(unit_text, alat)
    return np.array(vectors, dtype=float) * scale


def _convert_positions_to_angstrom(atoms, unit_text, lattice_matrix, alat):
    unit = unit_text.lower()
    converted = []
    for element, coords in atoms:
        if unit.startswith("crystal"):
            if lattice_matrix is None:
                raise ValueError("ATOMIC_POSITIONS crystal requires lattice vectors")
            position = np.dot(coords, lattice_matrix)
        else:
            position = coords * _qe_length_scale(unit, alat)
        converted.append((_clean_qe_element(element), position))
    return converted


def _qe_length_scale(unit_text, alat):
    unit = unit_text.lower()
    alat_match = re.search(r"alat\s*=\s*([-+0-9.Ee]+)", unit)
    if "angstrom" in unit:
        return 1.0
    if "bohr" in unit or "a.u." in unit:
        return BOHR_TO_ANGSTROM
    if "alat" in unit:
        if alat_match:
            return float(alat_match.group(1)) * BOHR_TO_ANGSTROM
        if alat is None:
            raise ValueError("QE output uses alat units but alat was not found")
        return alat
    raise ValueError(f"Unsupported QE length unit: {unit_text}")


def _clean_qe_element(label):
    match = re.match(r"([A-Z][a-z]?)", label)
    return match.group(1) if match else label


def _float_values(line):
    return [
        float(value)
        for value in re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?", line)
    ]


def _parenthesized_floats(line):
    matches = re.findall(r"\(([^)]*)\)", line)
    if not matches:
        return []
    return _float_values(matches[-1])
