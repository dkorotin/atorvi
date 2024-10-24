__author__ = "Dmitry Korotin"
__author_email__ = "dmitry@korotin.name"

import numpy as np
from .atomic_orbitals import *

class OrbitalFile:
    """
    A class for creating and managing files containing atomic orbital data.

    This class provides methods for adding atoms and orbitals, setting up the
    computational box, and writing the orbital data to a file.

    """

    def __init__(self, filename, grid_step=0.05) -> None:
        self.file = open(filename, "w")
        self.box_origin = [-2, -2, -2]
        self.box_size = [4, 4, 4]
        self.orbitals = []
        self.grid_step = grid_step
        self.atoms = []
        self.is_crystal = False

    def __del__(self):
        if not self.file.closed:
            self.file.close()

    def _calculate_orbital_grid(self, orbital, position=[0.0, 0.0, 0.0], znumber=8):

        origin = self.box_origin - position
        x = np.linspace(origin[0], origin[0] + self.box_size[0], self.grid[0])
        y = np.linspace(origin[1], origin[1] + self.box_size[1], self.grid[1])
        z = np.linspace(origin[2], origin[2] + self.box_size[2], self.grid[2])

        xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")

        sph_coords = cart2sph(xx, yy, zz)

        orbital_grid = get_orbital(orbital, sph_coords, znumber)

        return orbital_grid

    def _write_orbitals(self, squared):

        self._generate_box()

        for orb in self.orbitals:
            self.datagrid += orb[3] * self._calculate_orbital_grid(
                orb[0], orb[1], orb[2]
            )

        if squared:
            self.datagrid = self.datagrid**2

        self.file.write("\nBEGIN_BLOCK_DATAGRID_3D\nmy_datagrid\n")
        self.file.write("BEGIN_DATAGRID_3D\n")
        self.file.write(f"{self.grid[0]} {self.grid[1]} {self.grid[2]}\n")
        self.file.write(
            f"{self.box_origin[0]} {self.box_origin[1]} {self.box_origin[2]}\n"
        )
        self.file.write(f"{self.box_size[0]} 0 0\n")
        self.file.write(f"0 {self.box_size[1]} 0\n")
        self.file.write(f"0 0 {self.box_size[2]}\n")

        outdata = self.datagrid.T
        for vec in outdata:
            np.savetxt(self.file, vec, fmt="%12.9f", footer=" ", comments="")

        self.file.write("END_DATAGRID_3D\n")
        self.file.write("END_BLOCK_DATAGRID_3D\n")

    def _write_atoms(self):

        if self.is_crystal:
            self.file.write("CRYSTAL\nPRIMVEC\n")
            self.file.writelines(
                f"{row[0]:9.6f} {row[1]:9.6f} {row[2]:9.6f}\n"
                for row in self.lattice.matrix
            )

            self.file.write("CONVVEC\n")
            self.file.writelines(
                f"{row[0]:9.6f} {row[1]:9.6f} {row[2]:9.6f}\n"
                for row in self.lattice.matrix
            )

            self.file.write("PRIMCOORD\n")
            self.file.write(f"{len(self.atoms)} 1\n")
        else:

            self.file.write("ATOMS\n")

        for at in self.atoms:
            atomic_number = get_atomic_number(at[0])
            x, y, z = at[1]
            self.file.write(f"{atomic_number:3n}  {x:9.6f}  {y:9.6f}  {z:9.6f}\n")

    def add_atoms(self, atoms):
        """
        Add atoms to the system.

        This method adds atoms to the system if it's not a periodic system.

        Parameters:
        ----------
        atoms : list
            A list of atoms to add, where each atom is represented by a tuple
            containing the element symbol and its coordinates.

        Raises:
        ------
        ValueError
            If the system is periodic (crystal structure).
        """

        if not self.is_crystal:
            for a in atoms:
                self.atoms.append(a)
        else:
            raise ValueError("It is a periodic system. Do not add atoms manually")

    def write_data(self, squared=False):
        """
        Write the orbital and atom data to the file.

        This method writes the atom positions and orbital data to the file
        and closes it.
        """

        if len(self.atoms) > 0:
            self._write_atoms()

        if len(self.orbitals) > 0:
            self._write_orbitals( squared )

        print(f"\nFile {self.file.name} successfully written")

        self.file.close()

    def _generate_box(self):
        min_orbital_size = 5.0

        orbitals_positions = [orb[1] for orb in self.orbitals]

        box_center, box_size = get_bbox_center_and_size(orbitals_positions)

        box_size += [min_orbital_size, min_orbital_size, min_orbital_size]
        self.box_size = np.array(box_size)

        self.box_origin = box_center - self.box_size / 2.0

        self.grid = np.array(self.box_size // self.grid_step, dtype=int)
        self.datagrid = np.zeros(self.grid)

    def add_orbital(self, orbital, position=[0.0, 0.0, 0.0], znumber=8, coeff=1.0):
        """
        Add an orbital to the system.

        This method adds an orbital with specified parameters to the system.

        Parameters:
        ----------
        orbital : str
            The type of orbital (e.g., "s", "p_x", "d_{xy}", etc.).
        position : list, optional
            The position of the orbital (default is [0.0, 0.0, 0]).
        znumber : int, optional
            The atomic number of the element (default is 8).
        coeff : float, optional
            The coefficient of the orbital (default is 1.0).
        """
        self.orbitals.append([orbital, position, znumber, coeff])

    def add_orbital_at_atom(self, orbital, atom_index, coeff=1.0):
        """
        Add an orbital at a specific atom's position.

        This method adds an orbital at the position of a specified atom.

        Parameters:
        ----------
        orbital : str
            The type of orbital (e.g., "s", "p_x", "d_{xy}", etc.).
        atom_index : int
            The index of the atom in the system.
        coeff : float, optional
            The coefficient of the orbital (default is 1.0).

        Raises:
        ------
        KeyError
            If the atom_index is out of range.
        """
        if atom_index in range(0, len(self.atoms)):
            position = self.atoms[atom_index][1]
            znumber = get_atomic_number(self.atoms[atom_index][0])
            self.orbitals.append([orbital, position, znumber, coeff])
        else:
            raise KeyError("Wrong atom_index")

    def crystal_from_file(self, filename):
        """
        Create a crystal structure from a file.

        This method reads a crystal structure from a file and sets up the system.

        Parameters:
        ----------
        filename : str
            The path to the file containing the crystal structure.

        Returns:
        -------
        pymatgen.core.structure.IStructure
            The created crystal structure.
        """
        try:
            from pymatgen.core.structure import IStructure
        except ImportError:
            raise ImportError(
                "\nOptional dependency 'pymatgen' required for crystal_from_file is not installed. \n\
                You can install it by running: pip install pymatgen"
            )

        structure = IStructure.from_file(filename)
        self.crystal_from_pymatgen(structure)

        return structure

    def crystal_from_pymatgen(self, structure):
        """
        Set up the crystal structure from a pymatgen Structure object.

        This method sets up the crystal structure using a pymatgen Structure object.

        Parameters:
        ----------
        structure : pymatgen.core.structure.IStructure
            The pymatgen Structure object representing the crystal.
        """

        crystal = structure
        self.is_crystal = True

        self.lattice = crystal.lattice

        for atom in crystal.sites:
            self.atoms.append((atom.specie.symbol, atom.coords))

    def add_orbital_at_element(self, orbital, element, coeff=1.0):
        """
        Add an orbital to all atoms of a specific element.

        This method adds the specified orbital to all atoms of the given element.

        Parameters:
        ----------
        orbital : str
            The type of orbital (e.g., "s", "p_x", "d_{xy}", etc.).
        element : str
            The symbol of the element.
        coeff : float, optional
            The coefficient of the orbital (default is 1.0).
        """
        for iatom, atom in enumerate(self.atoms):
            if element == atom[0]:
                self.add_orbital_at_atom(orbital, iatom, coeff)

def get_bbox_center_and_size(points):
    """
    Calculate the center and size of a bounding box for a set of points.

    This function takes a list of points and calculates the center and size of
    the smallest bounding box that contains all the points.

    Parameters:
    ----------
    points : list of tuples or numpy.ndarray
        A list of points, where each point is represented by a tuple or array
        of coordinates.

    Returns:
    -------
    tuple
        A tuple containing two numpy arrays:
        - center: The coordinates of the center of the bounding box.
        - bounding_box_size: The size of the bounding box in each dimension.
    """
    points_array = np.array(points)

    center = np.mean(points_array, axis=0)
    distances = points_array - center
    max_distances = np.max(np.abs(distances), axis=0)
    bounding_box_size = 2 * max_distances

    return center, bounding_box_size

def main():
    """
    Interactive mode for atorvi orbital visualization.
    Allows users to specify orbitals and their positions interactively.
    """
    from importlib.metadata import version

    print("\nATORVI - ATomic ORbital VIsualization")
    print(f"Version: {version('atorvi')}")
    print(f"Author: {__author__}")
    print(f"Contact: {__author_email__}\n")

    # Get output filename
    filename = input("Enter output filename (e.g. 'orbitals.xsf'): ").strip()
    if not filename:
        filename = "orbitals.xsf"
    
    # Initialize orbital file
    orbital_file = OrbitalFile(filename)

    print("\nSupported orbitals are:")
    print("s, p_z, p_x, p_y")
    print(f"{' '.join(map(str, d_orbitals))}")
    print(f"{' '.join(map(str, f_orbitals))}")
    
    print("\nEnter orbital data in format: Element orbital [x y z]")
    print("\nExample: 'O p_z 0 0 1' or 'Fe d_{xy}'")
    print("Type 'done' to finish and generate visualization")
    
    while True:
        try:
            # Get user input
            user_input = input("\nEnter orbital data: ").strip()
            
            # Check for exit condition
            if user_input.lower() == 'done':
                break
            
            # Parse input
            parts = user_input.split()
            if len(parts) < 2:
                print("Error: Please provide at least element and orbital type")
                continue
            
            element = parts[0]
            orbital_type = parts[1]

            if orbital_type not in supported_orbitals:
                print(f"Error: Unsupported orbital type '{orbital_type}'")
                continue
            
            # Parse position if provided
            if len(parts) == 5:
                try:
                    position = [float(parts[2]), float(parts[3]), float(parts[4])]
                except ValueError:
                    print("Error: Position coordinates must be numbers")
                    continue
            else:
                position = [0.0, 0.0, 0.0]
            
            # Get atomic number
            try:
                z_number = get_atomic_number(element)
            except KeyError:
                print(f"Error: Unknown element '{element}'")
                continue
            
            # Add orbital
            orbital_file.add_orbital(
                orbital=orbital_type,
                position=position,
                znumber=z_number,
                coeff=1.0
            )
            print(f"Added {orbital_type} orbital for {element} at position {position}")
            
        except KeyboardInterrupt:
            print("\nOperation cancelled by user")
            return
        except Exception as e:
            print(f"Error: {str(e)}")
            continue
    
    try:
        print("\nGenerating visualization...")
        orbital_file.write_data(squared=False)
        print("\nDone! You can now visualize the orbitals using XCrySDen, VESTA or similar software.")
    except Exception as e:
        print(f"Error while generating visualization: {str(e)}")


if __name__ == "__main__":
    main()