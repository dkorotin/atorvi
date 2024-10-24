# atorvi - ATomic ORbitals VIsualization

![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)

**atorvi** is a Python package for visualizing individual atomic orbitals and their various linear combinations. The library generates requested atomic orbitals on a 3D mesh and exports them, along with the inputted crystal structure (molecule or periodic crystal), into a file in [XCrysDen .xsf format](http://www.xcrysden.org/doc/XSF.html).

The resulting `.xsf` file can be opened and visualized using your favorite visualization software, such as [XCrysDen](http://www.xcrysden.org/), [VESTA](https://jp-minerals.org/vesta/en/) or [VMD](https://www.ks.uiuc.edu/Research/vmd/).

## Installation

To install **atorvi**, you can use `pip`:

```bash
pip install atorvi
```

## Features

- **Atomic Orbitals Generation**: Create atomic orbitals at a specific position in space or at a designated atom/element within a crystal structure. 

- **Flexible Structure Input**: Input crystal structures either manually, atom-by-atom, or directly from standard file formats such as POSCAR, XSF, CIF, and more. The integration with the `pymatgen` package simplifies structure handling.

- **Orbital Hybridization**: Implement orbital hybridization by mixing different orbitals with custom coefficients, allowing for the exploration of complex bonding interactions.

- **Orbital Squared Moduli**: Generate squared moduli of orbitals, providing insight into their spatial distribution and probability densities.

Examples of orbitals generated with **atorvi**: could be found in [examples](./examples/) folder.

## Quick start

The full user manual is available in [docs/atorvi_manual.md](./docs/atorvi_manual.md).

**atorvi** can be used in two modes: package mode (for scripting and notebooks) and CLI interactive mode (only basic functionality is available).

### 1.Package Mode 
In your Python script or Jupyter notebook, you can generate and visualize atomic orbitals with the following example:

```python
import atorvi

# Example: visualize a d_{3z^2-r^2} orbital for an Ni atom (Z = 28)
outfile = atorvi.OrbitalFile("Ni_orbital.xsf")

outfile.add_orbital("d_{3z^2-r^2}", position=[0, 0, 0], znumber=28)

outfile.write_data()
```

Then just open the `Ni_orbital.xsf` file in your favorite visualization software.

The orbitals available for generation are:
$$s$$
$$p_z, p_x, p_y$$
$$d_{3z^2-r^2}, d_{xz}, d_{yz}, d_{xy}, d_{x^2-y^2} $$
$$f_{z^3}, f_{xz^2}, f_{yz^2}, f_{xyz}, f_{z(x^2-y^2)}, f_{x(x^2-3y^2)}, f_{y(3x^2-y^2)}$$

```python
print(atorvi.supported_orbitals)

['s', 
'p_z', 'p_x', 'p_y', 
'd_{3z^2-r^2}', 'd_{xz}', 'd_{yz}', 'd_{xy}', 'd_{x^2-y^2}', 
'f_{z^3}', 'f_{xz^2}', 'f_{yz^2}', 'f_{xyz}', 'f_{z(x^2-y^2)}', 'f_{x(x^2-3y^2)}', 'f_{y(3x^2-y^2)}']
```

### 2. CLI Interactive Mode
atorvi also offers an interactive mode via the command line interface (CLI). You can simply run the script directly in a terminal:

```bash
atorvi_cli
```

Follow the prompts to generate orbitals and export them to XCrysDen-compatible formats.

Once the file is generated, you can open it using visualization tools like XCrysDen or VESTA.

## Author

`atorvi` is developed and maintained by [Dmitry Korotin](https://www.researchgate.net/profile/Dmitry-Korotin). Contributions, suggestions, and feedback are welcome to help improve the project.

## License

`atorvi` is released under the MIT License. You are free to use, modify, and distribute the software, provided that the original copyright and permission notice are included in all copies or substantial portions of the software.

The author kindly asks that you cite this GitHub repository [github.com/dkorotin/atorvi](https://github.com/dkorotin/atorvi) and the related paper (link will be available soon) in any publications that use images or data generated with the `atorvi` package.

For more details, refer to the full [MIT License](./LICENSE).
