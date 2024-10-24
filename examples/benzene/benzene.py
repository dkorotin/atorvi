import atorvi

outfile = atorvi.OrbitalFile("benzene.xsf")

atoms = [
    ("C", [1.4 * 2.866, 1 * 1.4, 0]),
    ("C", [1.4 * 2.0, 0.5 * 1.4, 0]),
    ("C", [1.4 * 3.732, 0.5 * 1.4, 0]),
    ("C", [1.4 * 2.0, -0.5 * 1.4, 0]),
    ("C", [1.4 * 3.732, -0.5 * 1.4, 0]),
    ("C", [1.4 * 2.866, -1 * 1.4, 0]),
    ("H", [1.4 * 2.866, 1.62 * 1.4, 0]),
    ("H", [1.4 * 2.866, -1.62 * 1.4, 0]),
    ("H", [1.4 * 1.4631, 0.81 * 1.4, 0]),
    ("H", [1.4 * 1.4631, -0.81 * 1.4, 0]),
    ("H", [1.4 * 4.269, 0.81 * 1.4, 0]),
    ("H", [1.4 * 4.269, -0.81 * 1.4, 0]),
]

outfile.add_atoms(atoms)

outfile.add_orbital_at_element("p_z", "C")

outfile.write_data()