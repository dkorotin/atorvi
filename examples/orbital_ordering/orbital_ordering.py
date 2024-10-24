import atorvi

outfile = atorvi.OrbitalFile('orbital_ordering.xsf')

atoms = [
    ('Cu', [0.0,0.0,0.0]),
    ('Cu', [2.0,0.0,0.0]),
    ('Cu', [0.0,2.0,0.0]),
    ('Cu', [2.0,2.0,0.0]),

    ('Cu', [0.0,0.0,2.0]),
    ('Cu', [2.0,0.0,2.0]),
    ('Cu', [0.0,2.0,2.0]),
    ('Cu', [2.0,2.0,2.0]),
    ]

outfile.add_atoms(atoms)

outfile.add_orbital_at_atom('d_{xz}', 0)
outfile.add_orbital_at_atom('d_{yz}', 1)
outfile.add_orbital_at_atom('d_{yz}', 2)
outfile.add_orbital_at_atom('d_{xz}', 3)

outfile.add_orbital_at_atom('d_{xz}', 5)
outfile.add_orbital_at_atom('d_{yz}', 4)
outfile.add_orbital_at_atom('d_{yz}', 7)
outfile.add_orbital_at_atom('d_{xz}', 6)

outfile.write_data()
