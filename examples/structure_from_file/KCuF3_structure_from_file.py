import atorvi

outfile = atorvi.OrbitalFile('KCuF3_orbitals.xsf')

structure = outfile.crystal_from_file('./KCuF3_structure.xsf')

outfile.add_orbital_at_atom('d_{xz}', 0, coeff=0.5)
outfile.add_orbital_at_atom('d_{yz}', 0, coeff=-0.5)

outfile.add_orbital_at_atom('d_{xz}', 1, coeff=0.5)
outfile.add_orbital_at_atom('d_{yz}', 1, coeff=0.5)

outfile.add_orbital_at_atom('d_{xz}', 2, coeff=0.5)
outfile.add_orbital_at_atom('d_{yz}', 2, coeff=0.5)

outfile.add_orbital_at_atom('d_{xz}', 3, coeff=0.5)
outfile.add_orbital_at_atom('d_{yz}', 3, coeff=-0.5)

outfile.add_orbital('d_{xz}',position=(structure.sites[0].coords+(structure.lattice.matrix[0]+structure.lattice.matrix[1])/2), znumber=29, coeff=0.5)
outfile.add_orbital('d_{yz}',position=(structure.sites[0].coords+(structure.lattice.matrix[0]+structure.lattice.matrix[1])/2), znumber=29, coeff=0.5)

outfile.add_orbital('d_{xz}',position=(structure.sites[1].coords+(structure.lattice.matrix[0]+structure.lattice.matrix[1])/2), znumber=29, coeff=0.5)
outfile.add_orbital('d_{yz}',position=(structure.sites[1].coords+(structure.lattice.matrix[0]+structure.lattice.matrix[1])/2), znumber=29, coeff=-0.5)

outfile.add_orbital('d_{xz}',position=(structure.sites[2].coords+(structure.lattice.matrix[0]+structure.lattice.matrix[1])/2), znumber=29, coeff=0.5)
outfile.add_orbital('d_{yz}',position=(structure.sites[2].coords+(structure.lattice.matrix[0]+structure.lattice.matrix[1])/2), znumber=29, coeff=-0.5)

outfile.add_orbital('d_{xz}',position=(structure.sites[3].coords+(structure.lattice.matrix[0]+structure.lattice.matrix[1])/2), znumber=29, coeff=0.5)
outfile.add_orbital('d_{yz}',position=(structure.sites[3].coords+(structure.lattice.matrix[0]+structure.lattice.matrix[1])/2), znumber=29, coeff=0.5)

outfile.write_data()
