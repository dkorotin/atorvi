import atorvi

outfile = atorvi.OrbitalFile("sp3_orbital.xsf")

outfile.add_atoms([("C", [0, 0, 0])])

outfile.add_orbital_at_atom("s", 0, coeff=0.5)
outfile.add_orbital_at_atom("p_x", 0, coeff=0.5)
outfile.add_orbital_at_atom("p_y", 0, coeff=0.5)
outfile.add_orbital_at_atom("p_z", 0, coeff=0.5)

outfile.write_data()
