import atorvi

# d-orbitals of Niobium (Z=41)
outfile = atorvi.OrbitalFile("d_orbitals.xsf")

outfile.add_orbital("d_{3z^2-r^2}", position=[0.0, 0, 0.0], znumber=41)
outfile.add_orbital("d_{xz}", position=[5, 0, 0.0], znumber=41)
outfile.add_orbital("d_{yz}", position=[10, 0, 0.0], znumber=41)
outfile.add_orbital("d_{x^2-y^2}", position=[15, 0, 0.0], znumber=41)
outfile.add_orbital("d_{xy}", position=[20, 0, 0.0], znumber=41)

outfile.write_data()

# Create squared atomic orbitals

outfile = atorvi.OrbitalFile("d_orbitals_squared.xsf")

outfile.add_orbital("d_{3z^2-r^2}", position=[0.0, 0, 0.0], znumber=41)
outfile.add_orbital("d_{xz}", position=[5, 0, 0.0], znumber=41)
outfile.add_orbital("d_{yz}", position=[10, 0, 0.0], znumber=41)
outfile.add_orbital("d_{x^2-y^2}", position=[15, 0, 0.0], znumber=41)
outfile.add_orbital("d_{xy}", position=[20, 0, 0.0], znumber=41)

outfile.write_data(squared=True)
