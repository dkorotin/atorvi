import atorvi

print(f"All f-orbitals are: {atorvi.f_orbitals}")
# f-orbitals of Gadalinium (Z=64)
outfile = atorvi.OrbitalFile("f_orbitals.xsf")

for i, orbital in enumerate(atorvi.f_orbitals):
    outfile.add_orbital(orbital, position=[0.0, 3.0*i, 0.0], znumber=64)

outfile.write_data()