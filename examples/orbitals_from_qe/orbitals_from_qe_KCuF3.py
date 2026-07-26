import atorvi


outfile = atorvi.OrbitalFile("KCuF3_orbitals.xsf")

outfile.orbitals_from_qe(
    qe_outfile="KCuF3_scf.out",
    atoms = [1,2,3,4],
    mode="hole",
)

outfile.write_data(squared=False)
