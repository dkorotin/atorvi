import atorvi


outfile = atorvi.OrbitalFile("LaMnO3_orbitals.xsf")

outfile.orbitals_from_qe(
    qe_outfile="LaMnO3_scf.out",
    atoms = [5,6,7,8],
    spin="up",
    mode="hole",
    eigenstate="dominant",
)

outfile.write_data(squared=False)
