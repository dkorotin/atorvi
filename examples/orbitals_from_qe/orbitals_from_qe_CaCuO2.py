import atorvi


outfile = atorvi.OrbitalFile("CaCuO2_orbitals.xsf")

outfile.orbitals_from_qe(
    qe_outfile="CaCuO2_scf.out",
    mode="hole",
    eigenstate="dominant",
)

outfile.write_data(squared=False)
