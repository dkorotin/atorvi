import atorvi

bmo = atorvi.OrbitalFile('H2_molecule_BMO.xsf')
bmo.add_atoms([
    ('H', [0, 0, 0]),
    ('H', [0.737, 0, 0])
])

# bonding molecular orbital = 1/sqrt(2)(s1+s2)
bmo.add_orbital_at_atom('s', 0, coeff=0.707)
bmo.add_orbital_at_atom('s', 1, coeff=0.707)

bmo.write_data()


abmo = atorvi.OrbitalFile('H2_molecule_ABMO.xsf')
abmo.add_atoms([
    ('H', [0, 0, 0]),
    ('H', [0.737, 0, 0])
])

# antibonding molecular orbital = 1/sqrt(2)(s1-s2)
abmo.add_orbital_at_atom('s', 0, coeff =   0.707)
abmo.add_orbital_at_atom('s', 1, coeff = - 0.707)

abmo.write_data()
