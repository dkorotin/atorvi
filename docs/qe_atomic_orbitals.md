# Atomic orbitals in Quantum ESPRESSO

This note documents the convention used by Quantum ESPRESSO 7.5 for atomic
wavefunctions and how it differs from the orbitals implemented in `atorvi`.

The relevant QE source files are:

- `upflib/ylmr2.f90`: real spherical harmonics and their ordering.
- `Modules/atomic_wfc_mod.f90`: construction of atomic wavefunctions in
  reciprocal space.
- `upflib/atwfc_mod.f90`: radial Fourier transform of pseudo-atomic
  wavefunctions from UPF files.
- `PW/src/write_ns.f90`: printing of Hubbard occupation eigenvectors.

## QE atomic wavefunction

For a collinear calculation, QE constructs an atomic wavefunction in reciprocal
space as

```text
psi_{n l m}^{I}(k + G) =
    i^l S_I(k + G) Y_{l m}^{QE}(k + G) chi_{n l}^{I}(|k + G|)
```

where

```text
S_I(k + G) = exp[-i (k + G) . tau_I]
```

is the structure factor of atom `I`, `Y_{l m}^{QE}` is a real spherical
harmonic from `ylmr2`, and `chi_{n l}^{I}(q)` is a radial Fourier transform of
the pseudo-atomic wavefunction stored in the UPF pseudopotential.

In `Modules/atomic_wfc_mod.f90`, this is the code path:

```fortran
CALL ylmr2( (lmax_wfc+1)**2, npw, gk, qg, ylm )
CALL interp_atwfc ( npw, qg, nwfcm, chiq )
lphase = (0.d0,1.d0)**l
wfcatom(ig,1,n_starting_wfc) =
    lphase * sk(ig) * CMPLX(ylm(ig,lm) * chiq(ig,nb,nt), KIND=DP)
```

The comment in QE says that the `i^l` factor must be present so that the
`k = 0` wavefunctions are real in real space.

## QE radial part

QE does not use hydrogen-like analytic radial functions for these atomic
wavefunctions. It uses the pseudo-atomic radial functions from the UPF file:

```text
chi_{n l}(q) =
    4 pi / sqrt(Omega) * integral chi_{n l}(r) j_l(q r) r dr
```

This comes from `upflib/atwfc_mod.f90`:

```fortran
CALL sph_bes( msh(nt), rgrid(nt)%r, q, l, aux )
vchi(ir) = upf(nt)%chi(ir,nb) * aux(ir) * rgrid(nt)%r(ir)
CALL simpson( msh(nt), vchi, rgrid(nt)%rab, vqint )
tab_atwfc( iq, nb, nt ) = vqint * pref
```

In contrast, `atorvi` currently uses an analytic hydrogen-like radial function
with an effective nuclear charge:

```text
R_{n l}(r) ~ exp[-rho / 2] rho^l L_{n-l-1}^{2l+1}(rho)
rho = 2 Z_eff r / (n a_0)
```

So matching QE exactly would require reading the UPF radial wavefunctions and
using QE's radial transform or reconstructing the corresponding real-space
pseudo-atomic orbital.

## QE real spherical harmonics

QE's `ylmr2` generates real spherical harmonics in this order:

```text
l = 0: m = 0
l = 1: m = 0, +1, -1
l = 2: m = 0, +1, -1, +2, -2
l = 3: m = 0, +1, -1, +2, -2, +3, -3
```

The angular functions are generated from associated Legendre polynomials
`P_l^m(cos theta)`, then multiplied by `cos(m phi)` for `m > 0` and
`sin(m phi)` for `m < 0`.

The corresponding orbital labels are:

```text
p: p_z, p_x, p_y

d: d_{3z^2-r^2}, d_{xz}, d_{yz}, d_{x^2-y^2}, d_{xy}

f: f_{z^3}, f_{xz^2}, f_{yz^2}, f_{z(x^2-y^2)}, f_{xyz},
   f_{x(x^2-3y^2)}, f_{y(3x^2-y^2)}
```

This ordering matters because `PW/src/write_ns.f90` prints the occupation
matrix eigenvectors as unnamed rows. Those rows are just `m1 = 1..2l+1`, in
the same order as the atomic wavefunctions.

## Difference from `atorvi` angular functions

`atorvi` builds real orbitals from SciPy complex spherical harmonics:

```text
p_z = Y_1^0
p_x = (Y_1^{-1} - Y_1^{+1}) / sqrt(2)
p_y = i (Y_1^{-1} + Y_1^{+1}) / sqrt(2)

d_{3z^2-r^2} = Y_2^0
d_{xz}       = (Y_2^{-1} - Y_2^{+1}) / sqrt(2)
d_{yz}       = i (Y_2^{-1} + Y_2^{+1}) / sqrt(2)
d_{x^2-y^2}  = (Y_2^{-2} + Y_2^{+2}) / sqrt(2)
d_{xy}       = i (Y_2^{-2} - Y_2^{+2}) / sqrt(2)
```

With SciPy's Condon-Shortley phase convention, direct numerical comparison of
`ylmr2` against these `atorvi` definitions gives this angular sign map:

```text
s:  +1

p_z: +1
p_x: -1
p_y: -1

d_{3z^2-r^2}: +1
d_{xz}:       -1
d_{yz}:       -1
d_{x^2-y^2}:  +1
d_{xy}:       +1

f_{z^3}:          +1
f_{xz^2}:         -1
f_{yz^2}:         -1
f_{z(x^2-y^2)}:   +1
f_{xyz}:          +1
f_{x(x^2-3y^2)}:  -1
f_{y(3x^2-y^2)}:  -1
```

This is the origin of `QE_ORBITAL_FIX` in `atorvi/qe_output.py`.

## Practical consequences for Hubbard occupation eigenvectors

The `eigenvectors (columns)` printed by QE are coefficients in QE's real
spherical-harmonic basis. To visualize those vectors with `atorvi` orbitals, we
must:

1. Use QE's orbital order.
2. Multiply each coefficient by the sign map above.
3. Weight the vector by either `occupation` or `1 - occupation`, depending on
   whether electron or hole density is requested.

Changing the sign of one whole eigenvector only swaps isosurface colors and
does not change the orbital shape. Incorrect relative signs inside one shell
can rotate or mix orbitals. Incorrect ordering, especially swapping
`d_{x^2-y^2}` and `d_{xy}`, changes the visible orbital symmetry.
