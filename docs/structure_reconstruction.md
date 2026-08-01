# Candidate structure reconstruction

`reconstruct_backbone(sequence, phi, psi; omega=180.0)` constructs only the N,
Cα, C, and O atoms, in ångström, using idealized standard peptide lengths and
angles. Explicit torsions are degrees. The output reports Cα radius of gyration,
an 8 Å nonlocal Cα contact map, coarse torsion-region assignments, and obvious
nonlocal clashes. These checks are diagnostics, not structure validation.

The legacy `Trajectory.dihedral` array is sites × saved frames. Its values are
dimensionless reduced coordinates in a double-well model; they are not φ or ψ.
`approximate_backbone_torsions` deliberately returns a `TorsionMapping` with a
45-degree uncertainty and warning. For visualization it maps
`φ=-60+15θ`, `ψ=-45-15θ`, and `ω=180` degrees. This arbitrary α-like reference
does not add structural information to the reduced simulation.

`write_pdb` and `write_mmcif` include model-derived warnings. Exported
coordinates should be refined and tested with an independently selected
atomistic workflow; they are not validated native structures.
