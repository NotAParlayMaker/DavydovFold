# Migration note

This change is additive. `Trajectory`, `run_soliton`, saved HDF5 datasets,
spectral APIs, CLI `predict`, REST behavior, and database schemas are unchanged.
New explicit `BackboneTorsions` and `TorsionMapping` objects avoid changing or
silently reinterpreting `Trajectory.dihedral`. The CLI adds `reconstruct` and
`prepare-openmm`; existing invocations remain valid.

The new APIs require angles in degrees and coordinates in ångström, while
OpenMM configuration fields carry `_K`, `_bar`, `_molar`, `_nm`, `_fs`, or `_ps`
suffixes. No implicit unit conversion is performed.
