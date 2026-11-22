# Repository review and compatibility assessment

The public trajectory is `Trajectory`: time plus exciton density, displacement,
and reduced dihedral matrices arranged as residues × saved times, followed by
sequence and `DavydovParams`. Simulation time and spectral waiting time are in
ps; structure reconstruction adds explicitly documented degrees and ångström.
The reduced model itself is described as internally scaled ps/Å/amu with a
dimensionless reduced dihedral coordinate.

Public APIs are exported by `src/DavydovFoldon.jl`. Heavy integrations are
function-local imports: HDF5 serialization, SQLite persistence, HTTP/JSON REST,
and optional CairoMakie plotting. Existing serialization is HDF5 format 1.0;
the REST interface returns JSON and the database uses SQLite. OpenMM remains an
external optional Python dependency and its job description is JSON.

Changing `Trajectory` would break positional construction and saved-data
consumers. Consequently explicit torsions are an auxiliary result and all
changes are backward-compatible additions. Draft experimental-feasibility work
is preserved: clocks remain separate, units remain explicit, and simulated
features retain cautious interpretations.
