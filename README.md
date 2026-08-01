# DavydovFoldon v1.0

DavydovFoldon is a Julia 1.10+ research package for exploring a stochastic,
dihedrally coupled Davydov–Scott model and its predicted 2D-IR response. It
contains simulation, spectral analysis, experiment-planning, SQLite, REST, and
command-line interfaces. The hypothesis is speculative; parameter estimates
and spectra are predictions rather than validated biological conclusions.

## Install

```bash
git clone <repository-url> DavydovFoldon
cd DavydovFoldon
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Julia creates `Manifest.toml` during `Pkg.instantiate`; commit that generated
file when an application-grade, platform-pinned environment is required.

## Quick start

```julia
using DavydovFoldon
p = DavydovParams(output_stride=10, seed=42)
traj, final_theta, final_rho = run_soliton("MKFLILFN", p; tmax=20, dt=0.01, T=310)
save_trajectory("results.h5", traj)
s = compute_2DIR_spectrum(traj, SpectrumParams())
beat = extract_beat_frequency(s)
fig = animate_soliton(traj) # requires: import Pkg; Pkg.add("CairoMakie")
```

HDF5 output contains `time`, `exciton_density`, `displacement`, `dihedral`, and
`contact_map`. Arrays use Julia/HDF5 column-major layout: sites × saved steps.

## CLI and service

```bash
julia --project=. bin/davydovfoldon predict MKFLILFN results.h5
julia --project=. -e 'using DavydovFoldon; serve(port=8080)'
curl -X POST localhost:8080/predict -d '{"sequence":"MKFLILFN"}'
```

The REST server is intended for trusted local use; deploy behind authentication,
request limits, and a production proxy. The notebooks provide Python-driven
examples that invoke Julia through PyJulia and plot results with matplotlib.

## Numerical model

The engine applies velocity Verlet to lattice coordinates, a norm-preserving
Crank–Nicolson/Cayley step to the exciton Hamiltonian, and Euler–Maruyama to
overdamped torsions. Boundary sites use reflecting (zero-gradient) neighbors.
Randomness is reproducible through `DavydovParams.seed`.

## License

Provided for research and educational use. Add the license appropriate to your
distribution before redistribution.
