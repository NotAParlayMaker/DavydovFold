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

## Experimental feasibility gate

Before treating a simulated signal as measurable, calculate the transient-state
population and compare it with an instrument-specific concentration threshold:

```julia
using DavydovFoldon

assumptions = FeasibilityParams(
    rnc_concentration_M = 100e-9,
    elongation_rate_s = 20.0,
    state_lifetime_s = 200e-12,
    detection_limit_M = 1e-6,
)

population = estimate_transient_population(assumptions)
assessment = assess_2dir_feasibility(assumptions)
plan = staged_validation_plan()
```

The default assumptions deliberately return a `:no_go` result for direct bulk
active-translation detection. The package recommends beginning with a labeled
synthetic peptide, advancing to purified stalled RNCs, and attempting active
translation only after the earlier gates pass. `detection_limit_M` is not a
universal constant; replace it with a measured limit for the chosen reporter,
path length, pulse sequence, sample background, and acquisition protocol.

The biological reaction-time scan is distinct from the internal 2D-IR waiting
time:

```julia
t2 = design_pulse_sequence((0.0, 2.0); step=0.02)       # picosecond-scale axis
reaction = design_reaction_scan((0.0, 10.0); step=0.5) # seconds-scale axis
```

See `docs/experimental_feasibility.md` for the staged validation program and
pre-registered go/no-go criteria.

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
