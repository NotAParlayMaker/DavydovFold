"""DavydovFoldon: stochastic Davydov–Scott simulations and predicted 2D-IR analysis."""
module DavydovFoldon

using LinearAlgebra, Random, Statistics, DelimitedFiles
import FFTW

include("utils.jl")
include("SolitonSim.jl")
include("Sol2DIR.jl")
include("FoldonExperiment.jl")
include("FoldonDB.jl")

export DavydovParams, Trajectory, aa_to_params, run_soliton, save_trajectory,
       contact_map, animate_soliton, SpectrumParams, Spectra, build_hamiltonian,
       compute_2DIR_spectrum, compute_2DIR, extract_beat_frequency,
       predict_signal_amplitude, design_pulse_sequence, analyze_beats,
       codon_optimizer, open_database, store_simulation!, serve

end
