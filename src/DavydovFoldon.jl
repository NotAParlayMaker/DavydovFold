module DavydovFoldon

using LinearAlgebra, Random, Statistics, DelimitedFiles, Printf
import FFTW

include("utils.jl")
include("SolitonSim.jl")
include("Sol2DIR.jl")
include("FoldonExperiment.jl")
include("FoldonDB.jl")
include("StructurePrediction.jl")
include("AtomisticExport.jl")

export DavydovParams, Trajectory, aa_to_params, run_soliton, save_trajectory,
       contact_map, animate_soliton, SpectrumParams, Spectra, build_hamiltonian,
       compute_2DIR_spectrum, compute_2DIR, extract_beat_frequency,
       predict_signal_amplitude, design_pulse_sequence, analyze_beats,
       ExperimentTimeAxis, CoherenceTime, WaitingTime, DetectionTime, BiologicalReactionTime,
       PopulationBudget, population_budget, feasibility_scenario, recommend_stage,
       UniformRange, monte_carlo_feasibility, sensitivity_analysis,
       codon_optimizer, open_database, store_simulation!, serve
export AtomCoordinate, CandidateStructure, StructureDiagnostics, BackboneTorsions,
       TorsionMapping, reconstruct_backbone, reconstruct_structure,
       approximate_backbone_torsions, write_pdb, write_mmcif,
       OpenMMSpec, prepare_openmm_system, write_openmm_spec

end
