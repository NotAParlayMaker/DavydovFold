"""Immutable OpenMM job specification. Every physical quantity is unit-labeled."""
Base.@kwdef struct OpenMMSpec
    input_structure_path::String
    forcefield::String="amber14-all.xml"
    water_model::String="amber14/tip3p.xml"
    temperature_K::Float64=310.0
    pressure_bar::Float64=1.0
    ionic_strength_molar::Float64=0.15
    solvent_padding_nm::Float64=1.0
    integration_timestep_fs::Float64=2.0
    minimization_max_iterations::Int=1000
    equilibration_duration_ps::Float64=10.0
    production_duration_ps::Float64=20.0
    trajectory_output_path::String="trajectory.dcd"
    final_coordinates_path::String="final.pdb"
    run_summary_path::String="run_summary.json"
    random_seed::Int=42
end

function _validate_openmm(s::OpenMMSpec)
    isfile(s.input_structure_path) || throw(ArgumentError("input structure does not exist: $(s.input_structure_path)"))
    for (name,x,positive) in ((:temperature_K,s.temperature_K,true),(:pressure_bar,s.pressure_bar,false),
        (:ionic_strength_molar,s.ionic_strength_molar,false),(:solvent_padding_nm,s.solvent_padding_nm,true),
        (:integration_timestep_fs,s.integration_timestep_fs,true),(:equilibration_duration_ps,s.equilibration_duration_ps,false),
        (:production_duration_ps,s.production_duration_ps,false))
        isfinite(x) && (positive ? x>0 : x>=0) || throw(ArgumentError("$name has an invalid value"))
    end
    s.minimization_max_iterations>0 || throw(ArgumentError("minimization_max_iterations must be positive")); s
end

"""Prepare an external OpenMM specification for an already exported candidate structure.

OpenMM remains optional and is never imported or executed by this package. A
force field is approximate; short trajectory survival does not establish a
Davydov mechanism or thermodynamic stability.
"""
function prepare_openmm_system(structure::CandidateStructure; input_structure_path::AbstractString, kwargs...)
    s=OpenMMSpec(;input_structure_path=String(input_structure_path),kwargs...); _validate_openmm(s)
end

"""Prepare an OpenMM specification from a PDB or mmCIF path."""
function prepare_openmm_system(input_structure_path::AbstractString; kwargs...)
    s=OpenMMSpec(;input_structure_path=String(input_structure_path),kwargs...); _validate_openmm(s)
end

"""Serialize an `OpenMMSpec` to JSON for the optional Python runner."""
function write_openmm_spec(path::AbstractString, spec::OpenMMSpec)
    _validate_openmm(spec); import JSON3
    fields=propertynames(spec); payload=Dict(String(k)=>getproperty(spec,k) for k in fields)
    payload["scientific_limitations"]=["Force fields are approximations.","Short MD does not prove thermodynamic stability.","Conformation survival does not establish a Davydov mechanism.","Failure may reflect reconstruction or force-field limitations."]
    open(path,"w") do io; JSON3.pretty(io,payload); end
    String(path)
end
