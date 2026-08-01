const _AA3 = Dict('A'=>"ALA",'R'=>"ARG",'N'=>"ASN",'D'=>"ASP",'C'=>"CYS",
    'Q'=>"GLN",'E'=>"GLU",'G'=>"GLY",'H'=>"HIS",'I'=>"ILE",'L'=>"LEU",
    'K'=>"LYS",'M'=>"MET",'F'=>"PHE",'P'=>"PRO",'S'=>"SER",'T'=>"THR",
    'W'=>"TRP",'Y'=>"TYR",'V'=>"VAL")

"""One model-derived backbone atom. Coordinates are in ångström."""
struct AtomCoordinate
    residue_index::Int
    residue_name::String
    atom_name::String
    xyz_A::NTuple{3,Float64}
    chain_id::Char
    confidence::Symbol
    warnings::Vector{String}
end

"""Diagnostics for an approximate candidate backbone, not a native-structure validation."""
struct StructureDiagnostics
    radius_of_gyration_A::Float64
    contact_map::BitMatrix
    secondary_structure::Vector{Symbol}
    clashes::Vector{Tuple{Int,Int,Float64}}
    warnings::Vector{String}
end

"""A model-derived candidate structure and its reconstruction diagnostics."""
struct CandidateStructure
    sequence::String
    atoms::Vector{AtomCoordinate}
    diagnostics::StructureDiagnostics
    metadata::Dict{String,Any}
end

"""Explicit backbone torsions in degrees. This is separate from `Trajectory.dihedral`."""
struct BackboneTorsions
    phi::Vector{Float64}
    psi::Vector{Float64}
    omega::Vector{Float64}
    function BackboneTorsions(phi, psi, omega)
        p, q, o = collect(Float64,phi), collect(Float64,psi), collect(Float64,omega)
        length(p)==length(q)==length(o) || throw(DimensionMismatch("phi, psi, and omega lengths must match"))
        all(isfinite, vcat(p,q,o)) || throw(ArgumentError("torsions must be finite degrees"))
        all(x -> -360 <= x <= 360, vcat(p,q,o)) || throw(ArgumentError("torsions must lie in [-360, 360] degrees"))
        new(p,q,o)
    end
end

"""Approximate torsions plus warnings and a per-residue uncertainty estimate in degrees."""
struct TorsionMapping
    torsions::BackboneTorsions
    uncertainty_degrees::Vector{Float64}
    warnings::Vector{String}
end

_v(a) = Float64[a...]
_tuple(v) = (v[1],v[2],v[3])

# Place D from A-B-C and the B-C-D angle and A-B-C-D dihedral.
function _place(a,b,c,length_A,angle_deg,dihedral_deg)
    bc=(_v(b)-_v(c)); bc ./= norm(bc)
    normal=cross(_v(b)-_v(a), _v(c)-_v(b))
    if norm(normal) < 1e-12
        normal=[0.0,0.0,1.0]
    else
        normal ./= norm(normal)
    end
    perp=cross(normal,bc)
    θ=deg2rad(angle_deg); τ=deg2rad(dihedral_deg)
    direction=cos(θ).*bc + sin(θ).*(cos(τ).*perp + sin(τ).*normal)
    _tuple(_v(c) + length_A.*direction)
end

function _secondary(phi,psi)
    (-100 <= phi <= -30 && -80 <= psi <= -5) && return :helix
    (-180 <= phi <= -70 && (psi >= 80 || psi <= -150)) && return :strand
    :coil
end

"""Reconstruct N, Cα, C, and O atoms using standard idealized peptide geometry.

`phi`, `psi`, and `omega` are degrees and contain one value per residue. The
result is an underdetermined, idealized candidate structure and is not a
validated native protein structure.
"""
function reconstruct_backbone(sequence::AbstractString, phi, psi; omega=180.0, chain_id::Char='A', clash_cutoff_A::Real=1.2)
    seq=uppercase(String(sequence)); n=length(seq)
    n>0 || throw(ArgumentError("sequence must not be empty"))
    all(haskey(_AA3,x) for x in seq) || throw(ArgumentError("sequence contains unsupported amino-acid codes"))
    p,q=collect(Float64,phi),collect(Float64,psi)
    length(p)==n && length(q)==n || throw(DimensionMismatch("phi and psi must each match sequence length"))
    o=omega isa Real ? fill(Float64(omega),n) : collect(Float64,omega)
    length(o)==n || throw(DimensionMismatch("omega must be scalar or match sequence length"))
    torsions=BackboneTorsions(p,q,o)
    isfinite(clash_cutoff_A) && clash_cutoff_A>0 || throw(ArgumentError("clash_cutoff_A must be finite and positive"))
    Ns=Vector{NTuple{3,Float64}}(undef,n); CAs=similar(Ns); Cs=similar(Ns); Os=similar(Ns)
    Ns[1]=(0.0,0.0,0.0); CAs[1]=(1.458,0.0,0.0)
    # N-Cα-C angle 111.2°.
    Cs[1]=(1.458-1.525*cosd(68.8), 1.525*sind(68.8), 0.0)
    Os[1]=_place(Ns[1],CAs[1],Cs[1],1.231,120.8,q[1]+180)
    for i in 2:n
        Ns[i]=_place(Ns[i-1],CAs[i-1],Cs[i-1],1.329,116.2,q[i-1])
        CAs[i]=_place(CAs[i-1],Cs[i-1],Ns[i],1.458,121.7,o[i-1])
        Cs[i]=_place(Cs[i-1],Ns[i],CAs[i],1.525,111.2,p[i])
        Os[i]=_place(Ns[i],CAs[i],Cs[i],1.231,120.8,q[i]+180)
    end
    atoms=AtomCoordinate[]
    for i in 1:n, (name,xyz) in zip(("N","CA","C","O"),(Ns[i],CAs[i],Cs[i],Os[i]))
        push!(atoms,AtomCoordinate(i,_AA3[seq[i]],name,xyz,chain_id,:model_derived,String[]))
    end
    clashes=Tuple{Int,Int,Float64}[]
    for i in eachindex(atoms), j in i+1:length(atoms)
        abs(atoms[i].residue_index-atoms[j].residue_index)<=1 && continue
        d=norm(_v(atoms[i].xyz_A)-_v(atoms[j].xyz_A))
        d<clash_cutoff_A && push!(clashes,(i,j,d))
    end
    cm=reduce(+,(_v(x) for x in CAs))./n
    rg=sqrt(sum(sum(abs2,_v(x)-cm) for x in CAs)/n)
    contacts=falses(n,n)
    for i in 1:n, j in i+3:n
        norm(_v(CAs[i])-_v(CAs[j]))<8.0 && (contacts[i,j]=contacts[j,i]=true)
    end
    warnings=["Coordinates use idealized bond geometry and have not been atomistically validated."]
    !isempty(clashes) && push!(warnings,"Detected $(length(clashes)) nonlocal atom clash(es); refinement is required.")
    diagnostics=StructureDiagnostics(rg,contacts,[_secondary(p[i],q[i]) for i=1:n],clashes,warnings)
    CandidateStructure(seq,atoms,diagnostics,Dict("coordinate_units"=>"angstrom", "source"=>"explicit_backbone_torsions", "model_derived"=>true))
end

"""Map one reduced trajectory frame to explicit torsions with stated uncertainty.

The reduced, dimensionless `dihedral` state is *not* a φ or ψ measurement. For
visualization only, this mapping perturbs an α-like reference as
`φ=-60+15θ`, `ψ=-45-15θ`, and sets trans peptide `ω=180` degrees.
"""
function approximate_backbone_torsions(traj::Trajectory, frame=:last)
    f=frame===:last ? length(traj.time) : frame===:first ? 1 : Int(frame)
    checkbounds(traj.time,f); θ=traj.dihedral[:,f]
    t=BackboneTorsions(-60 .+ 15 .* θ, -45 .- 15 .* θ, fill(180.0,length(θ)))
    TorsionMapping(t,fill(45.0,length(θ)),["Reduced dihedral is dimensionless and does not uniquely determine φ and ψ; the mapping is approximate."])
end

"""Reconstruct a candidate structure from a `Trajectory` using an explicit approximate mapping."""
function reconstruct_structure(traj::Trajectory; frame=:last, torsion_mapping=:approximate)
    torsion_mapping===:approximate || throw(ArgumentError("only the explicit :approximate mapping is supported"))
    m=approximate_backbone_torsions(traj,frame); t=m.torsions
    s=reconstruct_backbone(traj.sequence,t.phi,t.psi;omega=t.omega)
    append!(s.diagnostics.warnings,m.warnings); s.metadata["source"]="reduced_trajectory_approximate_mapping"
    s.metadata["torsion_uncertainty_degrees"]=m.uncertainty_degrees; s
end

"""Write a model-derived candidate structure in fixed-column PDB format."""
function write_pdb(path::AbstractString, structure::CandidateStructure)
    open(path,"w") do io
        println(io,"REMARK 250 MODEL-DERIVED CANDIDATE; NOT A VALIDATED NATIVE STRUCTURE")
        println(io,"REMARK 250 IDEALIZED BACKBONE GEOMETRY; ATOMISTIC VALIDATION REQUIRED")
        for (serial,a) in enumerate(structure.atoms)
            x,y,z=a.xyz_A; element=first(a.atom_name)
            @printf(io,"ATOM  %5d %-4s %-3s %c%4d    %8.3f%8.3f%8.3f  1.00  0.00           %c\n",serial,a.atom_name,a.residue_name,a.chain_id,a.residue_index,x,y,z,element)
        end
        println(io,"TER"); println(io,"END")
    end
    String(path)
end

"""Write a model-derived candidate structure as an mmCIF atom-site table."""
function write_mmcif(path::AbstractString, structure::CandidateStructure)
    open(path,"w") do io
        println(io,"data_davydovfoldon_candidate")
        println(io,"_struct.title 'Model-derived candidate; not a validated native structure'")
        println(io,"loop_")
        for h in ("_atom_site.group_PDB","_atom_site.id","_atom_site.type_symbol","_atom_site.label_atom_id","_atom_site.label_comp_id","_atom_site.label_asym_id","_atom_site.label_seq_id","_atom_site.Cartn_x","_atom_site.Cartn_y","_atom_site.Cartn_z") println(io,h) end
        for (serial,a) in enumerate(structure.atoms)
            x,y,z=a.xyz_A
            @printf(io,"ATOM %d %c %s %s %c %d %.3f %.3f %.3f\n",serial,first(a.atom_name),a.atom_name,a.residue_name,a.chain_id,a.residue_index,x,y,z)
        end
        println(io,"#")
    end
    String(path)
end
