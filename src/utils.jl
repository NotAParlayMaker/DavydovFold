"""Physical and numerical parameters in scaled ps/Å/amu units.

All fields are concrete so a parameter object can be passed through the hot
integration loop without dynamic dispatch. `output_stride` controls the number
of integration steps between stored frames and `seed` makes thermal runs
reproducible.
"""
Base.@kwdef struct DavydovParams
    E0::Float64 = 0.0
    J::Float64 = 1.0
    chi::Float64 = 0.35
    xi::Float64 = 0.12
    mass::Float64 = 114.0
    spring::Float64 = 13.0
    gamma::Float64 = 0.08
    gamma_theta::Float64 = 0.15
    theta0::Float64 = 1.0
    epsilon::Float64 = 0.8
    K::Float64 = 0.25
    kB::Float64 = 0.008314462618
    seed::Int = 42
    output_stride::Int = 10
end

const AA_CLASS = Dict(
    'G'=>(1.08,0.82,0.72), 'P'=>(0.72,1.18,1.35), 'A'=>(1.04,0.94,0.86),
    'V'=>(0.96,1.08,1.02), 'I'=>(0.94,1.10,1.02), 'L'=>(0.95,1.09,1.01),
    'F'=>(0.90,1.16,1.12), 'W'=>(0.86,1.22,1.20), 'Y'=>(0.89,1.17,1.16),
    'S'=>(1.02,0.96,1.07), 'T'=>(0.99,1.00,1.08), 'C'=>(0.98,1.04,1.04),
    'M'=>(0.94,1.08,1.05), 'N'=>(1.00,1.02,1.12), 'Q'=>(0.98,1.04,1.10),
    'D'=>(1.03,0.91,1.18), 'E'=>(1.01,0.93,1.16), 'K'=>(0.97,1.01,1.20),
    'R'=>(0.95,1.04,1.24), 'H'=>(0.96,1.06,1.20))

"""Return `(J, χ, ξ)` multipliers for a one-letter amino-acid code."""
aa_to_params(aa::Char)::NTuple{3,Float64} = get(AA_CLASS, uppercase(aa), (1.0, 1.0, 1.0))

"""A saved simulation trajectory (sites × time points for each field)."""

struct Trajectory
    time::Vector{Float64}
    exciton_density::Matrix{Float64}
    displacement::Matrix{Float64}
    dihedral::Matrix{Float64}
    sequence::String
    params::DavydovParams
end

@inline left(v, n) = n == 1 ? v[1] : v[n-1]
@inline right(v, n) = n == length(v) ? v[end] : v[n+1]

"""Binary proximity map from a simple planar backbone embedding."""
function contact_map(theta::AbstractVector{<:Real}; cutoff::Real=2.2)
    cutoff > 0 || throw(ArgumentError("cutoff must be positive"))
    N = length(theta); xy = zeros(Float64, 2, N)
    angle = 0.0
    for n in 2:N
        angle += theta[n-1]
        xy[:,n] .= xy[:,n-1] .+ (cos(angle), sin(angle))
    end
    C = falses(N,N)
    for j in 1:N, i in 1:j-3
        C[i,j] = C[j,i] = sum(abs2, xy[:,i] - xy[:,j]) < cutoff^2
    end
    C
end

"""Write `traj` and its final contact map to an HDF5 file at `path`."""
function save_trajectory(path::AbstractString, traj::Trajectory)
    import HDF5
    HDF5.h5open(path, "w") do f
        f["time"] = traj.time; f["exciton_density"] = traj.exciton_density
        f["displacement"] = traj.displacement; f["dihedral"] = traj.dihedral
        f["contact_map"] = UInt8.(contact_map(traj.dihedral[:,end]))
        HDF5.attributes(f)["sequence"] = traj.sequence
    end
    path
end
