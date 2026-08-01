Base.@kwdef struct SpectrumParams
    t2_values::Vector{Float64}=collect(0.0:0.1:5.0)
    coherence_points::Int=64
    coherence_dt::Float64=0.025
    relaxation::Float64=0.5
end
"""A real absorptive 2D-IR signal sampled on `(omega1, omega3, t2)`."""
struct Spectra
    omega1::Vector{Float64}
    omega3::Vector{Float64}
    t2::Vector{Float64}
    signal::Array{Float64,3}
    function Spectra(omega1::Vector{Float64}, omega3::Vector{Float64},
                     t2::Vector{Float64}, signal::Array{Float64,3})
        size(signal) == (length(omega1), length(omega3), length(t2)) ||
            throw(DimensionMismatch("signal dimensions must match its coordinate axes"))
        new(omega1, omega3, t2, signal)
    end
end

Spectra(omega1::AbstractVector{<:Real}, omega3::AbstractVector{<:Real},
        t2::AbstractVector{<:Real}, signal::AbstractArray{<:Real,3}) =
    Spectra(collect(Float64, omega1), collect(Float64, omega3),
            collect(Float64, t2), Array{Float64,3}(signal))

"""Construct the instantaneous real symmetric one-exciton Hamiltonian."""
function build_hamiltonian(traj::Trajectory, frame::Integer)
    checkbounds(traj.time, frame)
    N=length(traj.sequence); p=traj.params; H=zeros(Float64,N,N)
    f=aa_to_params.(collect(traj.sequence))
    for n in 1:N
        H[n,n]=p.E0+p.chi*f[n][2]*(right(view(traj.displacement,:,frame),n)-left(view(traj.displacement,:,frame),n))+p.xi*f[n][3]*traj.dihedral[n,frame]
        n<N && (H[n,n+1]=H[n+1,n]=-p.J*(f[n][1]+f[n+1][1])/2)
    end
    H
end

"""Compute a deterministic semi-impulsive absorptive 2D-IR spectrum."""
function compute_2DIR_spectrum(traj::Trajectory, sp::SpectrumParams=SpectrumParams())
    sp.coherence_points >= 2 || throw(ArgumentError("coherence_points must be at least 2"))
    isfinite(sp.coherence_dt) && sp.coherence_dt > 0 ||
        throw(ArgumentError("coherence_dt must be finite and positive"))
    isfinite(sp.relaxation) && sp.relaxation >= 0 ||
        throw(ArgumentError("relaxation must be finite and nonnegative"))
    isempty(sp.t2_values) && throw(ArgumentError("t2_values must not be empty"))
    all(isfinite, sp.t2_values) || throw(ArgumentError("t2_values must contain only finite values"))
    all(>=(0), sp.t2_values) || throw(ArgumentError("t2_values must be nonnegative"))
    issorted(sp.t2_values) || throw(ArgumentError("t2_values must be sorted"))
    m=sp.coherence_points; dt=sp.coherence_dt
    corr=zeros(ComplexF64,m,m,length(sp.t2_values)); nframes=length(traj.time)
    for (q,t2) in pairs(sp.t2_values)
        frame=clamp(searchsortedfirst(traj.time,t2),1,nframes)
        vals=eigvals(Symmetric(build_hamiltonian(traj,frame)))
        for i in 1:m, j in 1:m
            corr[i,j,q]=sum(exp.(-im.*vals.*((i-j)*dt))) / length(vals) * exp(-sp.relaxation*(i+j-2)*dt)
        end
    end
    raw=FFTW.fftshift(FFTW.fft(corr,(1,2)),(1,2))
    omega=2pi .* FFTW.fftshift(FFTW.fftfreq(m, inv(dt)))
    Spectra(omega,omega,sp.t2_values,real.(raw))
end
"""Convenience wrapper for computing a spectrum at selected waiting times."""
compute_2DIR(traj::Trajectory, t2_values) =
    compute_2DIR_spectrum(traj, SpectrumParams(t2_values=collect(Float64, t2_values)))

"""Return the dominant non-zero waiting-time beat frequency in THz."""
function extract_beat_frequency(s::Spectra)
    length(s.t2)<3 && throw(ArgumentError("at least three waiting times are required"))
    dt = sampling_interval(s.t2; name="waiting times")
    trace=vec(sum(abs.(s.signal),dims=(1,2))); trace .-= mean(trace)
    power=abs.(FFTW.rfft(trace)); freqs=FFTW.rfftfreq(length(trace),inv(dt))
    freqs[argmax(power[2:end])+1]
end
