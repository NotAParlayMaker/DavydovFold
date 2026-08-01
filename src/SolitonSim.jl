"""Integrate the stochastic extended Davydov–Scott model.

The lattice uses velocity Verlet, the exciton a norm-preserving Cayley
(Crank–Nicolson) step, and the overdamped torsions Euler–Maruyama.
"""
function run_soliton(seq::AbstractString, params::DavydovParams=DavydovParams();
                     tmax::Real=20.0, dt::Real=0.01, T::Real=310.0)
    isempty(seq) && throw(ArgumentError("sequence must not be empty"))
    dt > 0 && tmax >= 0 && T >= 0 || throw(ArgumentError("dt must be positive; tmax and T must be nonnegative"))
    params.mass > 0 && params.output_stride > 0 || throw(ArgumentError("mass and output_stride must be positive"))
    tmax64, dt64, temperature = Float64(tmax), Float64(dt), Float64(T)
    sequence=uppercase(String(seq))
    all(isletter, sequence) || throw(ArgumentError("sequence must contain one-letter amino-acid codes"))
    N=length(sequence); rng=MersenneTwister(params.seed)
    factors=aa_to_params.(collect(sequence)); J=params.J .* first.(factors)
    chi=params.chi .* getindex.(factors,2); xi=params.xi .* last.(factors)
    B=ComplexF64[exp(-0.5*((n-min(5,N))/2)^2) for n in 1:N]; B ./= norm(B)
    u=zeros(N); p=sqrt(params.mass*params.kB*temperature).*randn(rng,N); p .-= mean(p)
    theta=zeros(N); stride=params.output_stride; steps=floor(Int,tmax64/dt64)
    saved=0:stride:steps; times=collect(saved).*dt64; rho=zeros(N,length(saved))
    us=similar(rho); th=similar(rho); k=1
    rho[:,k].=abs2.(B); us[:,k].=u; th[:,k].=theta
    force=zeros(N)
    function lattice_force!()
        @inbounds for n in 1:N
            force[n]=params.spring*(right(u,n)-2u[n]+left(u,n)) +
                     chi[n]*(abs2(right(B,n))-abs2(left(B,n))) - params.gamma*p[n]
        end
    end
    Iden=Matrix{ComplexF64}(I,N,N)
    for step in 1:steps
        lattice_force!(); p .+= 0.5dt64 .* force
        u .+= dt64/params.mass .* p
        H=zeros(ComplexF64,N,N)
        @inbounds for n in 1:N
            H[n,n]=params.E0+chi[n]*(right(u,n)-left(u,n))+xi[n]*theta[n]
            n<N && (H[n,n+1]=H[n+1,n]=-0.5*(J[n]+J[n+1]))
        end
        B .= (Iden + 0.5im*dt64*H) \ ((Iden - 0.5im*dt64*H)*B)
        lattice_force!()
        p .+= 0.5dt64 .* force .+ sqrt(2params.gamma*params.kB*temperature*params.mass*dt64).*randn(rng,N)
        @inbounds for n in 1:N
            grad=4params.epsilon*theta[n]*(theta[n]^2-params.theta0^2)
            lap=right(theta,n)-2theta[n]+left(theta,n)
            theta[n] += dt64*(-grad+params.K*lap+xi[n]*abs2(B[n])-params.gamma_theta*theta[n]) +
                        sqrt(2params.gamma_theta*params.kB*temperature*dt64)*randn(rng)
        end
        if step % stride == 0
            k+=1; rho[:,k].=abs2.(B); us[:,k].=u; th[:,k].=theta
        end
    end
    traj=Trajectory(times,rho,us,th,sequence,params)
    (traj, copy(theta), abs2.(B))
end

"""Create an interactive exciton kymograph (Makie is loaded on demand)."""
function animate_soliton(traj::Trajectory)
    try
        @eval import CairoMakie
    catch
        throw(ArgumentError("install CairoMakie to plot trajectories"))
    end
    fig=CairoMakie.Figure(); ax=CairoMakie.Axis(fig[1,1],xlabel="time (ps)",ylabel="site")
    CairoMakie.heatmap!(ax,traj.time,1:size(traj.exciton_density,1),permutedims(traj.exciton_density))
    fig
end
