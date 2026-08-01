@testset "soliton" begin
    p = DavydovParams(output_stride=2, seed=7, gamma=0.0, gamma_theta=0.0)
    traj, theta, rho = run_soliton("MAGL", p; tmax=0.04, dt=0.01, T=0.0)
    @test size(traj.exciton_density) == (4, 3)
    @test sum(rho) ≈ 1 atol=1e-10
    @test all(isfinite, theta)
    @test aa_to_params('X') == (1.0, 1.0, 1.0)
    @test contact_map(theta) isa BitMatrix
    @test_throws ArgumentError run_soliton("MAGL"; T=-1)
    @test_throws ArgumentError contact_map(theta; cutoff=0)
end
