@testset "2D-IR" begin
    traj, _, _ = run_soliton("MAGL"; tmax=0.04, dt=0.01)
    h = build_hamiltonian(traj, 1)
    @test h ≈ h'
    spectra = compute_2DIR_spectrum(traj, SpectrumParams(t2_values=collect(0:0.1:0.4), coherence_points=8))
    @test size(spectra.signal) == (8, 8, 5)
    @test extract_beat_frequency(spectra) >= 0
    @test_throws DimensionMismatch Spectra([0.0], [0.0], [0.0], zeros(2, 1, 1))
    irregular = Spectra([0.0], [0.0], [0.0, 0.1, 0.3], zeros(1, 1, 3))
    @test_throws ArgumentError extract_beat_frequency(irregular)
    @test_throws ArgumentError compute_2DIR_spectrum(
        traj, SpectrumParams(t2_values=[-0.1, 0.0], coherence_points=8))
    @test_throws ArgumentError compute_2DIR_spectrum(
        traj, SpectrumParams(coherence_dt=Inf, coherence_points=8))
end
