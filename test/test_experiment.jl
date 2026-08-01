@testset "experiment design" begin
    @test predict_signal_amplitude(2.0, 3.0).snr > 0
    @test length(design_pulse_sequence((0.0, 1.0); step=0.1).waiting_times) == 11
    @test codon_optimizer("MAG", [1, -1, 1]) == ["ATG", "GCG", "GGT"]
    @test_throws DimensionMismatch codon_optimizer("MA", [1])
    @test codon_optimizer("WR", [1, -1]) == ["TGG", "AGG"]
    @test_throws ArgumentError codon_optimizer("X", [1])
    @test_throws ArgumentError design_pulse_sequence((0.0, 1.0); step=Inf)

    mktemp() do path, io
        writedlm(io, [0.0 0.0; 0.1 1.0; 0.25 0.0; 0.3 -1.0], ',')
        close(io)
        @test_throws ArgumentError analyze_beats(path)
    end
    @test_throws ArgumentError analyze_beats("unused.csv"; predicted_beat_freq=Inf)
end

@testset "population budgets and stages" begin
    p = PopulationBudget(rnc_concentration_mol_L=1e-6, launch_rate_s_inv=10.0,
        state_lifetime_s=0.02, synchronized_fraction=0.5, competent_fraction=0.5,
        reporters_per_active_state=2, labeling_efficiency=0.8, detection_threshold_mol_L=1e-7)
    r = population_budget(p)
    @test r.active_fraction ≈ 0.2
    @test r.active_state_concentration_mol_L ≈ 5e-8
    @test r.effective_reporter_concentration_mol_L ≈ 8e-8
    @test r.required_state_lifetime_s ≈ 0.025
    @test population_budget(PopulationBudget(rnc_concentration_mol_L=1e-6,
        launch_rate_s_inv=10, state_lifetime_s=1, detection_threshold_mol_L=1e-7)).active_fraction == 1
    for rate_life in ((0.0,1.0),(1.0,0.0))
        z=PopulationBudget(rnc_concentration_mol_L=1e-6, launch_rate_s_inv=rate_life[1],
            state_lifetime_s=rate_life[2], detection_threshold_mol_L=1e-7)
        @test population_budget(z).active_state_concentration_mol_L == 0
        @test isinf(population_budget(z).required_state_lifetime_s)
    end
    @test_throws ArgumentError population_budget(PopulationBudget(rnc_concentration_mol_L=-1,
        launch_rate_s_inv=1, state_lifetime_s=1, detection_threshold_mol_L=1e-6))
    @test_throws ArgumentError population_budget(PopulationBudget(rnc_concentration_mol_L=1,
        launch_rate_s_inv=1, state_lifetime_s=1, synchronized_fraction=1.01, detection_threshold_mol_L=1e-6))
    @test population_budget(feasibility_scenario(:active_translation)).feasibility_status == :no_go
    @test population_budget(feasibility_scenario(:synthetic_peptide)).feasibility_status == :go
    @test recommend_stage(population_budget(feasibility_scenario(:active_translation))).stage == :A
    stalled=feasibility_scenario(:stalled_rnc)
    @test recommend_stage(population_budget(stalled); inputs=stalled).stage == :C
    @test PopulationBudget(rnc_concentration_mol_L=100e-9, launch_rate_s_inv=1,
        state_lifetime_s=1, detection_threshold_mol_L=1e-6).rnc_concentration_mol_L == 1e-7
end

@testset "time axes, uncertainty, and beat evidence" begin
    @test CoherenceTime([0,1e-15]).values_s[end] == 1e-15
    @test WaitingTime([0,1e-12]) isa ExperimentTimeAxis{:waiting}
    @test_throws ArgumentError DetectionTime([0.0,0.0])
    @test_throws ArgumentError BiologicalReactionTime([1.0,0.5])
    args=(rnc_concentration_mol_L=UniformRange(1e-7,2e-7), launch_rate_s_inv=10.0,
        state_lifetime_s=UniformRange(1e-3,2e-3), synchronized_fraction=0.5,
        competent_fraction=0.5, reporters_per_active_state=2.0, labeling_efficiency=0.8,
        detection_threshold_mol_L=1e-9)
    m1=monte_carlo_feasibility(; samples=100, seed=7, args...)
    m2=monte_carlo_feasibility(; samples=100, seed=7, args...)
    @test m1 == m2
    @test m1.credible_interval_mol_L[1] <= m1.median_active_concentration_mol_L <= m1.credible_interval_mol_L[2]
    @test 0 <= m1.probability_go <= 1
    @test length(sensitivity_analysis(feasibility_scenario(:active_translation);
        rnc_concentrations_mol_L=[1e-9,1e-6], number_of_averages=[1,4])) == 4
    mktemp() do path, io
        t=collect(0.0:0.1:1.0); writedlm(io,hcat(t,sin.(2pi.*t)),','); close(io)
        a=analyze_beats(path; predicted_beat_freq=1.0)
        @test a.candidate_frequency >= 0
        @test a.warning !== nothing
        @test_throws ArgumentError analyze_beats(path; predicted_beat_freq=6.0)
    end
end
