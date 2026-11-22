@testset "experiment design" begin
    @test predict_signal_amplitude(2.0, 3.0).snr > 0
    @test length(design_pulse_sequence((0.0, 1.0); step=0.1).waiting_times) == 11
    @test length(design_reaction_scan((0.0, 2.0); step=0.5).reaction_times) == 5
    @test codon_optimizer("MAG", [1, -1, 1]) == ["ATG", "GCG", "GGT"]
    @test_throws DimensionMismatch codon_optimizer("MA", [1])
    @test_throws ArgumentError design_pulse_sequence((0.0, 1.0); step=0)
    @test_throws ArgumentError predict_signal_amplitude(1, 1; noise_floor=0)
    @test length(codon_optimizer("ACDEFGHIKLMNPQRSTVWY", ones(20))) == 20
    @test_throws ArgumentError codon_optimizer("X", [1])
    @test codon_optimizer("WR", [1, -1]) == ["TGG", "AGG"]
    @test_throws ArgumentError codon_optimizer("X", [1])
    @test_throws ArgumentError design_pulse_sequence((0.0, 1.0); step=Inf)
    @test_throws ArgumentError design_reaction_scan((0.0, 1.0); step=0.0)

    baseline = FeasibilityParams()
    population = estimate_transient_population(baseline)
    @test population.occupancy_fraction < 1e-8
    @test population.reporter_concentration_M < 1e-12
    @test !population.feasible
    @test isinf(population.required_state_lifetime_s)

    assessment = assess_2dir_feasibility(baseline)
    @test assessment.status == :no_go
    @test assessment.recommended_stage == :synthetic_peptide
    @test !isempty(assessment.recommendations)

    optimistic = FeasibilityParams(rnc_concentration_M=1e-3,
                                   elongation_rate_s=20.0,
                                   state_lifetime_s=100e-6,
                                   detection_limit_M=1e-6,
                                   labels_per_rnc=1)
    @test assess_2dir_feasibility(optimistic).status == :go
    @test length(staged_validation_plan()) == 3
    @test staged_validation_plan()[1].stage == :synthetic_peptide

    scaled = predict_signal_amplitude(1e-6, 2.0; active_fraction=0.25, labels_per_rnc=2)
    @test scaled.reporter_concentration == 5e-7
    @test_throws ArgumentError predict_signal_amplitude(1.0, 1.0; active_fraction=1.5)
    @test_throws ArgumentError FeasibilityParams(labels_per_rnc=0) |> estimate_transient_population

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
