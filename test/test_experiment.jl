@testset "experiment design" begin
    @test predict_signal_amplitude(2.0, 3.0).snr > 0
    @test length(design_pulse_sequence((0.0, 1.0); step=0.1).waiting_times) == 11
    @test length(design_reaction_scan((0.0, 2.0); step=0.5).reaction_times) == 5
    @test codon_optimizer("MAG", [1, -1, 1]) == ["ATG", "GCG", "GGT"]
    @test_throws DimensionMismatch codon_optimizer("MA", [1])
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
